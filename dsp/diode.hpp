/*
    Copyright (C) Tim Blechmann

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



    adapted from julian parker's diode ringmod

 */

#ifndef DIODE_HPP
#define DIODE_HPP

#include <tuple>

#include <boost/simd/operator/operator.hpp>
#include <boost/simd/include/functions/fast_divides.hpp>
#include <boost/simd/include/functions/fast_rec.hpp>
#include <boost/simd/include/functions/if_else.hpp>
#include <boost/simd/include/functions/if_else_zero.hpp>

//#include <boost/simd/constant/constants/one.hpp>
//#include <boost/simd/constant/constants/zero.hpp>

namespace nova {

namespace detail {

template <typename T>
auto double_( T arg )
{
    return arg + arg;
}

template <typename T>
auto square_( T arg )
{
    return arg * arg;
}

}

template <typename T>
auto calcDiodeCoefficients( T vB, T vL, T gain )
{
    using namespace detail;

    const T squareSection = vL - vB;

    const T curveCoeff = boost::simd::fast_divides( gain, double_( squareSection ) );
    const T linCoeff   = curveCoeff * detail::square_( squareSection ) - gain * vL;
    return std::make_tuple( linCoeff, curveCoeff );
}

template <typename SampleType, typename ArgType >
SampleType diodeTransferFunction( SampleType input, ArgType linCoeff, ArgType curveCoeff,
                                  ArgType vB, ArgType vL, ArgType gain )
{
    using namespace detail;
    SampleType linearSection    = gain * input + linCoeff;
    SampleType quadraticSection = curveCoeff * square_( input - vB );

    SampleType nonZeroSection   = boost::simd::if_else( input < vL, quadraticSection, linearSection );
    SampleType ret              = boost::simd::if_else_zero( input > vB, nonZeroSection );
    return ret;
}

// implicit gain: 1
template <typename SampleType, typename ArgType >
SampleType diodeTransferFunction( SampleType input, ArgType linCoeff, ArgType curveCoeff,
                                  ArgType vB, ArgType vL )
{
    using namespace detail;
    SampleType linearSection    = input + linCoeff;
    SampleType quardaticSection = curveCoeff * square_( input - vB );

    SampleType nonZeroSection   = boost::simd::if_else( input < vL, quardaticSection, linearSection );
    SampleType ret              = boost::simd::if_else_zero( input > vB, nonZeroSection );
    return ret;
}


template <typename ArgType>
auto calcDiodeSaturationCoefficients( ArgType vB, ArgType vL )
{
    ArgType linCoeff, curveCoeff;
    std::tie( linCoeff, curveCoeff ) = calcDiodeCoefficients( vB, vL, boost::simd::One<ArgType>() );

    const ArgType scaleFactor        = boost::simd::fast_rec( diodeTransferFunction( ArgType{1.f}, linCoeff, curveCoeff, vB, vL ) );

    return std::make_tuple( linCoeff, curveCoeff, scaleFactor );
}


template <typename SampleType, typename SaturationParameter >
SampleType diodeSaturation( SampleType input, SaturationParameter saturationParameter )
{
    const auto vB          = std::get<0>( saturationParameter );
    const auto vL          = std::get<1>( saturationParameter );
    const auto linCoeff    = std::get<2>( saturationParameter );
    const auto curveCoeff  = std::get<3>( saturationParameter );
    const auto scaleFactor = std::get<4>( saturationParameter );

    const SampleType scaledArgument = SampleType{1} - boost::simd::abs( input );
    const SampleType diodeShaped       = diodeTransferFunction( scaledArgument,
                                                                linCoeff, curveCoeff, vB, vL );

    const SampleType scaled = SampleType{1} - scaleFactor * diodeShaped;

    return boost::simd::copysign( scaled, input );
}

template <typename SampleType, typename ArgType >
SampleType diodeSaturation( SampleType input, ArgType vB, ArgType vL )
{
    ArgType linCoeff, curveCoeff, scaleFactor;
    std::tie( linCoeff, curveCoeff, scaleFactor ) = calcDiodeSaturationCoefficients( vB, vL );

    return diodeSaturation( input, std::make_tuple( vB, vL, linCoeff, curveCoeff, scaleFactor ) );
}

}

#endif
