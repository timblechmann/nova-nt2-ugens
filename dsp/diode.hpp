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
#include <boost/simd/include/functions/if_else.hpp>
#include <boost/simd/include/functions/if_else_zero.hpp>


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
std::tuple<T, T> calcDiodeCoefficients( T vB, T vL, T gain )
{
    using namespace detail;

    const T squareSection = vL - vB;

    const T curveCoeff = boost::simd::fast_divides( vB, double_( squareSection ) );
    const T linCoeff   = curveCoeff * square( squareSection ) - gain * vL;
    return std::make_tuple( linCoeff, curveCoeff );
}

template <typename SampleType, typename ArgType >
SampleType diodeTransferFunction( SampleType input, ArgType curveCoeff, ArgType linCoeff,
                                  ArgType vB, ArgType vL, ArgType gain )
{
    using namespace detail;
    SampleType linearSection    = gain * input + linCoeff;
    SampleType quardaticSection = curveCoeff * square_( input - vB );

    SampleType nonZeroSection   = boost::simd::if_else( input < vL, quardaticSection, linearSection );
    SampleType ret              = boost::simd::if_else_zero( input > vB, nonZeroSection );
    return ret;
}

// implicit gain: 1
template <typename SampleType, typename ArgType >
SampleType diodeTransferFunction( SampleType input, ArgType curveCoeff, ArgType linCoeff,
                                  ArgType vB, ArgType vL )
{
    using namespace detail;
    SampleType linearSection    = input + linCoeff;
    SampleType quardaticSection = curveCoeff * square_( input - vB );

    SampleType nonZeroSection   = boost::simd::if_else( input < vL, quardaticSection, linearSection );
    SampleType ret              = boost::simd::if_else_zero( input > vB, nonZeroSection );
    return ret;
}

template <typename SampleType >
SampleType diodeSaturation( SampleType input )
{



}


#endif
