/*
    Copyright (C) 2015 Tim Blechmann

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
 */

#ifndef SATURATORS_HPP
#define SATURATORS_HPP



#include <boost/simd/include/constants/one.hpp>
#include <boost/simd/include/constants/two.hpp>
#include <boost/simd/include/constants/quarter.hpp>
#include <boost/simd/include/constants/zero.hpp>

#include <boost/simd/arithmetic/include/functions/abs.hpp>
#include <boost/simd/arithmetic/include/functions/raw_rec.hpp>
#include <boost/simd/arithmetic/include/functions/fast_rec.hpp>

#include <boost/simd/operator/include/functions/fast_divides.hpp>

#include <boost/simd/ieee/include/functions/copysign.hpp>

#include <nt2/include/functions/pow_abs.hpp>
#include <nt2/include/functions/tanh.hpp>

#include <dsp/utils.hpp>

#include <dsp/tanh_approximation.hpp>

namespace nova      {
namespace saturator {

template< typename Arg >
auto distort ( Arg sample )
{
    using namespace boost::simd;
    return fast_div( sample, One<Arg>() + abs(sample) );
}

template< typename Arg >
auto raw_distort ( Arg sample )
{
    using namespace boost::simd;
    return sample * raw_rec( One<Arg>() + abs(sample) );
}

template< typename Arg0, typename Arg1 >
auto pow ( Arg0 sig, Arg1 level )
{
	return boost::simd::copysign( nt2::pow_abs( sig, level ), sig );
}

template< typename Arg0, typename Arg1 >
auto hyperbol ( Arg0 sig, Arg1 level )
{
    using namespace boost::simd;

    auto saturated = level - (level * level * fast_rec( level + abs( sig ) ) );
    return copysign( saturated, sig );
}


template< typename Arg0, typename Arg1 >
auto parabol ( Arg0 sig, Arg1 level )
{
    using namespace boost::simd;

    // TODO: can boost.simd optimize this?
//    auto limit         = Two<SampleType>() * level;
    auto limit         = level + level;
    auto clippedSignal = nova::clip2( sig, limit );
    auto factor        = One<Arg0>() - ( abs( clippedSignal ) * Quarter<Arg0>() * fast_rec( level ) );

    return clippedSignal * factor;
}


template< typename Arg0, typename Arg1 >
auto tanh_saturator( Arg0 sig, Arg1 preGain, Arg1 postGain ) -> Arg0
{
    return nt2::tanh( sig * preGain ) * postGain;
}


template< typename Arg0, typename Arg1 >
auto fast_tanh_saturator( Arg0 sig, Arg1 preGain, Arg1 postGain ) -> Arg0
{
    return nova::approximations::fast_tanh<Arg0>( sig * preGain ) * postGain;
}

template< typename Arg0, typename Arg1 >
auto faster_tanh_saturator( Arg0 sig, Arg1 preGain, Arg1 postGain ) -> Arg0
{
    return nova::approximations::faster_tanh<Arg0>( sig * preGain ) * postGain;
}


}}

#endif // SATURATORS_HPP

