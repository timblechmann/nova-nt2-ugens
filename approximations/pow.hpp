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

    adapted from fastapprox:

 *=====================================================================*
 *                   Copyright (C) 2011 Paul Mineiro                   *
 * All rights reserved.                                                *
 *                                                                     *
 * Redistribution and use in source and binary forms, with             *
 * or without modification, are permitted provided that the            *
 * following conditions are met:                                       *
 *                                                                     *
 *     * Redistributions of source code must retain the                *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer.                                       *
 *                                                                     *
 *     * Redistributions in binary form must reproduce the             *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer in the documentation and/or            *
 *     other materials provided with the distribution.                 *
 *                                                                     *
 *     * Neither the name of Paul Mineiro nor the names                *
 *     of other contributors may be used to endorse or promote         *
 *     products derived from this software without specific            *
 *     prior written permission.                                       *
 *                                                                     *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND              *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,         *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES               *
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE             *
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER               *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,                 *
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES            *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE           *
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR                *
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF          *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT           *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY              *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE             *
 * POSSIBILITY OF SUCH DAMAGE.                                         *
 *                                                                     *
 * Contact: Paul Mineiro <paul@mineiro.com>                            *
 *=====================================================================*

 */

#ifndef POW_APPROXIMATION_HPP
#define POW_APPROXIMATION_HPP

#include <boost/simd/include/constants/one.hpp>
#include <boost/simd/include/constants/zero.hpp>

#include <boost/simd/include/functions/splat.hpp>
#include <boost/simd/operator/operator.hpp>

#include <boost/simd/arithmetic/include/functions/abs.hpp>
#include <boost/simd/arithmetic/include/functions/raw_rec.hpp>
#include <boost/simd/arithmetic/include/functions/fast_rec.hpp>

#include <boost/simd/ieee/include/functions/copysign.hpp>
#include <boost/simd/include/functions/if_one_else_zero.hpp>
#include <boost/simd/include/functions/if_else.hpp>
#include <boost/simd/include/functions/if_else_zero.hpp>
#include <boost/simd/include/functions/fast_toint.hpp>
#include <boost/simd/include/functions/tofloat.hpp>
#include <boost/simd/include/functions/trunc.hpp>
#include <boost/simd/include/functions/bitwise_cast.hpp>

#include <nt2/include/functions/pow2.hpp>
#include <nt2/include/functions/log2.hpp>

namespace nova {
namespace approximations {


struct PowPrecise {};
struct PowFast    {};
struct PowFaster  {};



// pow


template <typename Arg>
BOOST_FORCEINLINE auto pow2( Arg x, PowPrecise) -> Arg
{
    return nt2::pow2( x );
}


template <typename Arg>
BOOST_FORCEINLINE auto pow2( Arg x, PowFast) -> Arg
{
    using namespace boost::simd;
    typedef typename boost::dispatch::meta::as_integer<Arg>::type IntType;

//    const Arg offset = if_one_else_zero( x < Zero<Arg>() );
    const Arg offset = if_else_zero( x < Zero<Arg>(), One<Arg>() );
    const Arg clipp  = if_else( x < splat<Arg>(-126),
                                splat<Arg>(-126),
                                x );

    const Arg z = clipp - trunc( clipp ) + offset;

    const Arg f = splat<Arg>(1 << 23) * (clipp + 121.2740575f + 27.7280233f * fast_rec(4.84252568f - z) - 1.49012907f * z);
    const IntType i = fast_toint( f );
    return boost::simd::bitwise_cast<Arg>( i );
}


template <typename Arg>
BOOST_FORCEINLINE auto pow2( Arg x, PowFaster) -> Arg
{
    using namespace boost::simd;
    typedef typename boost::dispatch::meta::as_integer<Arg>::type IntType;

    const Arg offset = if_else_zero( x < Zero<Arg>(), One<Arg>() );
    const Arg clipp  = if_else( x < splat<Arg>(-126),
                                splat<Arg>(-126),
                                x );

    const Arg z = clipp - trunc( clipp ) + offset;

    const Arg f = splat<Arg>(1 << 23) * (clipp + 121.2740575f + 27.7280233f * fast_rec(4.84252568f - z) - 1.49012907f * z);
    const IntType i = fast_toint( f );
    return boost::simd::bitwise_cast<Arg>( i );
}


// log2

template <typename Arg>
BOOST_FORCEINLINE auto log2( Arg x, PowPrecise) -> Arg
{
    return nt2::log2( x );
}

template <typename Arg>
BOOST_FORCEINLINE auto log2( Arg x, PowFast) -> Arg
{
    using namespace boost::simd;
    typedef typename boost::dispatch::meta::as_integer<Arg>::type IntType;

    const IntType i     = bitwise_cast<IntType>( x );
    const IntType mask  = (i & 0x007FFFFF) | 0x3f000000;
    const Arg maskFloat = bitwise_cast<Arg>( mask );

    const Arg y = tofloat( i ) * 1.1920928955078125e-7f;

    return y - 124.22551499f
             - 1.498030302f * maskFloat
             - 1.72587999f  * fast_rec( 0.3520887068f + maskFloat );
}

template <typename Arg>
BOOST_FORCEINLINE auto log2( Arg x, PowFaster) -> Arg
{
    using namespace boost::simd;
    typedef typename boost::dispatch::meta::as_integer<Arg>::type IntType;

    const IntType i = bitwise_cast<IntType>( x );
    const Arg y     = tofloat( i ) * 1.1920928955078125e-7f;

    return y - 126.94269504f;
}


// pow

template <typename Arg0, typename Arg1, typename PrecisionTag>
BOOST_FORCEINLINE auto pow( Arg0 x, Arg1 p, PrecisionTag precision)
{
    Arg0 pow2Arg = p * log2( x, precision );
    return pow2( pow2Arg, precision );
}


}}

#endif // POW_APPROXIMATION_HPP
