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
 *=====================================================================* */

#ifndef TAN_APPROXIMATION_HPP
#define TAN_APPROXIMATION_HPP

#include <boost/simd/include/constants/four.hpp>
#include <boost/simd/include/constants/one.hpp>
#include <boost/simd/include/constants/pio_2.hpp>
#include <boost/simd/include/constants/pi.hpp>
#include <boost/simd/include/functions/splat.hpp>
#include <boost/simd/include/functions/bitwise_and.hpp>
#include <boost/simd/include/functions/bitwise_or.hpp>
#include <boost/simd/operator/operator.hpp>

#include <boost/simd/arithmetic/include/functions/abs.hpp>
#include <boost/simd/arithmetic/include/functions/raw_rec.hpp>
#include <boost/simd/arithmetic/include/functions/fast_rec.hpp>

#include <boost/simd/ieee/include/functions/copysign.hpp>
#include <boost/simd/ieee/include/functions/bitofsign.hpp>

#include <nt2/include/functions/sin.hpp>
#include <nt2/include/functions/tan.hpp>

namespace nova {
namespace approximations {

struct SinPrecise {};
struct SinFast    {};
struct SinFaster  {};


struct TanPrecise {};
struct TanFast    {};
struct TanFaster  {};


template <typename Arg>
BOOST_FORCEINLINE auto sin( Arg x, SinPrecise ) -> Arg
{
    return nt2::sin( x );
}

template <typename Arg>
BOOST_FORCEINLINE auto sin( Arg x, SinFast ) -> Arg
{
    using namespace boost::simd;

    const Arg fourOverPi   = boost::simd::Four<Arg>() / boost::simd::Pi<Arg>(); // 1.2732395447351627f
    const Arg fourOverPiSq = splat<Arg>( 0.40528473456935109f );

    const Arg q = splat<Arg>( 0.78444488374548933f);
    const Arg p = splat<Arg>( 0.20363937680730309f);
    const Arg r = splat<Arg>( 0.015124940802184233f);
    const Arg s = splat<Arg>(-0.0032225901625579573f);

    const Arg qpprox   = fourOverPi * x - fourOverPiSq * x * abs( x );
    const Arg qpproxsq = qpprox * qpprox;

    Arg y              = qpproxsq * (p + qpproxsq * (r + qpproxsq * s));

    y        = bitwise_xor( y, bitofsign( x ) );

    return q * qpprox + y;
}

template <typename Arg>
BOOST_FORCEINLINE auto sin( Arg x, SinFaster ) -> Arg
{

    using namespace boost::simd;

    const Arg fourOverPi   = boost::simd::Four<Arg>() / boost::simd::Pi<Arg>(); // 1.2732395447351627f
    const Arg fourOverPiSq = splat<Arg>( 0.40528473456935109f );

    const Arg q    = splat<Arg>(0.77633023248007499f);
    const Arg plit = splat<Arg>(0.22308510060189463f);

    const Arg qpprox = fourOverPi * x - fourOverPiSq * x * boost::simd::abs( x );

    const Arg p    = bitwise_or( bitofsign( x ), plit );
    // const Arg p    = copysign( plit, x );

    return qpprox * (q + p * qpprox);
}







template <typename Arg>
BOOST_FORCEINLINE auto tan( Arg x, TanPrecise ) -> Arg
{
    return nt2::tan( x );
}

template <typename Arg>
BOOST_FORCEINLINE auto tan( Arg x, TanFast ) -> Arg
{
    Arg piOver2 = boost::simd::Pio_2<Arg>();

    return sin( x, SinFast() ) * boost::simd::fast_rec( sin( x * piOver2, SinFast() ) );
}

template <typename Arg>
BOOST_FORCEINLINE auto tan( Arg x, TanFaster ) -> Arg
{
    Arg piOver2 = boost::simd::Pio_2<Arg>();

    return sin( x, SinFaster() ) * boost::simd::raw_rec( sin( x * piOver2, SinFaster() ) );
}



}}

#endif // TAN_APPROXIMATION_HPP
