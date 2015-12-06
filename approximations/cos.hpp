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

#ifndef COS_APPROXIMATION_HPP
#define COS_APPROXIMATION_HPP


#include <approximations/sin.hpp>

#include <nt2/include/functions/cos.hpp>

#include <boost/simd/include/constants/one.hpp>
#include <boost/simd/include/constants/pio_2.hpp>
#include <boost/simd/include/functions/abs.hpp>
#include <boost/simd/include/functions/if_else.hpp>
#include <boost/simd/include/functions/splat.hpp>

namespace nova {
namespace approximations {

struct CosPrecise {};
struct CosFast    {};
struct CosFaster  {};

template <typename Arg>
BOOST_FORCEINLINE auto cos( Arg x, CosPrecise ) -> Arg
{
    return nt2::cos( x );
}

template <typename Arg>
BOOST_FORCEINLINE auto cos( Arg x, CosFast ) -> Arg
{
    using namespace boost::simd;

    const Arg piOver2           = Pio_2<Arg>();
    const Arg piOver2MinusTwoPi = splat<Arg>(-4.7123889803846899f);

    const Arg offset = if_else( x > piOver2, piOver2MinusTwoPi, piOver2 );

    return sin( x + offset, SinFast() );
}

template <typename Arg>
BOOST_FORCEINLINE auto cos( Arg x, CosFaster ) -> Arg
{
    using namespace boost::simd;

    const Arg twoOverPi = splat<Arg>( 0.63661977236758134f );
    const Arg p         = splat<Arg>( 0.54641335845679634f );

    const Arg qpprox    = One<Arg>() - twoOverPi * abs( x );

    return qpprox + p * qpprox * ( One<Arg>() - qpprox * qpprox );
}

}}

#endif // COS_APPROXIMATION_HPP
