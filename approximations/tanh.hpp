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

#ifndef TANH_APPROXIMATION_HPP
#define TANH_APPROXIMATION_HPP

#include <boost/simd/include/constants/one.hpp>
#include <boost/simd/include/functions/splat.hpp>
#include <boost/simd/operator/operator.hpp>

#include <boost/simd/arithmetic/include/functions/abs.hpp>
#include <boost/simd/arithmetic/include/functions/raw_rec.hpp>
#include <boost/simd/arithmetic/include/functions/fast_rec.hpp>

#include <boost/simd/ieee/include/functions/copysign.hpp>

namespace nova {
namespace approximations {

template< typename Type >
Type faster_tanh( Type x )
{
    // adapted from http://www.kvraudio.com/forum/viewtopic.php?p=3778524
    namespace bs = boost::simd;

    const Type one = bs::One<Type>();
    const Type a   = bs::splat<Type>( 0.66422417311781 );
    const Type b   = bs::splat<Type>( 0.36483285408241 );

    const Type absX = bs::abs( x );
    const Type x2   = x * x;

    const Type absResult = one - bs::fast_rec(one + absX + x2 + a * x2*absX + b * x2 * x2);

    return bs::copysign( absResult, x );
}

template< typename Type >
Type fast_tanh( Type x )
{
    // adapted from http://www.kvraudio.com/forum/viewtopic.php?p=3778524
    namespace bs = boost::simd;

    const Type one = boost::simd::One<Type>();
    const Type a   = boost::simd::splat<Type>( 0.58576695 );
    const Type b   = boost::simd::splat<Type>( 0.55442112 );
    const Type c   = boost::simd::splat<Type>( 0.057481508 );

    const Type absX = bs::abs( x );
    const Type x2   =  x *  x;
    const Type x3   = x2 * absX;
    const Type x4   = x2 * x2;
    const Type x7   = x3 * x4;

    const Type absResult = one - bs::fast_rec( one + absX + x2 + a * x3 + b * x4 + c * x7 );

    return bs::copysign( absResult, x );
}


}}

#endif // TANH_APPROXIMATION_HPP
