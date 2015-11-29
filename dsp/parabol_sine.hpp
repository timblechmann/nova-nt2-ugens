/*
    Copyright (C) 2013 Tim Blechmann

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

#ifndef PARABOL_SINE_HPP
#define PARABOL_SINE_HPP

#include <boost/simd/arithmetic/functions.hpp>

#include <tuple>

#include <boost/simd/include/functions/if_else.hpp>
#include <boost/simd/include/functions/is_greater.hpp>

#include <boost/simd/include/functions/combine.hpp>
#include <boost/simd/include/functions/slice.hpp>

#include <boost/simd/include/constants/four.hpp>
#include <boost/simd/include/constants/half.hpp>
#include <boost/simd/include/constants/pi.hpp>
#include <boost/simd/include/constants/two.hpp>

#include <boost/dispatch/meta/scalar_of.hpp>
#include <boost/simd/sdk/meta/cardinal_of.hpp>

namespace nova {

/* fast sine approximation
 * adapted from http://www.devmaster.net/forums/showthread.php?t=5784
 *
 * ok between -pi..pi
 * */
template <typename Type>
Type parabol_sin( Type x )
{
    using namespace boost::simd;

    const Type B =  Four<Type>() / Pi<Type>();
    const Type C = -Four<Type>() / ( Pi<Type>() * Pi<Type>() );

    const Type y = B * x + C * x * abs(x);

    const Type P (0.225);

    return P * (y * abs(y) - y) + y;   // Q * y + P * y * abs(y)
}


template <typename Type>
Type parabol_cos( Type x )
{
    using namespace boost::simd;

    x += Pi<Type>() * Half<Type>();

    x = if_else( is_greater( x, Pi<Type>() ),
                 x - Two<Type>() * Pi<Type>(),
                 x );

    return parabol_sin( x );
}

template <typename Type>
std::tuple<Type, Type> parabol_sincos( Type x )
{
    using namespace boost::simd;

    Type sinArg = x;
    Type cosArg = x;
    cosArg += Pi<Type>() * Half<Type>();

    cosArg = if_else( is_greater( cosArg, Pi<Type>() ),
                      cosArg - Two<Type>() * Pi<Type>(),
                      cosArg );


    typedef typename boost::simd::pack< typename boost::dispatch::meta::scalar_of< Type >::type,
                                        boost::simd::meta::cardinal_of< Type >::value * 2 > CombinedType;

    CombinedType combinedArg = combine( sinArg, cosArg );

    CombinedType result = parabol_sin( combinedArg );

    Type cos;
    Type sin = slice( result, cos );

    return std::make_tuple( sin, cos );
}

}

#endif // PARABOL_SINE_HPP
