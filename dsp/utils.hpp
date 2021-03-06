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

#ifndef UTILS_HPP
#define UTILS_HPP

#include "boost/simd/include/pack.hpp"
#include "boost/simd/include/functions/max.hpp"
#include "boost/simd/include/functions/min.hpp"
#include "boost/simd/include/functions/negs.hpp"

namespace nova {

template < typename SampleType0, typename SampleType1, typename SampleType2 >
inline auto clip ( SampleType0 x, SampleType1 lo, SampleType2 hi)
{
    using namespace boost::simd;
    return max( lo, min (x, hi) );
}

template <typename SampleType0, typename SampleType1>
inline auto clip2 ( SampleType0 x, SampleType1 hi)
{
    using namespace boost::simd;
    auto lo = neg(hi);
    return clip(x, lo, hi);
}

template <typename ScalarType, size_t Size>
struct as_pack
{
    typedef boost::simd::pack<ScalarType, Size>  type;
};

template <typename ScalarType>
struct as_pack< ScalarType, 1 >
{
    typedef ScalarType type;
};


namespace impl {

template <typename ReturnType, typename ArgumentType>
BOOST_FORCEINLINE ReturnType castType( ArgumentType const & arg, std::true_type isSame )
{
    return arg;
}

template <typename ReturnType, typename ArgumentType>
BOOST_FORCEINLINE ReturnType castType( ArgumentType const & arg, std::false_type isSame )
{
    using namespace boost::simd;
    typedef typename meta::scalar_of<ArgumentType>::type ArgScalar;
    typedef typename meta::scalar_of<ReturnType>::type   ReturnScalar;

    const size_t size = meta::cardinal_of<ReturnType>::value;

    typedef typename as_pack< ArgScalar, size >::type EvaluatedArgType;

    const EvaluatedArgType evaluatedArg = arg;
    ReturnType ret;
    for (size_t i = 0; i != size; ++i) {
        auto scalar = extract(evaluatedArg, i);
        insert(scalar, ret, i);
    }

    return ret;
}

}

template <typename ReturnType, typename ArgumentType>
BOOST_FORCEINLINE ReturnType castType( ArgumentType const & arg )
{
    using namespace boost::simd;

    typedef typename meta::scalar_of<ArgumentType>::type ArgScalar;
    typedef typename meta::scalar_of<ReturnType>::type   ReturnScalar;

    return impl::castType<ReturnType>( arg, typename std::is_same< ArgScalar, ReturnScalar>::type() );
}

}

#endif // UTILS_HPP
