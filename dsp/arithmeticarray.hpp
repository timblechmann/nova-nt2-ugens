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

#ifndef ARITHMETICARRAY_HPP
#define ARITHMETICARRAY_HPP

#include <algorithm>
#include <array>
#include <initializer_list>

#include <boost/range/irange.hpp>


namespace nova   {
namespace detail {

template <typename Type, int N>
struct ArithmeticArray
{
    typedef Type value_type;

    ArithmeticArray()                                          = default;
    ArithmeticArray( ArithmeticArray const & rhs )             = default;
    ArithmeticArray & operator=( ArithmeticArray const & rhs ) = default;

    ArithmeticArray( std::initializer_list<Type> const & init )
    {
        std::copy( init.begin(), init.end(), data.begin() );
    }

    template <typename RhsType>
    ArithmeticArray( ArithmeticArray<RhsType, N> const & rhs )
    {
        for( int i : boost::irange(0, N) )
            data[i] = rhs[i];
    }

    Type const & get(size_t n) { return data[n]; }

    ArithmeticArray operator+(ArithmeticArray const & rhs) const
    {
        ArithmeticArray ret;
        for( int i : boost::irange(0, N) )
            ret[i] = data[i] + rhs[i];
        return ret;
    }

    template <typename RhsType>
    ArithmeticArray & operator+=(ArithmeticArray<RhsType, N> const & rhs)
    {
        ArithmeticArray ret;
        for( int i : boost::irange(0, N) )
            data[i] += rhs[i];
        return *this;
    }

    ArithmeticArray operator-(ArithmeticArray const & rhs) const
    {
        ArithmeticArray ret;
        for( int i : boost::irange(0, N) )
            ret[i] = data[i] - rhs[i];
        return ret;
    }

    template <typename ArgType>
    ArithmeticArray operator*(ArgType const & value) const
    {
        ArithmeticArray ret;
        for( int i : boost::irange(0, N) )
            ret[i] = data[i] * value;
        return ret;
    }

    Type & operator[](size_t n)             { return data[n]; }
    Type const & operator[](size_t n) const { return data[n]; }

private:
    std::array<Type, N> data;
};

}}

#endif // ARITHMETICARRAY_HPP
