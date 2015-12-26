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

#ifndef RANGE_HPP
#define RANGE_HPP

#include <type_traits>

#include <boost/range/irange.hpp>

namespace nova
{

template <typename T,
          typename = std::enable_if_t<std::is_integral<T>::value> >
inline auto range( T const & a )
{
	return boost::irange( T(0), a );
}

template <typename T>
inline auto range( T const & a, T const & b )
{
	return boost::irange( a, b );
}

}

#endif // RANGE_HPP
