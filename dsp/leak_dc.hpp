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


#ifndef LEAK_DC_HPP
#define LEAK_DC_HPP

#include <type_traits>

#include <boost/simd/constant/include/constants/zero.hpp>
#include <boost/simd/operator/include/functions/multiplies.hpp>
#include <boost/simd/operator/include/functions/plus.hpp>

namespace nova {

template <typename SampleType, typename ParameterType>
struct LeakDC
{
	LeakDC (ParameterType a = ParameterType(0.f)):
        _y_1(0.f),
		_x_1(0.f),
        _a(a)
    {}

    void set_a (ParameterType a)
    {
        _a = a;
    }

	typedef decltype(boost::simd::Zero<ParameterType>()) Zero;

	template < typename InputFunctor, typename OutputFunctor, typename ASlope = Zero >
	inline void run ( InputFunctor & in, OutputFunctor & out, size_t count, ASlope aSlope = ASlope() )
    {
        SampleType    y_1 = _y_1;
		SampleType    x_1 = _x_1;
        ParameterType a   = _a;

		const size_t unroll4 = count / 2;
		const size_t remain  = count & 1;

		for (size_t i = 0; i != unroll4; ++i) {
			SampleType x0 = in();
			SampleType y0 = tick(x0, x_1, y_1, a);

			if (!std::is_same<ASlope, Zero>::value)
				a += aSlope;

			SampleType x1 = in();
			SampleType y1 = tick(x1, x_1, y_1, a);
			a += aSlope;
			if (!std::is_same<ASlope, Zero>::value)
				a += aSlope;

//			SampleType x2 = in();
//			SampleType y2 = tick(x2, x_1, y_1, a);
//			a += aSlope;

//			SampleType x3 = in();
//			SampleType y3 = tick(x3, x_1, y_1, a);
//			a += aSlope;

			out(y0);
			out(y1);
//			out(y2);
//			out(y3);
        }

		for (size_t i = 0; i != remain; ++i) {
			SampleType x = in();
			SampleType y = tick(x, x_1, y_1, a);
			a += aSlope;

			out(y);
		}

		if (!std::is_same<ASlope, Zero>::value)
			_a = a;

        _y_1 = y_1;
		_x_1 = x_1;
    }

	static inline SampleType tick( SampleType input, SampleType & x_1, SampleType & y_1, ParameterType a )
    {
		SampleType output = input - x_1 + a * y_1;
        y_1 = output;
		x_1 = input;
        return output;
    }

	ParameterType _a;
	SampleType _x_1;
private:
	SampleType _y_1;
};

}

#endif // LEAK_DC_HPP
