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

#include <type_traits>

#include <boost/simd/include/pack.hpp>
#include <boost/simd/constant/include/constants/zero.hpp>
#include <boost/simd/include/functions/split.hpp>
#include <boost/simd/include/functions/group.hpp>
#include <boost/simd/include/functions/groups.hpp>
#include <boost/simd/operator/include/functions/multiplies.hpp>
#include <boost/simd/operator/include/functions/plus.hpp>

#ifndef BIQUAD_HPP
#define BIQUAD_HPP

namespace nova {

template <typename SampleType, typename ParameterStruct>
struct Biquad
{
    Biquad() :
        _y_1(0), _y_2(0)
    {}

    typedef typename ParameterStruct::type ParameterType;
    typedef decltype(boost::simd::Zero<ParameterType>()) Zero;

    template < typename InputFunctor, typename OutputFunctor, typename SlopeType>
    BOOST_FORCEINLINE void run ( InputFunctor & in, OutputFunctor & out, size_t count,
                                 SlopeType a0Slope, SlopeType a1Slope, SlopeType a2Slope,
                                 SlopeType b1Slope, SlopeType b2Slope)
    {
        SampleType    y_1 = _y_1;
        SampleType    y_2 = _y_2;

        ParameterType a0 = _parameters.a0();
        ParameterType a1 = _parameters.a1();
        ParameterType a2 = _parameters.a2();
        ParameterType b1 = _parameters.b1();
        ParameterType b2 = _parameters.b2();

        for (size_t i = 0; i != count; ++i) {
            SampleType input = in();

            SampleType y_0;
            SampleType output = tick( input, y_0, y_1, y_2, a0, a1, a2, b1, b2 );

            y_2 = y_1;
            y_1 = y_0;

            a0 += a0Slope;
            a1 += a1Slope;
            a2 += a2Slope;
            b1 += b1Slope;
            b2 += b2Slope;

            out(output);
        }

        _y_1 = y_1;
        _y_2 = y_2;
    }

    template < typename InputFunctor, typename OutputFunctor >
    BOOST_FORCEINLINE void run ( InputFunctor & in, OutputFunctor & out, size_t count)
    {
        SampleType    y_1 = _y_1;
        SampleType    y_2 = _y_2;

        auto a0 = _parameters.a0();
        auto a1 = _parameters.a1();
        auto a2 = _parameters.a2();
        auto b1 = _parameters.b1();
        auto b2 = _parameters.b2();

        for (size_t i = 0; i != count; ++i) {
            SampleType input = in();

            SampleType y_0;
            SampleType output = tick( input, y_0, y_1, y_2, a0, a1, a2, b1, b2 );

            y_2 = y_1;
            y_1 = y_0;

            out(output);
        }

        _y_1 = y_1;
        _y_2 = y_2;
    }

    template <typename ParameterTypeA0, typename ParameterTypeA1, typename ParameterTypeA2,
              typename ParameterTypeB1, typename ParameterTypeB2>
    BOOST_FORCEINLINE SampleType tick(SampleType input, SampleType & y_0, SampleType y_1, SampleType y_2,
                                      ParameterTypeA0 a0, ParameterTypeA1 a1, ParameterTypeA2 a2,
                                      ParameterTypeB1 b1, ParameterTypeB2 b2)
    {
        using namespace boost::simd;
#if 0 //def __AVX__
        // LATER: try to group/split into double-sized array
        pack<double, 4> par1 = group(SampleType(a1), SampleType(b1));
        pack<double, 4> par2 = group(SampleType(a2), SampleType(b2));

        auto res1 = par1 * group(y_1, y_1);
        auto res2 = par2 * group(y_2, y_2);

        SampleType a1y1, b1y1;
        boost::simd::split (res1, a1y1, b1y1);

        SampleType a2y2, b2y2;
        boost::simd::split (res2, a2y2, b2y2);

        y_0               = input      - b1y1 - b2y2;
        SampleType output = (a0 * y_0) + a1y1 + a2y2;
#else
        y_0               = input      - (b1 * y_1) - (b2 * y_2);
        SampleType output = (a0 * y_0) + (a1 * y_1) + (a2 * y_2);
#endif
        return output;
    }

    template < typename InputFunctor, typename OutputFunctor >
    BOOST_FORCEINLINE void run_unrolled ( InputFunctor & in, OutputFunctor & out, size_t count_3, size_t count)
    {
        SampleType    y_1 = _y_1;
        SampleType    y_2 = _y_2;

        auto a0 = _parameters.a0();
        auto a1 = _parameters.a1();
        auto a2 = _parameters.a2();
        auto b1 = _parameters.b1();
        auto b2 = _parameters.b2();

        // unroll by 3
        for (size_t i = 0; i != count_3; ++i) {
            SampleType input0 = in();
            SampleType input1 = in();
            SampleType input2 = in();

            SampleType y_0;
            SampleType output0 = tick( input0, y_0, y_1, y_2, a0, a1, a2, b1, b2 );
            SampleType output1 = tick( input1, y_2, y_0, y_1, a0, a1, a2, b1, b2 );
            SampleType output2 = tick( input2, y_1, y_2, y_0, a0, a1, a2, b1, b2 );

            out(output0);
            out(output1);
            out(output2);
        }

        for (size_t i = 0; i != count; ++i) {
            SampleType input = in();

            SampleType y_0;
            SampleType output = tick( input, y_0, y_1, y_2, a0, a1, a2, b1, b2 );
            y_2 = y_1;
            y_1 = y_0;

            out(output);
        }

        _y_1 = y_1;
        _y_2 = y_2;
    }

    template < typename InputFunctor, typename OutputFunctor, typename SlopeType = Zero >
    BOOST_FORCEINLINE void run_unrolled ( InputFunctor & in, OutputFunctor & out, size_t count_3, size_t count,
                                          SlopeType a0Slope, SlopeType a1Slope, SlopeType a2Slope,
                                          SlopeType b1Slope, SlopeType b2Slope)
    {
        SampleType    y_1 = _y_1;
        SampleType    y_2 = _y_2;

        ParameterType a0 = _parameters.a0();
        ParameterType a1 = _parameters.a1();
        ParameterType a2 = _parameters.a2();
        ParameterType b1 = _parameters.b1();
        ParameterType b2 = _parameters.b2();

        // unroll by 3
        for (size_t i = 0; i != count_3; ++i) {
            SampleType input0 = in();
            SampleType input1 = in();
            SampleType input2 = in();

            SampleType y_0;
            SampleType output0 = tick( input0, y_0, y_1, y_2, a0, a1, a2, b1, b2 );

            a0 += a0Slope;
            a1 += a1Slope;
            a2 += a2Slope;
            b1 += b1Slope;
            b2 += b2Slope;

            SampleType output1 = tick( input1, y_2, y_0, y_1, a0, a1, a2, b1, b2 );

            a0 += a0Slope;
            a1 += a1Slope;
            a2 += a2Slope;
            b1 += b1Slope;
            b2 += b2Slope;

            SampleType output2 = tick( input2, y_1, y_2, y_0, a0, a1, a2, b1, b2 );

            a0 += a0Slope;
            a1 += a1Slope;
            a2 += a2Slope;
            b1 += b1Slope;
            b2 += b2Slope;

            out(output0);
            out(output1);
            out(output2);
        }

        for (size_t i = 0; i != count; ++i) {
            SampleType input = in();

            SampleType y_0;
            SampleType output = tick( input, y_0, y_1, y_2, a0, a1, a2, b1, b2 );
            y_2 = y_1;
            y_1 = y_0;

            a0 += a0Slope;
            a1 += a1Slope;
            a2 += a2Slope;
            b1 += b1Slope;
            b2 += b2Slope;

            out(output);
        }

        _y_1 = y_1;
        _y_2 = y_2;
    }


    SampleType _y_1, _y_2;
    ParameterStruct _parameters;
};

}

#endif // BIQUAD_HPP
