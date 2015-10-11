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
#include <boost/simd/constant/include/constants/two.hpp>
#include <boost/simd/constant/include/constants/pi.hpp>

#include <boost/simd/operator/include/functions/multiplies.hpp>
#include <boost/simd/operator/include/functions/plus.hpp>

#include <nt2/include/functions/exp.hpp>

namespace nova {

template <typename SampleType, typename ParameterType>
struct LeakDC
{
    LeakDC (ParameterType a = ParameterType(0.f)):
        _a(a)
    {}

    // internal state
    void set_x1( SampleType x1 ) { _x_1 = x1; }

    // parameters
    void set_a (ParameterType a) {  _a = a; }

    auto getA()
    {
        return [state = _a] { return state; };
    }

    template <typename FreqType, typename SampleRateType>
    static auto computeAForCutoff( FreqType cutoff, SampleRateType sampleRate, SampleRateType sampleDur )
    {
        using namespace boost::simd;

        auto cutoffFreq = nova::clip(cutoff, SampleRateType(0.1f), sampleRate);
        auto parameter = nt2::exp( - Two<SampleRateType>() * Pi<SampleRateType>() * cutoffFreq * (SampleRateType)sampleDur );

        return parameter;
    }


    // DSP
    template < typename InputFunctor, typename OutputFunctor >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count )
    {
        run( in, out, count, getA() );
    }

    template < typename InputFunctor, typename OutputFunctor, typename A >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count, A && a )
    {
        SampleType    y_1 = _y_1;
        SampleType    x_1 = _x_1;

        const size_t unroll4 = count / 2;
        const size_t remain  = count & 1;

        for (size_t i = 0; i != unroll4; ++i) {
            SampleType x0   = in();
            SampleType y0 = tick(x0, x_1, y_1, a());


            SampleType x1 = in();
            SampleType y1 = tick(x1, x_1, y_1, a());

            out(y0);
            out(y1);
        }

        for (size_t i = 0; i != remain; ++i) {
            SampleType x = in();
            SampleType y = tick(x, x_1, y_1, a());

            out(y);
        }

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

private:
    ParameterType _a;
    SampleType _x_1 = {0};
    SampleType _y_1 = {0};
};

}

#endif // LEAK_DC_HPP
