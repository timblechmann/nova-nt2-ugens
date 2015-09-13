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


#ifndef FEEDBACK_AM_HPP
#define FEEDBACK_AM_HPP

#include <type_traits>

#include <boost/simd/constant/include/constants/one.hpp>
#include <boost/simd/constant/include/constants/zero.hpp>
#include <boost/simd/operator/include/functions/multiplies.hpp>
#include <boost/simd/operator/include/functions/plus.hpp>

#include "utils.hpp"

namespace nova {

template <typename SampleType, typename ParameterType>
struct FeedbackAM
{
    FeedbackAM (ParameterType fb = ParameterType(0.f)):
        _fb(fb)
    {}

    typedef decltype(boost::simd::Zero<ParameterType>()) Zero;

    template < typename InputFunctor, typename OutputFunctor, typename ASlope = Zero >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count, ASlope fbSlope = ASlope() )
    {
        SampleType    y_1 = _y_1;
        ParameterType fb   = _fb;

#if 0
        const size_t unroll2 = count / 2;
        const size_t remain  = count & 1;

        for (size_t i = 0; i != unroll2; ++i) {
            SampleType x0 = in();
            SampleType y0 = tick(x0, y_1, fb);

            if (!std::is_same<ASlope, Zero>::value)
                fb += fbSlope;

            SampleType x1 = in();
            SampleType y1 = tick(x1, y_1, fb);
            fb += fbSlope;
            if (!std::is_same<ASlope, Zero>::value)
                fb += fbSlope;
            out(y0);
            out(y1);
        }
#endif
        size_t remain = count;

        for (size_t i = 0; i != remain; ++i) {
            SampleType x = in();
            SampleType y = tick(x, y_1, fb);
            if (!std::is_same<ASlope, Zero>::value)
                fb += fbSlope;

            out(y);
        }

        if (!std::is_same<ASlope, Zero>::value)
            _fb = fb;

        _y_1 = y_1;
    }

    template < typename SigInputFunctor, typename ParamInputFunctor, typename OutputFunctor>
    inline void run_ar ( SigInputFunctor & sigIn, ParamInputFunctor & paramIn, OutputFunctor & out, size_t count)
    {
        SampleType    y_1 = _y_1;
        for (size_t i = 0; i != count; ++i)
        {
            auto x  = sigIn();
            auto fb = paramIn();

            auto clippedFB = checkParameter(fb);
            auto y  = tick(x, y_1, toDouble<ParameterType>( clippedFB ) );
            out(y);
        }
        _y_1 = y_1;
    }

    template <typename Arg>
    static auto checkParameter( Arg fb )
    {
        return clip( fb, Arg(-2.1f), Arg(2.1f) );
    }

    static inline SampleType tick( SampleType input, SampleType & y_1, ParameterType fb )
    {
        SampleType output = input * ( boost::simd::One<SampleType>() + y_1 * fb );
        y_1 = output;
        return output;
    }

    ParameterType _fb;

private:
    SampleType _y_1 = { 0.f };
};

}

#endif // FEEDBACK_AM_HPP
