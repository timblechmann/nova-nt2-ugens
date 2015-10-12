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
    typedef boost::simd::aligned_array<ParameterType, 1, 4> ParameterState;

    enum {
        stateFB
    };

    FeedbackAM (ParameterState fb = ParameterState{0.f}):
        _state(fb)
    {}

    // internal state
    void setState( ParameterState const & newState ) { _state = newState; }
    auto getState()
    {
        auto ret = [=] { return _state; };
        return ret;
    }

    ParameterState currentState()
    {
        return _state;
    }

    template <typename DSPContext>
    ParameterState computeState( ParameterType fb, DSPContext && )
    {
        auto clippedParameter = clip( fb, ParameterType(-2.1f), ParameterType(2.1f) );
        return ParameterState { clippedParameter } ;
    }

    template < typename InputFunctor, typename OutputFunctor, typename FeedbackFunctor >
    inline void run ( InputFunctor && in, OutputFunctor && out, size_t count, FeedbackFunctor && state )
    {
        SampleType    y_1 = _y_1;

        size_t remain = count;

        for (size_t i = 0; i != remain; ++i) {
            SampleType x = in();
            SampleType y = tick(x, y_1, state()[stateFB]);
            out(y);
        }
        _y_1 = y_1;
    }

    template < typename InputFunctor, typename OutputFunctor >
    inline void run ( InputFunctor && in, OutputFunctor && out, size_t count )
    {
        run( in, out, count, getState() );
    }

    template < typename SigInputFunctor, typename StateInputFunctor, typename OutputFunctor>
    inline void run_ar ( SigInputFunctor && in, StateInputFunctor && stateIn, OutputFunctor && out, size_t count)
    {
        run( in, out, count, stateIn );
    }

    static inline SampleType tick( SampleType input, SampleType & y_1, ParameterType fb )
    {
        SampleType output = input * ( boost::simd::One<SampleType>() + y_1 * fb );
        y_1 = output;
        return output;
    }


private:
    ParameterState _state;
    SampleType _y_1 = { 0.f };
};

}

#endif // FEEDBACK_AM_HPP
