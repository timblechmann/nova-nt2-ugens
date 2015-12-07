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


#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <type_traits>

#include <boost/simd/constant/include/constants/one.hpp>
#include <boost/simd/constant/include/constants/zero.hpp>
#include <boost/simd/operator/include/functions/multiplies.hpp>
#include <boost/simd/operator/include/functions/plus.hpp>

#include "utils.hpp"

namespace nova {

template <typename SampleType, typename ParameterType>
struct Integrator
{
    typedef boost::simd::aligned_array<ParameterType, 1, 4> ParameterState;

    enum {
        stateA
    };

    Integrator (ParameterState a = ParameterState{0.f}):
        _state(a)
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
    ParameterState computeState( ParameterType leak, DSPContext const & )
    {
        return computeState(leak);
    }

    ParameterState computeState( ParameterType leak )
    {
        auto clippedParameter = clip( leak, boost::simd::Zero<ParameterType>(), boost::simd::One<ParameterType>() );
        return ParameterState { clippedParameter } ;
    }

    // dsp
    template < typename InputFunctor, typename OutputFunctor, typename ParamInput >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count, ParamInput a )
    {
        SampleType    y_1 = _y_1;

        const size_t unroll2 = count / 2;
        const size_t remain  = count & 1;

        for (size_t i = 0; i != unroll2; ++i) {
            SampleType x0 = in();
            SampleType y0 = tick(x0, y_1, a()[stateA]);

            SampleType x1 = in();
            SampleType y1 = tick(x1, y_1, a()[stateA]);
            out(y0);
            out(y1);
        }

        for (size_t i = 0; i != remain; ++i) {
            SampleType x = in();
            SampleType y = tick(x, y_1, a()[stateA]);

            out(y);
        }

        _y_1 = y_1;
    }

    template < typename InputFunctor, typename OutputFunctor >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count )
    {
        run( in, out, count, [=]{ return _state; } );
    }

	template < typename SigInputFunctor, typename StateInputFunctor, typename OutputFunctor>
	inline void run_ar ( SigInputFunctor && sigIn, StateInputFunctor && stateIn, OutputFunctor && out, size_t count )
	{
		SampleType    y_1 = _y_1;
		for (size_t i = 0; i != count; ++i)
		{
			auto x  = sigIn();
			auto y  = tick( x, y_1, castType<ParameterType>( stateIn()[stateA] ) );
			out(y);
		}
		_y_1 = y_1;
	}

    static inline SampleType tick( SampleType input, SampleType & y_1, ParameterType a )
    {
        SampleType output = input + y_1 * a;
        y_1 = output;
        return output;
    }


private:
    ParameterState _state;
    SampleType _y_1 = {0};
};

}

#endif // INTEGRATOR_HPP
