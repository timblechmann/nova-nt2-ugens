/*
    Copyright (C) Tim Blechmann

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

#ifndef OSC_DPW_HPP
#define OSC_DPW_HPP

#include <cmath>

#include <boost/simd/include/pack.hpp>
#include <boost/simd/include/functions/fast_divides.hpp>
#include <boost/simd/include/functions/floor.hpp>
#include <boost/simd/include/functions/if_else.hpp>
#include <boost/simd/include/functions/is_greater_equal.hpp>

namespace nova {

namespace impl {

template <typename Derived>
struct DPWBase
{
protected:
    typedef float ParameterType;
    typedef boost::simd::aligned_array<ParameterType, 2, 4> ParameterState;

    enum {
        PhaseIncrement,
        Scale
    };

    explicit DPWBase( ParameterState state ):
        mState( state )
    {}

public:
    // internal state
    void setState( ParameterState const & newState ) { mState = newState;             }
    auto getState()                                  { return [=] { return mState; }; }
    ParameterState currentState()                    { return mState;                 }

    // DSP
    template < typename OutputFunctor >
    inline void run ( OutputFunctor & out, size_t count )
    {
        Derived::run( out, count, getState() );
    }

private:
    ParameterState mState;
};

}


struct SawDPW:
    public impl::DPWBase<SawDPW>
{
    typedef impl::DPWBase<SawDPW>                          DPWBase;
    typedef typename impl::DPWBase<SawDPW>::ParameterState ParameterState;
    typedef typename impl::DPWBase<SawDPW>::ParameterType  ParameterType;

    enum {
        PhaseIncrement,
        Scale
    };

    // parameters
    template <typename DSPContext>
    static auto computeState( ParameterType freq, DSPContext const & dspContext )
    {
        // we simulate negative frequencies by scaling positive frequencies by -1
        const ParameterType absScale = freq >= 0 ? 1 : -1;
        freq                         = std::abs( freq );

        const float sampleDuration = dspContext.sampleDur();
        ParameterType normalizedFrequency = freq * sampleDuration;
        ParameterType twoNFreq            = normalizedFrequency + normalizedFrequency;

        ParameterType scale = boost::simd::fast_divides( absScale,
                                                         (4.f * twoNFreq * ( 1.f - twoNFreq * sampleDuration) ) );

        ParameterType phaseIncrement = twoNFreq;

        return ParameterState{ phaseIncrement, scale };
    }


    SawDPW( ParameterState state, ParameterType initialPhase = 0.f ):
        DPWBase( state ), mPhase( initialPhase ), mLastValue( initialPhase )
    {}

    // DSP
    template < typename OutputFunctor >
    inline void run ( OutputFunctor & out, size_t count )
    {
        run( out, count, getState() );
    }

    template < typename OutputFunctor, typename State >
    inline void run ( OutputFunctor & out, size_t count, State && a )
    {
        ParameterType phase     = mPhase;
        ParameterType lastValue = mLastValue;

        const size_t unroll2 = count / 2;
        const size_t remain  = count & 1;

        for (size_t i = 0; i != unroll2; ++i) {
            ParameterState state = a();
            float y0 = tick( phase, lastValue, state[PhaseIncrement], state[Scale] );

            state = a();
            float y1 = tick( phase, lastValue, state[PhaseIncrement], state[Scale] );

            out(y0);
            out(y1);
        }

		if( BOOST_UNLIKELY( remain ) ) {
			for (size_t i = 0; i != remain; ++i) {
				ParameterState state = a();
				float y0 = tick( phase, lastValue, state[PhaseIncrement], state[Scale] );
				out(y0);
			}
		}

        mPhase     = phase;
        mLastValue = lastValue;
    }


    static inline float tick(ParameterType & phase, ParameterType & lastVal, float phaseIncrement, float scale)
    {
        phase = incrementPhase(phase, phaseIncrement);

        // squared saw
        float val = phase * phase;

        // differentiate parabolic wave
        float saw = (val - lastVal) * scale;

        lastVal = val;

        return saw;
    }

    static inline ParameterType incrementPhase(ParameterType phase, ParameterType phaseIncrement)
    {
        phase += phaseIncrement;

        return boost::simd::if_else( boost::simd::is_greater_equal( phase, ParameterType(1.f) ),
                                     phase - 2.f,
                                     phase );
    }

protected:
    ParameterType mPhase = 0.f, mLastValue = 0.f;
};


struct TriDPW:
    public impl::DPWBase<TriDPW>
{
    typedef impl::DPWBase<TriDPW>                          DPWBase;
    typedef typename impl::DPWBase<TriDPW>::ParameterState ParameterState;
    typedef typename impl::DPWBase<TriDPW>::ParameterType  ParameterType;

    // parameters
    template <typename DSPContext>
    static auto computeState( ParameterType freq, DSPContext const & dspContext )
    {
        // we simulate negative frequencies by scaling positive frequencies by -1
        const ParameterType absScale = freq >= 0 ? -1 : 1;
        freq                         = std::abs( freq );

        const float sampleDuration = dspContext.sampleDur();
        ParameterType normalizedFrequency = freq * sampleDuration;
        ParameterType twoNFreq            = normalizedFrequency + normalizedFrequency;

        ParameterType scale = boost::simd::fast_divides( absScale,
                                                         (1.f * twoNFreq * ( 1.f - twoNFreq * sampleDuration) ) );

        ParameterType phaseIncrement = twoNFreq;

        return ParameterState{ phaseIncrement, scale };
    }

    template <typename SampleType>
    inline SampleType wrap(SampleType x)
    {
        namespace bs = boost::simd;

        SampleType y = bs::floor( (x+1.f) * 0.5f );
        return x - y - y ;
    }

    TriDPW( ParameterState state, ParameterType initialPhase = 0.f ):
        DPWBase( state ), mPhase( wrap( initialPhase + 0.5 ) ), mLastValue( mPhase )
    {}

    // DSP
    template < typename OutputFunctor >
    inline void run ( OutputFunctor & out, size_t count )
    {
        run( out, count, getState() );
    }

    template < typename OutputFunctor, typename State >
    inline void run ( OutputFunctor & out, size_t count, State && a )
    {
        ParameterType phase     = mPhase;
        ParameterType lastValue = mLastValue;

        const size_t unroll2 = count / 2;
        const size_t remain  = count & 1;

        for (size_t i = 0; i != unroll2; ++i) {
            ParameterState state = a();
            float y0 = tick( phase, lastValue, state[PhaseIncrement], state[Scale] );

            state = a();
            float y1 = tick( phase, lastValue, state[PhaseIncrement], state[Scale] );

            out(y0);
            out(y1);
        }

		if( BOOST_UNLIKELY( remain ) ) {
			for (size_t i = 0; i != remain; ++i) {
				ParameterState state = a();
				float y0 = tick( phase, lastValue, state[PhaseIncrement], state[Scale] );
				out(y0);
			}
		}

        mPhase     = phase;
        mLastValue = lastValue;
    }


    static inline float tick(ParameterType & phase, ParameterType & lastVal, float phaseIncrement, float scale)
    {
        phase = incrementPhase(phase, phaseIncrement);

#if 0
		float val;
		if ( phase <  1.f )
			val = 1.f - phase*phase;
		else {
			float tmp = phase - 2.f;
			val = tmp * tmp - 1.f;
		}
#else
        float val = phase - phase * std::abs(phase);
#endif
        // differentiate
        float tri = (val - lastVal) * scale;

        lastVal = val;

        return tri;
    }

    static inline ParameterType incrementPhase(ParameterType phase, ParameterType phaseIncrement)
    {
        phase += phaseIncrement;

        return boost::simd::if_else( boost::simd::is_greater_equal( phase, 1.f),
                                     phase - 2.f,
                                     phase );
    }

protected:
    ParameterType mPhase = 0.f, mLastValue = 0.f;
};


struct PulseDPW:
    public impl::DPWBase<PulseDPW>
{
    typedef impl::DPWBase<PulseDPW>                          DPWBase;
    typedef typename impl::DPWBase<PulseDPW>::ParameterState ParameterState;
    typedef typename impl::DPWBase<PulseDPW>::ParameterType  ParameterType;

    // parameters
    template <typename DSPContext>
    static auto computeState( ParameterType freq, DSPContext const & dspContext )
    {
        // we simulate negative frequencies by scaling positive frequencies by -1
        const ParameterType absScale = freq >= 0 ? -1 : 1;
        freq                         = std::abs( freq );

        const float sampleDuration = dspContext.sampleDur();
        ParameterType normalizedFrequency = freq * sampleDuration;
        ParameterType twoNFreq            = normalizedFrequency + normalizedFrequency;

        ParameterType scale = boost::simd::fast_divides( absScale,
                                                         (4.f * twoNFreq * ( 1.f - twoNFreq * sampleDuration) ) );

        ParameterType phaseIncrement = twoNFreq;

        return ParameterState{ phaseIncrement, scale };
    }

    template <typename SampleType>
    inline SampleType wrap(SampleType x)
    {
        namespace bs = boost::simd;

        SampleType y = bs::floor( (x+1.f) * 0.5f );
        return x - y - y ;
    }

    PulseDPW( ParameterState state, ParameterType width = 0.5f, ParameterType initialPhase = 0.f ):
        DPWBase( state ),
        mPhase0( initialPhase ),                                      mLastValue0( mPhase0 ),
        mPhase1( wrap(initialPhase + width + width ) ), mLastValue1( mPhase1 )
    {}

    // DSP
    template < typename OutputFunctor >
    inline void run ( OutputFunctor & out, size_t count )
    {
        run( out, count, getState() );
    }


    template < typename OutputFunctor, typename State >
    inline void run ( OutputFunctor & out, size_t count, State && a )
    {
        ParameterType phase0     = mPhase0;
        ParameterType phase1     = mPhase1;
        ParameterType lastValue0 = mLastValue0;
        ParameterType lastValue1 = mLastValue1;

        const size_t unroll2 = count / 2;
        const size_t remain  = count & 1;

        for (size_t i = 0; i != unroll2; ++i) {
            ParameterState state = a();
            float y0 = tick( phase0, lastValue0, phase1, lastValue1, state[PhaseIncrement], state[Scale] );

            state = a();
            float y1 = tick( phase0, lastValue0, phase1, lastValue1, state[PhaseIncrement], state[Scale] );

            out(y0);
            out(y1);
        }

		if( BOOST_UNLIKELY( remain ) ) {
			for (size_t i = 0; i != remain; ++i) {
				ParameterState state = a();
				float y0 = tick( phase0, lastValue0, phase1, lastValue1, state[PhaseIncrement], state[Scale] );
				out(y0);
			}
		}

        mPhase0     = phase0;
        mPhase1     = phase1;
        mLastValue0 = lastValue0;
        mLastValue1 = lastValue1;
    }


    static inline float tick(ParameterType & phase0, ParameterType & lastVal0,
                             ParameterType & phase1, ParameterType & lastVal1,
                             float phaseIncrement, float scale)
    {
        phase0 = incrementPhase(phase0, phaseIncrement);
        phase1 = incrementPhase(phase1, phaseIncrement);

        // squared saw
        ParameterType val0 = phase0 * phase0;
        ParameterType val1 = phase1 * phase1;

        // differentiate parabolic wave
        ParameterType differentiatedPhase0 = val0 - lastVal0;
        ParameterType differentiatedPhase1 = val1 - lastVal1;

        // substracting saws to get a pulse
        ParameterType pulse = (differentiatedPhase0 - differentiatedPhase1) * scale;

        lastVal0 = val0;
        lastVal1 = val1;

        return pulse;
    }

    template <typename Type>
    static inline Type incrementPhase(Type phase, Type phaseIncrement)
    {
        phase += phaseIncrement;

        return boost::simd::if_else( boost::simd::is_greater_equal( phase, Type(1.f)),
                                     phase - Type(2.f),
                                     phase );
    }

protected:
    ParameterType mPhase0 = 0.f, mLastValue0 = 0.f;
    ParameterType mPhase1 = 0.f, mLastValue1 = 0.f;
};


}

#endif // OSC_DPW_HPP
