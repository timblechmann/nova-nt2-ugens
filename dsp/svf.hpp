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

#include <boost/simd/include/constants/pi.hpp>
#include <boost/simd/include/constants/two.hpp>

#include <boost/simd/sdk/simd/logical.hpp>
#include <boost/simd/include/functions/compare_equal.hpp>
#include <boost/simd/include/functions/eq.hpp>
#include <boost/simd/include/functions/negate.hpp>
#include <boost/simd/include/functions/fast_divides.hpp>

#include <approximations/tan.hpp>

#include <dsp/arithmeticarray.hpp>
#include <dsp/utils.hpp>


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace nova {

enum SVFType {
    svfLPF,
    svfHPF,
    svfBPF,
    svfBRF,
    svfPeak,
    svfLowShelf,
    svfHighShelf,
    svfEq
};


/* adapted from Andrew Simper, Cytomic, 2013, andy@cytomic.com
 * http://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf
 * */
template <size_t Size, size_t ParameterSize, SVFType FilterType>
struct SVFFilter
{
    using SampleType    = typename nova::as_pack<float, Size>::type;
    using ParameterType = typename nova::as_pack<float, ParameterSize>::type;

    typedef nova::detail::ArithmeticArray<ParameterType, 4> ParameterState;

    SVFFilter() {}

    void setState( ParameterState const & newState ) { _state = newState; }
    auto getState()
    {
        auto ret = [=] { return _state; };
        return ret;
    }
    ParameterState currentState() { return _state; }


    // parameters
    template <typename DSPContext>
    static auto computeState( ParameterType cutoff, ParameterType resonance, DSPContext const & dspContext )
    {
        using namespace boost::simd;

        typedef decltype(dspContext.sampleRate()) SampleRateType;
        auto cutoffFreq = nova::clip( (SampleRateType)cutoff, SampleRateType(0.01f), dspContext.sampleRate() * (SampleRateType)(0.5));

        ParameterType normalizedFrequency = cutoffFreq * dspContext.sampleDur();
        ParameterType g = nova::approximations::tan( normalizedFrequency * Pi<ParameterType>(),
                                                     nova::approximations::TanFast() );

        ParameterType two = boost::simd::Two<ParameterType>();
        ParameterType k   = two - two * resonance;

        ParameterType a1 = fast_rec(1 + g * (g + k));
        ParameterType a2 = g*a1;
        ParameterType a3 = g*a2;

        return ParameterState{ a1, a2, a3, k };
    }

    // DSP
    template < typename InputFunctor, typename OutputFunctor >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count )
    {
        doRun<true>( in, out, count, getState() );
    }

    template < typename InputFunctor, typename OutputFunctor, typename State >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count, State && a )
    {
        doRun<false>( in, out, count, getState() );
    }

private:
    template < bool StateImmutable, typename InputFunctor, typename OutputFunctor, typename State >
    inline void doRun ( InputFunctor & in, OutputFunctor & out, size_t count, State && state )
    {
        SampleType two = boost::simd::Two<SampleType>();

        ParameterType a1, a2, a3, k;

        if( StateImmutable ) {
            ParameterState a = state();
            a1 = a[0], a2 = a[1], a3 = a[2], k = a[3];
        }

        SampleType ic1eq = ic1eq_, ic2eq = ic2eq_;
        for( size_t i = 0; i != count; ++i ) {
            if( !StateImmutable ) {
                ParameterState a = state();
                a1 = a[0], a2 = a[1], a3 = a[2]; k = a[3];
            }

            SampleType v0 = in();
            SampleType v3 = v0 - ic2eq;
            SampleType v1 =         a1*ic1eq + a2*v3;
            SampleType v2 = ic2eq + a2*ic1eq + a3*v3;
            ic1eq         = two * v1 - ic1eq;
            ic2eq         = two * v2 - ic2eq;

            SampleType result;

            switch( FilterType ) {
            case svfLPF:  result = v2;                 break;
            case svfHPF:  result = v0 - k*v1 - v2;     break;
            case svfBPF:  result = v1;                 break;
            case svfBRF:  result = v0 - k*v1;          break;
            case svfPeak: result = v0 - k*v1 - two*v2; break;
            default: assert(false);
            }

            out( result );
        }

        ic1eq_ = ic1eq;
        ic2eq_ = ic2eq;
    }

    // state
    SampleType ic1eq_ = SampleType{0}, ic2eq_ =  SampleType{0};

    ParameterState _state; // a1_, a2_, a3_, k;
};


template <size_t Size, size_t ParameterSize>
struct SVFEQ
{
    using SampleType    = typename nova::as_pack<float, Size>::type;
    using ParameterType = typename nova::as_pack<float, ParameterSize>::type;

    typedef nova::detail::ArithmeticArray<ParameterType, 4> ParameterState;

    SVFEQ() {}

    void setState( ParameterState const & newState ) { _state = newState; }
    auto getState()
    {
        auto ret = [=] { return _state; };
        return ret;
    }
    ParameterState currentState() { return _state; }


    // parameters
    template <typename DSPContext>
    static auto computeState( ParameterType cutoff, ParameterType amp, ParameterType resonance, DSPContext const & dspContext )
    {
        using namespace boost::simd;

        typedef decltype(dspContext.sampleRate()) SampleRateType;
        auto cutoffFreq = nova::clip( (SampleRateType)cutoff, SampleRateType(0.01f), dspContext.sampleRate() * (SampleRateType)(0.5));

        ParameterType normalizedFrequency = cutoffFreq * dspContext.sampleDur();
        ParameterType g = nova::approximations::tan( normalizedFrequency * Pi<ParameterType>(),
                                                     nova::approximations::TanFast() );

        ParameterType A = amp;

        ParameterType two = boost::simd::Two<ParameterType>();
        ParameterType k   = (two - two * resonance) * fast_rec(A);

        ParameterType a1 = fast_rec(1 + g * (g + k));
        ParameterType a2 = g*a1;
        ParameterType a3 = g*a2;
        ParameterType m1 = k * ( A * A + One<ParameterType>() );

        return ParameterState{ a1, a2, a3, m1 };
    }

    // DSP
    template < typename InputFunctor, typename OutputFunctor >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count )
    {
        doRun<true>( in, out, count, getState() );
    }

    template < typename InputFunctor, typename OutputFunctor, typename State >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count, State && a )
    {
        doRun<false>( in, out, count, getState() );
    }

private:
    template < bool StateImmutable, typename InputFunctor, typename OutputFunctor, typename State >
    inline void doRun ( InputFunctor & in, OutputFunctor & out, size_t count, State && state )
    {
        SampleType two = boost::simd::Two<SampleType>();

        ParameterType a1, a2, a3, m1;

        if( StateImmutable ) {
            ParameterState a = state();
            a1 = a[0], a2 = a[1], a3 = a[2], m1 = a[3];
        }

        SampleType ic1eq = ic1eq_, ic2eq = ic2eq_;
        for( size_t i = 0; i != count; ++i ) {
            if( !StateImmutable ) {
                ParameterState a = state();
                a1 = a[0], a2 = a[1], a3 = a[2]; m1 = a[3];
            }

            SampleType v0 = in();
            SampleType v3 = v0 - ic2eq;
            SampleType v1 =         a1*ic1eq + a2*v3;
            SampleType v2 = ic2eq + a2*ic1eq + a3*v3;
            ic1eq         = two * v1 - ic1eq;
            ic2eq         = two * v2 - ic2eq;

            SampleType result = v0 + m1*v1;

            out( result );
        }

        ic1eq_ = ic1eq;
        ic2eq_ = ic2eq;
    }

    // state
    SampleType ic1eq_ = SampleType{0}, ic2eq_ =  SampleType{0};

    ParameterState _state; // a1_, a2_, a3_, m1;
};

template <size_t Size, size_t ParameterSize>
struct SVFLowShelf
{
    using SampleType    = typename nova::as_pack<float, Size>::type;
    using ParameterType = typename nova::as_pack<float, ParameterSize>::type;

    typedef nova::detail::ArithmeticArray<ParameterType, 5> ParameterState;

    SVFLowShelf() {}

    void setState( ParameterState const & newState ) { _state = newState; }
    auto getState()
    {
        auto ret = [=] { return _state; };
        return ret;
    }
    ParameterState currentState() { return _state; }


    // parameters
    template <typename DSPContext>
    static auto computeState( ParameterType cutoff, ParameterType amp, ParameterType resonance, DSPContext const & dspContext )
    {
        using namespace boost::simd;

        typedef decltype(dspContext.sampleRate()) SampleRateType;
        auto cutoffFreq = nova::clip( (SampleRateType)cutoff, SampleRateType(0.01f), dspContext.sampleRate() * (SampleRateType)(0.5));

        ParameterType normalizedFrequency = cutoffFreq * dspContext.sampleDur();
        ParameterType g = nova::approximations::tan( normalizedFrequency * Pi<ParameterType>(),
                                                     nova::approximations::TanFast() );

        ParameterType A = amp;

        ParameterType two = boost::simd::Two<ParameterType>();
        resonance = nova::clip(resonance, 0.f, 1.f );
        ParameterType k   = (two - two * resonance) * fast_rec(A);

        ParameterType a1 = fast_rec(1 + g * (g + k));
        ParameterType a2 = g*a1;
        ParameterType a3 = g*a2;
        ParameterType m1 = k * (     A + One<ParameterType>() );
        ParameterType m2 =     ( A * A + One<ParameterType>() );

        return ParameterState{ a1, a2, a3, m1, m2 };
    }

    // DSP
    template < typename InputFunctor, typename OutputFunctor >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count )
    {
        doRun<true>( in, out, count, getState() );
    }

    template < typename InputFunctor, typename OutputFunctor, typename State >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count, State && a )
    {
        doRun<false>( in, out, count, getState() );
    }

private:
    template < bool StateImmutable, typename InputFunctor, typename OutputFunctor, typename State >
    inline void doRun ( InputFunctor & in, OutputFunctor & out, size_t count, State && state )
    {
        SampleType two = boost::simd::Two<SampleType>();

        ParameterType a1, a2, a3, m1, m2;

        if( StateImmutable ) {
            ParameterState a = state();
            a1 = a[0], a2 = a[1], a3 = a[2], m1 = a[3], m2 = a[4];
        }

        SampleType ic1eq = ic1eq_, ic2eq = ic2eq_;
        for( size_t i = 0; i != count; ++i ) {
            if( !StateImmutable ) {
                ParameterState a = state();
                a1 = a[0], a2 = a[1], a3 = a[2]; m1 = a[3], m2 = a[4];
            }

            SampleType v0 = in();
            SampleType v3 = v0 - ic2eq;
            SampleType v1 =         a1*ic1eq + a2*v3;
            SampleType v2 = ic2eq + a2*ic1eq + a3*v3;
            ic1eq         = two * v1 - ic1eq;
            ic2eq         = two * v2 - ic2eq;

            SampleType result = v0 + m1*v1 + m2*v2;
            out( result );
        }

        ic1eq_ = ic1eq;
        ic2eq_ = ic2eq;
    }

    // state
    SampleType ic1eq_ = SampleType{0}, ic2eq_ =  SampleType{0};
    ParameterState _state; // a1_, a2_, a3_, m1, m2;
};

template <size_t Size, size_t ParameterSize>
struct SVFHighShelf
{
    using SampleType    = typename nova::as_pack<float, Size>::type;
    using ParameterType = typename nova::as_pack<float, ParameterSize>::type;

    typedef nova::detail::ArithmeticArray<ParameterType, 6> ParameterState;

    SVFHighShelf() {}

    void setState( ParameterState const & newState ) { _state = newState; }
    auto getState()
    {
        auto ret = [=] { return _state; };
        return ret;
    }
    ParameterState currentState() { return _state; }


    // parameters
    template <typename DSPContext>
    static auto computeState( ParameterType cutoff, ParameterType amp, ParameterType resonance, DSPContext const & dspContext )
    {
        using namespace boost::simd;

        typedef decltype(dspContext.sampleRate()) SampleRateType;
        auto cutoffFreq = nova::clip( (SampleRateType)cutoff, SampleRateType(0.01f), dspContext.sampleRate() * (SampleRateType)(0.5));

        ParameterType normalizedFrequency = cutoffFreq * dspContext.sampleDur();
        ParameterType g = nova::approximations::tan( normalizedFrequency * Pi<ParameterType>(),
                                                     nova::approximations::TanFast() );

        ParameterType A = amp;

        ParameterType two = boost::simd::Two<ParameterType>();
        ParameterType k   = (two - two * resonance) * fast_rec(A);

        ParameterType a1 = fast_rec(1 + g * (g + k));
        ParameterType a2 = g*a1;
        ParameterType a3 = g*a2;

        ParameterType m0 = A*A;
        ParameterType m1 = k * (One<ParameterType>() - A) * A;
        ParameterType m2 = ( One<ParameterType>() - A*A);

        return ParameterState{ a1, a2, a3, m0, m1, m2 };
    }

    // DSP
    template < typename InputFunctor, typename OutputFunctor >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count )
    {
        doRun<true>( in, out, count, getState() );
    }

    template < typename InputFunctor, typename OutputFunctor, typename State >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count, State && a )
    {
        doRun<false>( in, out, count, getState() );
    }

private:
    template < bool StateImmutable, typename InputFunctor, typename OutputFunctor, typename State >
    inline void doRun ( InputFunctor & in, OutputFunctor & out, size_t count, State && state )
    {
        SampleType two = boost::simd::Two<SampleType>();

        ParameterType a1, a2, a3, m0, m1, m2;

        if( StateImmutable ) {
            ParameterState a = state();
            a1 = a[0], a2 = a[1], a3 = a[2], m0 = a[3], m1 = a[4], m1 = a[5];
        }

        SampleType ic1eq = ic1eq_, ic2eq = ic2eq_;
        for( size_t i = 0; i != count; ++i ) {
            if( !StateImmutable ) {
                ParameterState a = state();
                a1 = a[0], a2 = a[1], a3 = a[2], m0 = a[3], m1 = a[4], m1 = a[5];
            }

            SampleType v0 = in();
            SampleType v3 = v0 - ic2eq;
            SampleType v1 =         a1*ic1eq + a2*v3;
            SampleType v2 = ic2eq + a2*ic1eq + a3*v3;
            ic1eq         = two * v1 - ic1eq;
            ic2eq         = two * v2 - ic2eq;

            SampleType result = m0*v0 + m1*v1 + m2*v2;
            out( result );
        }

        ic1eq_ = ic1eq;
        ic2eq_ = ic2eq;
    }

    // state
    SampleType ic1eq_ = SampleType{0}, ic2eq_ =  SampleType{0};
    ParameterState _state; // a1_, a2_, a3_, m0, m1, m2;
};

}
