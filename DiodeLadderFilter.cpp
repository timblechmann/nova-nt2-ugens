/*
 *
 *    Copyright (C) 2013 Tim Blechmann
 *
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "SC_PlugIn.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

#include <dsp/utils.hpp>
#include <dsp/saturators.hpp>

#include "producer_consumer_functors.hpp"
#include "NovaUGensCommon.hpp"

#include <boost/simd/include/functions/fast_divides.hpp>
#include <boost/simd/include/functions/fast_rec.hpp>
#include <boost/simd/include/functions/compare_equal.hpp>
#include <boost/simd/sdk/meta/as_logical.hpp>

static InterfaceTable *ft;

namespace nova {

// adapted from DiodeLadderFilter
//
// Copyright (c) 2012 Dominique Wurtz (www.blaukraut.info)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.



// TODO: add some oversampling
template <size_t Channels, bool hasScalarArguments>
struct DiodeLadderFilter:
    public NovaUnit,
    public nova::multichannel::SignalInput< DiodeLadderFilter<Channels, hasScalarArguments>, 0, Channels >,
    public nova::multichannel::OutputSink<  DiodeLadderFilter<Channels, hasScalarArguments>, 0, Channels >,


    // Controls:
    public nova::multichannel::ControlInput<  DiodeLadderFilter<Channels, hasScalarArguments>, Channels,                                           (hasScalarArguments ? 1 : Channels) >,
    public nova::multichannel::ControlInput<  DiodeLadderFilter<Channels, hasScalarArguments>, Channels + (hasScalarArguments ? 1 : Channels) * 1, (hasScalarArguments ? 1 : Channels) >,
    public nova::multichannel::ControlInput<  DiodeLadderFilter<Channels, hasScalarArguments>, Channels + (hasScalarArguments ? 1 : Channels) * 2, (hasScalarArguments ? 1 : Channels) >
{
    static const size_t numberOfChannels    = Channels;
    static const size_t parameterInputSize  = hasScalarArguments ? 1 : Channels;

    typedef typename nova::as_pack<float,  parameterInputSize>::type ParameterType;
#if 0
    typedef typename nova::as_pack<double, numberOfChannels>::type   InternalType;
    typedef typename nova::as_pack<double, parameterInputSize>::type InternalParameterType;
#else
    typedef typename nova::as_pack<float, numberOfChannels>::type   InternalType;
    typedef typename nova::as_pack<float, parameterInputSize>::type InternalParameterType;
#endif

    static const size_t freqInputIndex = numberOfChannels;
    static const size_t qInputIndex    = freqInputIndex   + parameterInputSize;
    static const size_t hpfInputIndex  = qInputIndex      + parameterInputSize;

    using SignalInput = nova::multichannel::SignalInput< DiodeLadderFilter<Channels, hasScalarArguments>, 0, Channels >;
    using OutputSink  = nova::multichannel::OutputSink<  DiodeLadderFilter<Channels, hasScalarArguments>, 0, Channels >;

    using FreqInput = nova::multichannel::ControlInput<  DiodeLadderFilter<Channels, hasScalarArguments>, Channels,                                           (hasScalarArguments ? 1 : Channels) >;
    using QInput    = nova::multichannel::ControlInput<  DiodeLadderFilter<Channels, hasScalarArguments>, Channels + (hasScalarArguments ? 1 : Channels) * 1, (hasScalarArguments ? 1 : Channels) >;
    using HPFInput  = nova::multichannel::ControlInput<  DiodeLadderFilter<Channels, hasScalarArguments>, Channels + (hasScalarArguments ? 1 : Channels) * 2, (hasScalarArguments ? 1 : Channels) >;


    // later: factor this out:
    typedef boost::simd::aligned_array<InternalParameterType, 2, 4> FeedbackHPFParameterState;
    enum {
        stateAH,
        stateBH
    };

    typedef boost::simd::aligned_array<InternalParameterType, 2, 4> QParameterState;
    enum {
        stateK,
        stateA
    };

public:
    DiodeLadderFilter()
    {
        _qParameterState  = calculateQControls( QInput::readInput() );
        _feedbackHPFState = calculateFeedbackHPFControls( HPFInput::readInput() * sampleDur() );

        const auto freqRate = inRate(freqInputIndex, freqInputIndex + parameterInputSize);
        const auto qRate    = inRate(qInputIndex,    qInputIndex    + parameterInputSize);
        const auto hpfRate  = inRate(hpfInputIndex,  hpfInputIndex  + parameterInputSize);

        if (freqRate == calc_FullRate) {
            switch ( qRate ) {
            case calc_ScalarRate: {
                switch (hpfRate) {
                case calc_ScalarRate:
                    set_calc_function<DiodeLadderFilter, &DiodeLadderFilter::next<true, calc_ScalarRate, calc_ScalarRate> >();
                    return;

                default:
                    set_calc_function<DiodeLadderFilter, &DiodeLadderFilter::next<true, calc_ScalarRate, calc_BufRate> >();
                    return;
                }
            }

            default: {
                switch ( hpfRate ) {
                case calc_ScalarRate:
                    set_calc_function<DiodeLadderFilter, &DiodeLadderFilter::next<true, calc_BufRate, calc_ScalarRate> >();
                    return;

                default:
                    set_calc_function<DiodeLadderFilter, &DiodeLadderFilter::next<true, calc_BufRate, calc_BufRate> >();
                    return;
                }
            }
            }
        } else {
            switch ( qRate ) {
            case calc_ScalarRate: {
                switch ( hpfRate ) {
                case calc_ScalarRate:
                    set_calc_function<DiodeLadderFilter, &DiodeLadderFilter::next<false, calc_ScalarRate, calc_ScalarRate> >();
                    return;

                default:
                    set_calc_function<DiodeLadderFilter, &DiodeLadderFilter::next<false, calc_ScalarRate, calc_BufRate> >();
                    return;
                }
            }

            default: {
                switch ( hpfRate ) {
                case calc_ScalarRate:
                    set_calc_function<DiodeLadderFilter, &DiodeLadderFilter::next<false, calc_BufRate, calc_ScalarRate> >();
                    return;

                default:
                    set_calc_function<DiodeLadderFilter, &DiodeLadderFilter::next<false, calc_BufRate, calc_BufRate> >();
                    return;
                }
            }
            }
        }
    }

private:
    template <bool AudioRateFrequency, typename QParam, typename HPFParam>
    void next_(int inNumSamples, QParam && __restrict__ qParameter, HPFParam && __restrict__ hpfParam)
    {
        if (AudioRateFrequency)
            next_a(inNumSamples, qParameter, hpfParam);
        else
            next_k(inNumSamples, qParameter, hpfParam);
    }

    template <bool AudioRateFrequency, int QRate, int HPFRate>
    void next(int inNumSamples)
    {
        using namespace boost::simd;
        switch ( QRate ) {
        case calc_ScalarRate: {
            next_selectHPF<AudioRateFrequency, HPFRate>( inNumSamples, makeQScalarParam() );
            return;
        }

        default: {
            if ( QInput::changed() ) {
                QParameterState currentState = _qParameterState;
                QParameterState nextState = calculateQControls( QInput::readAndUpdateInput() );
                _qParameterState = nextState;

                auto qParam = makeRamp( currentState, nextState );
                next_selectHPF<AudioRateFrequency, HPFRate>(inNumSamples, qParam);
            } else {
                next_selectHPF<AudioRateFrequency, HPFRate>(inNumSamples, makeQScalarParam() );
            }
            return;
        }
        }
    }

    template <bool AudioRateFrequency, int HPFRate, typename QParam>
    void next_selectHPF(int inNumSamples, QParam && __restrict__ qParam)
    {
        using namespace boost::simd;

        switch (HPFRate) {
        case calc_ScalarRate: {
            next_<AudioRateFrequency>(inNumSamples, qParam, makeHPFScalarParm());
            return;
        }

        default: {
            if ( HPFInput::changed() ) {
                FeedbackHPFParameterState currentState = _feedbackHPFState;
                InternalParameterType fc = HPFInput::readAndUpdateInput() * sampleDur();
                FeedbackHPFParameterState nextState = calculateFeedbackHPFControls( fc );
                _feedbackHPFState = nextState;

                auto hpfParam = makeRamp( currentState, nextState );
                next_<AudioRateFrequency>(inNumSamples, qParam, hpfParam);
            } else {
                next_<AudioRateFrequency>(inNumSamples, qParam, makeHPFScalarParm() );
            }
        }
        }
    }

    template <typename QParam, typename HPFParam>
    void next_k(int inNumSamples, QParam && __restrict__ qParameter, HPFParam && __restrict__ hpfParam)
    {
        InternalParameterType a, a_inv, a2, b, b2, c, g, g0;

        auto inSig   = SignalInput::template makeInputSignal<InternalType>();
        auto outSink = OutputSink ::template makeSink<InternalType>();

        if ( !FreqInput::changed() && !HPFInput::changed() ) {
            auto hpfState = hpfParam();
            calcFilterCoefficients( FreqInput::readInput(), a, a2, a_inv, b, b2, c, g, g0, hpfState[stateBH]);

            loop(inNumSamples, [&] {
                InternalType x = inSig();

                auto out = tick(x, a, a2, a_inv, b, b2, c, g, g0, qParameter, hpfState);
                outSink( out );
            });
        } else {
            auto freq = FreqInput::template makeRampSignal<ParameterType>();

            loop(inNumSamples, [&] {
                InternalType x = inSig();
                auto hpfState = hpfParam();
                calcFilterCoefficients(freq(), a, a2, a_inv, b, b2, c, g, g0, hpfState[stateBH]);
                auto out = tick(x, a, a2, a_inv, b, b2, c, g, g0, qParameter, hpfState);
                outSink( out );
            });
        }
    }

    template <typename QParam, typename HPFParam>
    void next_a(int inNumSamples, QParam && __restrict__ qParameter, HPFParam && __restrict__ hpfParam)
    {
        auto inSig   = SignalInput::template makeInputSignal<InternalType>();
        auto inFreq  = FreqInput  ::template makeAudioInputSignal<ParameterType>();
        auto outSink = OutputSink ::template makeSink<InternalType>();

        for (int i = 0; i != inNumSamples; ++i) {

            auto hpfState = hpfParam();
            ParameterType freq = inFreq();

            InternalParameterType a, a_inv, a2, b, b2, c, g, g0;
            calcFilterCoefficients(freq, a, a2, a_inv, b, b2, c, g, g0, hpfState[stateBH]);

            InternalType x = inSig();
            InternalType out = tick(x, a, a2, a_inv, b, b2, c, g, g0, qParameter, hpfState);
            outSink(out);
        }
    }

    template <typename QParam>
    BOOST_FORCEINLINE
    InternalType tick(InternalType x, InternalParameterType a, InternalParameterType a2, InternalParameterType ainv,
                      InternalParameterType b, InternalParameterType b2,
                      InternalParameterType c, InternalParameterType g, InternalParameterType g0,
                      QParam && __restrict__ qParameter, FeedbackHPFParameterState & __restrict__ hpfState)
    {
        using namespace boost::simd;

        auto qState = qParameter();
        InternalParameterType k = qState[stateK];
        InternalParameterType A = qState[stateA];

        InternalParameterType bh = hpfState[stateBH];

        auto two = boost::simd::Two<InternalParameterType>();

        // current state
#if 0
        const InternalType s0 = (a2*a*z[0] + a2*b*z[1] + z[2] * (b2 - two*a2) * a + z[3] * (b2 - 3.0 * a2) * b ) * c;
#else
        const InternalType s0 = (a2*a*z[0] + a2*b*z[1] + z[2] * (b2 - (a2+a2)) * a + z[3] * (b2 - 3.0 * a2) * b ) * c;
#endif
        const InternalType s = bh * s0 - z[4];

        // solve feedback loop (linear)
        InternalType y5 = fast_div( g*x + s, One<InternalType>() + g*k );

        // input clipping
        const InternalType y0 = saturator::distort<InternalType>( x - k*y5 );
        y5 = g*y0 + s;

        // compute integrator outputs
        const InternalType y4 = g0*y0 + s0;
        const InternalType y3 = (b*y4 - z[3]) * ainv;
        const InternalType y2 = (b*y3 - a*y4 - z[2]) * ainv;
        const InternalType y1 = (b*y2 - a*y3 - z[1]) * ainv;

#if 0
        const auto two_a = a * two;

        // update filter state
        z[0] += two * two_a * (y0 - y1 + y2);
        z[1] +=       two_a * (y1 - two*y2 + y3);
        z[2] +=       two_a * (y2 - two*y3 + y4);
        z[3] +=       two_a * (y3 - two*y4);
        z[4] = bh*y4 + ah*y5;

#else

        const auto two_a = a + a;

        // update filter state
        z[0] += two * two_a * (y0 -     y1  + y2);
        z[1] +=       two_a * (y1 - (y2+y2) + y3);
        z[2] +=       two_a * (y2 - (y3+y3) + y4);
        z[3] +=       two_a * (y3 - (y4+y4)     );


        InternalParameterType ah = hpfState[stateAH];
        z[4] = bh*y4 + ah*y5;
#endif

        InternalType result = A * y4;

        return result;
    }

    BOOST_FORCEINLINE
    void calcFilterCoefficients(const ParameterType newFreq, InternalParameterType & __restrict__ a,
                                InternalParameterType & __restrict__ a2, InternalParameterType & __restrict__ a_inv,
                                InternalParameterType & __restrict__ b, InternalParameterType & __restrict__ b2,
                                InternalParameterType & __restrict__ c, InternalParameterType & __restrict__ g,
                                InternalParameterType & __restrict__ g0, InternalParameterType & bh)
    {
        using namespace boost::simd;

        InternalParameterType fc = max(ParameterType(10), newFreq) * sampleDur();
        fc = min(fc, InternalParameterType(0.25) );
        a = Pi<InternalParameterType>() * fc; // PI is Nyquist frequency
        //		a = 2 * tan(0.5*a); // dewarping, not required with 2x oversampling

        a_inv = fast_rec( a );

        auto two = Two<InternalParameterType>();

        a2 = a*a;
        b  = two*a + One<InternalParameterType>();
        b2 = b*b;
        c  = fast_rec ( two*a2*a2 - 4.0*a2*b2 + b2*b2 );
        g0 = two*a2*a2*c;
        g  = g0 * bh;
    }

    static auto calculateFeedbackHPFControls( InternalParameterType fc )
    {
        using namespace boost::simd;
        const InternalParameterType K = fc * Pi<InternalParameterType>();

        const InternalParameterType two = Two<InternalParameterType>();

        const InternalParameterType rec_k_2 = fast_rec( K + two );
        const InternalParameterType ah = ( K - two ) * rec_k_2;
        const InternalParameterType bh =         two * rec_k_2;

        return FeedbackHPFParameterState{ ah, bh };
    }

    auto calculateQControls( InternalParameterType q )
    {
        using namespace boost::simd;
        q = nova::clip( q, Zero<InternalParameterType>(), One<InternalParameterType>() );
        auto k = splat<InternalParameterType>(20.0) * q;
        auto A = splat<InternalParameterType>(1.0) + splat<InternalParameterType>(0.5)*k; // resonance gain compensation
        return QParameterState{ k, A };
    }

    auto makeHPFScalarParm() { return [&] { return _feedbackHPFState; }; }
    auto makeQScalarParam()  { return [&] { return _qParameterState;  }; }


    InternalType z[5] = { {0}, {0}, {0}, {0}, {0} };

    QParameterState           _qParameterState;
    FeedbackHPFParameterState _feedbackHPFState;
};

}

PluginLoad(NovaFilters)
{
    using namespace nova;

    ft = inTable;

    registerUnit< nova::DiodeLadderFilter<1, true> >( ft, "DiodeLadderFilter"  );
    registerUnit< nova::DiodeLadderFilter<2, true> >( ft, "DiodeLadderFilter2" );
    registerUnit< nova::DiodeLadderFilter<4, true> >( ft, "DiodeLadderFilter4" );
}
