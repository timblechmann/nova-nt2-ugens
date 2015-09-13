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

#include "dsp/utils.hpp"

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


typedef enum {scalar, slope} ArgType;

template <int ArgType, typename InternalType>
struct QParameter
{};

template <typename InternalType>
struct QParameter<scalar, InternalType>
{
    QParameter(InternalType k, InternalType A):
        k_(k), A_(A)
    {}

    auto getParameters()
    {
        return std::make_tuple( k_, A_ );
    }

    InternalType k_, A_;
};

template <typename InternalType>
struct QParameter<slope, InternalType>
{
    QParameter(InternalType k, InternalType kSlope, InternalType A, InternalType ASlope):
        k_(k), kSlope_(kSlope), A_(A), ASlope_(ASlope)
    {}

    auto getParameters()
    {
        auto ret = std::make_tuple( k_, A_ );
        k_ += kSlope_;
        A_ += ASlope_;
        return ret;
    }

    InternalType k_, A_;
    InternalType kSlope_, ASlope_;
};


template <int ArgType, typename InternalType>
struct HPFParameter
{};


template <typename InternalType>
struct HPFParameter<scalar, InternalType>
{
    HPFParameter(InternalType a, InternalType b):
        a_(a), b_(b)
    {}

    auto getParameters()
    {
        return std::make_tuple( a_, b_ );
    }

    InternalType a_, b_;
};

template <typename InternalType>
struct HPFParameter<slope, InternalType>
{
    HPFParameter(InternalType a, InternalType aSlope, InternalType b, InternalType bSlope):
        a_(a), b_(b), aSlope(aSlope), bSlope(bSlope)
    {}

    auto getParameters()
    {
        auto ret = std::make_tuple( a_, b_ );
        a_ += aSlope;
        b_ += bSlope;
        return ret;
    }

    InternalType a_, b_;
    InternalType aSlope, bSlope;
};



// TODO: add some oversampling
template <size_t Channels, bool hasScalarArguments>
struct DiodeLadderFilter:
    public NovaUnit
{
    static const size_t numberOfChannels    = Channels;
    static const size_t parameterInputSize  = hasScalarArguments ? 1 : Channels;

    typedef typename nova::as_pack<float,  parameterInputSize>::type ParameterType;
    typedef typename nova::as_pack<double, numberOfChannels>::type   InternalType;
    typedef typename nova::as_pack<double, parameterInputSize>::type InternalParameterType;

    typedef nova::InputInterleaver<parameterInputSize> ParameterReader;


    static const size_t freqInputIndex = numberOfChannels;
    static const size_t qInputIndex    = freqInputIndex   + parameterInputSize;
    static const size_t hpfInputIndex  = qInputIndex      + parameterInputSize;

public:
    DiodeLadderFilter()
    {
        _freq            = nova::Interleaver< ParameterType, freqInputIndex >(this)();
        ParameterType newQ        = nova::Interleaver< ParameterType, qInputIndex    >(this)();
        ParameterType newHPCutoff = nova::Interleaver< ParameterType, hpfInputIndex  >(this)();
        set_q(newQ);

        setFeedbackHPF(newHPCutoff * sampleDur());
        _hpCutoff = newHPCutoff;

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
    void next_(int inNumSamples, QParam & __restrict__ qParameter, HPFParam & __restrict__ hpfParam)
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
            QParameter<scalar, InternalParameterType> qParam(m_k, m_A);
            next_selectHPF<AudioRateFrequency, HPFRate>(inNumSamples, qParam);
            return;
        }

        default: {
            ParameterType newQ = ParameterReader::read(this, qInputIndex );
            newQ = nova::clip(newQ, Zero<ParameterType>(), One<ParameterType>());

            auto qConstant = compare_equal( newQ, _q );

            if (!qConstant) {
                InternalParameterType oldA = m_A, oldk = m_k;
                InternalParameterType newA, newk;

                update_kA(newQ, newk, newA);
                m_A = newA;
                m_k = newk;
                _q = newQ;

                InternalParameterType kSlope = calcSlope(newk, oldk);
                InternalParameterType ASlope = calcSlope(newA, oldA);

                QParameter<slope, InternalParameterType> qParam(oldk, kSlope, oldA, ASlope);
                next_selectHPF<AudioRateFrequency, HPFRate>(inNumSamples, qParam);
            } else {
                QParameter<scalar, InternalParameterType> qParam(m_k, m_A);
                next_selectHPF<AudioRateFrequency, HPFRate>(inNumSamples, qParam);
            }
            return;
        }
        }
    }

    template <bool AudioRateFrequency, int HPFRate, typename QParam>
    void next_selectHPF(int inNumSamples, QParam & __restrict__ qParam)
    {
        using namespace boost::simd;

        switch (HPFRate) {
        case calc_ScalarRate: {
            HPFParameter<scalar, InternalParameterType> hpfParam(m_ah, m_bh);
            next_<AudioRateFrequency>(inNumSamples, qParam, hpfParam);
            return;
        }

        default: {
            ParameterType newHPCutoff = ParameterReader::read(this, hpfInputIndex );

            auto hpfConstant = compare_equal( newHPCutoff, _hpCutoff );

            if (!hpfConstant) {
                InternalParameterType oldA = m_ah, oldB = m_bh;
                InternalParameterType newA, newB;
                InternalParameterType fc = newHPCutoff * sampleDur();

                feedbackHPF(fc, newA, newB);
                m_ah = newA;
                m_bh = newB;
                _hpCutoff = newHPCutoff;

                InternalParameterType slopeA = calcSlope(newA, oldA);
                InternalParameterType slopeB = calcSlope(newB, oldB);

                HPFParameter<slope, InternalParameterType> hpfParam(oldA, slopeA, oldB, slopeB);
                next_<AudioRateFrequency>(inNumSamples, qParam, hpfParam);
            } else {
                HPFParameter<scalar, InternalParameterType> hpfParam(m_ah, m_bh);
                next_<AudioRateFrequency>(inNumSamples, qParam, hpfParam);
            }
        }
        }
    }

    template <typename QParam, typename HPFParam>
    void next_k(int inNumSamples, QParam & __restrict__ qParameter, HPFParam & __restrict__ hpfParam)
    {
        InternalParameterType a, a_inv, a2, b, b2, c, g, g0;

        nova::Interleaver<InternalType>  inSig( this );
        nova::Interleaver<ParameterType, freqInputIndex> inFreq( this );
        ParameterType newFreq = inFreq();

        nova::Deinterleaver<InternalType> outSig( this );

        auto freqConstant = boost::simd::compare_equal( newFreq, _freq );

        if (freqConstant) {
            calcFilterCoefficients(_freq, a, a2, a_inv, b, b2, c, g, g0);

            loop(inNumSamples, [&] {
                InternalType x = inSig();

                auto out = tick(x, a, a2, a_inv, b, b2, c, g, g0, qParameter, hpfParam);
                outSig( out );
            });
        } else {
            ParameterType oldfreq = _freq;
            ParameterType freqSlope = calcSlope(newFreq, oldfreq);
            _freq = newFreq;

            loop(inNumSamples, [&] {
                InternalType x = inSig();
                calcFilterCoefficients(oldfreq, a, a2, a_inv, b, b2, c, g, g0);
                oldfreq += freqSlope;

                auto out = tick(x, a, a2, a_inv, b, b2, c, g, g0, qParameter, hpfParam);
                outSig( out );
            });
        }
    }

    template <typename QParam, typename HPFParam>
    void next_a(int inNumSamples, QParam & __restrict__ qParameter, HPFParam & __restrict__ hpfParam)
    {
        nova::Interleaver<InternalType>                  inSig( this );
        nova::Interleaver<ParameterType, freqInputIndex> inFreq( this );
        nova::Deinterleaver<InternalType>                outSig( this );

        for (int i = 0; i != inNumSamples; ++i) {
            ParameterType freq = inFreq();

            InternalParameterType a, a_inv, a2, b, b2, c, g, g0;
            calcFilterCoefficients(freq, a, a2, a_inv, b, b2, c, g, g0);

            InternalType x = inSig();
            InternalType out = tick(x, a, a2, a_inv, b, b2, c, g, g0, qParameter, hpfParam);
            outSig(out);
        }
    }

    template <typename QParam, typename HPFParam>
    InternalType tick(InternalType x, InternalParameterType a, InternalParameterType a2, InternalParameterType ainv,
                      InternalParameterType b, InternalParameterType b2,
                      InternalParameterType c, InternalParameterType g, InternalParameterType g0,
                      QParam & __restrict__ qParameter, HPFParam & __restrict__ hpfParam)
    {
        using namespace boost::simd;

        InternalParameterType k, A;
        std::tie(k, A) = qParameter.getParameters();

        InternalParameterType ah, bh;
        std::tie(ah, bh) = hpfParam.getParameters();
        auto two = boost::simd::Two<InternalParameterType>();

        // current state
        const InternalType s0 = (a2*a*z[0] + a2*b*z[1] + z[2] * (b2 - two*a2) * a + z[3] * (b2 - 3.0 * a2) * b ) * c;
        const InternalType s = bh * s0 - z[4];

        // solve feedback loop (linear)
        InternalType y5 = fast_div( g*x + s, One<InternalType>() + g*k );

        // input clipping
        const InternalType y0 = saturate(x - k*y5);
        y5 = g*y0 + s;

        // compute integrator outputs
        const InternalType y4 = g0*y0 + s0;
        const InternalType y3 = (b*y4 - z[3]) * ainv;
        const InternalType y2 = (b*y3 - a*y4 - z[2]) * ainv;
        const InternalType y1 = (b*y2 - a*y3 - z[1]) * ainv;

        const auto two_a = a * two;

        // update filter state
        z[0] += two * two_a * (y0 - y1 + y2);
        z[1] +=       two_a * (y1 - two*y2 + y3);
        z[2] +=       two_a * (y2 - two*y3 + y4);
        z[3] +=       two_a * (y3 - two*y4);
        z[4] = bh*y4 + ah*y5;

        InternalType result = A * y4;

        return result;
    }

    BOOST_FORCEINLINE
    void calcFilterCoefficients(const ParameterType newFreq, InternalParameterType & __restrict__ a,
                                InternalParameterType & __restrict__ a2, InternalParameterType & __restrict__ a_inv,
                                InternalParameterType & __restrict__ b, InternalParameterType & __restrict__ b2,
                                InternalParameterType & __restrict__ c, InternalParameterType & __restrict__ g,
                                InternalParameterType & __restrict__ g0)
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
        g  = g0 * m_bh;
    }

    void setFeedbackHPF(InternalParameterType fc)
    {
        feedbackHPF(fc, m_ah, m_bh);
    }

    static void feedbackHPF(InternalParameterType fc, InternalParameterType & __restrict__ ah, InternalParameterType & __restrict__ bh)
    {
        using namespace boost::simd;
        const InternalParameterType K = fc * M_PI;

        auto two = Two<InternalParameterType>();

        const InternalParameterType rec_k_2 = fast_rec( K + two );
        ah = ( K - two ) * rec_k_2;
        bh =         two * rec_k_2;
    }

    void set_q(InternalParameterType q_)
    {
        using namespace boost::simd;
        _q = nova::clip(q_, Zero<InternalParameterType>(), One<InternalParameterType>());
        update_kA(q_, m_k, m_A);
    }

    static void update_kA(InternalParameterType q, InternalParameterType & __restrict__ k, InternalParameterType & __restrict__ A)
    {
        k = 20.0 * q;
        A = 1.0 + 0.5*k; // resonance gain compensation
    }

    static InternalType saturate (InternalType sample)
    {
        using namespace boost::simd;
        return fast_div (sample, One<InternalType>() + abs(sample) );
    }

    ParameterType _freq, _q, _hpCutoff;
    InternalParameterType m_k, m_A;
    InternalType z[5] = { {0}, {0}, {0}, {0}, {0} };
    InternalParameterType m_ah, m_bh;
};

}

typedef nova::DiodeLadderFilter<1, true> DiodeLadderFilter;
typedef nova::DiodeLadderFilter<2, true> DiodeLadderFilter2;
typedef nova::DiodeLadderFilter<4, true> DiodeLadderFilter4;
typedef nova::DiodeLadderFilter<4, false> DiodeLadderFilter4_4;

DEFINE_XTORS(DiodeLadderFilter)
DEFINE_XTORS(DiodeLadderFilter2)
DEFINE_XTORS(DiodeLadderFilter4)
//DEFINE_XTORS(DiodeLadderFilter4_4)

PluginLoad(NovaFilters)
{
    ft = inTable;
    DefineSimpleUnit(DiodeLadderFilter);
    DefineSimpleUnit(DiodeLadderFilter2);
    DefineSimpleUnit(DiodeLadderFilter4);
//    DefineSimpleUnit(DiodeLadderFilter4_4);
}
