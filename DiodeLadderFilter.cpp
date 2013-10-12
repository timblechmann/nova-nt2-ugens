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

#include "dsp/utils.hpp"

#include <boost/simd/include/functions/fast_divides.hpp>
#include <boost/simd/include/functions/fast_rec.hpp>

static InterfaceTable *ft;

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

	void getParameters(InternalType & k, InternalType & A)
	{
		k = k_;
		A = A_;
	}

	InternalType k_, A_;
};

template <typename InternalType>
struct QParameter<slope, InternalType>
{
	QParameter(InternalType k, InternalType kSlope, InternalType A, InternalType ASlope):
		k_(k), kSlope_(kSlope), A_(A), ASlope_(ASlope)
	{}

	void getParameters(InternalType & k, InternalType & A)
	{
		k = k_;
		A = A_;
		k_ += kSlope_;
		A_ += ASlope_;
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

	void getParameters(InternalType & a, InternalType & b)
	{
		a = a_;
		b = b_;
	}

	InternalType a_, b_;
};

template <typename InternalType>
struct HPFParameter<slope, InternalType>
{
	HPFParameter(InternalType a, InternalType aSlope, InternalType b, InternalType bSlope):
		a_(a), b_(b), aSlope(aSlope), bSlope(bSlope)
	{}

	void getParameters(InternalType & a, InternalType & b)
	{
		a = a_;
		b = b_;
		a_ += aSlope;
		b_ += bSlope;
	}

	InternalType a_, b_;
	InternalType aSlope, bSlope;
};



// TODO: add some oversampling
struct DiodeLadderFilter:
		public SCUnit
{
	typedef float  ParameterType;
	typedef double InternalType;

	static const size_t numberOfChannels    = 1;
	static const size_t parameterInputSize  = 1;

	static const size_t freqInputIndex = numberOfChannels + 1;
	static const size_t qInputIndex    = freqInputIndex   + parameterInputSize;
	static const size_t hpfInputIndex  = qInputIndex      + parameterInputSize;


	int inRate(int start, int end)
	{
		// dummy
		return SCUnit::inRate(start);
	}

public:
	DiodeLadderFilter()
	{
		std::fill(z, z + 5, boost::simd::Zero<InternalType>());

		_freq                      = in0(freqInputIndex);
		ParameterType newQ         = in0(qInputIndex);
		ParameterType newHPCutoff  = in0(hpfInputIndex);
		set_q(newQ);

		setFeedbackHPF(newHPCutoff * sampleDur());
		_hpCutoff = newHPCutoff;

		const auto freqRate = inRate(freqInputIndex, parameterInputSize);
		const auto qRate    = inRate(qInputIndex,    parameterInputSize);
		const auto hpfRate  = inRate(hpfInputIndex,  parameterInputSize);


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
	void next_(int inNumSamples, QParam & qParameter, HPFParam & hpfParam)
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
			QParameter<scalar, InternalType> qParam(m_k, m_A);
			next_selectHPF<AudioRateFrequency, HPFRate>(inNumSamples, qParam);
			return;
		}

		default: {
			ParameterType newQ = in0( qInputIndex );
			newQ = nova::clip(newQ, Zero<ParameterType>(), One<ParameterType>());

			if (newQ != _q) {
				InternalType oldA = m_A, oldk = m_k;
				InternalType newA, newk;

				update_kA(newQ, newk, newA);
				m_A = newA;
				m_k = newk;
				_q = newQ;

				InternalType kSlope = calcSlope(newk, oldk);
				InternalType ASlope = calcSlope(newA, oldA);

				QParameter<slope, InternalType> qParam(oldk, kSlope, oldA, ASlope);
				next_selectHPF<AudioRateFrequency, HPFRate>(inNumSamples, qParam);
			} else {
				QParameter<scalar, InternalType> qParam(m_k, m_A);
				next_selectHPF<AudioRateFrequency, HPFRate>(inNumSamples, qParam);
			}
			return;
		}
		}
	}

	template <bool AudioRateFrequency, int HPFRate, typename QParam>
	void next_selectHPF(int inNumSamples, QParam & qParam)
	{
		switch (HPFRate) {
		case calc_ScalarRate: {
			HPFParameter<scalar, InternalType> hpfParam(m_ah, m_bh);
			next_<AudioRateFrequency>(inNumSamples, qParam, hpfParam);
			return;
		}

		default: {
			ParameterType newHPCutoff    = in0( hpfInputIndex );
			if (newHPCutoff != _hpCutoff) {
				InternalType oldA = m_ah, oldB = m_bh;
				InternalType newA, newB;
				InternalType fc = newHPCutoff * sampleDur();

				feedbackHPF(fc, newA, newB);
				m_ah = newA;
				m_bh = newB;
				_hpCutoff = newHPCutoff;

				InternalType slopeA = calcSlope(newA, oldA);
				InternalType slopeB = calcSlope(newB, oldB);

				HPFParameter<slope, InternalType> hpfParam(oldA, slopeA, oldB, slopeB);
				next_<AudioRateFrequency>(inNumSamples, qParam, hpfParam);
			} else {
				HPFParameter<scalar, InternalType> hpfParam(m_ah, m_bh);
				next_<AudioRateFrequency>(inNumSamples, qParam, hpfParam);
			}
		}
		}
	}

	template <typename QParam, typename HPFParam>
	void next_k(int inNumSamples, QParam & qParameter, HPFParam & hpfParam)
	{
		InternalType a, a_inv, a2, b, b2, c, g, g0;

		InternalType z0 = z[0];
		InternalType z1 = z[1];
		InternalType z2 = z[2];
		InternalType z3 = z[3];
		InternalType z4 = z[4];

		const float * inSig = zin(0);
		float * outSig = zout(0);

		ParameterType newFreq = in0( freqInputIndex );

		if (newFreq == _freq) {
			calcFilterCoefficients(_freq, a, a2, a_inv, b, b2, c, g, g0);

			loop(inNumSamples, [&] {
				InternalType x = ZXP(inSig);

				ZXP(outSig) = tick(x, a, a2, a_inv, b, b2, c, g, g0, z0, z1, z2, z3, z4, qParameter, hpfParam);
			});
		} else {
			ParameterType oldfreq = _freq;
			ParameterType freqSlope = calcSlope(newFreq, oldfreq);
			_freq = newFreq;

			loop(inNumSamples, [&] {
				InternalType x = ZXP(inSig);
				calcFilterCoefficients(oldfreq, a, a2, a_inv, b, b2, c, g, g0);
				oldfreq += freqSlope;

				ZXP(outSig) = tick(x, a, a2, a_inv, b, b2, c, g, g0, z0, z1, z2, z3, z4, qParameter, hpfParam);
			});
		}

		z[0] = z0;
		z[1] = z1;
		z[2] = z2;
		z[3] = z3;
		z[4] = z4;
	}

	template <typename QParam, typename HPFParam>
	void next_a(int inNumSamples, QParam & qParameter, HPFParam & hpfParam)
	{
		const float * inSig = zin(0);
		const ParameterType * inFreq = zin( freqInputIndex );
		float * outSig = zout(0);

		InternalType z0 = z[0];
		InternalType z1 = z[1];
		InternalType z2 = z[2];
		InternalType z3 = z[3];
		InternalType z4 = z[4];
		for (int i = 0; i != inNumSamples; ++i) {
			ParameterType freq = ZXP(inFreq);

			InternalType a, a_inv, a2, b, b2, c, g, g0;
			calcFilterCoefficients(freq, a, a2, a_inv, b, b2, c, g, g0);

			InternalType x = ZXP(inSig);
			ZXP(outSig) = tick(x, a, a2, a_inv, b, b2, c, g, g0, z0, z1, z2, z3, z4, qParameter, hpfParam);
		}
		z[0] = z0;
		z[1] = z1;
		z[2] = z2;
		z[3] = z3;
		z[4] = z4;
	}

	template <typename QParam, typename HPFParam>
	static InternalType tick(InternalType x, InternalType a, InternalType a2, InternalType ainv, InternalType b, InternalType b2, InternalType c, InternalType g, InternalType g0,
							 InternalType & z0, InternalType & z1, InternalType & z2, InternalType & z3, InternalType & z4, QParam & qParameter, HPFParam & hpfParam)
	{
		using namespace boost::simd;

		InternalType k, A;
		qParameter.getParameters(k, A);

		InternalType ah, bh;
		hpfParam.getParameters(ah, bh);

		// current state
		const auto s0 = (a2*a*z0 + a2*b*z1 + z2 * (b2 - 2*a2) * a + z3 * (b2 - 3 * a2) * b ) * c;
		const auto s = bh * s0 - z4;

		// solve feedback loop (linear)
		InternalType y5 = fast_div( g*x + s, One<InternalType>() + g*k );

		// input clipping
		const auto y0 = saturate(x - k*y5);
		y5 = g*y0 + s;

		// compute integrator outputs
		const auto y4 = g0*y0 + s0;
		const auto y3 = (b*y4 - z3) * ainv;
		const auto y2 = (b*y3 - a*y4 - z2) * ainv;
		const auto y1 = (b*y2 - a*y3 - z1) * ainv;

		const auto two_a = a * 2;

		// update filter state
		z0 += 2 * two_a * (y0 - y1 + y2);
		z1 +=     two_a * (y1 - 2*y2 + y3);
		z2 +=     two_a * (y2 - 2*y3 + y4);
		z3 +=     two_a * (y3 - 2*y4);
		z4 = bh*y4 + ah*y5;

		InternalType result = A * y4;

		return result;
	}

	void calcFilterCoefficients(const ParameterType newFreq, InternalType & a, InternalType & a2, InternalType & a_inv,
								InternalType & b, InternalType & b2, InternalType & c, InternalType & g, InternalType & g0)
	{
		using namespace boost::simd;

		InternalType fc = max(ParameterType(10), newFreq) * sampleDur();
		fc = min(fc, InternalType(0.25) );
		a = Pi<InternalType>() * fc; // PI is Nyquist frequency
		//		a = 2 * tan(0.5*a); // dewarping, not required with 2x oversampling

		a_inv = fast_rec( a );

		a2 = a*a;
		b = 2*a + 1;
		b2 = b*b;
		c = fast_rec ( 2*a2*a2 - 4*a2*b2 + b2*b2 );
		g0 = 2*a2*a2*c;
		g = g0 * m_bh;
	}

	void setFeedbackHPF(InternalType fc)
	{
		feedbackHPF(fc, m_ah, m_bh);
	}

	static void feedbackHPF(InternalType fc, InternalType & ah, InternalType & bh)
	{
		using namespace boost::simd;
		const InternalType K = fc * M_PI;

		auto rec_k_2 = fast_rec( K + 2 );
		ah = ( K - 2 ) * rec_k_2;
		bh =         2 * rec_k_2;
	}

	void set_q(InternalType q_)
	{
		using namespace boost::simd;
		_q = nova::clip(q_, Zero<InternalType>(), One<InternalType>());
		update_kA(q_, m_k, m_A);
	}

	static void update_kA(InternalType q, InternalType & k, InternalType & A)
	{
		k = 20 * q;
		A = 1 + 0.5*k; // resonance gain compensation
	}

	static InternalType saturate (InternalType sample)
	{
		using namespace boost::simd;
		return fast_div (sample, 1 + abs(sample) );
	}

	ParameterType _freq, _q, _hpCutoff;
	InternalType m_k, m_A;
	InternalType z[5];
	InternalType m_ah, m_bh;
};

DEFINE_XTORS(DiodeLadderFilter)

PluginLoad(NovaFilters)
{
	ft = inTable;
	DefineSimpleUnit(DiodeLadderFilter);
}
