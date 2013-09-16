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

#include "SC_PlugIn.hpp"

#include "leak_dc.hpp"
#include "biquad.hpp"
#include "producer_consumer_functors.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/simd/constant/constants.hpp>

#include "utils.hpp"

#include <cmath>

namespace {

InterfaceTable *ft;

namespace constants = boost::math::constants;

struct NovaLeakDC2:
	public SCUnit
{
	typedef boost::simd::pack<double, 2> v2d;
	typedef boost::simd::pack<float, 2>  v2f;
	typedef nova::LeakDC<v2d, double> Filter;

	NovaLeakDC2()
	{
		initFilter(in0(2));

		auto inFn  = nova::Interleaver2<v2d>(this);
		_filter._x_1 = inFn();

		switch (inRate(2))
		{
		case calc_ScalarRate:
			set_calc_function<NovaLeakDC2, &NovaLeakDC2::next_i>();
			break;

		case calc_BufRate:
		default:
			_freq = std::numeric_limits<float>::quiet_NaN();
			set_calc_function<NovaLeakDC2, &NovaLeakDC2::next_k>();
		}
	}

	void initFilter(float cutoffFreq)
	{
		_filter.set_a( designFilter(cutoffFreq) );
	}

	float designFilter ( float cutoffFreq )
	{
		cutoffFreq = nova::clip(cutoffFreq, 0.1f, (float)sampleRate());

		float parameter = std::exp( -2.f * constants::pi<float>() * cutoffFreq * (float)sampleDur());
		return parameter;
	}

	void next_i(int inNumSamples)
	{
		auto inFn  = nova::Interleaver2<v2d>(this);
		auto outFn = nova::Deinterleaver2<v2d>(this);

		_filter.run(inFn, outFn, inNumSamples);
	}

	void next_k(int inNumSamples)
	{
		float newFreq = in0(2);
		if (newFreq != _freq) {
			float oldA = _filter._a;
			float newA = designFilter( newFreq );

			float slopeA = calcSlope(newA, oldA);

			_freq = newFreq;

			auto inFn  = nova::Interleaver2<v2d>(this);
			auto outFn = nova::Deinterleaver2<v2d>(this);

			_filter.run(inFn, outFn, inNumSamples, slopeA);
			_filter.set_a( newA );
		} else {
			next_i(inNumSamples);
		}
	}

	Filter _filter;
	float _freq;
};

template <typename ParameterType>
struct BiquadParameterStruct
{
	ParameterType a0, a1, a2, b1, b2;
};

template <typename FilterDesigner>
struct NovaBiquadBase:
	public SCUnit
{
	typedef boost::simd::pack<double, 2> v2d;

	typedef double ParameterType;
	typedef nova::Biquad<v2d, ParameterType> Filter;
	typedef BiquadParameterStruct<ParameterType> BiquadParameterStruct;

	NovaBiquadBase()
	{
		initFilter();

		if ((inRate(2) == calc_ScalarRate) && (inRate(3) == calc_ScalarRate))
			set_vector_calc_function<NovaBiquadBase, &NovaBiquadBase::next_i,
					&NovaBiquadBase::next_1>();
		else
			set_vector_calc_function<NovaBiquadBase, &NovaBiquadBase::next_k,
					&NovaBiquadBase::next_1>();
	}

	void initFilter()
	{
		float newFreq = in0(2);
		float newQ    = in0(3);

		_freqArg = newFreq;
		_qArg    = newQ;
		auto params = FilterDesigner::template designFilter<BiquadParameterStruct>( newFreq * (float)sampleDur(), newQ );
		storeFilterParameters( params );
	}

	void next_1(int)
	{
		auto inFn  = nova::Interleaver2<v2d>(this);
		auto outFn = nova::Deinterleaver2<v2d>(this);

		_filter.run(inFn, outFn, 1);
	}

	void next_i(int)
	{
		auto inFn  = nova::Interleaver2<v2d>(this);
		auto outFn = nova::Deinterleaver2<v2d>(this);

		_filter.run_unrolled(inFn, outFn, mRate->mFilterLoops, mRate->mFilterRemain);
	}

	void next_k(int inNumSamples)
	{
		float newFreq = in0(2);
		float newQ    = in0(3);

		if ( (newFreq == _freqArg) && (newQ == _qArg) ) {
			next_i(inNumSamples);
			return;
		}

		_freqArg = newFreq;
		_qArg    = newQ;
		BiquadParameterStruct newParameters = FilterDesigner::template designFilter<BiquadParameterStruct>( newFreq * (float)sampleDur(), newQ );

		auto a0Slope = calcSlope( newParameters.a0, _filter._a0 );
		auto a1Slope = calcSlope( newParameters.a1, _filter._a1 );
		auto a2Slope = calcSlope( newParameters.a2, _filter._a2 );
		auto b1Slope = calcSlope( newParameters.b1, _filter._b1 );
		auto b2Slope = calcSlope( newParameters.b2, _filter._b2 );

		auto inFn  = nova::Interleaver2<v2d>(this);
		auto outFn = nova::Deinterleaver2<v2d>(this);

		_filter.run_unrolled(inFn, outFn, mRate->mFilterLoops, mRate->mFilterRemain,
					a0Slope, a1Slope, a2Slope, b1Slope, b2Slope );
		storeFilterParameters( newParameters );
	}

	void storeFilterParameters(BiquadParameterStruct const & parameters)
	{
		_filter._a0 = parameters.a0;
		_filter._a1 = parameters.a1;
		_filter._a2 = parameters.a2;
		_filter._b1 = parameters.b1;
		_filter._b2 = parameters.b2;
	}

	float _freqArg, _qArg;
	Filter _filter;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct DesignLPF
{
	template < typename BiquadParameterStruct, typename ArgumentType >
	static BiquadParameterStruct designFilter ( ArgumentType cutoff, ArgumentType q )
	{
		using namespace boost::simd;

		auto K  = std::tan( boost::simd::Pi<ArgumentType>() * cutoff );
		auto K2 = K*K;

		BiquadParameterStruct ret;
		auto denom = (K2*q + K + q);
		ret.a0 =     (K2*q) / denom;
		ret.a1 = 2 * ret.a0;
		ret.a2 =     ret.a0;

		ret.b1 = 2 * q * (K2-1) / denom;
		ret.b2 = (K2*q - K + q) / denom;
		return ret;
	}
};

struct DesignHPF
{
	template < typename BiquadParameterStruct, typename ArgumentType >
	static BiquadParameterStruct designFilter ( ArgumentType cutoff, ArgumentType q )
	{
		auto K  = std::tan( boost::simd::Pi<ArgumentType>() * cutoff );
		auto K2 = K*K;

		BiquadParameterStruct ret;

		auto denom = (K2*q + K + q);
		ret.a0 = q      / denom;
		ret.a1 = -2 * q / denom;
		ret.a2 = ret.a0;

		ret.b1 = 2 * q * (K2-1) / denom;
		ret.b2 = (K2*q - K + q) / denom;
		return ret;
	}
};

struct DesignBPF
{
	template < typename BiquadParameterStruct, typename ArgumentType >
	static BiquadParameterStruct designFilter ( ArgumentType cutoff, ArgumentType q )
	{
		auto K  = std::tan( boost::simd::Pi<ArgumentType>() * cutoff );
		auto K2 = K*K;

		BiquadParameterStruct ret;

		auto denom = (K2*q + K + q);
		ret.a0 = K / denom;
		ret.a1 = 0;
		ret.a2 = -ret.a0;

		ret.b1 = 2 * q * (K2-1) / denom;
		ret.b2 = (K2*q - K + q) / denom;
		return ret;
	}
};

struct DesignBRF
{
	template < typename BiquadParameterStruct, typename ArgumentType >
	static BiquadParameterStruct designFilter ( ArgumentType cutoff, ArgumentType q )
	{
		auto K  = std::tan( boost::simd::Pi<ArgumentType>() * cutoff );
		auto K2 = K*K;

		BiquadParameterStruct ret;

		auto denom = (K2*q + K + q);
		ret.a0 = q * (1 + K2)     / denom;
		ret.a1 = 2 * q * (K2 - 1) / denom;
		ret.a2 = ret.a0;

		ret.b1 = 2 * q * (K2-1) / denom;
		ret.b2 = (K2*q - K + q) / denom;
		return ret;
	}
};

struct DesignAPF
{
	template < typename BiquadParameterStruct, typename ArgumentType >
	static BiquadParameterStruct designFilter ( ArgumentType cutoff, ArgumentType q )
	{
		auto K  = std::tan( boost::simd::Pi<ArgumentType>() * cutoff );
		auto K2 = K*K;

		BiquadParameterStruct ret;

		auto denom = (K2*q + K + q);
		ret.a0 = (K2*q - K + q)   / denom;
		ret.a1 = 2 * q * (K2 - 1) / denom;
		ret.a2 = 1;

		ret.b1 = 2 * q * (K2-1) / denom;
		ret.b2 = (K2*q - K + q) / denom;
		return ret;
	}
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct NovaLowPass2:
	NovaBiquadBase<DesignLPF>
{
	NovaLowPass2() {}
};

struct NovaHighPass2:
	NovaBiquadBase<DesignHPF>
{
	NovaHighPass2() {}
};

struct NovaBandPass2:
	NovaBiquadBase<DesignBPF>
{
	NovaBandPass2() {}
};

struct NovaBandReject2:
	NovaBiquadBase<DesignBRF>
{
	NovaBandReject2() {}
};

struct NovaAllPass2:
	NovaBiquadBase<DesignAPF>
{
	NovaAllPass2() {}
};


DEFINE_XTORS(NovaLeakDC2)
DEFINE_XTORS(NovaLowPass2)
DEFINE_XTORS(NovaHighPass2)
DEFINE_XTORS(NovaBandPass2)
DEFINE_XTORS(NovaBandReject2)
DEFINE_XTORS(NovaAllPass2)

}

PluginLoad(NovaFilters)
{
	ft = inTable;
	DefineSimpleUnit(NovaLeakDC2);
	DefineSimpleUnit(NovaLowPass2);
	DefineSimpleUnit(NovaHighPass2);
	DefineSimpleUnit(NovaBandPass2);
	DefineSimpleUnit(NovaBandReject2);
	DefineSimpleUnit(NovaAllPass2);
}
