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

#include "boost/simd/include/pack.hpp"
#include "boost/simd/include/native.hpp"
#include "boost/simd/operator/include/functions/minus.hpp"
#include "boost/simd/operator/include/functions/multiplies.hpp"
#include "boost/simd/operator/include/functions/plus.hpp"

#include "boost/simd/include/constants/one.hpp"
#include "boost/simd/include/constants/two.hpp"
#include "boost/simd/include/constants/quarter.hpp"

#include "producer_consumer_functors.hpp"

#include "boost/simd/arithmetic/include/functions/abs.hpp"
#include "boost/simd/arithmetic/include/functions/fast_rec.hpp"
#include "boost/simd/ieee/include/functions/copysign.hpp"

#include <nt2/include/functions/pow.hpp>
#include <nt2/include/functions/pow_abs.hpp>

#include "dsp/utils.hpp"

namespace {

InterfaceTable *ft;

template <typename Parent>
struct SaturationBase:
	public SCUnit
{
	static float verifyLevel ( float arg ) { return arg; }

	SaturationBase()
	{
		switch (inRate(1)) {
		case calc_FullRate:
			if (boost::simd::is_aligned( bufferSize(), 8 ) )
				set_calc_function<SaturationBase, &SaturationBase::run_a< boost::simd::pack<float, 8> > >();
			else
				set_calc_function<SaturationBase, &SaturationBase::run_a<float>>();
			break;

		case calc_BufRate:
			_level = Parent::verifyLevel( in0(1) );
			if (boost::simd::is_aligned( bufferSize(), 8 ) )
				set_calc_function<SaturationBase, &SaturationBase::run_k< boost::simd::pack<float, 8> > >();
			else
				set_calc_function<SaturationBase, &SaturationBase::run_k<float>>();
			break;

		case calc_ScalarRate:
		default:
			_level = Parent::verifyLevel( in0(1) );
			if (boost::simd::is_aligned( bufferSize(), 4 ) )
				set_calc_function<SaturationBase, &SaturationBase::run_i< boost::simd::pack<float, 8> > >();
			else
				set_calc_function<SaturationBase, &SaturationBase::run_i<float>>();
		}

	}

	template <typename Functor>
	inline void loop (int loops, Functor const & f)
	{
		for (int i = 0; i != loops; ++i)
			f();
	}

	template <typename SampleType,
			  typename Input0,
			  typename Input1,
			  typename Output0>
	inline void perform(int inNumSamples, Input0 & input0, Input1 & input1, Output0 & output)
	{
		using namespace boost::simd;

		const size_t unroll = meta::cardinal_of<SampleType>::value;
		loop( inNumSamples / unroll, [&] {
			auto in0 = input0();
			auto in1 = input1();
			auto result = Parent::doDistort( in0, in1 );
			output(result);
		});
	}


	template <typename SampleType>
	void run_a (int inNumSamples)
	{
		auto input0 = nova::Packer<SampleType, 0>(this);
		auto input1 = nova::Packer<SampleType, 1>(this);
		auto output = nova::Unpacker<SampleType, 0>(this);

		perform<SampleType>( inNumSamples, input0, input1, output );
	}

	template <typename SampleType>
	void run_k (int inNumSamples)
	{
		float newLevel = Parent::verifyLevel( in0(1) );
		if (newLevel != _level) {
			float slope = calcSlope(newLevel, _level);

			auto input0 = nova::Packer<SampleType, 0>(this);
			auto input1 = nova::makeRamp<SampleType>(_level, slope);
			auto output = nova::Unpacker<SampleType, 0>(this);
			_level = newLevel;

			perform<SampleType>( inNumSamples, input0, input1, output );
		} else {
			auto input0 = nova::Packer<SampleType, 0>(this);
			auto input1 = nova::Scalar<SampleType>(_level);
			auto output = nova::Unpacker<SampleType, 0>(this);

			perform<SampleType>( inNumSamples, input0, input1, output );
		}
	}

	template <typename SampleType>
	void run_i (int inNumSamples)
	{
		auto input0 = nova::Packer<SampleType, 0>(this);
		auto input1 = nova::Scalar<SampleType>(_level);
		auto output = nova::Unpacker<SampleType, 0>(this);

		perform<SampleType>( inNumSamples, input0, input1, output );
	}

	float _level;
};

class HyperbolSaturation : public SaturationBase<HyperbolSaturation>
{
public:
	HyperbolSaturation()
	{}

	template <typename SampleType>
	static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, SampleType level)
	{
		using namespace boost::simd;

		auto saturated = level - (level * level * fast_rec( level + abs(sig) ) );

		return copysign(saturated, sig);
	}

	static float verifyLevel(float arg)
	{
		return std::max( arg, 1e-10f );
	}
};

class ParabolSaturation : public SaturationBase<ParabolSaturation>
{
public:
	ParabolSaturation()
	{}

	template <typename SampleType>
	static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, SampleType level)
	{
		using namespace boost::simd;

		auto limit = Two<SampleType>() * level;
		auto clippedSignal = nova::clip2<SampleType>(sig, limit);

		auto factor = One<SampleType>() - ( abs(clippedSignal) * Quarter<SampleType>() * fast_rec(level));

		return sig * factor;
	}

	static float verifyLevel(float arg)
	{
		return std::max( arg, 1e-10f );
	}
};


class PowSaturation : public SaturationBase<PowSaturation>
{
public:
	PowSaturation()
	{}

	template <typename SampleType>
	static BOOST_FORCEINLINE FLATTEN SampleType doDistort(SampleType sig, SampleType level)
	{
		using namespace boost::simd;

//		auto pow = nt2::pow(abs(sig), level);
//		auto ret = boost::simd::copysign(pow, sig);
//		return pow;
		return boost::simd::copysign(nt2::pow_abs(sig, level), sig);
	}

	static float verifyLevel(float arg)
	{
		return std::max( arg, 1e-10f );
	}
};



DEFINE_XTORS(HyperbolSaturation)
DEFINE_XTORS(ParabolSaturation)
DEFINE_XTORS(PowSaturation)

}

PluginLoad(NovaSaturators)
{
	ft = inTable;
	DefineDtorUnit(HyperbolSaturation);
	DefineDtorUnit(ParabolSaturation);
	DefineDtorUnit(PowSaturation);
}

