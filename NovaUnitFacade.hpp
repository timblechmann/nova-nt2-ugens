/*
 *
 *    Copyright (C) Tim Blechmann
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

#ifndef NOVAUNITFACADE_HPP
#define NOVAUNITFACADE_HPP

#include "SC_PlugIn.hpp"

#include "NovaUGensCommon.hpp"
#include "dsp/utils.hpp"

#include "producer_consumer_functors.hpp"

#include "boost/simd/include/functions/compare_equal.hpp"

namespace nova {


// TODO:
// * nova::toDouble -> toParameterDSPType
//
template < class DerivedClass,
		   int NumberOfChannels,
		   typename ScalarSampleType = float,
		   bool ScalarArguments = false>
struct NovaUnitUnary:
	public NovaUnit
{
	static const size_t ParameterSize = ScalarArguments ? 1
														: NumberOfChannels;

	typedef float            FrontendSampleType;
	typedef ScalarSampleType BackendSampleType;

	typedef typename nova::as_pack<double, NumberOfChannels>::type SampleType;
	typedef typename nova::as_pack<float, ParameterSize>::type  ParameterType;
	typedef typename nova::as_pack<double, ParameterSize>::type ParameterDSPType;

	typedef InputInterleaver<ParameterSize> ParameterReader;

	static const size_t IndexOfParameter = NumberOfChannels;

	NovaUnitUnary()
	{
		initDSP();

		if (isScalarRate(IndexOfParameter, IndexOfParameter + ParameterSize)) {
			set_calc_function<NovaUnitUnary, &NovaUnitUnary::next_i>();
			return;
		}

		if (isBufRate(IndexOfParameter, IndexOfParameter + ParameterSize)) {
			set_calc_function<NovaUnitUnary, &NovaUnitUnary::next_k>();
			return;
		}

		set_calc_function<NovaUnitUnary, &NovaUnitUnary::next_a>();
	}

	void initDSP()
	{
		ParameterType newParameter = ParameterReader::read(this, IndexOfParameter);

		ParameterType checkedParameter = DerivedClass::checkParameter( newParameter );

		auto & dspEngine = getEngine();
		dspEngine.setParameter( nova::toDouble<ParameterDSPType>( checkedParameter ) );

		_cachedParameter = newParameter;
	}

	// hooks

	template <typename AType>
	static auto checkParameter(AType const & a) { return a; }

	// static void getDSPEngine(DerivedClass * self);


	void next_i(int inNumSamples)
	{
		auto inFn  = nova::Interleaver<SampleType>(this);
		auto outFn = nova::Deinterleaver<SampleType>(this);

		auto & dspEngine = getEngine();
		dspEngine.run(inFn, outFn, inNumSamples);
	}

	void next_k(int inNumSamples)
	{
		ParameterType newParameter = ParameterReader::read(this, IndexOfParameter);
		using namespace boost::simd;

		logical<float> coeffConstant = compare_equal(newParameter, _cachedParameter);

		if (coeffConstant) {
			next_i(inNumSamples);
		} else {
			ParameterType checkedParameter = DerivedClass::checkParameter( newParameter );

			auto & dspEngine = getEngine();

			auto slope = calcSlope( nova::toDouble<ParameterDSPType>( checkedParameter ), dspEngine.getParameter() );
			_cachedParameter = newParameter;

			auto inFn  = nova::Interleaver<SampleType>(this);
			auto outFn = nova::Deinterleaver<SampleType>(this);

			dspEngine.run(inFn, outFn, inNumSamples, slope);
			dspEngine.setParameter( nova::toDouble<ParameterDSPType>( checkedParameter ) );
		}
	}

	void next_a(int inNumSamples)
	{
		auto inFn      = nova::Interleaver<SampleType>(this);
		auto inParamFn = nova::Interleaver<ParameterDSPType, IndexOfParameter>(this);
		auto outFn = nova::Deinterleaver<SampleType>(this);

		auto & dspEngine = getEngine();

		dspEngine.run_ar(inFn, inParamFn, outFn, inNumSamples);
	}

	// helpers
	auto & getEngine()
	{
		return DerivedClass::getDSPEngine( static_cast<DerivedClass*>( this ) );
	}

	ParameterType _cachedParameter;
};

}

#endif // NOVAUNITFACADE_HPP
