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

#include "producer_consumer_functors.hpp"

#include "dsp/feedback_am.hpp"
#include "dsp/utils.hpp"

static InterfaceTable *ft;

namespace nova {

template <int NumberOfChannels>
struct NovaFeedbackAM:
	public SCUnit
{
	typedef typename nova::as_pack<double, 1>::type vDouble;
	typedef nova::FeedbackAM<vDouble, double> Filter;

	static const size_t IndexOfCoefficient = NumberOfChannels;

	NovaFeedbackAM()
	{
		initFilter(in0(IndexOfCoefficient));

		switch (inRate(IndexOfCoefficient))
		{
		case calc_ScalarRate:
			set_calc_function<NovaFeedbackAM, &NovaFeedbackAM::next_i>();
			break;

		case calc_BufRate:
		default:
			_coeff = std::numeric_limits<float>::quiet_NaN();
			set_calc_function<NovaFeedbackAM, &NovaFeedbackAM::next_k>();
		}
	}

	void initFilter(float leakFactor)
	{
		_filter.set_fb( leakFactor );
		_coeff = leakFactor;
	}

	void next_i(int inNumSamples)
	{
		auto inFn  = nova::Interleaver<vDouble>(this);
		auto outFn = nova::Deinterleaver<vDouble>(this);

		_filter.run(inFn, outFn, inNumSamples);
	}

	void next_k(int inNumSamples)
	{
		float newCoeff = in0(IndexOfCoefficient);
		if (newCoeff != _coeff) {
			auto slopeA = calcSlope(newCoeff, _coeff);
			_coeff = newCoeff;

			auto inFn  = nova::Interleaver<vDouble>(this);
			auto outFn = nova::Deinterleaver<vDouble>(this);

			_filter.run(inFn, outFn, inNumSamples, slopeA);
			_filter.set_fb( newCoeff );
		} else {
			next_i(inNumSamples);
		}
	}

	Filter _filter;
	float _coeff;
};

}

typedef nova::NovaFeedbackAM<1> NovaFeedbackAM;
typedef nova::NovaFeedbackAM<2> NovaFeedbackAM2;
typedef nova::NovaFeedbackAM<4> NovaFeedbackAM4;
typedef nova::NovaFeedbackAM<8> NovaFeedbackAM8;

DEFINE_XTORS(NovaFeedbackAM)
DEFINE_XTORS(NovaFeedbackAM2)
DEFINE_XTORS(NovaFeedbackAM4)
DEFINE_XTORS(NovaFeedbackAM8)

PluginLoad(FBAM)
{
	ft = inTable;
	DefineSimpleUnit(NovaFeedbackAM);
	DefineSimpleUnit(NovaFeedbackAM2);
	DefineSimpleUnit(NovaFeedbackAM4);
	DefineSimpleUnit(NovaFeedbackAM8);
}
