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

#include "NovaUnitFacade.hpp"

#include "dsp/feedback_am.hpp"

static InterfaceTable *ft;

namespace nova {
namespace      {


template <int NumberOfChannels, bool ScalarArguments = false>
struct NovaFeedbackAM:
	public NovaUnitUnary<NovaFeedbackAM<NumberOfChannels, ScalarArguments>, NumberOfChannels, double, ScalarArguments>
{
	typedef NovaUnitUnary<NovaFeedbackAM<NumberOfChannels, ScalarArguments>, NumberOfChannels, double, ScalarArguments> Base;

	typedef nova::FeedbackAM<typename Base::SampleType, typename Base::ParameterDSPType> Filter;

    NovaFeedbackAM()
    {}

	struct DSPEngine:
		Filter
	{
		template <typename T>
		void setParameter( T const & t) { Filter::_fb = t; }

		auto getParameter() { return Filter::_fb; }
	};

    template <typename AType>
	static auto checkParameter(AType const & a)
	{
		return Filter::checkParameter(a);
	}

	static DSPEngine & getDSPEngine( NovaFeedbackAM * self )
	{
		return self->_filter;
	}

	DSPEngine _filter;
};

}
}

typedef nova::NovaFeedbackAM<1, true> NovaFeedbackAM;
typedef nova::NovaFeedbackAM<2, true> NovaFeedbackAM2;
typedef nova::NovaFeedbackAM<4, true> NovaFeedbackAM4;
typedef nova::NovaFeedbackAM<8, true> NovaFeedbackAM8;

typedef nova::NovaFeedbackAM<2> NovaFeedbackAM2_2;
typedef nova::NovaFeedbackAM<4> NovaFeedbackAM4_4;
typedef nova::NovaFeedbackAM<8> NovaFeedbackAM8_8;

DEFINE_XTORS(NovaFeedbackAM)
DEFINE_XTORS(NovaFeedbackAM2)
DEFINE_XTORS(NovaFeedbackAM4)
DEFINE_XTORS(NovaFeedbackAM8)

DEFINE_XTORS(NovaFeedbackAM2_2)
DEFINE_XTORS(NovaFeedbackAM4_4)
DEFINE_XTORS(NovaFeedbackAM8_8)

PluginLoad(FBAM)
{
	ft = inTable;
	DefineSimpleUnit(NovaFeedbackAM);
	DefineSimpleUnit(NovaFeedbackAM2);
	DefineSimpleUnit(NovaFeedbackAM4);
	DefineSimpleUnit(NovaFeedbackAM8);

	DefineSimpleUnit(NovaFeedbackAM2_2);
	DefineSimpleUnit(NovaFeedbackAM4_4);
	DefineSimpleUnit(NovaFeedbackAM8_8);
}
