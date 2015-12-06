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

    typedef nova::FeedbackAM<typename Base::SampleType, typename Base::ParameterDSPType> DSPEngine;

    NovaFeedbackAM()
    {}

    static DSPEngine & getDSPEngine( NovaFeedbackAM * self )
    {
        return self->_filter;
    }

    DSPEngine _filter;
};

}
}

PluginLoad(FBAM)
{
    using namespace nova;

    ft = inTable;

    registerUnit< nova::NovaFeedbackAM<1, true> >( ft, "NovaFeedbackAM"  );
    registerUnit< nova::NovaFeedbackAM<2, true> >( ft, "NovaFeedbackAM2"  );
    registerUnit< nova::NovaFeedbackAM<4, true> >( ft, "NovaFeedbackAM4"  );
    registerUnit< nova::NovaFeedbackAM<8, true> >( ft, "NovaFeedbackAM8"  );

    registerUnit< nova::NovaFeedbackAM<2, false> >( ft, "NovaFeedbackAM2_2" );
    registerUnit< nova::NovaFeedbackAM<4, false> >( ft, "NovaFeedbackAM4_4" );
    registerUnit< nova::NovaFeedbackAM<8, false> >( ft, "NovaFeedbackAM8_8" );
}
