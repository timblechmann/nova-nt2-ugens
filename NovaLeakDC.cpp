/*
    Copyright (C) Tim Blechmann

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


#include "NovaUGensCommon.hpp"

#include "dsp/leak_dc.hpp"
#include "producer_consumer_functors.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/simd/include/constants/pi.hpp>

#include "dsp/utils.hpp"

#include <cmath>

using nova::NovaUnit;

InterfaceTable *ft;

namespace nova {

namespace constants = boost::math::constants;

template <int NumberOfChannels>
struct NovaLeakDC:
    public NovaUnit,
    public nova::multichannel::SignalInput< NovaLeakDC<NumberOfChannels>, 0, NumberOfChannels >,

    public nova::multichannel::OutputSink< NovaLeakDC<NumberOfChannels>, 0, NumberOfChannels >

{
    typedef typename nova::as_pack< double, NumberOfChannels >::type vDouble;
    typedef nova::LeakDC< vDouble, double > Filter;

    static const size_t IndexOfCoefficient = NumberOfChannels;

    typedef nova::multichannel::SignalInput< NovaLeakDC<NumberOfChannels>, 0, NumberOfChannels > InputSignal;
    typedef nova::multichannel::OutputSink<  NovaLeakDC<NumberOfChannels>, 0, NumberOfChannels > OutputSink;

    NovaLeakDC()
    {
        initFilter(in0(IndexOfCoefficient));

        auto inFn = InputSignal::template makeInputSignal<vDouble>();
        _filter._x_1 = inFn();

        switch (inRate(IndexOfCoefficient))
        {
        case calc_ScalarRate:
            set_calc_function<NovaLeakDC, &NovaLeakDC::next_i>();
            break;

        case calc_BufRate:
        default:
            _freq = std::numeric_limits<float>::quiet_NaN();
            set_calc_function<NovaLeakDC, &NovaLeakDC::next_k>();
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
        auto inFn  = InputSignal::template makeInputSignal<vDouble>();
        auto outFn = OutputSink:: template makeSink<vDouble>();

        _filter.run(inFn, outFn, inNumSamples);
    }

    void next_k(int inNumSamples)
    {
        float newFreq = in0(IndexOfCoefficient);
        if (newFreq != _freq) {
            float oldA = _filter._a;
            float newA = designFilter( newFreq );

            float slopeA = calcSlope(newA, oldA);

            _freq = newFreq;

            auto inFn  = InputSignal::template makeInputSignal<vDouble>();
            auto outFn = OutputSink:: template makeSink<vDouble>();

            _filter.run(inFn, outFn, inNumSamples, slopeA);
            _filter.set_a( newA );
        } else {
            next_i(inNumSamples);
        }
    }

    Filter _filter;
    float _freq;
};

}

typedef nova::NovaLeakDC<1> NovaLeakDC;
typedef nova::NovaLeakDC<2> NovaLeakDC2;
typedef nova::NovaLeakDC<4> NovaLeakDC4;
typedef nova::NovaLeakDC<8> NovaLeakDC8;




DEFINE_XTORS(NovaLeakDC )
DEFINE_XTORS(NovaLeakDC2)
DEFINE_XTORS(NovaLeakDC4)
DEFINE_XTORS(NovaLeakDC8)


PluginLoad(NovaLeakDC)
{
    ft = inTable;
	DefineSimpleUnit(NovaLeakDC );
    DefineSimpleUnit(NovaLeakDC2);
    DefineSimpleUnit(NovaLeakDC4);
    DefineSimpleUnit(NovaLeakDC8);
}
