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

#include "dsp/utils.hpp"

#include <cmath>

using nova::NovaUnit;

InterfaceTable *ft;

namespace nova {

template <int NumberOfChannels>
struct NovaLeakDC:
    public NovaUnit,
    public nova::multichannel::SignalInput< NovaLeakDC<NumberOfChannels>, 0, NumberOfChannels >,

    public nova::ControlInput< NovaLeakDC<NumberOfChannels>, NumberOfChannels >,

    public nova::multichannel::OutputSink< NovaLeakDC<NumberOfChannels>, 0, NumberOfChannels >

{
    typedef typename nova::as_pack< double, NumberOfChannels >::type vDouble;
    typedef nova::LeakDC< vDouble, double > Filter;

    static const size_t IndexOfCoefficient = NumberOfChannels;

    typedef nova::multichannel::SignalInput< NovaLeakDC<NumberOfChannels>, 0, NumberOfChannels > InputSignal;
    typedef nova::ControlInput< NovaLeakDC<NumberOfChannels>, NumberOfChannels >                  FreqInput;

    typedef nova::multichannel::OutputSink<  NovaLeakDC<NumberOfChannels>, 0, NumberOfChannels > OutputSink;

    NovaLeakDC():
        _filter( computeState(FreqInput::readInput()),
                 InputSignal::template readInputs<vDouble>() )
    {
        setCalcFunction< NovaLeakDC >( IndexOfCoefficient );
    }

    auto computeState ( float cutoffFreq )
    {
        return Filter::computeState( cutoffFreq, makeDSPContext() );
    }

    template <typename ControlSignature>
    void run(int inNumSamples)
    {
        run( inNumSamples, ControlSignature() );
    }

    void run(int inNumSamples, nova::control_signature_a)
    {
        auto inFn       = InputSignal::template makeInputSignal<vDouble>();
        auto outFn      = OutputSink:: template makeSink<vDouble>();
        auto freqInFunc = FreqInput::template makeAudioInputSignal<float>();

        auto parameterFunctor = [=] () mutable {
            auto freq = freqInFunc();
            return computeState( freq );
        };

        _filter.run(inFn, outFn, inNumSamples, parameterFunctor );
    }

    void run(int inNumSamples, nova::control_signature_i)
    {
        auto inFn  = InputSignal::template makeInputSignal<vDouble>();
        auto outFn = OutputSink:: template makeSink<vDouble>();

        _filter.run(inFn, outFn, inNumSamples );
    }

    void run(int inNumSamples, nova::control_signature_k)
    {
        if ( FreqInput::changed() ) {
            auto currentState = _filter.currentState();
            auto newState     = computeState( FreqInput::readInput() );
            _filter.setState( newState );

            auto state = parameter::makeRamp( currentState, newState, (float)this->mRate->mSlopeFactor );

            auto inFn  = InputSignal::template makeInputSignal<vDouble>();
            auto outFn = OutputSink:: template makeSink<vDouble>();
            _filter.run(inFn, outFn, inNumSamples, state );
        } else {
            run(inNumSamples, nova::control_signature_i());
        }
    }

    Filter _filter;
};

}

PluginLoad(NovaLeakDC)
{
    using namespace nova;

    ft = inTable;

    nova::registerUnit<NovaLeakDC<1>>( ft, "NovaLeakDC"  );
    nova::registerUnit<NovaLeakDC<2>>( ft, "NovaLeakDC2" );
    nova::registerUnit<NovaLeakDC<4>>( ft, "NovaLeakDC4" );
    nova::registerUnit<NovaLeakDC<8>>( ft, "NovaLeakDC8" );
}
