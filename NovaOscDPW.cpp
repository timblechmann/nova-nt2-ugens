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

#include <dsp/osc_dpw.hpp>


#include "NovaUGensCommon.hpp"


InterfaceTable *ft;

namespace nova {

template <typename Oscillator>
struct NovaOscDPW:
    public NovaUnit,
    public nova::ControlInput< NovaOscDPW<Oscillator>, 0 >, // Frequency
    public nova::ControlInput< NovaOscDPW<Oscillator>, 1 >, // Phase
    public nova::OutputSink<   NovaOscDPW<Oscillator>, 0 >
{
    typedef nova::ControlInput< NovaOscDPW, 0 > FreqInput;
    typedef nova::ControlInput< NovaOscDPW, 1 > PhaseInput;
    typedef nova::OutputSink<   NovaOscDPW, 0 > OutputSink;

    NovaOscDPW():
        mOsc( computeState( FreqInput::readInput() ),
              PhaseInput::readInput() )
    {
        setCalcFunction< NovaOscDPW >( 0 );
    }

	auto computeState ( float cutoffFreq ) -> typename Oscillator::ParameterState
    {
        return Oscillator::computeState( cutoffFreq, makeDSPContext() );
    }

    template <typename ControlSignature>
    void run(int inNumSamples)
    {
        run( inNumSamples, ControlSignature() );
    }

    void run(int inNumSamples, nova::control_signature_a)
    {
        auto outFn      = OutputSink:: template makeSink<float>();
        auto freqInFunc = FreqInput::template makeAudioInputSignal<float>();

        auto parameterFunctor = [=] () mutable {
            auto freq = freqInFunc();
            return computeState( freq );
        };

        mOsc.run( outFn, inNumSamples, parameterFunctor );
    }

    void run(int inNumSamples, nova::control_signature_i)
    {
        auto outFn = OutputSink:: template makeSink<float>();

        mOsc.run( outFn, inNumSamples );
    }

    void run(int inNumSamples, nova::control_signature_k)
    {
        if ( FreqInput::changed() ) {
            auto currentState = mOsc.currentState();
            auto newState     = computeState( FreqInput::readInput() );
            mOsc.setState( newState );

            auto state = parameter::makeRamp( currentState, newState, (float)this->mRate->mSlopeFactor );

            auto outFn = OutputSink:: template makeSink<float>();
            mOsc.run( outFn, inNumSamples, state );
        } else {
            run(inNumSamples, nova::control_signature_i());
        }
    }

    Oscillator mOsc;
};


struct NovaPulseDPW:
    public NovaUnit,
    public nova::ControlInput< NovaPulseDPW, 0 >, // Frequency
    public nova::ControlInput< NovaPulseDPW, 1 >, // Width
    public nova::ControlInput< NovaPulseDPW, 2 >, // Phase
    public nova::OutputSink<   NovaPulseDPW, 0 >
{
    typedef nova::ControlInput< NovaPulseDPW, 0 > FreqInput;
    typedef nova::ControlInput< NovaPulseDPW, 1 > WidthInput;
    typedef nova::ControlInput< NovaPulseDPW, 2 > PhaseInput;
    typedef nova::OutputSink<   NovaPulseDPW, 0 > OutputSink;

    NovaPulseDPW():
        mOsc( computeState( FreqInput::readInput() ),
              WidthInput::readInput(),
              PhaseInput::readInput()
            )
    {
        setCalcFunction< NovaPulseDPW >( 0 );
    }

	auto computeState ( float cutoffFreq ) -> typename PulseDPW::ParameterState
    {
        return PulseDPW::computeState( cutoffFreq, makeDSPContext() );
    }

    template <typename ControlSignature>
    void run(int inNumSamples)
    {
        run( inNumSamples, ControlSignature() );
    }

    void run(int inNumSamples, nova::control_signature_a)
    {
        auto outFn      = OutputSink:: template makeSink<float>();
        auto freqInFunc = FreqInput::template makeAudioInputSignal<float>();

        auto parameterFunctor = [=] () mutable {
            auto freq = freqInFunc();
            return computeState( freq );
        };

        mOsc.run( outFn, inNumSamples, parameterFunctor );
    }

    void run(int inNumSamples, nova::control_signature_i)
    {
        auto outFn = OutputSink:: template makeSink<float>();

        mOsc.run( outFn, inNumSamples );
    }

    void run(int inNumSamples, nova::control_signature_k)
    {
        if ( FreqInput::changed() ) {
            auto currentState = mOsc.currentState();
            auto newState     = computeState( FreqInput::readInput() );
            mOsc.setState( newState );

            auto state = parameter::makeRamp( currentState, newState, (float)this->mRate->mSlopeFactor );

            auto outFn = OutputSink:: template makeSink<float>();
            mOsc.run( outFn, inNumSamples, state );
        } else {
            run(inNumSamples, nova::control_signature_i());
        }
    }

    PulseDPW mOsc;
};

}

PluginLoad(NovaLeakDC)
{
    using namespace nova;

    ft = inTable;

	nova::registerUnit<NovaOscDPW<SawDPW>>( ft, "NovaSawDPW"   );
	nova::registerUnit<NovaOscDPW<TriDPW>>( ft, "NovaTriDPW"   );
    nova::registerUnit<NovaPulseDPW      >( ft, "NovaPulseDPW" );
}
