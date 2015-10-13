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
    public NovaUnit,
    public nova::multichannel::SignalInput< NovaUnitUnary<DerivedClass, NumberOfChannels, ScalarSampleType, ScalarArguments>, 0, NumberOfChannels >,
    public nova::multichannel::OutputSink<  NovaUnitUnary<DerivedClass, NumberOfChannels, ScalarSampleType, ScalarArguments>, 0, NumberOfChannels >,

    public nova::multichannel::ControlInput< NovaUnitUnary<DerivedClass, NumberOfChannels, ScalarSampleType, ScalarArguments>, NumberOfChannels, ScalarArguments ? 1 : NumberOfChannels >
{
    static const size_t ParameterSize = ScalarArguments ? 1
                                                        : NumberOfChannels;

    typedef float            FrontendSampleType;
    typedef ScalarSampleType BackendSampleType;

    typedef typename nova::as_pack<double, NumberOfChannels>::type SampleType;
    typedef typename nova::as_pack<float, ParameterSize>::type  ParameterType;
    typedef typename nova::as_pack<double, ParameterSize>::type ParameterDSPType;

    // inputs:
    using SignalInput = typename nova::multichannel::SignalInput< NovaUnitUnary<DerivedClass, NumberOfChannels, ScalarSampleType, ScalarArguments>, 0, NumberOfChannels >;
    using OutputSink  = typename nova::multichannel::OutputSink<  NovaUnitUnary<DerivedClass, NumberOfChannels, ScalarSampleType, ScalarArguments>, 0, NumberOfChannels >;

    using ParameterInput = nova::multichannel::ControlInput< NovaUnitUnary<DerivedClass, NumberOfChannels, ScalarSampleType, ScalarArguments>, NumberOfChannels, ScalarArguments ? 1 : NumberOfChannels >;

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
        auto & dspEngine = getEngine();
        dspEngine.setState( computeState( ParameterInput::readInput() ) );
    }

    auto computeState ( ParameterType parameter )
    {
        return computeState( nova::toDouble<ParameterDSPType>( parameter ) );
    }

    auto computeState ( ParameterDSPType parameter )
    {
        auto & dspEngine = getEngine();
        return dspEngine.computeState( parameter, makeDSPContext() );
    }


    // process function
    void next_i(int inNumSamples)
    {
        auto inSig   = SignalInput::template makeInputSignal<SampleType>();
        auto outSink = OutputSink ::template makeSink<SampleType>();

        auto & dspEngine = getEngine();
        dspEngine.run(inSig, outSink, inNumSamples);
    }

    void next_k(int inNumSamples)
    {
        if( !ParameterInput::changed() ) {
            next_i(inNumSamples);
        } else {
            auto & dspEngine = getEngine();

            auto currentState = dspEngine.currentState();
            auto nextState    = computeState( ParameterInput::readInput() );
            auto state        = parameter::makeRamp( currentState, nextState, mRate->mSlopeFactor );
            dspEngine.setState( nextState );

            auto inSig   = SignalInput::template makeInputSignal<SampleType>();
            auto outSink = OutputSink ::template makeSink<SampleType>();

            dspEngine.run(inSig, outSink, inNumSamples, state);
        }
    }

    void next_a(int inNumSamples)
    {
        auto inSig     = SignalInput::template makeInputSignal<SampleType>();
        auto inParamFn = ParameterInput::template makeAudioInputSignal<ParameterType>();
        auto outSink   = OutputSink ::template makeSink<SampleType>();

        auto & dspEngine = getEngine();

        auto stateIn = [&]() {
            return computeState( inParamFn() );
        };

        dspEngine.run_ar(inSig, stateIn, outSink, inNumSamples);
    }

    // helpers
    auto & getEngine()
    {
        return DerivedClass::getDSPEngine( static_cast<DerivedClass*>( this ) );
    }
};

}

#endif // NOVAUNITFACADE_HPP
