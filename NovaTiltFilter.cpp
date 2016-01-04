/*
    Copyright (C) 2016 Tim Blechmann

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

#include <dsp/tilt_filter.hpp>

#include "NovaUnitFacade.hpp"

#include <producer_consumer_functors.hpp>

#include <boost/mpl/if.hpp>
#include <boost/mpl/map.hpp>

namespace nova
{

InterfaceTable *ft;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <size_t Size, bool ScalarArguments = true>
struct NovaTiltFilter:
    public NovaUnit,
    public nova::multichannel::SignalInput< NovaTiltFilter<Size, ScalarArguments>, 0, Size >,
    public nova::multichannel::OutputSink<  NovaTiltFilter<Size, ScalarArguments>, 0, Size >,

    // Controls:
    public nova::multichannel::ControlInput<  NovaTiltFilter<Size, ScalarArguments>, Size,                                (ScalarArguments ? 1 : Size) >,
    public nova::multichannel::ControlInput<  NovaTiltFilter<Size, ScalarArguments>, Size + (ScalarArguments ? 1 : Size), (ScalarArguments ? 1 : Size) >
{
    static const size_t ParameterSize = ScalarArguments ? 1 : Size;
    static const int FreqInputIndex = Size;
    static const int QInputIndex    = Size + ParameterSize;

    using SignalInput = nova::multichannel::SignalInput< NovaTiltFilter, 0, Size >;
    using OutputSink  = nova::multichannel::OutputSink<  NovaTiltFilter, 0, Size >;

    using FreqInput = nova::multichannel::ControlInput<  NovaTiltFilter, Size,                                ScalarArguments ? 1 : Size >;
    using GainInput = nova::multichannel::ControlInput<  NovaTiltFilter, Size + (ScalarArguments ? 1 : Size), ScalarArguments ? 1 : Size >;

    typedef double InternalType;
    using SampleType    = typename nova::as_pack<InternalType, Size>::type;
    using ParameterType = typename nova::as_pack<float, ParameterSize>::type;
    using InternalParameterType = typename nova::as_pack<InternalType, ParameterSize>::type;


    typedef TiltFilter< Size, ParameterSize, InternalType > Filter;
    Filter _filter;

    NovaTiltFilter()
    {
        initFilter();
        setCalcFunction< NovaTiltFilter >( FreqInput::index, GainInput::index );
    }

    void initFilter()
    {
        auto newFreq = FreqInput::readInput();
        auto res    = GainInput::readInput();

        auto state = Filter::computeState( castType<InternalParameterType>(newFreq), castType<InternalParameterType>(res), makeDSPContext() );
        _filter.setState( state );
    }

    template <typename ControlSignature>
    void run(int inNumSamples)
    {
        next( inNumSamples, ControlSignature() );
    }

    void next(int inNumSamples, nova::control_signature_kk)
    {
        if ( !FreqInput::changed() && !GainInput::changed() ) {
            next( inNumSamples, nova::control_signature_ii() );
            return;
        }

        auto oldState = _filter.currentState();
        auto newState = Filter::computeState( castType<InternalParameterType>( FreqInput::readAndUpdateInput() ),
                                              castType<InternalParameterType>( GainInput::readAndUpdateInput()  ),
                                              makeDSPContext() );
        _filter.setState( newState );

        auto stateSlope = calcSlope( newState, oldState );
        auto stateRamp = nova::makeScalarRamp<decltype(stateSlope)>( oldState, stateSlope );

        auto inFn  = SignalInput::template makeInputSignal<SampleType>();
        auto outFn = OutputSink ::template makeSink<SampleType>();
        _filter.run(inFn, outFn, inNumSamples, stateRamp );
    }

    void next(int inNumSamples, nova::control_signature_ii)
    {
        auto inFn  = SignalInput::template makeInputSignal<SampleType>();
        auto outFn = OutputSink ::template makeSink<SampleType>();
        _filter.run( inFn, outFn, inNumSamples, _filter.getState() );
    }

    void next(int inNumSamples, nova::control_signature_aa)
    {
        auto context = makeDSPContext();

        auto makeState = [this, index = 0, context] () mutable {
            ParameterType freq = FreqInput::template readInputs<ParameterType>( index );
            ParameterType res  = GainInput:: template readInputs<ParameterType>( index );
            index += 1;

            return Filter::computeState( castType<InternalParameterType>(freq),
                                         castType<InternalParameterType>(res),
                                         context );
        };

        auto inFn  = SignalInput::template makeInputSignal<SampleType>();
        auto outFn = OutputSink ::template makeSink<SampleType>();
        _filter.run( inFn, outFn, inNumSamples, makeState );
    }
};

}

using namespace nova;

PluginLoad(NovaTiltFilter)
{
    ft = inTable;

    nova::registerUnit< NovaTiltFilter<1, true>  >( ft, "NovaTiltFilter"    );
    nova::registerUnit< NovaTiltFilter<2, true>  >( ft, "NovaTiltFilter2"   );
    nova::registerUnit< NovaTiltFilter<4, true>  >( ft, "NovaTiltFilter4"   );

    nova::registerUnit< NovaTiltFilter<2, false> >( ft, "NovaTiltFilter2_2" );
    nova::registerUnit< NovaTiltFilter<4, false> >( ft, "NovaTiltFilter4_4" );
}
