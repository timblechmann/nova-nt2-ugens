/*
    Copyright (C) 2015 Tim Blechmann

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

#include <dsp/svf.hpp>

#include "NovaUnitFacade.hpp"

#include <producer_consumer_functors.hpp>

#include <boost/mpl/if.hpp>
#include <boost/mpl/map.hpp>

namespace nova
{

InterfaceTable *ft;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <SVFType FilterType, size_t Size, bool ScalarArguments = true>
struct NovaSVFBase:
    public NovaUnit,
    public nova::multichannel::SignalInput< NovaSVFBase<FilterType, Size, ScalarArguments>, 0, Size >,
    public nova::multichannel::OutputSink<  NovaSVFBase<FilterType, Size, ScalarArguments>, 0, Size >,

    // Controls:
    public nova::multichannel::ControlInput<  NovaSVFBase<FilterType, Size, ScalarArguments>, Size,                                (ScalarArguments ? 1 : Size) >,
    public nova::multichannel::ControlInput<  NovaSVFBase<FilterType, Size, ScalarArguments>, Size + (ScalarArguments ? 1 : Size), (ScalarArguments ? 1 : Size) >

{
    static const size_t ParameterSize = ScalarArguments ? 1 : Size;
    static const int FreqInputIndex = Size;
    static const int QInputIndex    = Size + ParameterSize;

    using SignalInput = nova::multichannel::SignalInput< NovaSVFBase, 0, Size >;
    using OutputSink  = nova::multichannel::OutputSink<  NovaSVFBase, 0, Size >;

    using FreqInput = nova::multichannel::ControlInput<  NovaSVFBase, Size,                                ScalarArguments ? 1 : Size >;
    using ResInput  = nova::multichannel::ControlInput<  NovaSVFBase, Size + (ScalarArguments ? 1 : Size), ScalarArguments ? 1 : Size >;

    using SampleType    = typename nova::as_pack<float, Size>::type;
    using ParameterType = typename nova::as_pack<float, ParameterSize>::type;

    typedef SVFFilter< Size, ParameterSize, FilterType > Filter;
    Filter _filter;

    NovaSVFBase()
    {
        initFilter();
        setCalcFunction< NovaSVFBase >( FreqInput::index, ResInput::index );
    }

    void initFilter()
    {
        auto newFreq = FreqInput::readInput();
        auto res    = ResInput::readInput();

        auto state = Filter::computeState( newFreq, res, makeDSPContext() );
        _filter.setState( state );
    }

    template <typename ControlSignature>
    void run(int inNumSamples)
    {
        next( inNumSamples, ControlSignature() );
    }

    void next(int inNumSamples, nova::control_signature_kk)
    {
        if ( !FreqInput::changed() && !ResInput::changed() ) {
            next( inNumSamples, nova::control_signature_ii() );
            return;
        }

        auto oldState = _filter.currentState();
        auto newState = Filter::computeState( FreqInput::readAndUpdateInput(), ResInput::readAndUpdateInput(), makeDSPContext() );
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
        _filter.run( inFn, outFn, inNumSamples );
    }

    void next(int inNumSamples, nova::control_signature_aa)
    {
        auto context = makeDSPContext();

        auto makeState = [this, index = 0, context] () mutable {
            ParameterType freq = FreqInput::template readInputs<ParameterType>( index );
            ParameterType res  = ResInput:: template readInputs<ParameterType>( index );
            index += 1;

            return Filter::computeState( freq, res, context );
        };

        auto inFn  = SignalInput::template makeInputSignal<SampleType>();
        auto outFn = OutputSink ::template makeSink<SampleType>();
        _filter.run( inFn, outFn, inNumSamples, makeState );
    }
};

// TODO: ternary controls
template <SVFType FilterType, size_t Size, bool ScalarArguments = true>
struct NovaSVFShelfBase:
    public NovaUnit,
    public nova::multichannel::SignalInput< NovaSVFShelfBase<FilterType, Size, ScalarArguments>, 0, Size >,
    public nova::multichannel::OutputSink<  NovaSVFShelfBase<FilterType, Size, ScalarArguments>, 0, Size >,

    // Controls:
    public nova::multichannel::ControlInput<  NovaSVFShelfBase<FilterType, Size, ScalarArguments>, Size,                                    (ScalarArguments ? 1 : Size) >,
    public nova::multichannel::ControlInput<  NovaSVFShelfBase<FilterType, Size, ScalarArguments>, Size + (ScalarArguments ? 1 : Size),     (ScalarArguments ? 1 : Size) >,
    public nova::multichannel::ControlInput<  NovaSVFShelfBase<FilterType, Size, ScalarArguments>, Size + (ScalarArguments ? 1 : Size) * 2, (ScalarArguments ? 1 : Size) >

{
    static const size_t ParameterSize = ScalarArguments ? 1 : Size;
    static const int FreqInputIndex = Size;
    static const int AmpInputIndex  = Size + ParameterSize;
    static const int ResInputIndex  = Size + ParameterSize * 2;

    using SignalInput = nova::multichannel::SignalInput< NovaSVFShelfBase, 0, Size >;
    using OutputSink  = nova::multichannel::OutputSink<  NovaSVFShelfBase, 0, Size >;

    using FreqInput = nova::multichannel::ControlInput<  NovaSVFShelfBase, Size,                     ParameterSize >;
    using AmpInput  = nova::multichannel::ControlInput<  NovaSVFShelfBase, Size + ParameterSize,     ParameterSize >;
    using ResInput  = nova::multichannel::ControlInput<  NovaSVFShelfBase, Size + ParameterSize * 2, ParameterSize >;

    using SampleType    = typename nova::as_pack<float, Size>::type;
    using ParameterType = typename nova::as_pack<float, ParameterSize>::type;

    typedef boost::mpl::map<
            boost::mpl::pair<boost::mpl::long_<svfLowShelf>,  SVFLowShelf  <Size, ParameterSize> >,
            boost::mpl::pair<boost::mpl::long_<svfHighShelf>, SVFHighShelf <Size, ParameterSize> >,
            boost::mpl::pair<boost::mpl::long_<svfEq>,        SVFEQ        <Size, ParameterSize> >
            > FilterMap;


    using Filter = typename boost::mpl::at< FilterMap, boost::mpl::long_<FilterType> >::type;
    Filter _filter;

    NovaSVFShelfBase()
    {
        initFilter();
        setCalcFunction< NovaSVFShelfBase >( FreqInput::index, ResInput::index );
    }

    void initFilter()
    {
        auto newFreq = FreqInput::readInput();
        auto amp     = AmpInput::readInput();
        auto res     = ResInput::readInput();

        auto state = Filter::computeState( newFreq, amp, res, makeDSPContext() );
        _filter.setState( state );
    }

    template <typename ControlSignature>
    void run(int inNumSamples)
    {
        next( inNumSamples, ControlSignature() );
    }

    void next(int inNumSamples, nova::control_signature_kk)
    {
        if ( !FreqInput::changed() && !AmpInput::changed() && !ResInput::changed() ) {
            next( inNumSamples, nova::control_signature_ii() );
            return;
        }

        auto oldState = _filter.currentState();
        auto newState = Filter::computeState( FreqInput::readAndUpdateInput(), AmpInput::readAndUpdateInput(),
                                              ResInput::readAndUpdateInput(), makeDSPContext() );
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
        _filter.run( inFn, outFn, inNumSamples );
    }

    void next(int inNumSamples, nova::control_signature_aa)
    {
        auto context = makeDSPContext();

        auto makeState = [this, index = 0, context] () mutable {
            ParameterType freq = FreqInput::template readInputs<ParameterType>( index );
            ParameterType amp  = AmpInput:: template readInputs<ParameterType>( index );
            ParameterType res  = ResInput:: template readInputs<ParameterType>( index );
            index += 1;

            return Filter::computeState( freq, amp, res, context );
        };

        auto inFn  = SignalInput::template makeInputSignal<SampleType>();
        auto outFn = OutputSink ::template makeSink<SampleType>();
        _filter.run( inFn, outFn, inNumSamples, makeState );
    }
};



typedef NovaSVFBase<svfLPF, 1, true> NovaLowPassSVF;
typedef NovaSVFBase<svfLPF, 2, true> NovaLowPassSVF2;
typedef NovaSVFBase<svfLPF, 4, true> NovaLowPassSVF4;
typedef NovaSVFBase<svfLPF, 8, true> NovaLowPassSVF8;
typedef NovaSVFBase<svfHPF, 1, true> NovaHighPassSVF;

typedef NovaSVFShelfBase<svfEq, 1, true> NovaEqSVF;
typedef NovaSVFShelfBase<svfEq, 2, true> NovaEqSVF2;

DEFINE_XTORS( NovaLowPassSVF )
DEFINE_XTORS( NovaLowPassSVF2 )
DEFINE_XTORS( NovaLowPassSVF4 )
DEFINE_XTORS( NovaLowPassSVF8 )

DEFINE_XTORS( NovaHighPassSVF )

DEFINE_XTORS( NovaEqSVF )
DEFINE_XTORS( NovaEqSVF2 )

}

using namespace nova;

PluginLoad(NovaSVF)
{
    ft = inTable;

    NovaDefineUnit( NovaLowPassSVF  );
    NovaDefineUnit( NovaLowPassSVF2 );
    NovaDefineUnit( NovaLowPassSVF4 );
    NovaDefineUnit( NovaLowPassSVF8 );

    NovaDefineUnit( NovaHighPassSVF  );

    NovaDefineUnit( NovaEqSVF  );
    NovaDefineUnit( NovaEqSVF2 );
}
