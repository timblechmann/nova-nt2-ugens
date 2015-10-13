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


#include "SC_PlugIn.hpp"

#include "boost/simd/include/pack.hpp"
#include "producer_consumer_functors.hpp"

#include "dsp/saturators.hpp"
#include "dsp/utils.hpp"

#include <NovaUGensCommon.hpp>

#include <boost/simd/include/functions/max.hpp>

namespace {

InterfaceTable *ft;

struct VerifyLevelInput
{
    template <typename Arg>
    auto operator()(Arg && arg) const
    {
        return boost::simd::max( arg, 1e-10f );
    }
};

template <typename Parent>
struct SaturationBase:
    public nova::NovaUnit,
    public nova::SignalInput<SaturationBase<Parent>, 0>,
    public nova::ControlInput<SaturationBase<Parent>, 1, VerifyLevelInput>,

    public nova::OutputSink<SaturationBase<Parent>, 0>
{
    typedef nova::SignalInput<SaturationBase<Parent>, 0>                    SignalInput;
    typedef nova::ControlInput<SaturationBase<Parent>, 1, VerifyLevelInput> LevelInput;

    typedef nova::OutputSink<SaturationBase<Parent>, 0>                     SignalOutput;


    SaturationBase()
    {
        switch (SCUnit::inRate(1)) {
        case calc_FullRate:
            if (boost::simd::is_aligned( bufferSize(), 8 ) )
                set_calc_function<SaturationBase, &SaturationBase::run_a< boost::simd::pack<float, 8> > >();
            else
                set_calc_function<SaturationBase, &SaturationBase::run_a<float>>();
            break;

        case calc_BufRate:
            if (boost::simd::is_aligned( bufferSize(), 8 ) )
                set_calc_function<SaturationBase, &SaturationBase::run_k< boost::simd::pack<float, 8> > >();
            else
                set_calc_function<SaturationBase, &SaturationBase::run_k<float>>();
            break;

        case calc_ScalarRate:
        default:
            if (boost::simd::is_aligned( bufferSize(), 8 ) )
                set_calc_function<SaturationBase, &SaturationBase::run_i< boost::simd::pack<float, 8> > >();
            else
                set_calc_function<SaturationBase, &SaturationBase::run_i<float>>();
        }

    }

    template <typename Functor>
    inline void loop (int loops, Functor const & f)
    {
#ifdef __GCC__
        _Pragma( GCC ivdep )
#endif
        for (int i = 0; i != loops; ++i)
            f();
    }

    template <typename SampleType,
              typename Input0,
              typename Input1,
              typename Output0>
    inline void perform(int inNumSamples, Input0 & input0, Input1 & input1, Output0 & output)
    {
        using namespace boost::simd;

        const size_t unroll = meta::cardinal_of<SampleType>::value;
        loop( inNumSamples / unroll, [&] {
            auto in0 = input0();
            auto in1 = input1();
            auto result = Parent::doDistort( in0, in1 );
            output(result);
        });
    }


    template <typename SampleType>
    void run_a (int inNumSamples)
    {
        auto input0 = SignalInput ::template makeInputSignal<SampleType>();
        auto input1 = LevelInput  ::template makeAudioInputSignal<SampleType>();
        auto output = SignalOutput::template makeSink<SampleType>();

        perform<SampleType>( inNumSamples, input0, input1, output );
    }

    template <typename SampleType>
    void run_k (int inNumSamples)
    {
        auto input0 = SignalInput::template makeInputSignal<SampleType>();
        if ( LevelInput::changed() ) {
            auto input1 = LevelInput  ::template makeRampSignal<SampleType>();
            auto output = SignalOutput::template makeSink<SampleType>();

            perform<SampleType>( inNumSamples, input0, input1, output );
        } else {
            auto input1 = LevelInput  ::template makeScalarInputSignal<SampleType>();
            auto output = SignalOutput::template makeSink<SampleType>();

            perform<SampleType>( inNumSamples, input0, input1, output );
        }
    }

    template <typename SampleType>
    void run_i (int inNumSamples)
    {
        auto input0 = SignalInput ::template makeInputSignal<SampleType>();
        auto input1 = LevelInput  ::template makeScalarInputSignal<SampleType>();
        auto output = SignalOutput::template makeSink<SampleType>();

        perform<SampleType>( inNumSamples, input0, input1, output );
    }
};

class HyperbolSaturation : public SaturationBase<HyperbolSaturation>
{
public:
    HyperbolSaturation() {}

    template <typename SampleType>
    static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, SampleType level)
    {
        return nova::saturator::hyperbol( sig, level );
    }
};

class ParabolSaturation : public SaturationBase<ParabolSaturation>
{
public:
    ParabolSaturation()  {}

    template <typename SampleType>
    static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, SampleType level)
    {
        return nova::saturator::parabol( sig, level );
    }
};


class PowSaturation : public SaturationBase<PowSaturation>
{
public:
    PowSaturation()  {}

    template <typename SampleType>
    static BOOST_FORCEINLINE FLATTEN SampleType doDistort(SampleType sig, SampleType level)
    {
        return nova::saturator::pow( sig, level );
    }
};



DEFINE_XTORS(HyperbolSaturation)
DEFINE_XTORS(ParabolSaturation)
DEFINE_XTORS(PowSaturation)

}

PluginLoad(NovaSaturators)
{
    ft = inTable;
    DefineDtorUnit(HyperbolSaturation);
    DefineDtorUnit(ParabolSaturation);
    DefineDtorUnit(PowSaturation);
}
