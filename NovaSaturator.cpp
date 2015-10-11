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


namespace {

InterfaceTable *ft;

template <typename Parent>
struct SaturationBase:
    public SCUnit
{
    static float verifyLevel ( float arg ) { return arg; }

    SaturationBase()
    {
        switch (inRate(1)) {
        case calc_FullRate:
            if (boost::simd::is_aligned( bufferSize(), 8 ) )
                set_calc_function<SaturationBase, &SaturationBase::run_a< boost::simd::pack<float, 8> > >();
            else
                set_calc_function<SaturationBase, &SaturationBase::run_a<float>>();
            break;

        case calc_BufRate:
            _level = Parent::verifyLevel( in0(1) );
            if (boost::simd::is_aligned( bufferSize(), 8 ) )
                set_calc_function<SaturationBase, &SaturationBase::run_k< boost::simd::pack<float, 8> > >();
            else
                set_calc_function<SaturationBase, &SaturationBase::run_k<float>>();
            break;

        case calc_ScalarRate:
        default:
            _level = Parent::verifyLevel( in0(1) );
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
        auto input0 = nova::Packer<SampleType, 0>(this);
        auto input1 = nova::Packer<SampleType, 1>(this);
        auto output = nova::Unpacker<SampleType, 0>(this);

        perform<SampleType>( inNumSamples, input0, input1, output );
    }

    template <typename SampleType>
    void run_k (int inNumSamples)
    {
        float newLevel = Parent::verifyLevel( in0(1) );
        if (newLevel != _level) {
            float slope = calcSlope(newLevel, _level);

            auto input0 = nova::Packer<SampleType, 0>(this);
            auto input1 = nova::makeRamp<SampleType>(_level, slope);
            auto output = nova::Unpacker<SampleType, 0>(this);
            _level = newLevel;

            perform<SampleType>( inNumSamples, input0, input1, output );
        } else {
            auto input0 = nova::Packer<SampleType, 0>(this);
            auto input1 = nova::Scalar<SampleType>(_level);
            auto output = nova::Unpacker<SampleType, 0>(this);

            perform<SampleType>( inNumSamples, input0, input1, output );
        }
    }

    template <typename SampleType>
    void run_i (int inNumSamples)
    {
        auto input0 = nova::Packer<SampleType, 0>(this);
        auto input1 = nova::Scalar<SampleType>(_level);
        auto output = nova::Unpacker<SampleType, 0>(this);

        perform<SampleType>( inNumSamples, input0, input1, output );
    }

    float _level;
};

class HyperbolSaturation : public SaturationBase<HyperbolSaturation>
{
public:
    HyperbolSaturation() = default;

    template <typename SampleType>
    static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, SampleType level)
    {
        return nova::saturator::hyperbol( sig, level );
    }

    static float verifyLevel(float arg)
    {
        return std::max( arg, 1e-10f );
    }
};

class ParabolSaturation : public SaturationBase<ParabolSaturation>
{
public:
    ParabolSaturation() = default;

    template <typename SampleType>
    static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, SampleType level)
    {
        return nova::saturator::parabol( sig, level );
    }

    static float verifyLevel(float arg)
    {
        return std::max( arg, 1e-10f );
    }
};


class PowSaturation : public SaturationBase<PowSaturation>
{
public:
    PowSaturation() = default;

    template <typename SampleType>
    static BOOST_FORCEINLINE FLATTEN SampleType doDistort(SampleType sig, SampleType level)
    {
        return nova::saturator::pow( sig, level );
    }

    static float verifyLevel(float arg)
    {
        return std::max( arg, 1e-10f );
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
