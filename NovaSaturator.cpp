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
#include <boost/simd/include/functions/floor.hpp>
#include <boost/simd/include/functions/abs.hpp>

#include <approximations/tanh.hpp>
#include <approximations/pow.hpp>
#include <approximations/sin.hpp>

namespace {

InterfaceTable *ft;

enum {
    computePrecise,
    computeFast,
    computeRaw
};


typedef nova::detail::Identity TrivialInputFunctor;

struct VerifyLevelInput :
    TrivialInputFunctor
{
    template <typename Arg>
    auto operator()(Arg && arg) const
    {
        return boost::simd::max( arg, 1e-10f );
    }
};

template <int Precision>
struct PowDistortion
{
    typedef VerifyLevelInput InputFunctor;

    template <typename SampleType, typename LevelType>
    static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, LevelType level)
    {
        switch( Precision ) {
        case computePrecise:
            return nova::saturator::pow( sig, level );

        case computeFast: {
            SampleType absSig = boost::simd::abs(sig);
            return boost::simd::copysign( nova::approximations::pow( absSig, level,
                                                                     nova::approximations::PowFast()
                                                                     ), sig );
        }
        case computeRaw: {
            SampleType absSig = boost::simd::abs(sig);
            return boost::simd::copysign( nova::approximations::pow( absSig, level,
                                                                     nova::approximations::PowFaster()
                                                                     ), sig );
        }
        }
    }
};

struct HyperbolDistortion
{
    typedef VerifyLevelInput InputFunctor;

    template <typename SampleType, typename LevelType>
    static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, LevelType level)
    {
        return nova::saturator::hyperbol( sig, level );
    }
};

struct ParabolDistortion
{
    typedef VerifyLevelInput InputFunctor;

    template <typename SampleType, typename LevelType>
    static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, LevelType level)
    {
        return nova::saturator::parabol( sig, level );
    }
};

template <int Precision>
struct TanhDistortion
{
    struct InputFunctor
    {
        template <typename Arg>
        auto operator()(Arg const & arg) const
        {
            Arg preGain  = boost::simd::max( arg, 1e-10f );
            Arg postGain = boost::simd::fast_rec( computeTanh( preGain ) );

            using Scalar                 = typename boost::simd::meta::scalar_of<Arg>::type;
            static const size_t Cardinal = boost::simd::meta::cardinal_of<Arg>::value;

            using ResultType = typename nova::as_pack< Scalar, Cardinal >::type;

            typedef nova::detail::ArithmeticArray<ResultType, 2> Result;

            return Result { preGain, postGain };
        }

        template <typename Type>
        struct State {
            typedef nova::detail::ArithmeticArray<Type, 2> type;
        };
    };

    template <typename SampleType, typename StateType>
    static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, StateType state)
    {
        SampleType preGained = sig * state[0];
        return computeTanh( preGained ) * state[1];
    }

    template <typename SampleType>
    static BOOST_FORCEINLINE SampleType computeTanh(SampleType x)
    {
        switch( Precision ) {
        case computePrecise: return nt2::tanh( x );

        case computeFast:    return nova::approximations::fast_tanh( x );

        case computeRaw:     return nova::approximations::faster_tanh( x );
        default:             assert(false); return x;
        }
    }
};


template <typename SampleType>
inline SampleType fold(SampleType x)
{
    namespace bs = boost::simd;

    SampleType one = bs::One<SampleType>();
    return one - bs::abs( x - 4.f * bs::floor( ( x + bs::One<SampleType>() ) * 0.25f ) - one );
}


struct FoldDistortion
{
    typedef TrivialInputFunctor InputFunctor;

    template <typename SampleType, typename ScaleType>
    static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, ScaleType scale)
    {
        return fold< SampleType >( SampleType(sig * scale) );
    }
};

template <int Precision>
struct SinFoldDistortion
{
    typedef TrivialInputFunctor InputFunctor;

    template <typename SampleType, typename ScaleType>
    static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, ScaleType scale)
    {
        SampleType folded = fold< SampleType >( SampleType(sig * scale) ) * (boost::simd::Pio_2<float>() );

        using namespace nova::approximations;

        switch( Precision ) {
        case computeFast:    return sin( folded, SinFast()    );
//        case computePrecise: return sin( folded, SinPrecise() );
        case computeRaw:     return sin( folded, SinFaster()  );

        default: assert(false);
        }
    }
};

template <typename SampleType>
inline SampleType wrap(SampleType x)
{
    namespace bs = boost::simd;

    SampleType y = bs::floor( (x+1.f) * 0.5f );
    return x - y - y ;
}

struct WrapDistortion
{
    typedef TrivialInputFunctor InputFunctor;

    template <typename SampleType, typename ScaleType>
    static BOOST_FORCEINLINE SampleType doDistort(SampleType sig, ScaleType scale)
    {
        return wrap< SampleType >( sig * scale );
    }
};


}

namespace {

using namespace nova;

template <typename Distortion>
struct SaturationBase:
    public nova::NovaUnit,
    public nova::SignalInput<  SaturationBase<Distortion>, 0 >,
    public nova::ControlInput< SaturationBase<Distortion>, 1, typename Distortion::InputFunctor >,

    public nova::OutputSink<   SaturationBase<Distortion>, 0 >
{
    typedef typename Distortion::InputFunctor InputFunctor;

    typedef nova::SignalInput<  SaturationBase, 0 >               SignalInput;
    typedef nova::ControlInput< SaturationBase, 1, InputFunctor > LevelInput;

    typedef nova::OutputSink<   SaturationBase, 0 >               SignalOutput;


    using vector_type     = boost::simd::pack<float, 8>;
    const int vector_size = boost::simd::meta::cardinal_of<vector_type>::value;

    SaturationBase()
    {
        setVectorCalcFunction< SaturationBase, vector_type >( 1 );
    }

    template <typename VectorType, typename ControlSignature>
    void run(int inNumSamples)
    {
        run<VectorType>( inNumSamples, ControlSignature() );
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

        const size_t unroll = boost::simd::meta::cardinal_of<SampleType>::value;
        loop( inNumSamples / unroll, [&] {
            auto in0 = input0();
            auto in1 = input1();
            auto result = Distortion::doDistort( in0, in1 );
            output(result);
        });
    }


    template <typename SampleType>
    void run (int inNumSamples, control_signature_a)
    {
        auto input0 = SignalInput ::template makeInputSignal<SampleType>();
        auto input1 = LevelInput  ::template makeAudioInputSignal<SampleType>();
        auto output = SignalOutput::template makeSink<SampleType>();

        perform<SampleType>( inNumSamples, input0, input1, output );
    }

    template <typename SampleType>
    void run (int inNumSamples, control_signature_k)
    {
        if ( LevelInput::changed() ) {
            auto input0 = SignalInput ::template makeInputSignal<SampleType>();
            auto input1 = LevelInput  ::template makeRampSignal<SampleType>();
            auto output = SignalOutput::template makeSink<SampleType>();

            perform<SampleType>( inNumSamples, input0, input1, output );
        } else {
            run<SampleType>( inNumSamples, control_signature_i() );
        }
    }

    template <typename SampleType>
    void run (int inNumSamples, control_signature_i)
    {
        auto input0 = SignalInput ::template makeInputSignal<SampleType>();
        auto input1 = LevelInput  ::template makeScalarInputSignal<SampleType>();
        auto output = SignalOutput::template makeSink<SampleType>();

        perform<SampleType>( inNumSamples, input0, input1, output );
    }

    template <typename SampleType>
    void run (int inNumSamples, control_signature_1)
    {
        auto input0 = SignalInput ::template makeInputSignal<SampleType>();
        auto input1 = LevelInput  ::template makeScalarInputSignal<SampleType>();
        auto output = SignalOutput::template makeSink<SampleType>();

        perform<SampleType>( 1, input0, input1, output );
    }
};

}

PluginLoad(NovaSaturators)
{
    using namespace nova;

    ft = inTable;

    nova::registerUnit<SaturationBase<PowDistortion<computeFast>>>( ft, "PowSaturation" );
    nova::registerUnit<SaturationBase<HyperbolDistortion>>        ( ft, "HyperbolSaturation" );
    nova::registerUnit<SaturationBase<ParabolDistortion>>         ( ft, "ParabolSaturation" );

    nova::registerUnit<SaturationBase<TanhDistortion<computePrecise>>>( ft, "TanhSaturation" );
    nova::registerUnit<SaturationBase<TanhDistortion<computeFast>>>   ( ft, "FastTanhSaturation" );
    nova::registerUnit<SaturationBase<TanhDistortion<computeRaw>>>    ( ft, "FasterTanhSaturation" );

    nova::registerUnit<SaturationBase<FoldDistortion>>                 ( ft, "NovaLinearWaveFolder" );
    nova::registerUnit<SaturationBase<SinFoldDistortion<computeFast>>> ( ft, "NovaSinWaveFolder" );
    nova::registerUnit<SaturationBase<SinFoldDistortion<computeRaw>>>  ( ft, "NovaSinWaveFolderRaw" );
    nova::registerUnit<SaturationBase<WrapDistortion>>                 ( ft, "NovaLinearWaveWrapper" );
}
