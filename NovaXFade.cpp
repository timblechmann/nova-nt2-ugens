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
#include "producer_consumer_functors.hpp"
#include "dsp/utils.hpp"
#include "approximations/parabol_sine.hpp"

#include <cmath>

#include <boost/simd/include/pack.hpp>
#include <boost/simd/constant/constants/half.hpp>
#include <boost/simd/constant/constants/one.hpp>
#include <boost/simd/constant/constants/pi.hpp>

#include <approximations/tan.hpp>

using nova::NovaUnit;

InterfaceTable *ft;

namespace nova {
namespace PanningLaws {

struct EqualPowerPanning
{
    template <typename SampleType>
    auto operator() (SampleType input)
    {
        SampleType piOver4 = boost::simd::Pi<SampleType>() * 0.25f;
        SampleType arg = input * piOver4;
        SampleType leftGain, rightGain;

        leftGain  = approximations::sin<SampleType>( piOver4 - input*piOver4, approximations::SinFast() );
        rightGain = approximations::sin<SampleType>( piOver4 + input*piOver4, approximations::SinFast() );

        typedef detail::ArithmeticArray<SampleType, 2> Result;

        Result ret;
        ret[0] = leftGain;
        ret[1] = rightGain;
        return ret;
    }

    template <typename Type>
    struct State {
        typedef detail::ArithmeticArray<Type, 2> type;
    };
};


struct LinearPanning
{
    template <typename SampleType>
    auto operator() (SampleType input)
    {
        SampleType half      = boost::simd::Half<SampleType>();
        SampleType rightGain = input * half + half;
        SampleType leftGain  = boost::simd::One<SampleType>() - rightGain;

        typedef detail::ArithmeticArray<SampleType, 2> Result;

        Result ret;
        ret[0] = leftGain;
        ret[1] = rightGain;
        return ret;
    }

    template <typename Type>
    struct State {
        typedef detail::ArithmeticArray<Type, 2> type;
    };
};

}


template <typename PanningLawFunctor>
struct NovaXFade:
    public NovaUnit,
    public nova::SignalInput<  NovaXFade<PanningLawFunctor>, 0 >,
    public nova::SignalInput<  NovaXFade<PanningLawFunctor>, 1 >,
    public nova::ControlInput< NovaXFade<PanningLawFunctor>, 2, PanningLawFunctor >,
    public nova::ControlInput< NovaXFade<PanningLawFunctor>, 3 >,

    public nova::OutputSink<   NovaXFade<PanningLawFunctor>, 0 >
{
    typedef nova::SignalInput<  NovaXFade<PanningLawFunctor>, 0 > LeftInput;
    typedef nova::SignalInput<  NovaXFade<PanningLawFunctor>, 1 > RightInput;
    typedef nova::ControlInput< NovaXFade<PanningLawFunctor>, 2, PanningLawFunctor > PanInput;
    typedef nova::ControlInput< NovaXFade<PanningLawFunctor>, 3 > LevelInput;


    typedef nova::OutputSink<   NovaXFade<PanningLawFunctor>, 0 > OutputSink;


    using vector_type     = boost::simd::pack<float, 8>;
    const int vector_size = boost::simd::meta::cardinal_of<vector_type>::value;

    NovaXFade()
    {
        if( LevelInput::scalarRate() && LevelInput::readInput() == 1.f )
            setVectorCalcFunction< NovaXFade, vector_type >( 2 );
        else
            setVectorCalcFunction< NovaXFade, vector_type >( 2, 3 );
    }

    template <typename VectorType, typename ControlSignature>
    void run(int inNumSamples)
    {
        next<VectorType>( inNumSamples, ControlSignature() );
    }


    template <typename VectorType>
    void next(int inNumSamples, control_signature_ak)
    {
        auto panFn      = PanInput::template makeAudioInputSignal<VectorType>();

        if( LevelInput::changed() ) {
            auto levelFn = LevelInput::template makeRampSignal<VectorType>();
            return next<VectorType, calc_FullRate>( panFn, levelFn, inNumSamples );
        } else {
            auto levelFn = LevelInput::template makeScalarInputSignal<VectorType>();
            next<VectorType, calc_BufRate>( panFn, levelFn, inNumSamples );
        }
    }

    template <typename VectorType>
    void next(int inNumSamples, control_signature_ka)
    {
        auto levelFn     = LevelInput::template makeAudioInputSignal<VectorType>();

        if( PanInput::changed() ) {
            auto panFn   = PanInput::  template makeRampSignal<VectorType>();
            next<VectorType, calc_FullRate>( panFn, levelFn, inNumSamples );
        } else {
            auto panFn   = PanInput::  template makeScalarInputSignal<VectorType>();
            next<VectorType, calc_BufRate>( panFn, levelFn, inNumSamples );
        }
    }

    template <typename VectorType>
    void next(int inNumSamples, control_signature_aa)
    {
        auto panFn      = PanInput::  template makeAudioInputSignal<VectorType>();
        auto levelFn    = LevelInput::template makeAudioInputSignal<VectorType>();
        next<VectorType, calc_FullRate>( panFn, levelFn, inNumSamples );
    }

    template <typename VectorType>
    void next(int inNumSamples, control_signature_kk)
    {
        if( PanInput::changed() || LevelInput::changed() ) {
            auto panFn   = PanInput::  template makeRampSignal<VectorType>();
            auto levelFn = LevelInput::template makeRampSignal<VectorType>();
            next<VectorType, calc_FullRate>( panFn, levelFn, inNumSamples );
        } else {
            auto panFn   = PanInput::  template makeScalarInputSignal<VectorType>();
            auto levelFn = LevelInput::template makeScalarInputSignal<VectorType>();
            next<VectorType, calc_BufRate>( panFn, levelFn, inNumSamples );
        }
    }

    template <typename VectorType>
    void next(int inNumSamples, control_signature_11)
    {
        auto panFn   = PanInput::  template makeScalarInputSignal<VectorType>();
        auto levelFn = LevelInput::template makeScalarInputSignal<VectorType>();
        next<VectorType, calc_BufRate>( panFn, levelFn, inNumSamples );
    }

    template <typename VectorType, int parameterRate, typename PanFn, typename LevelFn>
    void next( PanFn && panFn, LevelFn && levelFn, int inNumSamples )
    {
        const int vector_size = boost::simd::meta::cardinal_of<VectorType>::value;

        auto leftFn     = LeftInput:: template makeInputSignal<VectorType>();
        auto rightFn    = RightInput::template makeInputSignal<VectorType>();
        auto sink       = OutputSink::template makeSink<VectorType>();

        if( parameterRate == calc_FullRate ) {
            loop( inNumSamples / vector_size, [&] {
                auto panGains  = panFn();
                auto levelGain = levelFn();
                auto leftGain  = panGains[0] * levelGain;
                auto rightGain = panGains[1] * levelGain;

                sink( leftGain * leftFn() + rightGain * rightFn() );
            });
        } else {
            auto panGains  = panFn();
            auto levelGain = levelFn();
            auto leftGain  = panGains[0] * levelGain;
            auto rightGain = panGains[1] * levelGain;

            loop( inNumSamples / vector_size, [&] {
                sink( leftGain * leftFn() + rightGain * rightFn() );
            });
        }
    }

    // level = 1
    template <typename VectorType>
    void next(int inNumSamples, control_signature_k)
    {
        if( PanInput::changed() ) {
            auto panFn = PanInput::template makeRampSignal<VectorType>();
            next<VectorType, calc_FullRate>( panFn, inNumSamples );
        } else {
            next<VectorType>( inNumSamples, control_signature_i() );
        }
    }

    template <typename VectorType>
    void next(int inNumSamples, control_signature_i)
    {
        auto panFn = PanInput::template makeScalarInputSignal<VectorType>();
        next<VectorType, calc_BufRate>( panFn, inNumSamples );
    }


    template <typename VectorType>
    void next(int inNumSamples, control_signature_a)
    {
        auto panFn = PanInput::template makeAudioInputSignal<VectorType>();
        next<VectorType, calc_BufRate>( panFn, inNumSamples );
    }

    template <typename VectorType, int parameterRate, typename PanFn>
    void next( PanFn && panFn, int inNumSamples )
    {
        const int vector_size = boost::simd::meta::cardinal_of<VectorType>::value;

        auto leftFn     = LeftInput:: template makeInputSignal<VectorType>();
        auto rightFn    = RightInput::template makeInputSignal<VectorType>();
        auto sink       = OutputSink::template makeSink<VectorType>();

        if( parameterRate == calc_FullRate ) {
            loop( inNumSamples / vector_size, [&] {
                auto panGains  = panFn();
                auto leftGain  = panGains[0];
                auto rightGain = panGains[1];

                sink( leftGain * leftFn() + rightGain * rightFn() );
            });
        } else {
            auto panGains  = panFn();
            auto leftGain  = panGains[0];
            auto rightGain = panGains[1];

            loop( inNumSamples / vector_size, [&] {
                sink( leftGain * leftFn() + rightGain * rightFn() );
            });
        }
    }
};

}

typedef nova::NovaXFade<nova::PanningLaws::EqualPowerPanning> NovaXFade;
typedef nova::NovaXFade<nova::PanningLaws::LinearPanning>     NovaLinXFade;

DEFINE_XTORS(NovaXFade)
DEFINE_XTORS(NovaLinXFade)

PluginLoad(NovaXFade)
{
    ft = inTable;
    DefineSimpleUnit(NovaXFade);
    DefineSimpleUnit(NovaLinXFade);
}
