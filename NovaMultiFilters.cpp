/*
    Copyright (C) 2013 Tim Blechmann

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


#include "NovaUnitFacade.hpp"

#include <boost/mpl/if.hpp>

#include "dsp/biquad.hpp"
#include "dsp/integrator.hpp"
#include "dsp/leak_dc.hpp"
#include "producer_consumer_functors.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/simd/include/constants/pi.hpp>

#include <boost/simd/sdk/simd/logical.hpp>
#include <boost/simd/include/functions/compare_equal.hpp>
#include <boost/simd/include/functions/eq.hpp>
#include <boost/simd/include/functions/negate.hpp>
#include <boost/simd/include/functions/fast_divides.hpp>

#include <nt2/include/functions/tan.hpp>


#include <dsp/utils.hpp>

#include <cmath>

using nova::NovaUnit;

InterfaceTable *ft;

namespace constants = boost::math::constants;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace nova {
namespace      {

template <size_t NumberOfChannels, bool ScalarArguments = true>
struct NovaIntegrator:
    public NovaUnitUnary<NovaIntegrator<NumberOfChannels, ScalarArguments>, NumberOfChannels, double, ScalarArguments>
{
    typedef NovaUnitUnary<NovaIntegrator<NumberOfChannels, ScalarArguments>, NumberOfChannels, double, ScalarArguments> Base;

    typedef Integrator<typename Base::SampleType, typename Base::ParameterDSPType> DSPEngine;

    NovaIntegrator() = default;

    static DSPEngine & getDSPEngine( NovaIntegrator * self )
    {
        return self->_filter;
    }

    DSPEngine _filter;
};

}
}

typedef nova::NovaIntegrator<1> NovaIntegrator1;
typedef nova::NovaIntegrator<2> NovaIntegrator2;
typedef nova::NovaIntegrator<4> NovaIntegrator4;
//typedef nova::NovaIntegrator<8> NovaIntegrator8;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename ParameterType>
struct BiquadParameterStruct
{
    ParameterType a0() { return _a0; }
    ParameterType a1() { return _a1; }
    ParameterType a2() { return _a2; }
    ParameterType b1() { return _b1; }
    ParameterType b2() { return _b2; }

    typedef ParameterType type;
    ParameterType _a0, _a1, _a2, _b1, _b2;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <typename FilterDesigner, size_t Size, bool ScalarArguments = true, bool Fast = true>
struct NovaBiquadBase:
    public NovaUnit,
    public nova::multichannel::SignalInput< NovaBiquadBase<FilterDesigner, Size, ScalarArguments, Fast>, 0, Size >,
    public nova::multichannel::OutputSink<  NovaBiquadBase<FilterDesigner, Size, ScalarArguments, Fast>, 0, Size >,

    // Controls:
    public nova::multichannel::ControlInput<  NovaBiquadBase<FilterDesigner, Size, ScalarArguments, Fast>, Size,                                (ScalarArguments ? 1 : Size) >,
    public nova::multichannel::ControlInput<  NovaBiquadBase<FilterDesigner, Size, ScalarArguments, Fast>, Size + (ScalarArguments ? 1 : Size), (ScalarArguments ? 1 : Size) >
{
    static const size_t ParameterSize = ScalarArguments ? 1 : Size;
    static const int FreqInputIndex = Size;
    static const int QInputIndex    = Size + ParameterSize;

    typedef typename boost::mpl::if_c<Fast, float, double>::type InternalType;

    typedef typename nova::as_pack<float,  Size>::type vFloat;
    typedef typename nova::as_pack<double, Size>::type vInternal;

    typedef typename boost::mpl::if_c<ScalarArguments, InternalType, vInternal>::type ParameterType;
    typedef typename boost::mpl::if_c<ScalarArguments, float,  vFloat>::type HostParameterType;

    typedef typename FilterDesigner::template makeParameter<ParameterType>::type BiquadParameterStruct;
    typedef nova::Biquad<vInternal, BiquadParameterStruct> Filter;

    using SignalInput = nova::multichannel::SignalInput< NovaBiquadBase, 0, Size >;
    using OutputSink  = nova::multichannel::OutputSink<  NovaBiquadBase, 0, Size >;

    using FreqInput = nova::multichannel::ControlInput<  NovaBiquadBase, Size,                                ScalarArguments ? 1 : Size >;
    using QInput    = nova::multichannel::ControlInput<  NovaBiquadBase, Size + (ScalarArguments ? 1 : Size), ScalarArguments ? 1 : Size >;


    NovaBiquadBase()
    {
        initFilter();

        if (isScalarRate (FreqInputIndex, FreqInputIndex + ParameterSize) )
            setCalcFunction<NovaBiquadBase, nova::control_signature_i>();
        else
            setCalcFunction<NovaBiquadBase, nova::control_signature_k>();
    }

    void initFilter()
    {
        auto newFreq = FreqInput::readInput();
        auto newQ    = QInput::readInput();

        HostParameterType normalizedFreq = newFreq * (float)sampleDur();
        auto params = FilterDesigner::template designFilter<BiquadParameterStruct>( normalizedFreq, newQ );

        storeFilterParameters( params );
    }

    template <typename ControlSignature>
    void run(int inNumSamples)
    {
        next( inNumSamples, ControlSignature() );
    }

    void next(int, nova::control_signature_1)
    {
        auto inFn  = SignalInput::template makeInputSignal<vInternal>();
        auto outFn = OutputSink ::template makeSink<vInternal>();

        _filter.run(inFn, outFn, 1);
    }

    void next(int, nova::control_signature_i)
    {
        auto inFn  = SignalInput::template makeInputSignal<vInternal>();
        auto outFn = OutputSink ::template makeSink<vInternal>();

        _filter.run_unrolled(inFn, outFn, mRate->mFilterLoops, mRate->mFilterRemain);
    }

    void next(int inNumSamples, nova::control_signature_k)
    {
        auto newFreq = FreqInput::readInput();
        auto newQ    = QInput   ::readInput();

        if ( !FreqInput::changed() && !QInput::changed() ) {
            next(inNumSamples, nova::control_signature_i());
            return;
        }

        HostParameterType normalizedFreq = newFreq * (float)sampleDur();
        auto newParameters = FilterDesigner::template designFilter<BiquadParameterStruct>( normalizedFreq, newQ );

        ParameterType a0Slope = calcSlope( newParameters.a0(), _filter._parameters.a0() );
        ParameterType a1Slope = calcSlope( newParameters.a1(), _filter._parameters.a1() );
        ParameterType a2Slope = calcSlope( newParameters.a2(), _filter._parameters.a2() );
        ParameterType b1Slope = calcSlope( newParameters.b1(), _filter._parameters.b1() );
        ParameterType b2Slope = calcSlope( newParameters.b2(), _filter._parameters.b2() );

        auto inFn  = SignalInput::template makeInputSignal<vInternal>();
        auto outFn = OutputSink ::template makeSink<vInternal>();

        _filter.run_unrolled(inFn, outFn, mRate->mFilterLoops, mRate->mFilterRemain,
                             a0Slope, a1Slope, a2Slope, b1Slope, b2Slope );
        storeFilterParameters( newParameters );
    }

    void storeFilterParameters(BiquadParameterStruct const & parameters)
    {
        _filter._parameters = parameters;
    }

    Filter _filter;
};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using nova::castType;

struct DesignLPF
{
    template <typename ParameterType>
    struct ParameterStruct
    {
        typedef ParameterType type;

        ParameterType a0() { return _a0_a1_a2; }
        ParameterType a1() { return _a0_a1_a2 + _a0_a1_a2; }
        ParameterType a2() { return _a0_a1_a2; }
        ParameterType b1() { return _b1;       }
        ParameterType b2() { return _b2;       }

        ParameterType _a0_a1_a2;
        ParameterType _b1, _b2;
    };

    template <typename ParameterType>
    struct makeParameter
    {
        typedef ParameterStruct<ParameterType> type;
    };

    template < typename BiquadParameterStruct, typename ArgumentType >
    static BiquadParameterStruct designFilter ( ArgumentType cutoff, ArgumentType q )
    {
        using namespace boost::simd;

        typedef typename meta::scalar_of<ArgumentType>::type ScalarArg1;

        typedef typename BiquadParameterStruct::type ParameterType;

        auto KArg = Pi<ScalarArg1>() * cutoff;
        auto K  = nt2::tan( KArg );
        auto K2 = K*K;

        BiquadParameterStruct ret;
        auto denom = (K2*q + K + q);

        ret._a0_a1_a2 = castType<ParameterType>( fast_div( (K2*q), denom ) );

        ret._b1 = castType<ParameterType>( fast_div( 2.f * q * (K2-1.f),  denom) );
        ret._b2 = castType<ParameterType>( fast_div( K2*q - K + q,  denom)       );
        return ret;
    }
};

struct DesignHPF
{
    template <typename ParameterType>
    struct ParameterStruct
    {
        typedef ParameterType type;

        ParameterType a0() { return _a0_a2; }
        ParameterType a1() { return _a1;    }
        ParameterType a2() { return _a0_a2; }
        ParameterType b1() { return _b1;    }
        ParameterType b2() { return _b2;    }

        ParameterType _a0_a2, _a1;
        ParameterType _b1, _b2;
    };

    template <typename ParameterType>
    struct makeParameter
    {
        typedef ParameterStruct<ParameterType> type;
    };

    template < typename BiquadParameterStruct, typename ArgumentType >
    static BiquadParameterStruct designFilter ( ArgumentType cutoff, ArgumentType q )
    {
        using namespace boost::simd;
        auto K  = nt2::tan( boost::simd::Pi<ArgumentType>() * cutoff );
        auto K2 = K*K;

        BiquadParameterStruct ret;
        typedef typename BiquadParameterStruct::type ParameterType;

        auto denom = (K2*q + K + q);
        ret._a0_a2 = castType<ParameterType>( fast_div( q,  denom) );
        ret._a1    = castType<ParameterType>( fast_div( -2.f * q,  denom) );

        ret._b1 = castType<ParameterType>( fast_div( 2.f * q * (K2-1.f), denom) );
        ret._b2 = castType<ParameterType>( fast_div( (K2*q - K + q),  denom) );
        return ret;
    }
};

struct DesignBPF
{
    template <typename ParameterType>
    struct ParameterStruct
    {
        typedef ParameterType type;

        ParameterType a0() { return _a0; }
        auto          a1() { return boost::simd::Zero<type>(); }
        auto          a2() { return a1() - _a0; }

        ParameterType b1() { return _b1;    }
        ParameterType b2() { return _b2;    }

        ParameterType _a0;
        ParameterType _b1, _b2;
    };

    template <typename ParameterType>
    struct makeParameter
    {
        typedef ParameterStruct<ParameterType> type;
    };

    template < typename BiquadParameterStruct, typename ArgumentType >
    static BiquadParameterStruct designFilter ( ArgumentType cutoff, ArgumentType q )
    {
        using namespace boost::simd;
        auto K  = nt2::tan( boost::simd::Pi<ArgumentType>() * cutoff );
        auto K2 = K*K;

        BiquadParameterStruct ret;
        typedef typename BiquadParameterStruct::type ParameterType;

        auto denom = (K2*q + K + q);
        ret._a0 = castType<ParameterType>( fast_div( K,  denom) );

        ret._b1 = castType<ParameterType>( fast_div( 2.f * q * (K2-1.f), denom ) );
        ret._b2 = castType<ParameterType>( fast_div( K2*q - K + q, denom ) );
        return ret;
    }
};

struct DesignBRF
{
    template <typename ParameterType>
    struct ParameterStruct
    {
        typedef ParameterType type;

        ParameterType a0() { return _a0_a2; }
        ParameterType a1() { return _a1_b1; }
        ParameterType a2() { return _a0_a2; }
        ParameterType b1() { return _a1_b1; }
        ParameterType b2() { return _b2;    }

        ParameterType _a0_a2, _a1_b1;
        ParameterType _b2;
    };

    template <typename ParameterType>
    struct makeParameter
    {
        typedef ParameterStruct<ParameterType> type;
    };

    template < typename BiquadParameterStruct, typename ArgumentType >
    static BiquadParameterStruct designFilter ( ArgumentType cutoff, ArgumentType q )
    {
        using namespace boost::simd;
        auto K  = nt2::tan( boost::simd::Pi<ArgumentType>() * cutoff );
        auto K2 = K*K;

        BiquadParameterStruct ret;
        typedef typename BiquadParameterStruct::type ParameterType;

        auto denom = (K2*q + K + q);
        ret._a0_a2 = castType<ParameterType>( fast_div( q * (1.f + K2), denom ) );
        ret._a1_b1 = castType<ParameterType>( fast_div( 2.f * q * (K2 - 1.f), denom ) );

        ret._b2 = castType<ParameterType>( fast_div( K2*q - K + q,  denom     ) );
        return ret;
    }
};

struct DesignAPF
{
    template <typename ParameterType>
    struct ParameterStruct
    {
        typedef ParameterType type;

        ParameterType a0() { return _a0_b2; }
        ParameterType a1() { return _a1_b1; }
        auto          a2() { return boost::simd::One<type>(); }
        ParameterType b1() { return _a1_b1; }
        ParameterType b2() { return _a0_b2; }

        ParameterType _a0_b2, _a1_b1;
    };

    template <typename ParameterType>
    struct makeParameter
    {
        typedef ParameterStruct<ParameterType> type;
    };

    template < typename BiquadParameterStruct, typename ArgumentType >
    static BiquadParameterStruct designFilter ( ArgumentType cutoff, ArgumentType q )
    {
        using namespace boost::simd;
        auto K  = nt2::tan( boost::simd::Pi<ArgumentType>() * cutoff );
        auto K2 = K*K;

        BiquadParameterStruct ret;
        typedef typename BiquadParameterStruct::type ParameterType;

        auto denom = (K2*q + K + q);
        ret._a0_b2 = castType<ParameterType>( fast_div( K2*q - K + q, denom ) );
        ret._a1_b1 = castType<ParameterType>( fast_div( 2.f * q * (K2 - 1.f), denom ) );
        return ret;
    }
};

PluginLoad(NovaFilters)
{
    ft = inTable;
    using namespace nova;

    nova::registerUnit<nova::NovaIntegrator<1>> ( ft, "NovaIntegrator"  );
    nova::registerUnit<nova::NovaIntegrator<2>> ( ft, "NovaIntegrator2" );
    nova::registerUnit<nova::NovaIntegrator<4>> ( ft, "NovaIntegrator4" );

    nova::registerUnit<NovaBiquadBase<DesignLPF, 1>> (ft, "NovaLowPass"  );
    nova::registerUnit<NovaBiquadBase<DesignLPF, 2>> (ft, "NovaLowPass2" );
    nova::registerUnit<NovaBiquadBase<DesignLPF, 4>> (ft, "NovaLowPass4" );

    nova::registerUnit<NovaBiquadBase<DesignHPF, 1>> (ft, "NovaHighPass"  );
    nova::registerUnit<NovaBiquadBase<DesignHPF, 2>> (ft, "NovaHighPass2" );
    nova::registerUnit<NovaBiquadBase<DesignHPF, 4>> (ft, "NovaHighPass4" );

    nova::registerUnit<NovaBiquadBase<DesignBPF, 1>> (ft, "NovaBandPass"  );
    nova::registerUnit<NovaBiquadBase<DesignBPF, 2>> (ft, "NovaBandPass2" );
    nova::registerUnit<NovaBiquadBase<DesignBPF, 4>> (ft, "NovaBandPass4" );

    nova::registerUnit<NovaBiquadBase<DesignBRF, 1>> (ft, "NovaBandReject"  );
    nova::registerUnit<NovaBiquadBase<DesignBRF, 2>> (ft, "NovaBandReject2" );
    nova::registerUnit<NovaBiquadBase<DesignBRF, 4>> (ft, "NovaBandReject4" );

    nova::registerUnit<NovaBiquadBase<DesignAPF, 1>> (ft, "NovaAllPass"  );
    nova::registerUnit<NovaBiquadBase<DesignAPF, 2>> (ft, "NovaAllPass2" );
    nova::registerUnit<NovaBiquadBase<DesignAPF, 4>> (ft, "NovaAllPass4" );
}
