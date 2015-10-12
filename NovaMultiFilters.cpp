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
typedef nova::NovaIntegrator<8> NovaIntegrator8;

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


template <typename FilterDesigner, size_t Size, bool ScalarArguments = true>
struct NovaBiquadBase:
	public NovaUnit,
    public nova::multichannel::SignalInput< NovaBiquadBase<FilterDesigner, Size, ScalarArguments>, 0, Size >,
    public nova::multichannel::OutputSink<  NovaBiquadBase<FilterDesigner, Size, ScalarArguments>, 0, Size >,

    // Controls:
    public nova::multichannel::ControlInput<  NovaBiquadBase<FilterDesigner, Size, ScalarArguments>, Size,                                (ScalarArguments ? 1 : Size) >,
    public nova::multichannel::ControlInput<  NovaBiquadBase<FilterDesigner, Size, ScalarArguments>, Size + (ScalarArguments ? 1 : Size), (ScalarArguments ? 1 : Size) >
{
    static const size_t ParameterSize = ScalarArguments ? 1 : Size;
    static const int FreqInputIndex = Size;
    static const int QInputIndex    = Size + ParameterSize;

	typedef typename nova::as_pack<float,  Size>::type vFloat;
	typedef typename nova::as_pack<double, Size>::type vDouble;

    typedef typename boost::mpl::if_c<ScalarArguments, double, vDouble>::type ParameterType;
    typedef typename boost::mpl::if_c<ScalarArguments, float,  vFloat>::type HostParameterType;

    typedef typename FilterDesigner::template makeParameter<ParameterType>::type BiquadParameterStruct;
    typedef nova::Biquad<vDouble, BiquadParameterStruct> Filter;

    using SignalInput = nova::multichannel::SignalInput< NovaBiquadBase<FilterDesigner, Size, ScalarArguments>, 0, Size >;
    using OutputSink  = nova::multichannel::OutputSink<  NovaBiquadBase<FilterDesigner, Size, ScalarArguments>, 0, Size >;

    using FreqInput = nova::multichannel::ControlInput<  NovaBiquadBase<FilterDesigner, Size, ScalarArguments>, Size,                                ScalarArguments ? 1 : Size >;
    using QInput    = nova::multichannel::ControlInput<  NovaBiquadBase<FilterDesigner, Size, ScalarArguments>, Size + (ScalarArguments ? 1 : Size), ScalarArguments ? 1 : Size >;


    NovaBiquadBase()
    {
        initFilter();

		if (isScalarRate (FreqInputIndex, FreqInputIndex + ParameterSize) )
            set_vector_calc_function<NovaBiquadBase, &NovaBiquadBase::next_i, &NovaBiquadBase::next_1>();
        else
            set_vector_calc_function<NovaBiquadBase, &NovaBiquadBase::next_k, &NovaBiquadBase::next_1>();
    }

    void initFilter()
    {
        auto newFreq = FreqInput::readInput();
        auto newQ    = QInput::readInput();

        HostParameterType normalizedFreq = newFreq * (float)sampleDur();
        auto params = FilterDesigner::template designFilter<BiquadParameterStruct>( normalizedFreq, newQ );

        storeFilterParameters( params );
    }

    void next_1(int)
    {
        auto inFn  = SignalInput::template makeInputSignal<vDouble>();
        auto outFn = OutputSink ::template makeSink<vDouble>();

        _filter.run(inFn, outFn, 1);
    }

    void next_i(int)
    {
        auto inFn  = SignalInput::template makeInputSignal<vDouble>();
        auto outFn = OutputSink ::template makeSink<vDouble>();

        _filter.run_unrolled(inFn, outFn, mRate->mFilterLoops, mRate->mFilterRemain);
    }

    void next_k(int inNumSamples)
    {
        auto newFreq = FreqInput::readInput();
        auto newQ    = QInput   ::readInput();

        if ( !FreqInput::changed() && !QInput::changed() ) {
            next_i(inNumSamples);
            return;
        }

        HostParameterType normalizedFreq = newFreq * (float)sampleDur();
        auto newParameters = FilterDesigner::template designFilter<BiquadParameterStruct>( normalizedFreq, newQ );

        ParameterType a0Slope = calcSlope( newParameters.a0(), _filter._parameters.a0() );
        ParameterType a1Slope = calcSlope( newParameters.a1(), _filter._parameters.a1() );
        ParameterType a2Slope = calcSlope( newParameters.a2(), _filter._parameters.a2() );
        ParameterType b1Slope = calcSlope( newParameters.b1(), _filter._parameters.b1() );
        ParameterType b2Slope = calcSlope( newParameters.b2(), _filter._parameters.b2() );

        auto inFn  = SignalInput::template makeInputSignal<vDouble>();
        auto outFn = OutputSink ::template makeSink<vDouble>();

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

using nova::toDouble;

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

        ret._a0_a1_a2 = toDouble<ParameterType>( fast_div( (K2*q), denom ) );

        ret._b1 = toDouble<ParameterType>( fast_div( 2.f * q * (K2-1.f),  denom) );
        ret._b2 = toDouble<ParameterType>( fast_div( K2*q - K + q,  denom)       );
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
        ret._a0_a2 = toDouble<ParameterType>( fast_div( q,  denom) );
        ret._a1    = toDouble<ParameterType>( fast_div( -2.f * q,  denom) );

        ret._b1 = toDouble<ParameterType>( fast_div( 2.f * q * (K2-1.f), denom) );
        ret._b2 = toDouble<ParameterType>( fast_div( (K2*q - K + q),  denom) );
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
        ret._a0 = toDouble<ParameterType>( fast_div( K,  denom) );

        ret._b1 = toDouble<ParameterType>( fast_div( 2.f * q * (K2-1.f), denom ) );
        ret._b2 = toDouble<ParameterType>( fast_div( K2*q - K + q, denom ) );
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
        ret._a0_a2 = toDouble<ParameterType>( fast_div( q * (1.f + K2), denom ) );
        ret._a1_b1 = toDouble<ParameterType>( fast_div( 2.f * q * (K2 - 1.f), denom ) );

        ret._b2 = toDouble<ParameterType>( fast_div( K2*q - K + q,  denom     ) );
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
        ret._a0_b2 = toDouble<ParameterType>( fast_div( K2*q - K + q, denom ) );
        ret._a1_b1 = toDouble<ParameterType>( fast_div( 2.f * q * (K2 - 1.f), denom ) );
        return ret;
    }
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2nd order, 2-channel

#define DEFINED_2ND_ORDER_FILTERS(Chans)            \
    struct NovaLowPass##Chans:                      \
        NovaBiquadBase<DesignLPF, Chans>            \
    {                                               \
        NovaLowPass##Chans() {}                     \
    };                                              \
                                                    \
    struct NovaHighPass##Chans:                     \
        NovaBiquadBase<DesignHPF, Chans>            \
    {                                               \
        NovaHighPass##Chans() {}                    \
    };                                              \
                                                    \
    struct NovaBandPass##Chans:                     \
        NovaBiquadBase<DesignBPF, Chans>            \
    {                                               \
        NovaBandPass##Chans() {}                    \
    };                                              \
                                                    \
    struct NovaBandReject##Chans:                   \
        NovaBiquadBase<DesignBRF, Chans>            \
    {                                               \
        NovaBandReject##Chans() {}                  \
    };                                              \
                                                    \
    struct NovaAllPass##Chans:                      \
        NovaBiquadBase<DesignAPF, Chans>            \
    {                                               \
        NovaAllPass##Chans() {}                     \
    };                                              \
                                                    \
    struct NovaLowPass##Chans##_##Chans:            \
        NovaBiquadBase<DesignLPF, Chans, false>     \
    {                                               \
        NovaLowPass##Chans##_##Chans() {}           \
    };                                              \
                                                    \
    struct NovaHighPass##Chans##_##Chans:           \
        NovaBiquadBase<DesignHPF, Chans, false>     \
    {                                               \
        NovaHighPass##Chans##_##Chans() {}          \
    };                                              \
                                                    \
    struct NovaBandPass##Chans##_##Chans:           \
        NovaBiquadBase<DesignBPF, Chans, false>     \
    {                                               \
        NovaBandPass##Chans##_##Chans() {}          \
    };                                              \
                                                    \
    struct NovaBandReject##Chans##_##Chans:         \
        NovaBiquadBase<DesignBRF, Chans, false>     \
    {                                               \
        NovaBandReject##Chans##_##Chans() {}        \
    };                                              \
                                                    \
    struct NovaAllPass##Chans##_##Chans:            \
        NovaBiquadBase<DesignAPF, Chans, false>     \
    {                                               \
        NovaAllPass##Chans##_##Chans() {}           \
    };                                              \
                                                    \
    DEFINE_XTORS(NovaLowPass##Chans)                \
    DEFINE_XTORS(NovaHighPass##Chans)               \
    DEFINE_XTORS(NovaBandPass##Chans)               \
    DEFINE_XTORS(NovaBandReject##Chans)             \
    DEFINE_XTORS(NovaAllPass##Chans)                \
                                                    \
    DEFINE_XTORS(NovaLowPass##Chans##_##Chans)      \
    DEFINE_XTORS(NovaHighPass##Chans##_##Chans)     \
    DEFINE_XTORS(NovaBandPass##Chans##_##Chans)     \
    DEFINE_XTORS(NovaBandReject##Chans##_##Chans)   \
    DEFINE_XTORS(NovaAllPass##Chans##_##Chans)



DEFINED_2ND_ORDER_FILTERS(1)
DEFINED_2ND_ORDER_FILTERS(2)
DEFINED_2ND_ORDER_FILTERS(4)
//DEFINED_2ND_ORDER_FILTERS(8)



DEFINE_XTORS(NovaIntegrator1)
DEFINE_XTORS(NovaIntegrator2)
DEFINE_XTORS(NovaIntegrator4)
//DEFINE_XTORS(NovaIntegrator8)


PluginLoad(NovaFilters)
{
    ft = inTable;

    DefineSimpleUnit(NovaIntegrator1);
    DefineSimpleUnit(NovaIntegrator2);
    DefineSimpleUnit(NovaIntegrator4);
//    DefineSimpleUnit(NovaIntegrator8);

#define REGISTER_BIQUADS(Chans)                         \
    DefineSimpleUnit(NovaLowPass##Chans);               \
    DefineSimpleUnit(NovaHighPass##Chans);              \
    DefineSimpleUnit(NovaBandPass##Chans);              \
    DefineSimpleUnit(NovaBandReject##Chans);            \
    DefineSimpleUnit(NovaAllPass##Chans);               \
    \
    DefineSimpleUnit(NovaLowPass##Chans##_##Chans);     \
    DefineSimpleUnit(NovaHighPass##Chans##_##Chans);    \
    DefineSimpleUnit(NovaBandPass##Chans##_##Chans);    \
    DefineSimpleUnit(NovaBandReject##Chans##_##Chans);  \
    DefineSimpleUnit(NovaAllPass##Chans##_##Chans);

    REGISTER_BIQUADS(1);
    REGISTER_BIQUADS(2);
    REGISTER_BIQUADS(4);
//    REGISTER_BIQUADS(8);
}
