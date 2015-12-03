/*
 *
 *    Copyright (C) 2013 Tim Blechmann
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

#include "SC_PlugIn.hpp"

#ifndef NOVAUGENSCOMMON_HPP
#define NOVAUGENSCOMMON_HPP

#include <algorithm>
#include <type_traits>
#include <boost/range/irange.hpp>

#include <boost/simd/sdk/meta/cardinal_of.hpp>

#include <producer_consumer_functors.hpp>

namespace nova {

// dispatch tags

struct control_signature_k {};
struct control_signature_i : control_signature_k {};
struct control_signature_1 : control_signature_i {};
struct control_signature_a : control_signature_k {};


struct control_signature_kk {};
struct control_signature_ki : control_signature_kk {};
struct control_signature_ka : control_signature_kk {};

struct control_signature_ii : control_signature_kk {};
struct control_signature_11 : control_signature_ii {};
struct control_signature_ik : control_signature_kk {};
struct control_signature_ia : control_signature_ka {};

struct control_signature_ak : control_signature_kk {};
struct control_signature_ai : control_signature_ak {};
struct control_signature_aa : control_signature_ak {};

namespace meta {

template <typename ControlSignature>
struct InitControlSignature
{
    typedef typename std::is_base_of< control_signature_k,  ControlSignature >::type is_unary_control;
    typedef typename std::is_base_of< control_signature_kk, ControlSignature >::type is_binary_control;

    typedef typename boost::mpl::if_< is_unary_control,
                                      control_signature_1,
                                      typename boost::mpl::if_< is_binary_control, control_signature_11, void >::type
                                    >::type type;
};

}

struct NovaUnit:
    public SCUnit
{
private:
    struct DSPContext {
        explicit DSPContext( SCUnit * unit ):
            _rate(unit->mRate)
        {}

        auto sampleRate() const { return _rate->mSampleRate; }
        auto sampleDur()  const { return _rate->mSampleDur;  }

        struct Rate * __restrict__ _rate;
    };

public:
    auto makeDSPContext() { return DSPContext(this); }

    // input functors
    template< typename ResultType, int input >
    inline auto makeScalarInput() const
    {
        ResultType inputSignal = boost::simd::splat<ResultType>( SCUnit::in0( input ) );
        return [=] { return inputSignal; };
    }


    // rate checks
    int inRate(size_t index) const
    {
        return SCUnit::inRate(index);
    }

    int inRate(size_t beginIndex, size_t endIndex) const
    {
        assert( beginIndex <= endIndex );
        if (isScalarRate(beginIndex, endIndex))
            return calc_ScalarRate;
        if (isBufRate(beginIndex, endIndex))
            return calc_BufRate;

        return calc_FullRate;
    }

    bool isScalarRate(size_t beginIndex, size_t endIndex) const
    {
        auto range = boost::irange(beginIndex, endIndex);

        return std::all_of( range.begin(), range.end(), [this] (auto index) {
            return this->SCUnit::inRate(index) == calc_ScalarRate;
        });
    }

    bool isBufRate(size_t beginIndex, size_t endIndex) const
    {
        auto range = boost::irange(beginIndex, endIndex);

        return std::all_of( range.begin(), range.end(), [this] (auto index) {
            return this->SCUnit::inRate(index) <= calc_BufRate;
        });
    }

    /// calculate slope value
    template <typename FloatTypeA, typename FloatTypeB>
    auto calcSlope(FloatTypeA next, FloatTypeB current) const
    {
        const Unit * unit = this;
        return ((next - current) * unit->mRate->mSlopeFactor);
    }

    template <typename FloatTypeA, typename FloatTypeB>
    auto makeRamp(FloatTypeA next, FloatTypeB current) const
    {
        return nova::parameter::makeRamp( current, next, mRate->mSlopeFactor );
    }


	int numInputs() const
	{
		return mNumInputs;
	}

	int numOutputs() const
	{
		return mNumOutputs;
	}

    enum class ArgumentSignature {
        unknownInputSignature,

        i,
        k,
        a,

        ii,
        ik,
        ia,
        ki,
        kk,
        ka,
        ai,
        ak,
        aa,
    };

    ArgumentSignature argumentSignature( int index ) const
    {
        if( isScalarRateIn( index )  )
            return ArgumentSignature::i;
        if( isControlRateIn( index ) )
            return ArgumentSignature::k;
        if( isAudioRateIn( index )   )
            return ArgumentSignature::a;
        return ArgumentSignature::unknownInputSignature;
    }

    ArgumentSignature combineArgumentSignatures( ArgumentSignature s1, ArgumentSignature s2 ) const
    {
#define combineSignatures( left, right )                                      \
        if( s1 == ArgumentSignature::left && s2 == ArgumentSignature::right ) \
            return ArgumentSignature::left##right

        combineSignatures( i, i );
        combineSignatures( i, k );
        combineSignatures( i, a );

        combineSignatures( k, i );
        combineSignatures( k, k );
        combineSignatures( k, a );

        combineSignatures( a, i );
        combineSignatures( a, k );
        combineSignatures( a, a );

        assert(false);
        return ArgumentSignature::unknownInputSignature;
#undef combineArgumentSignatures
    }

    // expected signature:
    //
    // template <typename VectorType, typename ControlSignature>
    // void run(int inNumSamples)
    //
    template < typename DerivedUnit, typename DispatchTag, typename VectorType >
    void setVectorCalcFunction()
    {
        static const size_t vectorSize = boost::simd::meta::cardinal_of<VectorType>::value;

        typedef typename meta::InitControlSignature<DispatchTag>::type InitDispatchTag;

        if (boost::simd::is_aligned( bufferSize(), vectorSize ) )
            set_vector_calc_function< DerivedUnit,
                                     &DerivedUnit::template run< VectorType, DispatchTag >,
                                     &DerivedUnit::template run< float     , InitDispatchTag > >();
        else if( bufferSize() == 1)
            set_calc_function< DerivedUnit,
                              &DerivedUnit::template run< float     , InitDispatchTag > >();
        else
            set_vector_calc_function< DerivedUnit,
                                     &DerivedUnit::template run< float , DispatchTag >,
                                     &DerivedUnit::template run< float , InitDispatchTag > >();
    }

    template < typename DerivedUnit, typename VectorType >
    void setVectorCalcFunction( int controlInputIndex )
    {
        ArgumentSignature signature = argumentSignature( controlInputIndex );

        switch( signature ) {
        case ArgumentSignature::i: setVectorCalcFunction<DerivedUnit, control_signature_i, VectorType>(); return;
        case ArgumentSignature::k: setVectorCalcFunction<DerivedUnit, control_signature_k, VectorType>(); return;
        case ArgumentSignature::a: setVectorCalcFunction<DerivedUnit, control_signature_a, VectorType>(); return;

        default:
            assert(false);
        }
    }

    template < typename DerivedUnit, typename VectorType >
    void setVectorCalcFunction( int controlInputIndex1, int controlInputIndex2 )
    {
        ArgumentSignature signature = combineArgumentSignatures( argumentSignature( controlInputIndex1 ),
                                                                 argumentSignature( controlInputIndex2 ) );


        switch( signature ) {
        case ArgumentSignature::ii: setVectorCalcFunction<DerivedUnit, control_signature_ii, VectorType>(); return;
        case ArgumentSignature::ik: setVectorCalcFunction<DerivedUnit, control_signature_ik, VectorType>(); return;
        case ArgumentSignature::ia: setVectorCalcFunction<DerivedUnit, control_signature_ia, VectorType>(); return;

        case ArgumentSignature::ki: setVectorCalcFunction<DerivedUnit, control_signature_ki, VectorType>(); return;
        case ArgumentSignature::kk: setVectorCalcFunction<DerivedUnit, control_signature_kk, VectorType>(); return;
        case ArgumentSignature::ka: setVectorCalcFunction<DerivedUnit, control_signature_ka, VectorType>(); return;

        case ArgumentSignature::ai: setVectorCalcFunction<DerivedUnit, control_signature_ai, VectorType>(); return;
        case ArgumentSignature::ak: setVectorCalcFunction<DerivedUnit, control_signature_ak, VectorType>(); return;
        case ArgumentSignature::aa: setVectorCalcFunction<DerivedUnit, control_signature_aa, VectorType>(); return;

        default:
            assert(false);
        }
    }


    // expected signature:
    //
    // template <typename ControlSignature>
    // void run(int inNumSamples)
    //
    template < typename DerivedUnit, typename DispatchTag >
    void setCalcFunction()
    {
        typedef typename meta::InitControlSignature<DispatchTag>::type InitDispatchTag;

        if( bufferSize() == 1)
            set_calc_function< DerivedUnit,
                              &DerivedUnit::template run< InitDispatchTag > >();
        else
            set_vector_calc_function< DerivedUnit,
                                     &DerivedUnit::template run< DispatchTag >,
                                     &DerivedUnit::template run< InitDispatchTag > >();
    }

    template < typename DerivedUnit >
    void setCalcFunction( int controlInputIndex )
    {
        ArgumentSignature signature = argumentSignature( controlInputIndex );

        switch( signature ) {
        case ArgumentSignature::i: setCalcFunction<DerivedUnit, control_signature_i>(); return;
        case ArgumentSignature::k: setCalcFunction<DerivedUnit, control_signature_k>(); return;
        case ArgumentSignature::a: setCalcFunction<DerivedUnit, control_signature_a>(); return;

        default:
            assert(false);
        }
    }

    template < typename DerivedUnit >
    void setCalcFunction( int controlInputIndex1, int controlInputIndex2 )
    {
        ArgumentSignature signature = combineArgumentSignatures( argumentSignature( controlInputIndex1 ),
                                                                 argumentSignature( controlInputIndex2 ) );


        switch( signature ) {
        case ArgumentSignature::ii: setCalcFunction<DerivedUnit, control_signature_ii>(); return;
        case ArgumentSignature::ik: setCalcFunction<DerivedUnit, control_signature_ik>(); return;
        case ArgumentSignature::ia: setCalcFunction<DerivedUnit, control_signature_ia>(); return;

        case ArgumentSignature::ki: setCalcFunction<DerivedUnit, control_signature_ki>(); return;
        case ArgumentSignature::kk: setCalcFunction<DerivedUnit, control_signature_kk>(); return;
        case ArgumentSignature::ka: setCalcFunction<DerivedUnit, control_signature_ka>(); return;

        case ArgumentSignature::ai: setCalcFunction<DerivedUnit, control_signature_ai>(); return;
        case ArgumentSignature::ak: setCalcFunction<DerivedUnit, control_signature_ak>(); return;
        case ArgumentSignature::aa: setCalcFunction<DerivedUnit, control_signature_aa>(); return;

        default:
            assert(false);
        }
    }
};

} // namespace nova

#define NovaDefineUnit(name) \
    do {                                                            \
        if( std::is_trivially_destructible<name>::value ) {         \
            DefineSimpleUnit(name);                                 \
        } else  {                                                   \
            DefineDtorUnit(name);                                   \
        }                                                           \
    } while(0)

#endif // NOVAUGENSCOMMON_HPP
