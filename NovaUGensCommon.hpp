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
struct control_signature_k1 : control_signature_kk {};

struct control_signature_ii : control_signature_kk {};
struct control_signature_11 : control_signature_ii {};
struct control_signature_ik : control_signature_kk {};
struct control_signature_ia : control_signature_ka {};

struct control_signature_ak : control_signature_kk {};
struct control_signature_ai : control_signature_ak {};
struct control_signature_aa : control_signature_ak {};


struct control_signature_kkk {};
struct control_signature_kki : control_signature_kkk {};
struct control_signature_kka : control_signature_kkk {};
struct control_signature_kik : control_signature_kkk {};
struct control_signature_kia : control_signature_kik {};
struct control_signature_kii : control_signature_kik {};
struct control_signature_kak : control_signature_kkk {};
struct control_signature_kai : control_signature_kak {};
struct control_signature_kaa : control_signature_kak {};

struct control_signature_akk : control_signature_kkk {};
struct control_signature_aki : control_signature_akk {};
struct control_signature_aka : control_signature_akk {};
struct control_signature_aik : control_signature_akk {};
struct control_signature_aii : control_signature_aik {};
struct control_signature_aia : control_signature_aka {};
struct control_signature_aak : control_signature_akk {};
struct control_signature_aai : control_signature_aak {};
struct control_signature_aaa : control_signature_aak {};

struct control_signature_ikk : control_signature_kkk {};
struct control_signature_iki : control_signature_ikk {};
struct control_signature_ika : control_signature_ikk {};
struct control_signature_iik : control_signature_ikk {};
struct control_signature_iii : control_signature_iik {};
struct control_signature_iia : control_signature_ika {};
struct control_signature_iak : control_signature_ikk {};
struct control_signature_iai : control_signature_iak {};
struct control_signature_iaa : control_signature_iak {};

struct control_signature_111 : control_signature_iii {};


namespace meta {

template <typename ControlSignature>
struct InitControlSignature
{
    typedef typename std::is_base_of< control_signature_k,   ControlSignature >::type is_unary_control;
    typedef typename std::is_base_of< control_signature_kk,  ControlSignature >::type is_binary_control;
    typedef typename std::is_base_of< control_signature_kkk, ControlSignature >::type is_ternary_control;

    typedef typename boost::mpl::if_< is_unary_control,
                                      control_signature_1,
                                      typename boost::mpl::if_< is_binary_control,
                                                                control_signature_11,
                                                                typename boost::mpl::if_< is_ternary_control, control_signature_111,void >::type
                                      >::type
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


        aii,
        aik,
        aia,
        aki,
        akk,
        aka,
        aai,
        aak,
        aaa,

        kii,
        kik,
        kia,
        kki,
        kkk,
        kka,
        kai,
        kak,
        kaa,

        iii,
        iik,
        iia,
        iki,
        ikk,
        ika,
        iai,
        iak,
        iaa,
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
#undef combineSignatures
    }

    ArgumentSignature combineArgumentSignatures( ArgumentSignature s1, ArgumentSignature s2, ArgumentSignature s3 ) const
    {
#define combineSignatures( S1, S2, S3 )   \
        if( s1 == ArgumentSignature::S1   \
         && s2 == ArgumentSignature::S2   \
         && s3 == ArgumentSignature::S3 ) \
            return ArgumentSignature::S1##S2##S3

        combineSignatures( a, i, i );
        combineSignatures( a, i, k );
        combineSignatures( a, i, a );

        combineSignatures( a, k, i );
        combineSignatures( a, k, k );
        combineSignatures( a, k, a );

        combineSignatures( a, a, i );
        combineSignatures( a, a, k );
        combineSignatures( a, a, a );


        combineSignatures( k, i, i );
        combineSignatures( k, i, k );
        combineSignatures( k, i, a );

        combineSignatures( k, k, i );
        combineSignatures( k, k, k );
        combineSignatures( k, k, a );

        combineSignatures( k, a, i );
        combineSignatures( k, a, k );
        combineSignatures( k, a, a );


        combineSignatures( i, i, i );
        combineSignatures( i, i, k );
        combineSignatures( i, i, a );

        combineSignatures( i, k, i );
        combineSignatures( i, k, k );
        combineSignatures( i, k, a );

        combineSignatures( i, a, i );
        combineSignatures( i, a, k );
        combineSignatures( i, a, a );

        assert(false);
        return ArgumentSignature::unknownInputSignature;
#undef combineSignatures
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

    template < typename DerivedUnit, typename VectorType >
    void setVectorCalcFunction( int controlInputIndex1, int controlInputIndex2, int controlInputIndex3 )
    {
        ArgumentSignature signature = combineArgumentSignatures( argumentSignature( controlInputIndex1 ),
                                                                 argumentSignature( controlInputIndex2 ),
                                                                 argumentSignature( controlInputIndex3 ));


        switch( signature ) {
        case ArgumentSignature::aii: setVectorCalcFunction<DerivedUnit, control_signature_aii, VectorType>(); return;
        case ArgumentSignature::aik: setVectorCalcFunction<DerivedUnit, control_signature_aik, VectorType>(); return;
        case ArgumentSignature::aia: setVectorCalcFunction<DerivedUnit, control_signature_aia, VectorType>(); return;

        case ArgumentSignature::aki: setVectorCalcFunction<DerivedUnit, control_signature_aki, VectorType>(); return;
        case ArgumentSignature::akk: setVectorCalcFunction<DerivedUnit, control_signature_akk, VectorType>(); return;
        case ArgumentSignature::aka: setVectorCalcFunction<DerivedUnit, control_signature_aka, VectorType>(); return;

        case ArgumentSignature::aai: setVectorCalcFunction<DerivedUnit, control_signature_aai, VectorType>(); return;
        case ArgumentSignature::aak: setVectorCalcFunction<DerivedUnit, control_signature_aak, VectorType>(); return;
        case ArgumentSignature::aaa: setVectorCalcFunction<DerivedUnit, control_signature_aaa, VectorType>(); return;


        case ArgumentSignature::kii: setVectorCalcFunction<DerivedUnit, control_signature_kii, VectorType>(); return;
        case ArgumentSignature::kik: setVectorCalcFunction<DerivedUnit, control_signature_kik, VectorType>(); return;
        case ArgumentSignature::kia: setVectorCalcFunction<DerivedUnit, control_signature_kia, VectorType>(); return;

        case ArgumentSignature::kki: setVectorCalcFunction<DerivedUnit, control_signature_kki, VectorType>(); return;
        case ArgumentSignature::kkk: setVectorCalcFunction<DerivedUnit, control_signature_kkk, VectorType>(); return;
        case ArgumentSignature::kka: setVectorCalcFunction<DerivedUnit, control_signature_kka, VectorType>(); return;

        case ArgumentSignature::kai: setVectorCalcFunction<DerivedUnit, control_signature_kai, VectorType>(); return;
        case ArgumentSignature::kak: setVectorCalcFunction<DerivedUnit, control_signature_kak, VectorType>(); return;
        case ArgumentSignature::kaa: setVectorCalcFunction<DerivedUnit, control_signature_kaa, VectorType>(); return;


        case ArgumentSignature::iii: setVectorCalcFunction<DerivedUnit, control_signature_iii, VectorType>(); return;
        case ArgumentSignature::iik: setVectorCalcFunction<DerivedUnit, control_signature_iik, VectorType>(); return;
        case ArgumentSignature::iia: setVectorCalcFunction<DerivedUnit, control_signature_iia, VectorType>(); return;

        case ArgumentSignature::iki: setVectorCalcFunction<DerivedUnit, control_signature_iki, VectorType>(); return;
        case ArgumentSignature::ikk: setVectorCalcFunction<DerivedUnit, control_signature_ikk, VectorType>(); return;
        case ArgumentSignature::ika: setVectorCalcFunction<DerivedUnit, control_signature_ika, VectorType>(); return;

        case ArgumentSignature::iai: setVectorCalcFunction<DerivedUnit, control_signature_iai, VectorType>(); return;
        case ArgumentSignature::iak: setVectorCalcFunction<DerivedUnit, control_signature_iak, VectorType>(); return;
        case ArgumentSignature::iaa: setVectorCalcFunction<DerivedUnit, control_signature_iaa, VectorType>(); return;

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
                                                                 argumentSignature( controlInputIndex2 ));


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
    template < typename DerivedUnit >
    void setCalcFunction( int controlInputIndex1, int controlInputIndex2, int controlInputIndex3 )
    {
        ArgumentSignature signature = combineArgumentSignatures( argumentSignature( controlInputIndex1 ),
                                                                 argumentSignature( controlInputIndex2 ),
                                                                 argumentSignature( controlInputIndex3 ));


        switch( signature ) {
        case ArgumentSignature::aii: setCalcFunction<DerivedUnit, control_signature_aii>(); return;
        case ArgumentSignature::aik: setCalcFunction<DerivedUnit, control_signature_aik>(); return;
        case ArgumentSignature::aia: setCalcFunction<DerivedUnit, control_signature_aia>(); return;

        case ArgumentSignature::aki: setCalcFunction<DerivedUnit, control_signature_aki>(); return;
        case ArgumentSignature::akk: setCalcFunction<DerivedUnit, control_signature_akk>(); return;
        case ArgumentSignature::aka: setCalcFunction<DerivedUnit, control_signature_aka>(); return;

        case ArgumentSignature::aai: setCalcFunction<DerivedUnit, control_signature_aai>(); return;
        case ArgumentSignature::aak: setCalcFunction<DerivedUnit, control_signature_aak>(); return;
        case ArgumentSignature::aaa: setCalcFunction<DerivedUnit, control_signature_aaa>(); return;


        case ArgumentSignature::kii: setCalcFunction<DerivedUnit, control_signature_kii>(); return;
        case ArgumentSignature::kik: setCalcFunction<DerivedUnit, control_signature_kik>(); return;
        case ArgumentSignature::kia: setCalcFunction<DerivedUnit, control_signature_kia>(); return;

        case ArgumentSignature::kki: setCalcFunction<DerivedUnit, control_signature_kki>(); return;
        case ArgumentSignature::kkk: setCalcFunction<DerivedUnit, control_signature_kkk>(); return;
        case ArgumentSignature::kka: setCalcFunction<DerivedUnit, control_signature_kka>(); return;

        case ArgumentSignature::kai: setCalcFunction<DerivedUnit, control_signature_kai>(); return;
        case ArgumentSignature::kak: setCalcFunction<DerivedUnit, control_signature_kak>(); return;
        case ArgumentSignature::kaa: setCalcFunction<DerivedUnit, control_signature_kaa>(); return;


        case ArgumentSignature::iii: setCalcFunction<DerivedUnit, control_signature_iii>(); return;
        case ArgumentSignature::iik: setCalcFunction<DerivedUnit, control_signature_iik>(); return;
        case ArgumentSignature::iia: setCalcFunction<DerivedUnit, control_signature_iia>(); return;

        case ArgumentSignature::iki: setCalcFunction<DerivedUnit, control_signature_iki>(); return;
        case ArgumentSignature::ikk: setCalcFunction<DerivedUnit, control_signature_ikk>(); return;
        case ArgumentSignature::ika: setCalcFunction<DerivedUnit, control_signature_ika>(); return;

        case ArgumentSignature::iai: setCalcFunction<DerivedUnit, control_signature_iai>(); return;
        case ArgumentSignature::iak: setCalcFunction<DerivedUnit, control_signature_iak>(); return;
        case ArgumentSignature::iaa: setCalcFunction<DerivedUnit, control_signature_iaa>(); return;

        default:
            assert(false);
        }
    }
};

namespace detail {

template <class UGenClass>
void constructClass( Unit * unit ) { new( static_cast<UGenClass*>(unit) ) UGenClass(); }
template <class UGenClass>
void destroyClass( Unit * unit )   { static_cast<UGenClass*>(unit)->~UGenClass();    }

}

template <class Unit>
void registerUnit( InterfaceTable * ft, const char * name )
{
    UnitCtorFunc ctor = detail::constructClass<Unit>;
    UnitDtorFunc dtor = std::is_trivially_destructible<Unit>::value ? nullptr
                                                                    : detail::destroyClass<Unit>;

    (*ft->fDefineUnit)( name, sizeof(Unit), ctor, dtor, 0 );
}

} // namespace nova


#endif // NOVAUGENSCOMMON_HPP
