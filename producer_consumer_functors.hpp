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

#ifndef PRODUCER_CONSUMER_FUNCTORS_HPP
#define PRODUCER_CONSUMER_FUNCTORS_HPP

#include "SC_PlugIn.hpp"

#include "dsp/utils.hpp"

#include <array>
#include <tuple>

#include <boost/range/irange.hpp>

#include "boost/simd/include/functions/aligned_load.hpp"
#include "boost/simd/include/functions/aligned_store.hpp"
#include "boost/simd/include/functions/extract.hpp"
#include "boost/simd/include/functions/groups.hpp"
#include "boost/simd/include/functions/insert.hpp"
#include <boost/simd/include/functions/splat.hpp>
#include "boost/simd/include/functions/split.hpp"

#include "boost/simd/include/functions/compare_equal.hpp"


#include "boost/simd/operator/include/functions/plus.hpp"
#include "boost/simd/operator/include/functions/multiplies.hpp"
#include "boost/simd/include/functions/assign.hpp"

#include <boost/dispatch/meta/is_scalar.hpp>
#include "boost/simd/sdk/meta/cardinal_of.hpp"

namespace nova {

namespace detail {

template <int N>
struct packGenerator
{};

template <>
struct packGenerator<1>
{
    template <typename Result, typename Functor>
    static Result generate(Functor && f)
    {
        return Result( f() );
    }
};

template <>
struct packGenerator<2>
{
    template <typename Result, typename Functor>
    static Result generate(Functor && f)
    {
        auto a = f();
        auto b = f();

        return Result(a, b);
    }
};

template <>
struct packGenerator<4>
{
    template <typename Result, typename Functor>
    static Result generate(Functor && f)
    {
        auto a = f();
        auto b = f();
        auto c = f();
        auto d = f();

        return Result(a, b, c, d);
    }
};

template <>
struct packGenerator<8>
{
    template <typename Result, typename Functor>
    static Result generate(Functor && f)
    {
        auto a = f();
        auto b = f();
        auto c = f();
        auto d = f();
        auto e = f();
        auto f_ = f();
        auto g = f();
        auto h = f();

        return Result(a, b, c, d, e, f_, g, h);
    }
};

}

template <typename OutputType, typename Scalar>
auto makeScalarRamp( Scalar base, Scalar slope )
{
    OutputType state = base;

    return [=] () mutable {
        state += slope;
        return state;
    };
}


template <typename OutputType, typename Scalar>
auto makeVectorRamp( Scalar base, Scalar slope )
{
    using namespace boost::simd;

    const size_t size = meta::cardinal_of<OutputType>::value;
    OutputType state = detail::packGenerator<size>::template generate<OutputType>( makeScalarRamp<Scalar>(base, slope) );

    OutputType increment = slope * size;

    return [=] () mutable {
        state += increment;
        return state;
    };
}

template <typename SampleType>
struct Wire
{
    inline SampleType operator() ()
    {
        return _data;
    }

    inline void operator() (SampleType arg)
    {
        _data = arg;
    }

    SampleType _data;
};



namespace detail {

struct Identity {
    template <typename Arg>
    auto operator()(Arg && arg) const { return arg; }

    template <typename Type>
    struct State {
        typedef Type type;
    };
};

template <typename Type, int N>
struct ArithmeticArray
{
    typedef Type value_type;

    ArithmeticArray()                                          = default;
    ArithmeticArray( ArithmeticArray const & rhs )             = default;
    ArithmeticArray & operator=( ArithmeticArray const & rhs ) = default;

    ArithmeticArray( std::initializer_list<Type> const & init )
    {
        std::copy( init.begin(), init.end(), data.begin() );
    }

    template <typename RhsType>
    ArithmeticArray( ArithmeticArray<RhsType, N> const & rhs )
    {
        for( int i : boost::irange(0, N) )
            data[i] = rhs[i];
    }

    Type const & get(size_t n) { return data[n]; }

    ArithmeticArray operator+(ArithmeticArray const & rhs) const
    {
        ArithmeticArray ret;
        for( int i : boost::irange(0, N) )
            ret[i] = data[i] + rhs[i];
        return ret;
    }

    template <typename RhsType>
    ArithmeticArray & operator+=(ArithmeticArray<RhsType, N> const & rhs)
    {
        ArithmeticArray ret;
        for( int i : boost::irange(0, N) )
            data[i] += rhs[i];
        return *this;
    }

    ArithmeticArray operator-(ArithmeticArray const & rhs) const
    {
        ArithmeticArray ret;
        for( int i : boost::irange(0, N) )
            ret[i] = data[i] - rhs[i];
        return ret;
    }

    template <typename ArgType>
    ArithmeticArray operator*(ArgType const & value) const
    {
        ArithmeticArray ret;
        for( int i : boost::irange(0, N) )
            ret[i] = data[i] * value;
        return ret;
    }

    Type & operator[](size_t n)             { return data[n]; }
    Type const & operator[](size_t n) const { return data[n]; }

private:
    std::array<Type, N> data;
};

}


template <typename OutputType, typename Type, int N>
auto makeVectorRamp( detail::ArithmeticArray<Type, N> base, detail::ArithmeticArray<Type, N> const & slope )
{
    using namespace boost::simd;

    typedef typename OutputType::value_type OutputVector;
    const size_t cardinal = meta::cardinal_of<OutputVector>::value;

    OutputType state;
    detail::ArithmeticArray<Type, N> increment = slope * cardinal;

    for( int n : boost::irange(0, N) )
        state[n] = detail::packGenerator<cardinal>::template generate<OutputVector>( makeScalarRamp<Type>( base[n], slope[n] ) );

    return [=] () mutable {
        state += increment;
        return state;
    };
}

template <typename UGenClass, size_t InputIndex, typename InputFunctor = detail::Identity>
struct ScalarInput:
    protected InputFunctor
{
protected:
    UGenClass * asUGen()             { return static_cast<UGenClass*>(this);       }
    const UGenClass * asUGen() const { return static_cast<const UGenClass*>(this); }

public:
    bool audioRate()   const         { return asUGen()->isAudioRateIn(   InputIndex );  }
    bool controlRate() const         { return asUGen()->isControlRateIn( InputIndex );  }
    bool scalarRate()  const         { return asUGen()->isScalarRateIn(  InputIndex );  }

    float slopeFactor() const        { return asUGen()->mRate->mSlopeFactor;            }

    float readRawInput() const
    {
        return asUGen()->in0( InputIndex );
    }

    auto readInput()
    {
        return InputFunctor::operator()( readRawInput() );
    }

    auto readRawAndMappedInput()
    {
        auto rawInputValue = readRawInput();
        auto mappedInput   = InputFunctor::operator()( rawInputValue );
        return std::make_tuple( rawInputValue, mappedInput );
    }

    template< typename OutputType >
    auto makeInputSignal()
    {
        OutputType inputSignal = boost::simd::splat<OutputType>( readInput() );
        return [=] { return inputSignal; };
    }
};


template <typename UGenClass, size_t InputIndex, typename InputFunctor = detail::Identity >
struct SlopedInput:
    ScalarInput< UGenClass, InputIndex, InputFunctor >
{
private:
    struct vector_slope {};
    struct scalar_slope {};

public:
    typedef typename InputFunctor::template State<float>::type ScalarState;

    typedef ScalarInput< UGenClass, InputIndex, InputFunctor > Base;

    SlopedInput():
        mState( Base::readRawInput() )
    {}

    auto readInput()
    {
        return Base::readInput();
    }

    template <typename Type>
    auto makeRampSignal()
    {
        using ScalarType             = typename boost::simd::meta::scalar_of<Type>::type;
        static const size_t cardinal =          boost::simd::meta::cardinal_of<Type>::value;

        ScalarState current          = InputFunctor::operator ()(mState);
        std::tie( mState, mXState )  = Base::readRawAndMappedInput();
        ScalarState slope = (mXState - current) * ScalarType( Base::slopeFactor() );

        typedef typename boost::mpl::if_c< cardinal == 1, scalar_slope, vector_slope >::type slope_tag;

        typedef typename InputFunctor::template State<Type>::type VectorState;


        return makeRamp<VectorState>( current, slope, slope_tag() );
    }

    template <typename ScalarType>
    auto makeMultiRampSignal()
    {
        ScalarState current = InputFunctor::operator ()(mState);
        std::tie( mState, mXState ) = Base::readRawAndMappedInput();
        ScalarState slope = (mXState - current) * Base::slopeFactor();
        return makeScalarRamp<ScalarState>( current, slope );
    }


    template <typename OutputType>
    auto makeScalarInputSignal()
    {
        std::tie( mState, mXState ) = Base::readRawAndMappedInput();
        return [=] { return mXState; };
    }

    bool changed() const
    {
        float next = Base::readRawInput();
        return next != mState;
    }

private:
    template <typename Type>
    auto makeRamp(ScalarState const & base, ScalarState const & slope, vector_slope )
    {
        return makeVectorRamp<Type>( base, slope );
    }

    template <typename Type>
    auto makeRamp(ScalarState const & base, ScalarState const & slope, scalar_slope )
    {
        return makeScalarRamp<Type>( base, slope );
    }

    float mState;
    ScalarState mXState;
};


template <typename UGenClass, size_t InputIndex >
struct SlopedInput< UGenClass, InputIndex, detail::Identity >:
    ScalarInput< UGenClass, InputIndex, detail::Identity >
{
private:
    struct vector_slope {};
    struct scalar_slope {};

public:
    typedef ScalarInput< UGenClass, InputIndex, detail::Identity > Base;

    SlopedInput():
        mState( Base::readRawInput() )
    {}

    auto readInput()
    {
        return Base::readInput();
    }

    template <typename Type>
    auto makeRampSignal()
    {
        static const size_t cardinal =          boost::simd::meta::cardinal_of<Type>::value;

        float current  = mState;
        float newState = Base::readRawInput();
        float slope    = (newState - current) * Base::slopeFactor();

        mState = newState;
        typedef typename boost::mpl::if_c< cardinal == 1, scalar_slope, vector_slope >::type slope_tag;

        return makeRamp<Type>( current, slope, slope_tag() );
    }

    template <typename OutputType>
    auto makeScalarInputSignal()
    {
        mState = Base::readRawInput();
        return [=] { return mState; };
    }

    bool changed() const
    {
        float next = Base::readRawInput();
        return next != mState;
    }

private:
    template <typename Type>
    auto makeRamp(float const & base, float const & slope, vector_slope )
    {
        return makeVectorRamp<Type>( base, slope );
    }

    template <typename Type>
    auto makeRamp(float const & base, float const & slope, scalar_slope )
    {
        return makeScalarRamp<Type>( base, slope );
    }

    float mState;
};



struct vector_tag {};
struct scalar_tag {};

template <typename UGenClass, size_t InputIndex, typename InputFunctor = detail::Identity>
struct SignalInput:
    InputFunctor
{
    const float * inputVector()
    {
        return static_cast<UGenClass*>(this)->in( InputIndex );
    }

    template <typename OutputType>
    auto makeInputSignal()
    {
        typedef typename boost::mpl::if_c< boost::dispatch::meta::is_scalar< OutputType >::value, scalar_tag, vector_tag >::type dispatch_tag;
        return makeInputSignal<OutputType>( dispatch_tag() );
    }

private:
    template <typename SIMDType>
    auto makeInputSignal( vector_tag )
    {
        const float * input = inputVector();
        return [=] () mutable {
            SIMDType in = boost::simd::aligned_load<SIMDType>( input );
            auto ret = InputFunctor::operator()( in );
            input += boost::simd::meta::cardinal_of<SIMDType>::value;
            return ret;
        };
    }

    template <typename ScalarType>
    auto makeInputSignal( scalar_tag )
    {
        const float * input = inputVector();
        return [=] () mutable {
            ScalarType in = *input;
            auto ret = InputFunctor::operator()( in );
            input += 1;
            return ret;
        };
    }
};


template <typename UGenClass, size_t InputIndex, typename InputFunctor = detail::Identity>
struct ControlInput:
    SlopedInput<UGenClass, InputIndex, InputFunctor>,
    SignalInput<UGenClass, InputIndex, InputFunctor>
{
    /* audio rate input */
    template <typename OutputType>
    auto makeAudioInputSignal()
    {
        return SignalInput< UGenClass, InputIndex, InputFunctor>::template makeInputSignal<OutputType>();
    }

    /* control rate input */
    template <typename OutputType>
    auto makeRampSignal()
    {
        return SlopedInput< UGenClass, InputIndex, InputFunctor >::template makeRampSignal<OutputType>();
    }

    // multichannel ramp
    template <typename OutputType>
    auto makeMultiRampSignal()
    {
        return SlopedInput< UGenClass, InputIndex, InputFunctor >::template makeMultiRampSignal<OutputType>();
    }

    /* scalar input */
    template <typename OutputType>
    auto makeScalarInputSignal()
    {
        return SlopedInput< UGenClass, InputIndex, InputFunctor >::template makeScalarInputSignal<OutputType>();
    }
};


template <typename UGenClass, size_t OutputIndex>
struct OutputSink
{
    float * outputVector()
    {
        return static_cast<UGenClass*>(this)->out( OutputIndex );
    }

    template <typename SIMDType,
              typename std::enable_if< !boost::dispatch::meta::is_scalar<SIMDType>::value
                                       >::type * = nullptr>
    auto makeSink()
    {
        float * output = outputVector();
        return [=] (SIMDType arg) mutable {
            boost::simd::aligned_store<SIMDType>( arg, output );
            output += boost::simd::meta::cardinal_of<SIMDType>::value;
        };
    }

    template <typename ScalarType,
              typename std::enable_if< boost::dispatch::meta::is_scalar<ScalarType>::value
                                       >::type * = nullptr>
    auto makeSink()
    {
        float * output = outputVector();
        return [=] (ScalarType arg) mutable {
            *output = arg;
            output += 1;
        };
    }
};

namespace multichannel {

template <typename UGenClass, size_t InputIndex, size_t NumberOfChannels, typename InputFunctor = detail::Identity>
struct ScalarInput:
    private InputFunctor
{
    template< typename OutputType >
    auto readInputs()
    {
        return detail::packGenerator<NumberOfChannels>::template generate<OutputType>( [this, index = 0] () mutable {
             return InputFunctor::operator ()( static_cast<UGenClass*>(this)->in0( InputIndex + index++ ) );
        });
    }


    template< typename OutputType >
    auto makeInputSignal()
    {
        OutputType inputSignal = readInputs<OutputType>();
        return [=] { return inputSignal; };
    }
};

template <typename UGenClass, size_t InputIndex,  size_t NumberOfChannels, typename InputFunctor = detail::Identity >
struct SlopedInput:
    ScalarInput< UGenClass, InputIndex, NumberOfChannels,  InputFunctor >
{
    using Vector = typename nova::as_pack<float, NumberOfChannels>::type;

    SlopedInput():
        mState( readInput() )
    {}

    auto readInput()
    {
        return ScalarInput< UGenClass, InputIndex, NumberOfChannels, InputFunctor >::template readInputs<Vector>();
    }

    auto readAndUpdateInput()
    {
        updateState();
        return mState;
    }

    auto readSlopeFactor()
    {
        return boost::simd::splat<Vector>( static_cast<UGenClass*>(this)->mRate->mSlopeFactor );
    }

    auto currentValue()
    {
        return mState;
    }

    template< typename OutputType >
    auto makeRampSignal( )
    {
        Vector current = mState;
        Vector next  = readInput();
        Vector slope = (next - current) * readSlopeFactor();
        mState = next;
        return makeScalarRamp<OutputType>( current, slope );
    }

    void updateState()
    {
        mState = readInput();
    }

    template< typename OutputType>
    auto makeScalarInputSignal()
    {
        OutputType state { mState };
        return [=] { return state; };
    }

    bool changed()
    {
        auto inputConstant = boost::simd::compare_equal( mState, readInput() );
        return !bool(inputConstant);
    }

    Vector mState;
};


template <typename UGenClass, size_t InputIndex, size_t NumberOfChannels, typename InputFunctor = detail::Identity>
struct SignalInput:
    private InputFunctor
{
    template< typename OutputType >
    auto readInputs( int sampleIndex )
    {
        static_assert( NumberOfChannels == boost::simd::meta::cardinal_of<OutputType>::value, "failed" );

        return detail::packGenerator<NumberOfChannels>::template generate<OutputType>( [=, channelIndex = 0] () mutable {
            auto * input = static_cast<UGenClass*>(this)->SCUnit::in( InputIndex + channelIndex++ );
            return InputFunctor::operator ()( input[sampleIndex] );
        });
    }

    template< typename OutputType >
    auto readInputs()
    {
        return readInputs<OutputType>( 0 );
    }

    template< typename OutputType >
    auto makeInputSignal()
    {
        return [=, sampleIndex = 0] () mutable {
            return readInputs<OutputType>( sampleIndex++ );
        };
    }
};

template <typename UGenClass, size_t InputIndex, size_t NumberOfChannels>
struct SignalInput< UGenClass, InputIndex, NumberOfChannels, detail::Identity>
{
    template< typename OutputType >
    auto readInputs( int sampleIndex )
    {
        static_assert( NumberOfChannels == boost::simd::meta::cardinal_of<OutputType>::value, "failed" );

        return detail::packGenerator<NumberOfChannels>::template generate<OutputType>( [=, channelIndex = 0] () mutable {
            auto * input = static_cast<UGenClass*>(this)->SCUnit::in( InputIndex + channelIndex++ );
            return input[sampleIndex];
        });
    }

    template< typename OutputType >
    auto readInputs()
    {
        return readInputs<OutputType>( 0 );
    }

    template< typename OutputType >
    auto makeInputSignal()
    {
        return [=, sampleIndex = 0] () mutable {
            return readInputs<OutputType>( sampleIndex++ );
        };
    }
};



template <typename UGenClass, size_t InputIndex, size_t NumberOfChannels, typename InputFunctor = detail::Identity>
struct ControlInput:
    SlopedInput<UGenClass, InputIndex, NumberOfChannels, InputFunctor>,
    SignalInput<UGenClass, InputIndex, NumberOfChannels, InputFunctor>
{
    static const size_t index = InputIndex;

    /* audio rate input */
    template <typename OutputType>
    auto makeAudioInputSignal()
    {
        return SignalInput<UGenClass, InputIndex, NumberOfChannels, InputFunctor>::template makeInputSignal<OutputType>();
    }

    /* control rate input */
    template <typename OutputType>
    auto makeRampSignal()
    {
        return SlopedInput< UGenClass, InputIndex, NumberOfChannels, InputFunctor >::template makeRampSignal<OutputType>();
    }

    /* scalar input */
    template <typename OutputType>
    auto makeScalarInputSignal()
    {
        return SlopedInput< UGenClass, InputIndex, NumberOfChannels, InputFunctor >::template makeScalarInputSignal<OutputType>();
    }

    template< typename OutputType >
    auto readInputs( int sampleIndex )
    {
        return SignalInput<UGenClass, InputIndex, NumberOfChannels, InputFunctor>::template readInputs<OutputType>( sampleIndex );
    }

};



template <typename UGenClass, size_t OutputIndex, size_t NumberOfChannels>
struct OutputSink
{
    float * outputVector( unsigned index )
    {
        return static_cast<UGenClass*>(this)->out( OutputIndex + index );
    }

    template <typename OutputType>
    auto makeSink()
    {
        return [=, sampleIndex = 0] (OutputType const & arg) mutable {

            for (int channelIndex = 0; channelIndex != NumberOfChannels; ++channelIndex) {
                float * output = outputVector( channelIndex );
                output[sampleIndex] = boost::simd::extract(arg, channelIndex);
            }
            sampleIndex += 1;
        };
    }
};



}


///////////////////////////
///
/// multichannel parameters

namespace parameter {

template< typename ParameterType, typename SlopeFactor >
auto makeRamp( ParameterType const & current, ParameterType const & next, SlopeFactor slopeFactor )
{
    // this would be nice:
#if 0
    ParameterType slopeVector = (next - current) * slopeFactor;

    return [state = current, slope = slopeVector] () mutable {
        state += slope;
        return state;
    };
#else
    ParameterType slopeVector;
    std::transform( next.begin(), next.end(), current.begin(), slopeVector.begin(), [=](auto nextValue, auto currentValue) {
        return (nextValue - currentValue) * slopeFactor;
    } );

    return [state = current, slope = slopeVector] () mutable {

        std::transform( state.begin(), state.end(), slope.begin(), state.begin(), [=](auto state, auto slope) {
            return state + slope;
        } );

        return state;
    };

#endif
}

} // namespace parameter
} // namespace nova

#endif // PRODUCER_CONSUMER_FUNCTORS_HPP
