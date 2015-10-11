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

#include "boost/simd/include/functions/aligned_load.hpp"
#include "boost/simd/include/functions/aligned_store.hpp"
#include "boost/simd/include/functions/extract.hpp"
#include "boost/simd/include/functions/groups.hpp"
#include "boost/simd/include/functions/insert.hpp"
#include <boost/simd/include/functions/splat.hpp>
#include "boost/simd/include/functions/split.hpp"


#include "boost/simd/operator/include/functions/plus.hpp"
#include "boost/simd/operator/include/functions/multiplies.hpp"
#include "boost/simd/include/functions/assign.hpp"

#include <boost/dispatch/meta/is_scalar.hpp>
#include "boost/simd/sdk/meta/cardinal_of.hpp"

namespace nova {

template <int N>
struct packGenerator
{};

template <>
struct packGenerator<1>
{
    template <typename Result, typename Functor>
    static Result generate(Functor const & f)
    {
        Result ret = f();
        return ret;
    }
};

template <>
struct packGenerator<2>
{
    template <typename Result, typename Functor>
    static Result generate(Functor const & f)
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
    static Result generate(Functor const & f)
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
    static Result generate(Functor const & f)
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


template <typename OutputType, size_t StartIndex = 0, class Enable = void>
struct Interleaver
{
    static const size_t N = boost::simd::meta::cardinal_of<OutputType>::value;

    Interleaver(SCUnit * unit):
        unit(unit)
    {}

    BOOST_FORCEINLINE OutputType operator() ()
    {
        //		OutputType ret;
        //		for (size_t i = 0; i != N; ++i)
        //			boost::simd::insert(unit->in( StartIndex + i )[cnt], ret, i);

        size_t i = 0;
        OutputType ret = packGenerator<N>::template generate<OutputType>([&](){
            return unit->in( StartIndex + i++ )[cnt];
        });

        cnt += 1;
        return ret;
    }


    SCUnit * unit;
    size_t cnt = 0;
};

template <typename OutputType, size_t StartIndex>
struct Interleaver<OutputType, StartIndex, std::enable_if<std::is_scalar<OutputType>::value>>
{
    Interleaver(SCUnit * unit):
        unit(unit)
    {}

    BOOST_FORCEINLINE OutputType operator() ()
    {
        return unit->in( StartIndex )[cnt++];
    }


    SCUnit * unit;
    size_t cnt = 0;
};

template <typename OutputType>
struct Deinterleaver
{
    static const size_t N = boost::simd::meta::cardinal_of<OutputType>::value;

    Deinterleaver(SCUnit * unit):
        unit(unit)
    {}

    BOOST_FORCEINLINE void operator() (OutputType arg)
    {
        for (int i = 0; i != N; ++i) {
            float * out = unit->out(i);
            out[cnt] = boost::simd::extract(arg, i);
        }
        cnt += 1;
    }

    SCUnit * unit;
    size_t cnt = 0;
};

template <typename OutputType, int Channel>
struct Packer
{
    Packer(SCUnit * unit):
        ptr(unit->in(Channel))
    {}

    inline OutputType operator() ()
    {
        OutputType ret = boost::simd::aligned_load<OutputType>( ptr + cnt);
        cnt += boost::simd::meta::cardinal_of<OutputType>::value;
        return ret;
    }

    const float * ptr;
    size_t cnt = 0;
};

template <typename OutputType, int Channel>
struct Unpacker
{
    Unpacker(SCUnit * unit):
        ptr(unit->out(Channel))
    {}

    inline void operator() (OutputType arg)
    {
        boost::simd::aligned_store<OutputType>( arg, ptr + cnt );
        cnt += boost::simd::meta::cardinal_of<OutputType>::value;
    }

    float * ptr;
    size_t cnt = 0;
};


template <typename OutputType>
struct Scalar
{
    Scalar(OutputType val):
        _val(val)
    {}

    inline OutputType operator() () const
    {
        return _val;
    }

    OutputType _val;
};


template <typename OutputType>
struct Ramp
{
    Ramp(OutputType val, OutputType increment):
        _val(val), _increment(increment)
    {}

    inline OutputType operator() ()
    {
        OutputType ret = _val;
        _val += _increment;
        return _val;
    }

    OutputType _val;
    OutputType _increment;
};


template <typename OutputType, typename Scalar>
auto makeRamp( Scalar base, Scalar slope )
{
    using namespace boost::simd;

    const size_t size = meta::cardinal_of<OutputType>::value;

    OutputType val;

    insert(base, val, 0);

    for (size_t i = 1; i != size; ++i) {
        val += slope;
        insert(base, val, i);
    }

    OutputType increment = slope * size;

    return Ramp<OutputType>(val, increment);
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


template <size_t N>
struct InputInterleaver
{
    typedef typename as_pack<float, N>::type HostParameterType;

    inline static HostParameterType read( SCUnit * unit, size_t index )
    {
        HostParameterType ret;

        for (size_t i = 0; i != N; ++i) {
            float input = unit->in0(index + i);
            boost::simd::insert(input, ret, i);
        }
        return ret;
    }
};

namespace detail {

struct Identity {
    template <typename Arg>
    auto operator()(Arg && arg) const { return arg; }
};

}


template <typename UGenClass, size_t InputIndex, typename InputFunctor = detail::Identity>
struct ScalarInput:
    private InputFunctor
{
    auto readInput()
    {
        return InputFunctor::operator()( static_cast<UGenClass*>(this)->in0( InputIndex ) );
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
    SlopedInput():
        mState( readInput() )
    {}

    auto readInput()
    {
        return ScalarInput< UGenClass, InputIndex, InputFunctor >::readInput();
    }

    template <typename SIMDType,
              typename std::enable_if< !boost::dispatch::meta::is_scalar<SIMDType>::value
                                       >::type * = nullptr>
    auto makeInputSignal( int numberOfSamples )
    {
        using ScalarType = typename boost::simd::meta::scalar_of<SIMDType>::type;

        ScalarType current = mState;
        ScalarType next  = readInput();
        ScalarType slope = (next - current) / ScalarType( numberOfSamples );
        mState = next;
        return makeRamp<SIMDType>( current, slope );
    }

    template <typename ScalarType,
              typename std::enable_if< boost::dispatch::meta::is_scalar<ScalarType>::value
                                       >::type * = nullptr>
    auto makeInputSignal( int numberOfSamples )
    {
        ScalarType current = mState;
        ScalarType next  = readInput();
        ScalarType slope = (next - current) / ScalarType( numberOfSamples );
        mState = next;
        return makeRamp<ScalarType>( current, slope );
    }

    template <typename OutputType>
    auto makeScalarInputSignal()
    {
        return ScalarInput< UGenClass, InputIndex, InputFunctor >::template makeInputSignal<OutputType>();
    }

    bool changed()
    {
        float next = readInput();
        return next != mState;
    }

    float mState;
};


template <typename UGenClass, size_t InputIndex, typename InputFunctor = detail::Identity>
struct SignalInput:
    InputFunctor
{
    const float * inputVector()
    {
        return static_cast<UGenClass*>(this)->in( InputIndex );
    }

    template <typename SIMDType,
              typename std::enable_if< !boost::dispatch::meta::is_scalar<SIMDType>::value
                                       >::type * = nullptr>
    auto makeInputSignal()
    {
        const float * input = inputVector();
        return [=] () mutable {
            SIMDType ret = InputFunctor::operator()( boost::simd::aligned_load<SIMDType>( input ) );
            input += boost::simd::meta::cardinal_of<SIMDType>::value;
            return ret;
        };
    }

    template <typename ScalarType,
              typename std::enable_if< boost::dispatch::meta::is_scalar<ScalarType>::value
                                       >::type * = nullptr>
    auto makeInputSignal()
    {
        const float * input = inputVector();
        return [=] () mutable {
            ScalarType ret = InputFunctor::operator()( *input );
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
        return SignalInput<UGenClass, InputIndex, InputFunctor>::template makeInputSignal<OutputType>();
    }

    /* control rate input */
    template <typename OutputType>
    auto makeSmoothedInputSignal( int numberOfSamples )
    {
        return SlopedInput< UGenClass, InputIndex, InputFunctor >::template makeInputSignal<OutputType>( numberOfSamples );
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


}

#endif // PRODUCER_CONSUMER_FUNCTORS_HPP
