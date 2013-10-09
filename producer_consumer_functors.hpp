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

#include "boost/simd/include/pack.hpp"
#include "boost/simd/include/functions/aligned_load.hpp"
#include "boost/simd/include/functions/aligned_store.hpp"
#include "boost/simd/include/functions/extract.hpp"
#include "boost/simd/include/functions/groups.hpp"
#include "boost/simd/include/functions/insert.hpp"
#include "boost/simd/include/functions/split.hpp"

#include "boost/simd/sdk/meta/cardinal_of.hpp"

namespace nova {

template <int N>
struct packGenerator
{
};

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


template <typename OutputType, size_t StartIndex = 0>
struct Interleaver
{
	static const size_t N = boost::simd::meta::cardinal_of<OutputType>::value;

	Interleaver(SCUnit * unit):
		unit(unit), cnt(0)
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
	size_t cnt;
};

template <typename OutputType>
struct Deinterleaver
{
	static const size_t N = boost::simd::meta::cardinal_of<OutputType>::value;

	Deinterleaver(SCUnit * unit):
		unit(unit), cnt(0)
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
	size_t cnt;
};

template <typename OutputType, int Channel>
struct Packer
{
	Packer(SCUnit * unit):
		ptr(unit->in(Channel)), cnt(0)
	{}

	inline OutputType operator() ()
	{
		OutputType ret = boost::simd::aligned_load<OutputType>( ptr + cnt);
		cnt += boost::simd::meta::cardinal_of<OutputType>::value;
		return ret;
	}

	const float * ptr;
	size_t cnt;
};

template <typename OutputType, int Channel>
struct Unpacker
{
	Unpacker(SCUnit * unit):
		ptr(unit->out(Channel)), cnt(0)
	{}

	inline void operator() (OutputType arg)
	{
		boost::simd::aligned_store<OutputType>( arg, ptr + cnt );
		cnt += boost::simd::meta::cardinal_of<OutputType>::value;
	}

	float * ptr;
	size_t cnt;
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
Ramp<OutputType> makeRamp( Scalar base, Scalar slope )
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



}

#endif // PRODUCER_CONSUMER_FUNCTORS_HPP
