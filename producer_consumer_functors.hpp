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

#include "boost/simd/include/pack.hpp"
#include "boost/simd/memory/include/functions/aligned_load.hpp"
#include "boost/simd/memory/include/functions/aligned_store.hpp"
#include "boost/simd/swar/include/functions/groups.hpp"
#include "boost/simd/swar/include/functions/split.hpp"

namespace nova {

template <typename OutputType>
struct Interleaver2
{
	typedef boost::simd::pack<double, 2> v2d;
	typedef boost::simd::pack<float,  4> v4f;

	Interleaver2(SCUnit * unit):
		unit(unit), cnt(0)
	{}

	BOOST_FORCEINLINE OutputType operator() ()
	{
		OutputType ret(unit->in(0)[cnt], unit->in(1)[cnt]);
//		v4f retFloat(unit->in(0)[cnt], unit->in(1)[cnt], 0.f, 0.f);
//		v2d ret, dummy;
//		boost::simd::split(retFloat, dummy, ret);
		cnt += 1;
		return ret;
	}

	SCUnit * unit;
	size_t cnt;
};

template <typename OutputType>
struct Deinterleaver2
{
	typedef boost::simd::pack<double, 2> v2d;
	typedef boost::simd::pack<float,  4> v4f;

	Deinterleaver2(SCUnit * unit):
		unit(unit), cnt(0)
	{}

	BOOST_FORCEINLINE void operator() (OutputType arg)
	{
		v4f asFloat = boost::simd::group(arg, arg);
		for (int i = 0; i != 2; ++i) {
			float * out = unit->out(i);
			out[cnt] = boost::simd::extract(asFloat, i);
//			out[cnt] = boost::simd::extract(arg, i);
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


}

#endif // PRODUCER_CONSUMER_FUNCTORS_HPP
