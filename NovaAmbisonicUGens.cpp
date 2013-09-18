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

#include "SC_PlugIn.hpp"

#include "boost/simd/include/pack.hpp"
#include "boost/simd/include/native.hpp"
#include "boost/simd/operator/include/functions/minus.hpp"
#include "boost/simd/operator/include/functions/multiplies.hpp"
#include "boost/simd/operator/include/functions/plus.hpp"
#include "boost/simd/constant/constants.hpp"

#include "producer_consumer_functors.hpp"

#include "boost/simd/arithmetic/include/functions/abs.hpp"
#include "boost/simd/arithmetic/include/functions/fast_rec.hpp"
#include "boost/simd/ieee/include/functions/copysign.hpp"

#include <boost/math/constants/constants.hpp>

#include <nt2/include/functions/sincos.hpp>

#include "utils.hpp"

namespace {

InterfaceTable *ft;

struct NovaPanB2D:
	public SCUnit
{
	static float verifyLevel ( float arg ) { return arg; }

	NovaPanB2D()
	{
		calcAmps(in0(1));

		if (boost::simd::is_aligned( bufferSize(), 8 ) )
			set_vector_calc_function<NovaPanB2D,
					&NovaPanB2D::run_i< boost::simd::pack<float, 8> >,
					&NovaPanB2D::run_i<float> >();
		else
			set_calc_function<NovaPanB2D, &NovaPanB2D::run_i<float>>();
	}

	template <typename SampleType>
	void run_i (int inNumSamples)
	{
		using namespace boost::simd;

		SampleType wamp(1/std::sqrt(2.f)), xamp(_xamp), yamp(_yamp);

		const size_t unroll = meta::cardinal_of<SampleType>::value;

		int loops = inNumSamples / unroll;

		const float * input = in(0);
		float * w = out(0);
		float * x = out(1);
		float * y = out(2);
		for (int i = 0; i != loops; ++i)
		{
			SampleType sig  = aligned_load<SampleType>(input + i * unroll);
			auto wOut       = sig * wamp;
			auto xOut       = sig * xamp;
			auto yOut       = sig * yamp;

			aligned_store(wOut, w + i * unroll);
			aligned_store(xOut, x + i * unroll);
			aligned_store(yOut, y + i * unroll);
		}
	}

	void calcAmps(float azimuth)
	{
		std::tie(_yamp, _xamp) = nt2::sincos(azimuth * boost::math::float_constants::pi);
	}

	float _xamp, _yamp;
};



DEFINE_XTORS(NovaPanB2D)
//DEFINE_XTORS(NovaPanB3D)

}

PluginLoad(NovaSaturators)
{
	ft = inTable;
	DefineDtorUnit(NovaPanB2D);
}


