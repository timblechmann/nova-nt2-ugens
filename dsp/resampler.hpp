/*
    Copyright (C) 2015 Tim Blechmann

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

#include <type_traits>

#include <boost/simd/include/pack.hpp>
#include <boost/simd/constant/include/constants/zero.hpp>
#include <boost/simd/include/functions/split.hpp>
#include <boost/simd/include/functions/group.hpp>
#include <boost/simd/include/functions/groups.hpp>
#include <boost/simd/operator/include/functions/multiplies.hpp>
#include <boost/simd/operator/include/functions/plus.hpp>

#ifndef RESAMPLER_HPP
#define RESAMPLER_HPP

namespace nova {

struct MovingAverageFilter
{
    float filter( float x0 )
    {
        float x1 = x1_;
        float y0 = (x1 + x0) * 0.5f;
        x1_ = x0;
        return y0;
    }

    std::tuple<float, float> polyphase_filter( float x0 )
    {
        float x1 = x1_;
        float y0 = (x1 + x0) * 0.5f;
        x1_ = x0;
        return y0;
    }

    float x1_ = {0.f};
};


////////////////////////
//
// sample/hold resampler

struct Upsampler_2_sh
{
    std::tuple<float, float> upsample( float in ) { return std::make_pair( in, in ); }
};

struct Downsampler_2_sh
{
    float downsample( float in0, float in1 )      {  return in0; }
};


////////////////////////
//
// linear interpolation resampler

struct Upsampler_2_lin
{
    std::tuple<float, float> upsample( float in )
    {
        float sample1 = (in - last_) * 0.5f;
        last_ = in;
        return std::make_pair( sample1, in );
    }

    float last_ = {0.f};
};

struct Downsampler_2_lin
{
    float downsample( float in0, float in1 )
    {
        return (in0 + in1) * 0.5f;
    }
};


////////////////////////
//
// halfband resampler


struct Upsampler_2_halfband
{
    std::tuple<float, float> upsample( float in )
    {
        float out0 = filter.filter( in    );
        float out1 = filter.filter( 0.f   );
        return std::make_pair( out0, out1 );
    }

    MovingAverageFilter filter;
};

struct Downsampler_2_halfband
{
    float downsample( float in0, float in1 )
    {
        std::ignore = filter.filter( in0 );
        return filter.filter( in1 );
    }

    MovingAverageFilter filter;
};


}

#endif // RESAMPLER_HPP
