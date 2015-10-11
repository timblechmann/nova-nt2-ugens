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
#include <boost/range/irange.hpp>

#include <producer_consumer_functors.hpp>

namespace nova {

struct NovaUnit:
    public SCUnit
{
#if 0
    template <typename FloatType, typename Functor>
    struct ScalarSignal
    {
        ScalarSignal(FloatType value, Functor const & f):
            value(f(value))
        {}

        FloatType consume() const
        {
            return value;
        }

        FloatType value;
    };

    template <typename FloatType, typename Functor>
    struct SlopeSignal:
            Functor
    {
        SlopeSignal(FloatType value, FloatType slope, Functor const & f):
            Functor(f), value(value), slope(slope)
        {}

        FloatType consume()
        {
            FloatType ret = value;
            value += slope;
            return Functor::operator()(ret);
        }

        FloatType value, slope;
    };

    template <typename FloatType, typename Functor>
    struct AudioSignal:
            Functor
    {
        AudioSignal(const FloatType * pointer, Functor const & f):
            Functor(f), pointer(pointer)
        {}

        FloatType consume()
        {
            return Functor::operator()(*pointer++);
        }

        const FloatType * pointer;
    };

    template <typename FloatType, typename Functor>
    inline auto makeScalar(FloatType value, Functor const & f) const
    {
        return ScalarSignal<FloatType, Functor>(value, f);
    }

    template <typename FloatType>
    inline auto makeScalar(FloatType value) const
    {
        return SCUnit::makeScalar(value);
    }

    template <typename FloatType, typename Functor>
    inline auto makeSlope(FloatType next, FloatType last, Functor const & f) const
    {
        return SlopeSignal<FloatType, Functor>(last, calcSlope(next, last), f);
    }

    template <typename FloatType>
    inline auto makeSlope(FloatType next, FloatType last) const
    {
        return SCUnit::makeSlope(next, last);
    }

    template <size_t Size, typename Functor>
    inline auto makeSignal(int index, Functor const & f) const
    {
        const float * input = in(index);
        return AudioSignal<float, Functor>(input, f);
    }

    template <typename Functor>
    inline auto makeSignal(int index) const
    {
        return SCUnit::makeSignal(index);
    }
#endif

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
    auto makeSlope(FloatTypeA next, FloatTypeB current) const
    {
        FloatTypeA increment = calcSlope( next, current );
        return [=, state = current] () mutable {
            auto ret = state;
            state += increment;
            return ret;
        };
    }
};

}



#endif // NOVAUGENSCOMMON_HPP
