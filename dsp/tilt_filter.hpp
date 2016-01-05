/*
    Copyright (C) 2016 Tim Blechmann

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

#include <boost/simd/include/constants/pi.hpp>
#include <boost/simd/include/constants/zero.hpp>

#include <boost/simd/include/pack.hpp>
#include <boost/simd/include/functions/if_else.hpp>
#include <boost/simd/include/functions/fast_divides.hpp>
#include <boost/simd/include/functions/fast_rec.hpp>
#include <nt2/include/functions/exp.hpp>

#include <approximations/exp.hpp>
#include <dsp/arithmeticarray.hpp>
#include <dsp/range.hpp>
#include <dsp/utils.hpp>


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace nova {

/* adapted from Lubomir I. Ivanov
 * http://www.musicdsp.org/showArchiveComment.php?ArchiveID=267
 * */
template <size_t Size, size_t ParameterSize, typename InternalType = float>
struct TiltFilter
{
    using SampleType    = typename nova::as_pack<InternalType, Size>::type;
    using ParameterType = typename nova::as_pack<InternalType, ParameterSize>::type;

    typedef nova::detail::ArithmeticArray<ParameterType, 4> ParameterState;

    enum {
        A0,
        B1,
        lGain,
        hGain
    };

    TiltFilter() {}

    void setState( ParameterState const & newState ) { _state = newState; }
    auto getState()
    {
        auto ret = [=] { return _state; };
        return ret;
    }
    ParameterState currentState() { return _state; }


    // parameters
    template <typename DSPContext>
    static auto computeState( ParameterType frequency, ParameterType gain, DSPContext const & dspContext )
    {
        using namespace boost::simd;

        typedef typename meta::scalar_of<ParameterType>::type ParameterScalar;
        typedef decltype(dspContext.sampleRate()) SampleRateType;
        frequency = nova::clip( frequency, ParameterScalar(0.01f), ParameterScalar( dspContext.sampleRate() * (SampleRateType)(0.5)) );

        constexpr InternalType ampFactor = std::log( 2.0 ) / 6.0;

        InternalType gfactor = 5; // gfactor is the proportional gain

        ParameterType g1, g2;

        g1 = if_else( gain > Zero<ParameterType>(), -gfactor*gain,        -gain);
        g2 = if_else( gain > Zero<ParameterType>(),          gain, gfactor*gain);

        //two separate gains
        ParameterType lowGain = nova::approximations::exp<ParameterType>( g1 * ampFactor, nova::approximations::ExpFast() ) - 1.f;
        ParameterType hiGain  = nova::approximations::exp<ParameterType>( g2 * ampFactor, nova::approximations::ExpFast() ) - 1.f;

        //filter
        ParameterType twopi = Pi<ParameterType>() + Pi<ParameterType>();
        ParameterType omega = twopi * frequency;
        ParameterType   sr3 = dspContext.sampleRate() + dspContext.sampleRate() + dspContext.sampleRate();
        ParameterType     n = fast_rec( sr3 + omega);
        ParameterType    a0 = (omega+omega)*n;
        ParameterType    b1 = (sr3 - omega)*n;

        return ParameterState{ a0, b1, lowGain, hiGain };
    }

    // DSP
    template < typename InputFunctor, typename OutputFunctor >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count )
    {
        doRun<true>( in, out, count, getState() );
    }

    template < typename InputFunctor, typename OutputFunctor, typename State >
    inline void run ( InputFunctor & in, OutputFunctor & out, size_t count, State && a )
    {
        doRun<false>( in, out, count, a );
    }

private:
    template < bool StateImmutable, typename InputFunctor, typename OutputFunctor, typename State >
    inline void doRun ( InputFunctor & in, OutputFunctor & out, size_t count, State && state )
    {
        ParameterType a0, b1, lgain, hgain;

        if( StateImmutable ) {
            ParameterState a = state();
            a0 = a[A0], b1 = a[B1], lgain = a[lGain], hgain = a[hGain];
        }

        SampleType lp_out = lp_out_;

        for( size_t i : nova::range( count ) ) {
            if( !StateImmutable ) {
                ParameterState a = state();
                a0 = a[A0], b1 = a[B1], lgain = a[lGain], hgain = a[hGain];
            }

            SampleType inSample = in();

            lp_out            = a0*inSample +    b1*lp_out;
            SampleType result =    inSample + lgain*lp_out + hgain*(inSample - lp_out);

            out( result );
        }

        lp_out_ = lp_out;
    }

    // state
    SampleType lp_out_ = SampleType{0};

    ParameterState _state; // a0, b1, lgain, hgain;
};

}
