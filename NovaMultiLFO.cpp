/*
    Copyright (C) Tim Blechmann

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


#include "NovaUGensCommon.hpp"

#include "approximations/sin.hpp"
#include "dsp/range.hpp"

#include <cmath>

using nova::NovaUnit;

InterfaceTable *ft;

namespace nova {

enum lfoMode {
    minVal,
    maxVal,
    randVal,
    sinVal
};


struct MultiLFO:
    public NovaUnit,

    public nova::ScalarInput< MultiLFO, 0 >, // mode
    public nova::ScalarInput< MultiLFO, 1 >, // freq
    public nova::ScalarInput< MultiLFO, 2 >, // phase (0..1)
    public nova::ScalarInput< MultiLFO, 3 >, // min
    public nova::ScalarInput< MultiLFO, 4 >  // max
{

    typedef nova::ScalarInput< MultiLFO, 0 > ModeInput;
    typedef nova::ScalarInput< MultiLFO, 1 > FreqInput;
    typedef nova::ScalarInput< MultiLFO, 2 > PhaseInput;
    typedef nova::ScalarInput< MultiLFO, 3 > MinInput;
    typedef nova::ScalarInput< MultiLFO, 4 > MaxInput;

public:
    MultiLFO()
    {
        const int mode = int( ModeInput::readInput() );

        float minIn       = MinInput::readInput();
        float maxIn       = MaxInput::readInput();
        const float range  = calcRange( minIn, maxIn );

        switch( mode ) {
        case minVal: mValue = minIn; break;
        case maxVal: mValue = maxIn; break;

        case randVal: {
            mValue = mParent->mRGen->frand() * range + minIn;
            break;
        }

        case sinVal: {
            mPhase = PhaseInput::readInput();
            mPhase = sc_wrap( mPhase + 0.5f, 0.f, 1.f ) * twopi - pi;

            set_calc_function< MultiLFO, &MultiLFO::sin >();
            return;
        }

        default:
            ;
        }

        if( bufferSize() == 1)
            set_calc_function< MultiLFO, &MultiLFO::setValue_1 >();
        else
            set_calc_function< MultiLFO, &MultiLFO::setValue >();
    }

private:
    void sin( int numSamples )
    {
        float minVal      = MinInput::readInput();
        float maxVal      = MaxInput::readInput();
        const float range = calcRange( minVal, maxVal );

        float newFreq = FreqInput::readInput();
        if( BOOST_UNLIKELY( newFreq != mFreq ) ) {
            newFreq = nova::clip2( newFreq, float( sampleRate() ) * 0.5f );

            mPhaseIncrement = newFreq * float( sampleDur() ) * float( twopi );
            mFreq           = newFreq;
        }

        float phaseIncrement = mPhaseIncrement;

        float phase = mPhase;
        for( int sampleIndex : nova::range( numSamples ) ) {
            float outSin = nova::approximations::sin( phase, nova::approximations::SinFast() );
            out( 0 )[  sampleIndex ] = outSin * range + minVal;

            phase += phaseIncrement;

            phase = (phase >= pi) ? phase - twopi
                                  : phase;
        }
        mPhase = phase;
    }

    void setValue( int numSamples )
    {
        for( int sampleIndex : nova::range( numSamples ) )
            out( 0 )[ sampleIndex ] = mValue;
    }

    void setValue_1( int numSamples )
    {
        out(0)[0] = mValue;
    }

    float calcRange( float & a, float & b )
    {
        float low = std::min( a, b );
        float hi  = std::max( a, b );

        a = low, b = hi;
        return hi - low;
    }

    float mValue = 0;
    float mPhase = 0.f;
    float mFreq  = std::numeric_limits<float>::quiet_NaN();
    float mPhaseIncrement = 0.f;
};

}

PluginLoad(NovaMultiLFO)
{
    using namespace nova;

    ft = inTable;

    nova::registerUnit<MultiLFO>( ft, "NovaMultiLFO" );
}
