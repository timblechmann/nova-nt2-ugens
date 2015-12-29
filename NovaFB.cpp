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

#include <NovaUGensCommon.hpp>

#include <dsp/range.hpp>

#include <boost/simd/include/functions/aligned_load.hpp>
#include <boost/simd/include/functions/aligned_store.hpp>
#include <boost/simd/include/functions/load.hpp>
#include <boost/simd/include/functions/store.hpp>
#include <boost/simd/include/functions/stream.hpp>
#include <boost/simd/memory/align_on.hpp>
#include <boost/simd/memory/is_aligned.hpp>

#include <iostream>

static InterfaceTable *ft;

struct aligned {};
struct unaligned {};
struct stream : aligned {};

template <typename T, typename U>
inline auto load( const U * ptr, aligned )
{
    return boost::simd::aligned_load<T>(ptr);
}

template <typename T, typename U>
inline auto load( const U * ptr, unaligned )
{
    return boost::simd::load<T>(ptr);
}

template <typename T, typename U>
inline void store( T * ptr, U value, aligned )
{
    boost::simd::aligned_store(value, ptr);
}

template <typename T, typename U>
inline void store( T * ptr, U value, unaligned )
{
    boost::simd::store(value, ptr);
}

template <typename T, typename U>
inline void store( T * ptr, U value, stream )
{
    boost::simd::stream(value, ptr);
}


template <typename T, typename InputTag = aligned, typename OutputTag = aligned>
inline void transformCacheAware( T * dest, const T * src, unsigned count, InputTag input = InputTag(), OutputTag output = OutputTag() )
{
    typedef boost::simd::pack<T> pack;
    constexpr size_t framesPerPack      = pack::static_size;
    constexpr size_t packsPerCacheline  = 4;
    constexpr size_t framesPerCacheline = packsPerCacheline * framesPerPack;

    const     size_t unrolledLoops   = count / framesPerCacheline;
    const     size_t remainingFrames = count - unrolledLoops * framesPerCacheline;

    for( int index : nova::range(unrolledLoops) ) {
        pack p0      = load<pack>( src + index * framesPerCacheline,                   input );
        pack p1      = load<pack>( src + index * framesPerCacheline + 1*framesPerPack, input );
        pack p2      = load<pack>( src + index * framesPerCacheline + 2*framesPerPack, input );
        pack p3      = load<pack>( src + index * framesPerCacheline + 3*framesPerPack, input );
        store( dest + index * framesPerCacheline,                   p0, output );
        store( dest + index * framesPerCacheline + 1*framesPerPack, p1, output );
        store( dest + index * framesPerCacheline + 2*framesPerPack, p2, output );
        store( dest + index * framesPerCacheline + 3*framesPerPack, p3, output );
    }

    if( BOOST_UNLIKELY( remainingFrames > 0 ) ) {
        src  += unrolledLoops * framesPerCacheline;
        dest += unrolledLoops * framesPerCacheline;
        for( int index : nova::range(remainingFrames) ) {
            pack p0      = load<T>( src + index, unaligned() );
            store( dest + index, p0, unaligned() );
        }
    }
}


struct NovaFBIn:
    public SCUnit
{
public:
    NovaFBIn()
    {
        mChannelCount = in0(0);

        size_t allocSize = mBufLength * sizeof(float) * mChannelCount + 32;
        realBuf = (float*)RTAlloc( mWorld, allocSize );
        if( !realBuf ) {
            std::cout << "alloc failure" << std::endl;
            mCalcFunc = ft->fClearUnitOutputs;
            return;
        }

        mBuff = boost::simd::align_on( realBuf, 32 ); // force alignment

        memset(realBuf, 0, allocSize);

        if( boost::simd::is_aligned( bufferSize() ) ) {
            auto channelRange = nova::range( mChannelCount );
            bool outputsAligned = std::all_of( channelRange.begin(), channelRange.end(), [this] (size_t channel) {
                return boost::simd::is_aligned( out( channel ) + 1 );
            });

            if( outputsAligned )
                set_vector_calc_function<NovaFBIn, &NovaFBIn::next_vec<aligned>,   &NovaFBIn::next_i>();
            else
                set_vector_calc_function<NovaFBIn, &NovaFBIn::next_vec<unaligned>, &NovaFBIn::next_i>();
        } else
            set_calc_function<NovaFBIn, &NovaFBIn::next_i>();
    }

    ~NovaFBIn()
    {
        RTFree(mWorld, realBuf);
    }

private:
    template< typename OutputAlignemnt>
    void next_vec(int inNumSamples)
    {
        for( int index : nova::range( mChannelCount ) ) {
            const float * buf = mBuff + mBufLength * index;
            float * dest = out(index + 1);

            transformCacheAware( dest, buf, inNumSamples, aligned(), OutputAlignemnt() );
        };
    }

    void next_i(int inNumSamples)
    {
        for( int index : nova::range( mChannelCount ) ) {
            const float * buf = mBuff + mBufLength * index;
            float * dest = out(index + 1);

            memcpy(dest, buf, inNumSamples * sizeof(float));
        };
    }

    friend class NovaFBOut;
    float * mBuff;
    float * realBuf;
    int mChannelCount;
};

struct NovaFBOut:
    public SCUnit
{
public:
    NovaFBOut()
    {
        mChannelCount = in0(1);

        NovaFBIn * inputNode = static_cast<NovaFBIn*>(mInput[0]->mFromUnit);
        mBuff = inputNode->mBuff;

        if( boost::simd::is_aligned( bufferSize() ) ) {
            auto channelRange = nova::range( mChannelCount );
            bool inputsAligned = std::all_of( channelRange.begin(), channelRange.end(), [this] (size_t channel) {
                return boost::simd::is_aligned( in( 2 + channel ) );
            });

            if( inputsAligned )
                set_vector_calc_function<NovaFBOut, &NovaFBOut::next_vec<aligned>,   &NovaFBOut::next_i>();
            else
                set_vector_calc_function<NovaFBOut, &NovaFBOut::next_vec<unaligned>, &NovaFBOut::next_i>();
        } else
            set_calc_function<NovaFBOut, &NovaFBOut::next_i>();
    }

private:
    template< typename InputAlignemnt>
    void next_vec(int inNumSamples)
    {
        for( int index : nova::range( mChannelCount ) ) {
            const float * source = in( 2 + index );
            float * buf = mBuff + mBufLength * index;

            transformCacheAware( buf, source, inNumSamples, InputAlignemnt(), aligned() );
        };
    }

    void next_i(int inNumSamples)
    {
        for( int index : nova::range( inNumSamples ) ) {
            const float * source = in( 2 + index );
            float * buf = mBuff + mBufLength * index;

            memcpy(buf, source, inNumSamples * sizeof(float));
        };
    }

    float * mBuff;
    int mChannelCount;
};


PluginLoad(NovaFB)
{
    ft = inTable;
    nova::registerUnit< NovaFBIn  >( ft, "NovaFBIn"  );
    nova::registerUnit< NovaFBOut >( ft, "NovaFBOut" );
}
