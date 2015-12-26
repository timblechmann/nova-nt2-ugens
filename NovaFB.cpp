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
#include <boost/simd/memory/is_aligned.hpp>

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


template <typename T, typename Functor = nova::detail::Identity, typename InputTag = aligned, typename OutputTag = aligned>
inline void transformCacheAware( T * dest, const T * src, unsigned count, Functor && f = Functor(), InputTag input = InputTag(), OutputTag output = OutputTag() )
{
    typedef boost::simd::pack<T> pack;
    constexpr size_t framesPerPack      = pack::static_size;
    constexpr size_t packsPerCacheline  = 4;

    const     size_t unrolledLoops  = count / packsPerCacheline;
    const     size_t remainingPacks = count - unrolledLoops * packsPerCacheline * framesPerPack;

    for( int index : nova::range(unrolledLoops) ) {
        pack p0      = load<pack>( src + index * packsPerCacheline,                   input );
        pack result0 = f( p0 );
        pack p1      = load<pack>( src + index * packsPerCacheline + 1*framesPerPack, input );
        pack result1 = f( p1 );
        pack p2      = load<pack>( src + index * packsPerCacheline + 2*framesPerPack, input );
        pack result2 = f( p2 );
        pack p3      = load<pack>( src + index * packsPerCacheline + 3*framesPerPack, input );
        pack result3 = f( p3 );
        store( dest + index * packsPerCacheline,                   result0, output );
        store( dest + index * packsPerCacheline + 1*framesPerPack, result1, output );
        store( dest + index * packsPerCacheline + 2*framesPerPack, result2, output );
        store( dest + index * packsPerCacheline + 3*framesPerPack, result3, output );
    }

    if( BOOST_UNLIKELY( remainingPacks > 0 ) ) {
        src  += unrolledLoops * packsPerCacheline * framesPerPack;
        dest += unrolledLoops * packsPerCacheline * framesPerPack;
        for( int index : nova::range(remainingPacks) ) {
            pack p0      = load<T>( src + index, unaligned() );
            pack result0 = f( p0 );
            store( dest + index, result0, unaligned() );
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

        size_t allocSize = mBufLength * sizeof(float) * mChannelCount;
        mBuff = (float*)RTAlloc( mWorld, allocSize );
        memset(mBuff, 0, allocSize);

        if( boost::simd::is_aligned( bufferSize() ) ) {
            auto channelRange = nova::range( mChannelCount );
            bool outputsAligned = std::any_of( channelRange.begin(), channelRange.end(), [this] (size_t channel) {
                return boost::simd::is_aligned( out( channel ) );
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
        RTFree(mWorld, mBuff);
    }

private:
    template< typename OutputAlignemnt>
    void next_vec(int inNumSamples)
    {
        for( int index : nova::range( mChannelCount ) ) {
            const float * buf = mBuff + mBufLength * index;
            float * dest = out(index + 1);

            transformCacheAware( dest, buf, inNumSamples, nova::detail::Identity(), aligned(), OutputAlignemnt() );
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
            bool inputsAligned = std::any_of( channelRange.begin(), channelRange.end(), [this] (size_t channel) {
                return boost::simd::is_aligned( in( channel ) );
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

            transformCacheAware( buf, source, inNumSamples, nova::detail::Identity(), InputAlignemnt(), aligned() );
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
