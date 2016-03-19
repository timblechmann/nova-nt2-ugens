#include <iostream>
#include <boost/range/irange.hpp>

#include <dsp/saturators.hpp>
#include <dsp/diode.hpp>


#include <boost/simd/include/pack.hpp>
#include <boost/simd/sdk/simd/meta/vector_of.hpp>

#include <boost/simd/include/functions/splat.hpp>

namespace bs = boost::simd;

template <typename Functor>
void run_test(Functor && f)
{
    for ( auto i : boost::irange( -200, 200 ) ) {
        float x = float(i) / 100 * 5;

#if 0
        float out = f( x, float(1) );
#else
        using pack = boost::simd::pack<float, 8>;
//        using pack = boost::simd::meta::vector_of<float, 4>::type;
        pack out = f( pack{x}, pack{1.f} );
#endif

        std::cout << x << " " << out << std::endl;
    }
}


int main()
{
//    run_test( [](auto a, auto b) { return nova::saturator::parabol( a, b ); } );
//    run_test( [](auto a, auto b) { return nova::saturator::hyperbol( a, b ); } );
//    run_test( [](auto a, auto b) { return nova::saturator::pow( a, b ); } );
    run_test( [](auto a, auto b) { return nova::diodeSaturation( a, 0.2, 0.8 ); } );
}
