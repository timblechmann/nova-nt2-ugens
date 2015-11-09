#include <iostream>

#include <cassert>
#include <cmath>

#include <dsp/parabol_sine.hpp>

#include <boost/simd/include/pack.hpp>

void testSine()
{
    for( float f = -3.143; f < 3.143; f += 0.01 ) {
        float parabol = nova::parabol_sin( f );
        float math    = std::sin( f );

        float absError = parabol - math;

        assert( absError < 0.002 );
    }
}

#if 0
void testSinCos()
{
    for( float f = -3.143; f < 3.143; f += 0.01 ) {
        auto parabol = nova::parabol_sincos( boost::simd::pack<float, 4>(f) );
    }
}
#endif


int main()
{
    testSine();
//    testSinCos();
}
