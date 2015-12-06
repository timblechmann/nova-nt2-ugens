#include <boost/math/constants/constants.hpp>

#include <approximations/sin.hpp>
#include <approximations/cos.hpp>
#include <approximations/tan.hpp>

#include <approximations/parabol_sine.hpp>

#include <boost/simd/include/pack.hpp>
#include <boost/simd/include/functions/splat.hpp>

using boost::math::constants::pi;

template <typename Test>
void runTestOnFloatRange( float min, float max, float steps, Test && functor ) __attribute__ ((noinline,flatten));

template <typename Test>
void runTestOnFloatRange( float min, float max, float steps, Test && functor )
{
    float inc = (max - min) / steps;
    for( float f = min; f <= max; f += inc) {
        functor( boost::simd::pack<float, 8>{f} );
    }
}

//float dummy;
boost::simd::pack<float, 8> dummy;

int main()
{
    float range = boost::math::float_constants::pi;
    namespace na = nova::approximations;

    for( int i = 0; i != 32; ++i) {
        runTestOnFloatRange( -range, range, 10000000, [](auto x) {
            dummy += na::sin( x, na::SinFaster()  );
//            dummy += na::parabol_sin( x );
        });

//        runTestOnFloatRange( -range, range, 10000000, [](auto x) {
//            dummy += na::sin( x, na::SinFaster()  );
////            dummy += na::parabol_sin( x );
//        });
    }

    std::printf("%g\n", dummy[0]);
    return 0;
}

