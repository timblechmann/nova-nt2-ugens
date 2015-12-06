#include <gtest/gtest.h>

#include <boost/math/constants/constants.hpp>

#include <approximations/sin.hpp>
#include <approximations/cos.hpp>
#include <approximations/tan.hpp>
#include <approximations/parabol_sine.hpp>

using boost::math::constants::pi;

template <typename Test>
void runTestOnFloatRange( float min, float max, float steps, Test && functor )
{
    float inc = (max - min) / steps;
    for( float f = min; f <= max; f += inc) {
        functor(f);
    }
}

//////////////
// sin

TEST(sin_approximations, sin_precise)
{
    float pi = boost::math::float_constants::pi;
    namespace na = nova::approximations;

    runTestOnFloatRange( -pi, pi, 1000, [](float x) {  ASSERT_NEAR( na::sin( x, na::SinPrecise() ), std::sin(x), 0.000001 ); });
}

TEST(sin_approximations, sin_fast)
{
    float pi = boost::math::float_constants::pi;
    namespace na = nova::approximations;

    runTestOnFloatRange( -pi, pi, 1000, [](float x) {  ASSERT_NEAR( na::sin( x, na::SinFast()    ), std::sin(x), 0.0001 ); });
}

TEST(sin_approximations, sin_faster)
{
    float pi = boost::math::float_constants::pi;
    namespace na = nova::approximations;

    runTestOnFloatRange( -pi, pi, 1000, [](float x) {  ASSERT_NEAR( na::sin( x, na::SinFaster()  ), std::sin(x), 0.001 ); });
}

TEST(sin_approximations, sin_parabol)
{
    float pi = boost::math::float_constants::pi;
    namespace na = nova::approximations;

    runTestOnFloatRange( -pi, pi, 1000, [](float x) {  ASSERT_NEAR( na::parabol_sin( x ), std::sin(x), 0.0011 ); });
}

//////////////
// cos

TEST(cos_approximations, cos_precise)
{
    float pi = boost::math::float_constants::pi;
    namespace na = nova::approximations;

    runTestOnFloatRange( -pi, pi, 1000, [](float x) {  ASSERT_NEAR( na::cos( x, na::CosPrecise() ), std::cos(x), 0.000001 ); });
}

TEST(cos_approximations, cos_fast)
{
    float pi = boost::math::float_constants::pi;
    namespace na = nova::approximations;

    runTestOnFloatRange( -pi, pi, 1000, [](float x) {  ASSERT_NEAR( na::cos( x, na::CosFast()    ), std::cos(x), 0.0001 ); });
}

TEST(cos_approximations, cos_faster)
{
    float pi = boost::math::float_constants::pi;
    namespace na = nova::approximations;

    runTestOnFloatRange( -pi, pi, 1000, [](float x) {  ASSERT_NEAR( na::cos( x, na::CosFaster()  ), std::cos(x), 0.01 ); });
}

//////////////
// tan

TEST(tan_approximations, tan_precise)
{
    float range = boost::math::float_constants::pi/2 * 0.95;
    namespace na = nova::approximations;

    runTestOnFloatRange( -range, range, 1000, [](float x) {  ASSERT_NEAR( na::tan( x, na::TanPrecise() ), std::tan(x), 0.000001 ); });
}

TEST(tan_approximations, tan_fast)
{
    float range = boost::math::float_constants::pi/4;
    namespace na = nova::approximations;

    runTestOnFloatRange( -range, range, 1000, [](float x) {  ASSERT_NEAR( na::tan( x, na::TanFast()    ), std::tan(x), 0.0001 ); });
}

TEST(tan_approximations, tan_faster)
{
    float range = boost::math::float_constants::pi/4;
    namespace na = nova::approximations;

    runTestOnFloatRange( -range, range, 1000, [](float x) {  ASSERT_NEAR( na::tan( x, na::TanFaster()  ), std::tan(x), 0.001 ); });
}


