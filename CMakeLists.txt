cmake_minimum_required (VERSION 3.3)
project (nova-nt2-ugens)

enable_testing()

include(sc-plugin-builder.cmake)

file(GLOB headers dsp/*hpp)
file(GLOB headers_2 *hpp)
file(GLOB headers_3 approximations/*hpp)

#########
# nt2

include(ExternalProject)

ExternalProject_Add( NT2Project
  SOURCE_DIR         ${CMAKE_CURRENT_SOURCE_DIR}/nt2
  BINARY_DIR         ${CMAKE_CURRENT_BINARY_DIR}/nt2/build
  CONFIGURE_COMMAND  ${CMAKE_COMMAND} -DNT2_VERSION_STRING="git" -DBOOST_ROOT=${SC_PATH}/external_libraries/boost ${CMAKE_CURRENT_SOURCE_DIR}/nt2
  BUILD_COMMAND      ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/nt2/build --target stage
  INSTALL_COMMAND    ""
)

add_library( NT2 INTERFACE )
target_include_directories( NT2 INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/nt2/build/include )
add_dependencies( NT2 NT2Project )
target_compile_definitions( NT2 INTERFACE BOOST_SIMD_NO_STRICT_ALIASING BOOST_SIMD_NO_NANS )


########################
# we need a recent boost
#
include_directories( ~/boost_1_59_0/ )
include_directories( ~/boost_1_60_0/ )
include_directories( ~/boost_1_61_0/ )
include_directories( ~/boost_1_62_0/ )

if(APPLE)
  # hopefully homebrew provides boost
  include_directories( /usr/local/include )
endif()



#########
# plugins

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

add_scplugin(DiodeLadderFilter  DiodeLadderFilter.cpp)
add_scplugin(NovaAmbisonicUGens NovaAmbisonicUGens.cpp ${headers} ${headers_2} ${headers_3})
add_scplugin(NovaFB             NovaFB.cpp)
add_scplugin(NovaFeedbackAM     NovaFeedbackAM.cpp)
add_scplugin(NovaLeakDC         NovaLeakDC.cpp)
add_scplugin(NovaMultiFilters   NovaMultiFilters.cpp)
add_scplugin(NovaMultiLFO       NovaMultiLFO.cpp)
add_scplugin(NovaOscDPW         NovaOscDPW.cpp)
add_scplugin(NovaSaturator      NovaSaturator.cpp)
add_scplugin(NovaSVF            NovaSVF.cpp)
add_scplugin(NovaTiltFilter     NovaTiltFilter.cpp)
add_scplugin(NovaXFade          NovaXFade.cpp)

foreach( plugin NovaAmbisonicUGens NovaLeakDC NovaMultiFilters NovaSaturator NovaFeedbackAM DiodeLadderFilter NovaXFade NovaSVF NovaFB NovaMultiLFO NovaOscDPW NovaTiltFilter )
  target_link_libraries( ${plugin}           NT2)
  target_link_libraries( ${plugin}_supernova NT2)
endforeach()

add_scclass( NovaNT2UGens
  sc/NovaFeedbackAM.sc
  sc/NovaFB.sc
  sc/NovaMultiFilters.sc
  sc/NovaMultiLFO.sc
  sc/NovaNT2UGens.sc
  sc/NovaOscDPW.sc
  sc/NovaSaturators.sc
  sc/NovaSVF.sc
  sc/NovaTiltFilter.sc
  sc/NovaXFade.sc )

##########
# tests

add_executable( saturator_test tests/saturator_tests.cpp )
target_link_libraries( saturator_test NT2 )
set_target_properties( saturator_test PROPERTIES CXX_STANDARD 14)
target_compile_definitions( saturator_test PUBLIC BOOST_NO_AUTO_PTR )


if( NATIVE )
  target_compile_options( saturator_test PUBLIC -march=native )
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  if(APPLE)
    target_compile_options(saturator_test PUBLIC -Wa,-q) # workaround for homebrew's gcc
  endif()
endif()


add_library( gtest gtest/googletest/src/gtest-all.cc gtest/googletest/src/gtest_main.cc )
target_include_directories( gtest PUBLIC  gtest/googletest/include
                                  PRIVATE gtest/googletest/        )

if(UNIX)
  target_link_libraries( gtest pthread )
endif()

add_executable( tests tests/approximation_tests.cpp )
target_link_libraries( tests gtest NT2 )
set_target_properties( tests PROPERTIES CXX_STANDARD 14)
target_compile_definitions( tests PUBLIC BOOST_NO_AUTO_PTR )


###########
# benchmark

add_executable( benchmark tests/approximation_benchmark.cpp )
target_link_libraries( benchmark NT2 )
set_target_properties( benchmark PROPERTIES CXX_STANDARD 14)
target_compile_definitions( benchmark PUBLIC BOOST_NO_AUTO_PTR )

if( NATIVE )
  target_compile_options( benchmark PUBLIC -march=native )
endif()
