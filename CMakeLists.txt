cmake_minimum_required (VERSION 3.3)
project (nova-nt2-ugens)

enable_testing()

include(sc-plugin-builder.cmake)

file(GLOB headers dsp/*hpp)
file(GLOB headers_2 *hpp)

#########
# nt2

include(ExternalProject)

ExternalProject_Add( NT2Project
  SOURCE_DIR         ${CMAKE_CURRENT_SOURCE_DIR}/nt2
  BINARY_DIR         ${CMAKE_CURRENT_BINARY_DIR}/nt2/build
  CONFIGURE_COMMAND  cmake -DNT2_VERSION_STRING="git" ${CMAKE_CURRENT_SOURCE_DIR}/nt2
  BUILD_COMMAND      cmake --build ${CMAKE_CURRENT_BINARY_DIR}/nt2/build --target stage
  INSTALL_COMMAND    ""
)

add_library( NT2 INTERFACE )
target_include_directories( NT2 INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/nt2/build/include )
add_dependencies( NT2 NT2Project )
target_compile_definitions( NT2 INTERFACE BOOST_SIMD_NO_STRICT_ALIASING )




#########
# plugins


add_scplugin(NovaAmbisonicUGens NovaAmbisonicUGens.cpp ${headers} ${headers_2})
add_scplugin(NovaLeakDC         NovaLeakDC.cpp)
add_scplugin(NovaMultiFilters   NovaMultiFilters.cpp)
add_scplugin(NovaSaturator      NovaSaturator.cpp)
add_scplugin(NovaFeedbackAM     NovaFeedbackAM.cpp)
add_scplugin(DiodeLadderFilter  DiodeLadderFilter.cpp)

include_directories( ~/boost_1_59_0/ )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

if(APPLE)
  # hopefully homebrew provides boost
  include_directories( /usr/local/include )
endif()
foreach( plugin NovaAmbisonicUGens NovaLeakDC NovaMultiFilters NovaSaturator NovaFeedbackAM DiodeLadderFilter )
  target_link_libraries( ${plugin}           NT2)
  target_link_libraries( ${plugin}_supernova NT2)
endforeach()

add_scclass(NovaNT2UGens sc/NovaNT2UGens.sc)


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
