cmake_minimum_required (VERSION 3.0)
project (nova-nt2-ugens)

enable_testing()

include(sc-plugin-builder.cmake)

file(GLOB headers dsp/*hpp)
file(GLOB headers_2 *hpp)

set(sources
  NovaAmbisonicUGens.cpp
  NovaMultiFilters.cpp
  NovaSaturator.cpp
  NovaFeedbackAM.cpp
  DiodeLadderFilter
)

include_directories(nt2/build/include)
add_definitions(-DBOOST_SIMD_NO_STRICT_ALIASING)

add_scplugin(NovaAmbisonicUGens NovaAmbisonicUGens.cpp ${headers} ${headers_2})
add_scplugin(NovaMultiFilters   NovaMultiFilters.cpp)
add_scplugin(NovaSaturator      NovaSaturator.cpp)
add_scplugin(NovaFeedbackAM     NovaFeedbackAM.cpp)
add_scplugin(DiodeLadderFilter  DiodeLadderFilter.cpp)

add_scclass(NovaNT2UGens sc/NovaNT2UGens.sc)
