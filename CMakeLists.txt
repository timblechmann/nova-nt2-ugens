cmake_minimum_required (VERSION 2.6)
project (nova-nt2-ugens)

enable_testing()

GET_FILENAME_COMPONENT(SOURCEPARENT "${CMAKE_CURRENT_SOURCE_DIR}" PATH)
find_path(SC_PATH NAMES include/plugin_interface/SC_PlugIn.h
	PATHS "${SOURCEPARENT}"
	PATH_SUFFIXES SuperCollider)

set(SC_FOUND FALSE)
if(IS_DIRECTORY ${SC_PATH})
	set(SC_FOUND TRUE)
endif()

include_directories(${SC_PATH}/include/plugin_interface)
include_directories(${SC_PATH}/include/common)
include_directories(${SC_PATH}/external_libraries/nova-tt)
include_directories(${SC_PATH}/external_libraries/boost)

include_directories(nt2/build/include)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    add_definitions(-stdlib=libc++)
    add_definitions(-std=c++1y)
    set(CMAKE_SHARED_LINKER_FLAGS "-stdlib=libc++ ${CMAKE_SHARED_LINKER_FLAGS}")
endif()

if(${CMAKE_COMPILER_IS_GNUCXX})
    if(APPLE)
      add_definitions(-Wa,-q) # workaround for homebrew's gcc
    endif()

    add_definitions(-std=c++1y)
    add_definitions(-DBOOST_SIMD_NO_STRICT_ALIASING)
    add_definitions(-fvisibility=hidden -ffast-math -fno-finite-math-only)
endif()

SET(CMAKE_SHARED_LIBRARY_PREFIX "")

file(GLOB headers dsp/*hpp)
file(GLOB headers_2 *hpp)

set (sources
  NovaAmbisonicUGens.cpp
  NovaMultiFilters.cpp
  NovaSaturator.cpp
  NovaFeedbackAM.cpp
  DiodeLadderFilter
)

foreach(src ${sources})
  string(REPLACE .cpp "" libname ${src} )

  add_library(${libname} SHARED ${src} ${headers} ${headers_2})

  set (snlib ${libname}_supernova)

  add_library(${snlib} SHARED ${src})

  set_property(TARGET ${snlib}
    APPEND PROPERTY COMPILE_DEFINITIONS SUPERNOVA)

  set(libs ${libname} ${snlib} ${libs})
endforeach()


install(TARGETS ${libs}
        DESTINATION "lib/SuperCollider/plugins/")

install(FILES sc/NovaNT2UGens.sc
        DESTINATION "share/SuperCollider/Extensions/Nova/")
