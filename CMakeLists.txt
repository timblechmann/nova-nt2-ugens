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

#include_directories(nt2/modules/boost/dispatch/include/)
#include_directories(nt2/modules/boost/simd/base/include/)
#include_directories(nt2/modules/boost/simd/config/include/)
#include_directories(nt2/modules/boost/simd/constant/include/)
#include_directories(nt2/modules/boost/simd/reduction/include/)
#include_directories(nt2/modules/boost/simd/sdk/include/)
#include_directories(nt2/modules/boost/simd/swar/include/)

include_directories(nt2/modules/boost/boostification)

add_definitions(-std=c++11)
    add_definitions(-DBOOST_SIMD_NO_STRICT_ALIASING)
if(${CMAKE_COMPILER_IS_GNUCXX})
    add_definitions(-fvisibility=hidden)
endif()

SET(CMAKE_SHARED_LIBRARY_PREFIX "")

file (GLOB headers *hpp)

add_library(NovaMultiFilters SHARED NovaMultiFilters.cpp ${headers})
add_library(NovaSaturator SHARED NovaSaturator.cpp)

add_library(NovaMultiFilters_supernova SHARED NovaMultiFilters.cpp ${headers})
add_library(NovaSaturator_supernova SHARED NovaSaturator.cpp)


set_property(TARGET NovaMultiFilters_supernova NovaSaturator_supernova
             APPEND PROPERTY COMPILE_DEFINITIONS SUPERNOVA)


install(TARGETS NovaMultiFilters NovaSaturator NovaMultiFilters_supernova NovaSaturator_supernova
        DESTINATION "lib/SuperCollider/plugins/")

install(FILES sc/NovaNT2UGens.sc
        DESTINATION "share/SuperCollider/Extensions/Nova/")
