# $Id$
#
# @SAKURA_LICENSE_HEADER_START@
# Copyright (C) 2013-2014
# National Astronomical Observatory of Japan
# 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
# 
# This file is part of Sakura.
# 
# Sakura is free software: you can redistribute it and/or modify it under 
# the terms of the GNU Lesser General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or (at your 
# option) any later version.
# 
# Sakura is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License 
# along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
# @SAKURA_LICENSE_HEADER_END@
#
cmake_minimum_required(VERSION 2.6)
find_package(Threads)
find_package(Eigen3 REQUIRED)
find_package(Log4Cxx)
find_package(FFTW3)

if(EIGEN3_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${EIGEN3_INCLUDE_DIR}")
endif(EIGEN3_FOUND)
if(FFTW3_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${FFTW3_INCLUDE_DIR}")
endif(FFTW3_FOUND)
if(LOG4CXX_FOUND)
  set(HAS_LOG4CXX 1)
else(LOG4CXX_FOUND)
  set(HAS_LOG4CXX 0)
endif(LOG4CXX_FOUND)

project(libsakura)
set(libsakura_VERSION_MAJOR 0)
set(libsakura_VERSION_MINOR 1)
set(libsakura_PREFIX sakura)
set(CMAKE_INSTALL_PREFIX /usr/local)

option(SCALAR "Disable auto-vectorization by compiler" OFF)
option(BUILD_DOC "Use Doxygen to create the HTML based API documentation" OFF)
set(DOC_LANG "Japanese" CACHE STRING "Language for the document generation")
message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
#set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(VECT_VERB "" CACHE STRING "vectorization verbose compile option, such as -fopt-info-vec-optimized -fopt-info-vec-missed -ftree-vectorizer-verbose=2")
set(PROF_MODE "" CACHE STRING "-fprofile-generate or -fprofile-use -fprofile-correction for GCC")

set(DefaultArch "")
set(SandyBridgeArch "")
set(HaswellArch "")

message("CMAKE_CXX_COMPILER_ID is ${CMAKE_CXX_COMPILER_ID}")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	set(DefaultArch "-mtune=generic")
	set(SandyBridgeArch "-march=corei7-avx")
	set(HaswellArch "-march=core-avx2 -ffast-math -mfma") # -ffast-math is required to enable FMA
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(DefaultArch "-mtune=generic")
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9")
		set(SandyBridgeArch "-march=corei7-avx")
		set(HaswellArch "-march=core-avx2")
	else(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9")
		set(SandyBridgeArch "-march=sandybridge")
		set(HaswellArch "-march=haswell")
	endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	set(DefaultArch "-mkl=sequential -xSSE4.2")
	set(SandyBridgeArch "-mkl=sequential -xAVX ")
	set(HaswellArch "-mkl=sequential -xCORE-AVX2 -mtune=core-avx2")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

set(CMAKE_CXX_FLAGS "-std=c++11 -m64 -Wall -Wno-deprecated-register")
message("default CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
message("default CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "-funroll-all-loops -O3 -DNDEBUG ${VECT_VERB} ${PROF_MODE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -funroll-all-loops -O3 -DNDEBUG ${VECT_VERB} ${PROF_MODE}")
message("CMAKE_CXX_FLAGS_RELEASE was set to ${CMAKE_CXX_FLAGS_RELEASE}")

#set(ADDITIONAL_CXXFLAGS_Default "${DefaultArch} -DARCH_SUFFIX=Default -DARCH_DEFAULT=1")
#set(ADDITIONAL_CXXFLAGS_Default "-march=core2 ${DefaultArch} -DARCH_SUFFIX=Default -DARCH_DEFAULT=1")
if(SCALAR)
	set(ADDITIONAL_CXXFLAGS_Default "${DefaultArch} -DARCH_SUFFIX=Default -DARCH_SCALAR=1 -fno-tree-vectorize")
	set(ADDITIONAL_CXXFLAGS_AfterSandyBridge "${SandyBridgeArch} -DARCH_SUFFIX=AfterSandyBridge -DARCH_SCALAR=1 -fno-tree-vectorize")
	set(ADDITIONAL_CXXFLAGS_AfterHaswell "${HaswellArch} -DARCH_SUFFIX=AfterHaswell -DARCH_SCALAR=1 -fno-tree-vectorize")
else(SCALAR)
	set(ADDITIONAL_CXXFLAGS_Default "${DefaultArch} -DARCH_SUFFIX=Default -DARCH_DEFAULT=1")
	set(ADDITIONAL_CXXFLAGS_AfterSandyBridge "${SandyBridgeArch} -DARCH_SUFFIX=AfterSandyBridge -DARCH_AFTER_SANDY_BRIDGE=1")
	set(ADDITIONAL_CXXFLAGS_AfterHaswell "${HaswellArch} -DARCH_SUFFIX=AfterHaswell -DARCH_AFTER_SANDY_BRIDGE=1 -DARCH_AFTER_HASWELL=1")
endif(SCALAR)

set(MULTI_SOURCE apply_calibration.cc baseline.cc bit_operation.cc 
	bool_filter_collection.cc convolution.cc gridding.cc interpolation.cc 
	numeric_operation.cc statistics.cc)
set(UNI_SOURCE gen_util.cc concurrent.cc optimized_implementation_factory.cc
    gen_apply_calibration.cc gen_bit_operation.cc gen_bool_filter_collection.cc 
    gen_convolution.cc gen_gridding.cc gen_interpolation.cc 
    gen_numeric_operation.cc gen_baseline.cc gen_statistics.cc)

configure_file (
  "${PROJECT_SOURCE_DIR}/libsakura/config.h.in"
  "${PROJECT_BINARY_DIR}/libsakura/config.h"
  @ONLY)

if(BUILD_DOC)
	find_package(Doxygen)
	if (NOT DOXYGEN_FOUND)
		message(FATAL_ERROR 
			"Doxygen is needed to build the API documentation.")
	endif(NOT DOXYGEN_FOUND)

	configure_file (
	  "${PROJECT_SOURCE_DIR}/Doxyfile.in"
	  "${PROJECT_BINARY_DIR}/Doxyfile"
	  @ONLY)
	add_custom_target(apidoc
		mkdir -p "${PROJECT_BINARY_DIR}/docs"
		COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
		COMMENT "Generating the API documentation with Doxygen")
endif(BUILD_DOC)

if(HAS_LOG4CXX)
    include_directories(${LOG4CXX_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR} "${PROJECT_BINARY_DIR}" ${PROJECT_SOURCE_DIR})
else(HAS_LOG4CXX)
    include_directories(${EIGEN3_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR} "${PROJECT_BINARY_DIR}" ${PROJECT_SOURCE_DIR})
endif(HAS_LOG4CXX)
    
foreach(_file ${MULTI_SOURCE})
	get_filename_component(_base ${_file} NAME_WE)
	get_filename_component(_ext ${_file} EXT)

	# default
	add_custom_command(
		OUTPUT "${PROJECT_BINARY_DIR}/${_base}Default${_ext}"
		COMMAND ln -sf "${PROJECT_SOURCE_DIR}/${_file}" "${PROJECT_BINARY_DIR}/${_base}Default${_ext}"
	)
	set_property(SOURCE ${_base}Default${_ext} PROPERTY
		COMPILE_FLAGS "${ADDITIONAL_CXXFLAGS_Default} \${VECT_VERB} ${PROF_MODE}"
	)
	set(OPTIMIZED_SOURCES ${OPTIMIZED_SOURCES} ${_base}Default${_ext})

	# after sandybridge
	add_custom_command(
		OUTPUT "${PROJECT_BINARY_DIR}/${_base}AfterSandyBridge${_ext}"
		COMMAND ln -sf "${PROJECT_SOURCE_DIR}/${_file}" "${PROJECT_BINARY_DIR}/${_base}AfterSandyBridge${_ext}"
	)
	set_property(SOURCE ${_base}AfterSandyBridge${_ext} PROPERTY
		COMPILE_FLAGS "${ADDITIONAL_CXXFLAGS_AfterSandyBridge} \${VECT_VERB} ${PROF_MODE}"
	)
	set(OPTIMIZED_SOURCES ${OPTIMIZED_SOURCES} ${_base}AfterSandyBridge${_ext})

	# after haswell
	add_custom_command(
		OUTPUT "${PROJECT_BINARY_DIR}/${_base}AfterHaswell${_ext}"
		COMMAND ln -sf "${PROJECT_SOURCE_DIR}/${_file}" "${PROJECT_BINARY_DIR}/${_base}AfterHaswell${_ext}"
	)
	set_property(SOURCE ${_base}AfterHaswell${_ext} PROPERTY
		COMPILE_FLAGS "${ADDITIONAL_CXXFLAGS_AfterHaswell} \${VECT_VERB} ${PROF_MODE}"
	)
	set(OPTIMIZED_SOURCES ${OPTIMIZED_SOURCES} ${_base}AfterHaswell${_ext})
endforeach()
message(${OPTIMIZED_SOURCES})

#find_library(pthread pthread)
add_library(sakura SHARED ${UNI_SOURCE} ${OPTIMIZED_SOURCES})
set(public_headers
	libsakura/sakura.h
	${PROJECT_BINARY_DIR}/libsakura/config.h
	)

set_target_properties(sakura PROPERTIES
	VERSION ${libsakura_VERSION_MAJOR}.${libsakura_VERSION_MINOR}
	SOVERSION ${libsakura_VERSION_MAJOR}
	PUBLIC_HEADER "${public_headers}"
)

if(HAS_LOG4CXX)
    target_link_libraries(sakura log4cxx fftw3f ${CMAKE_THREAD_LIBS_INIT})
else(HAS_LOG4CXX)
    target_link_libraries(sakura fftw3f ${CMAKE_THREAD_LIBS_INIT})
endif(HAS_LOG4CXX)

install(TARGETS sakura
	LIBRARY DESTINATION lib
	PUBLIC_HEADER DESTINATION include/libsakura)
install(DIRECTORY "${PROJECT_BINARY_DIR}/docs/api" DESTINATION share/doc/libsakura)
