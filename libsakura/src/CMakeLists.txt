# $Id$
#
# @SAKURA_LICENSE_HEADER_START@
# Copyright (C) 2013-2016
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

project(libsakura)

set(CMAKE_MODULE_PATH ../cmake-modules CACHE STRING "List of directories to search for CMake modules")

find_package(Threads)
find_package(Eigen3 3.2.0 REQUIRED)
find_package(Log4Cxx)
find_package(FFTW3 3.3.2 REQUIRED)
include(SetArchFlags)

set(CMAKE_INSTALL_PREFIX /usr/local)

set(libsakura_PREFIX sakura)

set(CMAKE_CXX_FLAGS "-I${PROJECT_BINARY_DIR} -I${PROJECT_SOURCE_DIR} ${CMAKE_CXX_FLAGS}")
if(EIGEN3_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${EIGEN3_INCLUDE_DIR}")
endif(EIGEN3_FOUND)
if(FFTW3_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${FFTW3_INCLUDE_DIR}")
endif(FFTW3_FOUND)

set(HAS_LOG4CXX 0)
if(LOG4CXX_FOUND AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(HAS_LOG4CXX 1)
endif(LOG4CXX_FOUND AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

option(SCALAR "Disable auto-vectorization by compiler" OFF)
option(BUILD_DOC "Enable/disable Doxygen generation of Sakura API HTML documentation" ON)

set(SIMD_ARCH "NATIVE" CACHE STRING "SIMD architecture: one of NATIVE SSE4 AVX AVX2" )

message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

set(VECT_VERB "" CACHE STRING "vectorization verbose compile option, such as -fopt-info-vec-optimized -fopt-info-vec-missed -ftree-vectorizer-verbose=2")
set(PROF_MODE "" CACHE STRING "-fprofile-generate or -fprofile-use -fprofile-correction for GCC")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -Wall -Wno-deprecated-register -fPIC -pipe")

message("CMAKE_CXX_COMPILER_ID is ${CMAKE_CXX_COMPILER_ID}")

message("SIMD_ARCH is ${SIMD_ARCH}")
set_cxx_flags_from_arch(${SIMD_ARCH})

if(SCALAR)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARCH_SCALAR=1 -fno-tree-vectorize")
endif(SCALAR)

message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")
message("default CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
message("default CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(UNROLL_OPTION "-funroll-loops")
else()
    set(UNROLL_OPTION "-funroll-all-loops")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "${UNROLL_OPTION} -O3 -DNDEBUG ${VECT_VERB} ${PROF_MODE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${UNROLL_OPTION} -O3 -DNDEBUG ${VECT_VERB} ${PROF_MODE}")
message("CMAKE_CXX_FLAGS_RELEASE was set to ${CMAKE_CXX_FLAGS_RELEASE}")

set(SOURCES baseline.cc bit_operation.cc bool_filter_collection.cc 
	convolution.cc gridding.cc interpolation.cc 
	normalization.cc numeric_operation.cc statistics.cc fft.cc
	gen_util.cc concurrent.cc mask_edge.cc
	)
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
	if (${DOXYGEN_VERSION} MATCHES "^1\\.6.*$")
	    message(STATUS "use template for version < 1.8")
	    set(DOXYFILE_TEMPLATE Doxyfile.in)
	elseif (${DOXYGEN_VERSION} MATCHES "^1\\.8.*$")
	    message(STATUS "use template for version 1.8")
	    set(DOXYFILE_TEMPLATE Doxyfile.in.1_8)
	else()
	    message(STATUS "use template for version < 1.8")
        set(DOXYFILE_TEMPLATE Doxyfile.in)
	endif()

	configure_file (
	  "${PROJECT_SOURCE_DIR}/${DOXYFILE_TEMPLATE}"
	  "${PROJECT_BINARY_DIR}/Doxyfile"
	  @ONLY)
	add_custom_target(apidoc
		mkdir -p "${PROJECT_BINARY_DIR}/docs"
		COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
		COMMAND cp ${PROJECT_SOURCE_DIR}/../INSTALL.txt ${PROJECT_BINARY_DIR}/docs/api/html/
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
		COMMENT "Generating the API documentation with Doxygen")
endif(BUILD_DOC)

if(HAS_LOG4CXX)
    include_directories("${PROJECT_BINARY_DIR}" ${PROJECT_SOURCE_DIR} ${LOG4CXX_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR})
else(HAS_LOG4CXX)
    include_directories("${PROJECT_BINARY_DIR}" ${PROJECT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR})
endif(HAS_LOG4CXX)

add_library(sakura SHARED ${SOURCES})
set(public_headers
	libsakura/sakura.h
	${PROJECT_BINARY_DIR}/libsakura/config.h
	)

set_target_properties(sakura PROPERTIES
	VERSION ${libsakura_VERSION_MAJOR}.${libsakura_VERSION_MINOR}
	SOVERSION ${libsakura_VERSION_MAJOR}
	PUBLIC_HEADER "${public_headers}"
)

set(LIBS ${FFTW3_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

if(HAS_LOG4CXX)
	set(LIBS -llog4cxx ${LIBS})
endif(HAS_LOG4CXX)
target_link_libraries(sakura ${LIBS})

install(TARGETS sakura
	LIBRARY DESTINATION lib
	PUBLIC_HEADER DESTINATION include/libsakura)

if (BUILD_DOC)
    install(DIRECTORY "${PROJECT_BINARY_DIR}/docs/api" DESTINATION share/doc/libsakura)
endif(BUILD_DOC)
