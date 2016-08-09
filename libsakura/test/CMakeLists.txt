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
find_package(Threads)
find_package(OpenMP)
find_package(FFTW3)
find_package(Log4Cxx)
find_library(pthread pthread)
include(SetArchFlags)

message("GTEST_INCLUDE_DIRS was set to ${GTEST_INCLUDE_DIRS}")
message("GTEST_LIBRARIES was set to ${GTEST_LIBRARIES}")

# Embed into test executables directory paths implicitely searched by the C++ compiler 
foreach(dir ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES})
  set(CMAKE_C_LINK_FLAGS   "${CMAKE_C_LINK_FLAGS} ${CMAKE_EXECUTABLE_RUNTIME_C_FLAG}${dir}")
  set(CMAKE_CXX_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} ${CMAKE_EXECUTABLE_RUNTIME_CXX_FLAG}${dir}")
endforeach(dir ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES})

message("CMAKE_C_LINK_FLAGS  : ${CMAKE_C_LINK_FLAGS}")
message("CMAKE_CXX_LINK_FLAGS: ${CMAKE_CXX_LINK_FLAGS}")

if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS          "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

set(LOG4CXX_LIB "")
if(LOG4CXX_FOUND AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(LOG4CXX_LIB "-llog4cxx")
endif(LOG4CXX_FOUND AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

project(libsakura_test)
enable_testing()

#set(CMAKE_BUILD_TYPE Debug)

set(PROF_MODE "" CACHE STRING "-fprofile-generate or -fprofile-use -fprofile-correction for GCC")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -Wall ${PROF_MODE} -pipe")
set_cxx_flags_from_arch(${SIMD_ARCH})
message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

set(LIBS "-lsakura ${LOG4CXX_LIB} -lpthread")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	set(LIBS "${LIBS} -lc++abi")
endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG   "-g")

message("PROJECT_SOURCE_DIR is ${PROJECT_SOURCE_DIR}")
message("PROJECT_BINARY_DIR is ${PROJECT_BINARY_DIR}")

include_directories(${PROJECT_SOURCE_DIR}/../src ${PROJECT_BINARY_DIR}/../bin "${GTEST_INCLUDE_DIRS}" ${PROJECT_SOURCE_DIR} ${FFTW3_INCLUDE_DIR})

add_custom_target(logConfig
	COMMAND cp -f "${PROJECT_SOURCE_DIR}/libsakura.log4j" "${PROJECT_BINARY_DIR}"
	DEPENDS "${PROJECT_SOURCE_DIR}/libsakura.log4j"
)

add_executable(testBaseline              baseline.cc)
add_executable(testBaselineKS            baseline_ks.cc)
add_executable(testBaselineSK            baseline_sk.cc)
add_executable(testBaselineWK            baseline_wk.cc)
add_executable(testBaselineAll           baseline_all.cc)
add_executable(testBitOperation          bit_operation.cc)
add_executable(testBoolFilterCollection  bool_filter_collection.cc)
add_executable(testConvolution           convolution.cc)
add_executable(testGridding              gridding.cc)
add_executable(testInit                  init.cc)
add_executable(testInterpolationX        interpolation_x.cc)
add_executable(testInterpolationY        interpolation_y.cc)
add_executable(testNormalization         normalization.cc)
add_executable(testNumericOperation      numeric_operation.cc)
add_executable(testStatistics            statistics.cc)
add_executable(testFFT fft.cc)
add_executable(testCInterface c_interface.c)
add_executable(testCreateMaskNearEdge    mask_edge.cc)

add_library(testutil SHARED testutil.cc)

target_link_libraries (testBaseline              gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testBaselineKS            gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testBaselineSK            gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testBaselineWK            gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testBaselineAll           gtest_main gtest sakura          -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testBitOperation          gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testBoolFilterCollection  gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testConvolution           gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testGridding              gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testInit                  gtest_main gtest sakura          -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testInterpolationX        gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testInterpolationY        gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testNormalization         gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testNumericOperation      gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testStatistics            gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testFFT 				     gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testCInterface                             sakura          -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})
target_link_libraries (testCreateMaskNearEdge    gtest_main gtest sakura testutil -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin ${LIBS})

add_dependencies(testInit logConfig)

add_custom_target(testBaselineRun              COMMAND ./testBaseline              DEPENDS testBaseline)
add_custom_target(testBaselineKSRun            COMMAND ./testBaselineKS            DEPENDS testBaselineKS)
add_custom_target(testBaselineSKRun            COMMAND ./testBaselineSK            DEPENDS testBaselineSK)
add_custom_target(testBaselineWKRun            COMMAND ./testBaselineWK            DEPENDS testBaselineWK)
add_custom_target(testBaselineAllRun           COMMAND ./testBaselineAll           DEPENDS testBaselineAll)
add_custom_target(testBitOperationRun          COMMAND ./testBitOperation          DEPENDS testBitOperation)
add_custom_target(testBoolFilterCollectionRun  COMMAND ./testBoolFilterCollection  DEPENDS testBoolFilterCollection)
add_custom_target(testConvolutionRun           COMMAND ./testConvolution           DEPENDS testConvolution)
add_custom_target(testGriddingRun              COMMAND ./testGridding              DEPENDS testGridding)
add_custom_target(testInitRun                  COMMAND ./testInit                  DEPENDS testInit)
add_custom_target(testInterpolationXRun        COMMAND ./testInterpolationX        DEPENDS testInterpolationX)
add_custom_target(testInterpolationYRun        COMMAND ./testInterpolationY        DEPENDS testInterpolationY)
add_custom_target(testNormalizationRun         COMMAND ./testNormalization         DEPENDS testNormalization)
add_custom_target(testNumericOperationRun      COMMAND ./testNumericOperation      DEPENDS testNumericOperation)
add_custom_target(testStatisticsRun            COMMAND ./testStatistics            DEPENDS testStatistics)
add_custom_target(testFFTRun                   COMMAND ./testFFT                   DEPENDS testFFT)
add_custom_target(testNormalizationRunShort    COMMAND ./testNormalization --gtest_filter="-*Performance*" DEPENDS testNormalization)
add_custom_target(testInterpolationXRunShort   COMMAND ./testInterpolationX   --gtest_filter="-*Performance*" DEPENDS testInterpolationX)
add_custom_target(testInterpolationYRunShort   COMMAND ./testInterpolationY   --gtest_filter="-*Performance*" DEPENDS testInterpolationY)
add_custom_target(testConvolutionRunShort      COMMAND ./testConvolution --gtest_filter="-*Performance*" DEPENDS testConvolution)
add_custom_target(testCInterfaceRun            COMMAND ./testCInterface            DEPENDS testCInterface)
add_custom_target(testCreateMaskNearEdgeRun    COMMAND ./testCreateMaskNearEdge    DEPENDS testCreateMaskNearEdge)

# Tests order
# 1 Low-level tests (prerequisistes to other tests, initialization, log, memory management)
# 2 Stable tests (debugging done, no fail)
# 3 Short tests
# 4 New stable tests 
add_custom_target(test
	COMMAND ./testInit
	COMMAND ./testStatistics
	COMMAND ./testBitOperation
	COMMAND ./testBoolFilterCollection
	COMMAND ./testNumericOperation
	COMMAND ./testBaseline
	COMMAND ./testBaselineKS
	COMMAND ./testBaselineSK
	COMMAND ./testBaselineWK
	COMMAND ./testBaselineAll
	COMMAND ./testGridding
	COMMAND ./testInterpolationX
	COMMAND ./testInterpolationY
	COMMAND ./testNormalization
	COMMAND ./testConvolution
	COMMAND ./testFFT
	COMMAND ./testCInterface
	COMMAND ./testCreateMaskNearEdge
    )
