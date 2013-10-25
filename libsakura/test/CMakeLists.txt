# $Id$
#
cmake_minimum_required(VERSION 2.6)
set(CMAKE_Fortran_COMPILER_NAMES gfortran)
set(CMAKE_Fortran_COMPILER gfortran)
enable_language(Fortran)
find_package(Threads)
find_package(OpenMP)
find_package(FFTW3)
find_library(pthread pthread)

message("GTEST_INCLUDE_DIRS was set to ${GTEST_INCLUDE_DIRS}")
message("GTEST_LIBRARIES was set to ${GTEST_LIBRARIES}")
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS          "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

project(libsakura_test)
enable_testing()

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall")
message("default CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -Ofast")
message("CMAKE_CXX_FLAGS_RELEASE was set to ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -Ofast -DNDEBUG -fopt-info-vec-optimized")
set(CMAKE_Fortran_FLAGS_DEBUG   "-g -funroll-all-loops -Ofast -fopt-info-vec-optimized")

message("PROJECT_SOURCE_DIR is ${PROJECT_SOURCE_DIR}")
message("PROJECT_BINARY_DIR is ${PROJECT_BINARY_DIR}")

include_directories(${PROJECT_SOURCE_DIR}/../src ${PROJECT_BINARY_DIR}/../bin "${GTEST_INCLUDE_DIRS}" ${PROJECT_SOURCE_DIR})

add_custom_target(logConfig
	COMMAND cp -f "${PROJECT_SOURCE_DIR}/libsakura.log4j" "${PROJECT_BINARY_DIR}"
	DEPENDS "${PROJECT_SOURCE_DIR}/libsakura.log4j"
)

add_executable(testApplyCalibration      apply_calibration.cc)
add_executable(testBitOperation          bit_operation.cc)
add_executable(testBoolFilterCollection  bool_filter_collection.cc)
add_executable(testConvolution           convolution.cc)
add_executable(testGridding              gridding.cc fgridsd.f)
add_executable(testInit                  init.cc)
add_executable(testInterpolation1D       interpolation_1d.cc)
add_executable(testInterpolationX        interpolation_x.cc)
add_executable(testInterpolationY        interpolation_y.cc)
add_executable(testLogicalOperation      logical_operation.cc)
add_executable(testNumericOperation      numeric_operation.cc)
add_executable(testStatistics            statistics.cc)

target_link_libraries (testApplyCalibration      ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)
target_link_libraries (testBitOperation          ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)
target_link_libraries (testBoolFilterCollection  ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)
target_link_libraries (testConvolution           ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)
target_link_libraries (testGridding              ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread -lgfortran)
target_link_libraries (testInit                  ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)
target_link_libraries (testInterpolation1D       ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)
target_link_libraries (testInterpolationX        ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)
target_link_libraries (testInterpolationY        ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)
target_link_libraries (testLogicalOperation      ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)
target_link_libraries (testNumericOperation      ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)
target_link_libraries (testStatistics            ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)

add_dependencies(testInit logConfig)

add_custom_target(testApplyCalibrationRun      COMMAND ./testApplyCalibration      DEPENDS testApplyCalibration)
add_custom_target(testBitOperationRun          COMMAND ./testBitOperation          DEPENDS testBitOperation)
add_custom_target(testBoolFilterCollectionRun  COMMAND ./testBoolFilterCollection  DEPENDS testBoolFilterCollection)
add_custom_target(testConvolutionRun           COMMAND ./testConvolution           DEPENDS testConvolution)
add_custom_target(testGriddingRun              COMMAND ./testGridding              DEPENDS testGridding)
add_custom_target(testInitRun                  COMMAND ./testInit                  DEPENDS testInit)
add_custom_target(testInterpolation1DRun       COMMAND ./testInterpolation1D       DEPENDS testInterpolation1D)
add_custom_target(testInterpolationXRun        COMMAND ./testInterpolationX        DEPENDS testInterpolationX)
add_custom_target(testInterpolationYRun        COMMAND ./testInterpolationY        DEPENDS testInterpolationY)
add_custom_target(testLogicalOperationRun      COMMAND ./testLogicalOperation      DEPENDS testLogicalOperation)
add_custom_target(testNumericOperationRun      COMMAND ./testNumericOperation      DEPENDS testNumericOperation)
add_custom_target(testStatisticsRun            COMMAND ./testStatistics            DEPENDS testStatistics)
add_custom_target(testApplyCalibrationRunShort COMMAND ./testApplyCalibration --gtest_filter="-*Performance*" DEPENDS testApplyCalibration)

# テストの順序
# 1 より低レベルのもの(他のテストの前提となる、初期化とかログとかメモリー管理とか)
# 2 ステーブル(デバッグ済み)なもの(テストに失敗しないこと)
# 3 短時間で済むもの
# 4 先にステーブルな状態になったもの(早い者勝ち)
add_custom_target(test
	COMMAND ./testInit
	COMMAND ./testStatistics
	COMMAND ./testBitOperation
	COMMAND ./testBoolFilterCollection
	COMMAND ./testLogicalOperation
	COMMAND ./testGridding
	COMMAND ./testInterpolation1D
	COMMAND ./testInterpolationX
	COMMAND ./testInterpolationY
	COMMAND ./testApplyCalibration
	COMMAND ./testNumericOperation
	COMMAND ./testConvolution
	)
