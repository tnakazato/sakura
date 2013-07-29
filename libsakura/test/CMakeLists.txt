# $Id$
#
cmake_minimum_required(VERSION 2.6)
set(CMAKE_Fortran_COMPILER_NAMES gfortran)
set(CMAKE_Fortran_COMPILER gfortran)
enable_language(Fortran)
find_package(Threads)
find_package(OpenMP)
find_library(pthread pthread)

message("GTEST_INCLUDE_DIRS was set to ${GTEST_INCLUDE_DIRS}")
message("GTEST_LIBRARIES was set to ${GTEST_LIBRARIES}")
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

project(libsakura_test)
enable_testing()

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "-std=c++11")
message("default CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -ftree-vectorizer-verbose=3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast -ftree-vectorizer-verbose=3")
message("CMAKE_CXX_FLAGS_RELEASE was set to ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -Ofast -DNDEBUG -ftree-vectorizer-verbose=3")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -funroll-all-loops -Ofast -ftree-vectorizer-verbose=3")

message("PROJECT_SOURCE_DIR is ${PROJECT_SOURCE_DIR}")
message("PROJECT_BINARY_DIR is ${PROJECT_BINARY_DIR}")

include_directories(${PROJECT_SOURCE_DIR}/../src ${PROJECT_BINARY_DIR}/../bin "${GTEST_INCLUDE_DIRS}" ${PROJECT_SOURCE_DIR})

add_executable(testInit init.cc)
add_executable(testGridding gridding.cc fgridsd.f)
add_executable(testBitOperation bit_operation.cc)
target_link_libraries (testInit ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)
target_link_libraries (testGridding ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread -lgfortran)
target_link_libraries (testBitOperation ${GTEST_MAIN_LIBRARIES} ${GTEST_LIBRARIES} -L${PROJECT_BINARY_DIR}/../bin -Wl,-rpath,${PROJECT_BINARY_DIR}/../bin -lsakura -lpthread)
add_test(NAME testInit COMMAND ./testInit)
add_test(NAME testGridding COMMAND ./testGridding)
add_test(NAME testBitOperation COMMAND ./testBitOperation)

add_custom_target(test
	COMMAND ./testInit
	COMMAND ./testGridding
	COMMAND ./testBitOperation
	)
