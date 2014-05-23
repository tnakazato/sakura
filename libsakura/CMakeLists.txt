cmake_minimum_required(VERSION 2.6)

#project(libsakura)
set(libsakura_VERSION_MAJOR 0)
set(libsakura_VERSION_MINOR 1)

# put Google Test as "gtest" under the directory which includes this file.

find_path(GTEST_INCLUDE_DIRS gtest/gtest.h PATHS ${PROJECT_SOURCE_DIR}/gtest/include)
set(GTEST_LIBRARIES ${PROJECT_BINARY_DIR}/bingtest/libgtest.a)
set(GTEST_MAIN_LIBRARIES ${PROJECT_BINARY_DIR}/bingtest/libgtest_main.a)

# sakura root directory
set(SAKURA_ROOT_DIR ${PROJECT_SOURCE_DIR})

# set CMAKE_BUILD_TYPE with cmake -D option
#set(CMAKE_BUILD_TYPE Release)

add_subdirectory(src bin)
add_subdirectory(gtest bingtest)
add_subdirectory(test bintest)
add_subdirectory(python-binding python-bin)
