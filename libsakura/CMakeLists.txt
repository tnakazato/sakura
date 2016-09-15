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

cmake_minimum_required(VERSION 2.8)

set(libsakura_VERSION_MAJOR 4)
set(libsakura_VERSION_MINOR 0)

# option for C++ library
# TODO: need to reconsider a way to select standard C++ library 
#       by default, we should follow OS standard on which library 
#       (libstdc++ or libc++) to use while currently we assume 
#       that clang prefer libc++ 
option(USE_LIBCXX
       "set ON to use libc++ for standard C++ library (only effective for clang)"
       ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if(USE_LIBCXX)
    	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    else(USE_LIBCXX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
    endif(USE_LIBCXX)
endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")


find_path(GTEST_INCLUDE_DIRS gtest/gtest.h PATHS ${PROJECT_SOURCE_DIR}/gtest/include)
set(GTEST_LIBRARIES ${PROJECT_BINARY_DIR}/bingtest/libgtest.a)
set(GTEST_MAIN_LIBRARIES ${PROJECT_BINARY_DIR}/bingtest/libgtest_main.a)

# sakura root directory
set(SAKURA_ROOT_DIR ${PROJECT_SOURCE_DIR})

#set(CMAKE_BUILD_TYPE Release)

add_subdirectory(src bin)
add_subdirectory(gtest bingtest)
add_subdirectory(test bintest)
