# $Id$
#
# @SAKURA_LICENSE_HEADER_START@
# Copyright (C) 2013-2022
# Inter-University Research Institute Corporation, National Institutes of Natural Sciences
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
cmake_minimum_required(VERSION 2.8.12)
find_package(Threads)
find_package(OpenMP)
find_package(Eigen3 REQUIRED)
find_package(Log4Cxx)

option(BUNDLE_CASA_APP "Bundle python-binding to CASA release package (OSX only)" OFF)
if(APPLE)
  if(BUNDLE_CASA_APP)
    set(WITH_CASA TRUE)
    if(NOT CASA_APP_PATH)
      set(CASA_APP_PATH /Applications/CASA.app)
    endif(NOT CASA_APP_PATH)
    if(NOT EXISTS ${CASA_APP_PATH})
      message(FATAL_ERROR "CASA package was not found at ${CASA_APP_PATH}")
    endif(NOT EXISTS ${CASA_APP_PATH})
    set(PYTHON_FRAMEWORK_PATH ${CASA_APP_PATH}/Contents/Frameworks/Python.framework)
    execute_process(COMMAND ls -F ${PYTHON_FRAMEWORK_PATH}/Versions
                    COMMAND grep -E "^[0-9].[0-9]"
                    COMMAND sed -e "s|\\(.*\\)/.*$|\\1|"
                    COMMAND tr -d '\n'
                    OUTPUT_VARIABLE CASA_APP_PYTHON_VERSION)
    set(PYTHON_ROOTDIR ${PYTHON_FRAMEWORK_PATH}/Versions/${CASA_APP_PYTHON_VERSION})
    message(STATUS TN_DEBUG START)
    message(STATUS PYTHON_FRAMEWORK_PATH=${PYTHON_FRAMEWORK_PATH})
    message(STATUS CASA_APP_PYTHON_VERSION=${CASA_APP_PYTHON_VERSION})
    message(STATUS PYTHON_ROOTDIR=${PYTHON_ROOTDIR})
    message(STATUS TN_DEBUG END)
  endif(BUNDLE_CASA_APP)
endif(APPLE)

# use python3 module to obtain system information
EXECUTE_PROCESS(COMMAND which python3
                COMMAND tr -d '\n'
                OUTPUT_VARIABLE PYTHON3_EXEC)
MESSAGE(STATUS "python3 exec: ${PYTHON3_EXEC}")
if(NOT PYTHON_INCLUDE_PATH)
  # first try python3 to get include path
  EXECUTE_PROCESS(COMMAND ${PYTHON3_EXEC} -c "import sysconfig; print(sysconfig.get_path('include'))"
                  COMMAND tr -d '\n'
                  OUTPUT_VARIABLE _PYTHON_INCLUDE_PATH)
  set(PYTHON_INCLUDE_PATH ${_PYTHON_INCLUDE_PATH} CACHE PATH "Python include path")
  MESSAGE(STATUS "python3 include path ${PYTHON_INCLUDE_PATH}")
endif()

if (NOT NUMPY_INCLUDE_DIR)
  # first try python3 to get include path
  EXECUTE_PROCESS(COMMAND ${PYTHON3_EXEC} -c "import numpy; print(numpy.get_include())"
                  COMMAND tr -d '\n'
                  OUTPUT_VARIABLE _NUMPY_INCLUDE_DIR)
  set(NUMPY_INCLUDE_DIR ${_NUMPY_INCLUDE_DIR} CACHE PATH "Numpy include path")
  MESSAGE(STATUS "numpy include path ${NUMPY_INCLUDE_DIR}")
endif()

if (PYTHON_INCLUDE_PATH AND NUMPY_INCLUDE_DIR)
  EXECUTE_PROCESS(COMMAND ${PYTHON3_EXEC} -c "import sysconfig; print(sysconfig.get_python_version())"
                  COMMAND tr -d '\n'
                  OUTPUT_VARIABLE VER)
  set(PYTHON_VERSION "${VER}" CACHE STRING "Python version")

  EXECUTE_PROCESS(COMMAND echo ${PYTHON_INCLUDE_PATH}
                  COMMAND sed -e "s|\\(.*\\)/include.*|\\1|g"
                  COMMAND tr -d '\n'
                  OUTPUT_VARIABLE _PYTHON_ROOTDIR)
  FIND_LIBRARY(PYTHON_LIBRARY
               NAMES python${VER} python${VER}m
               PATHS ${_PYTHON_ROOTDIR}/lib ${_PYTHON_ROOTDIR}/lib64
               NO_DEFAULT_PATH)
  MESSAGE(STATUS "python library path is ${PYTHON_LIBRARY}")
endif()

if(FALSE)
  find_package(PythonLibs REQUIRED)
else(FALSE)
 message(STATUS "PYTHON and NUMPY search")
 foreach(VER 2.7 2.6)
  if (PYTHON_ROOTDIR)
   message("PYTHON_ROOTDIR is specified: ${PYTHON_ROOTDIR}")
   find_path(PYTHON_INCLUDE_PATH
             NAMES Python.h
             PATHS ${PYTHON_ROOTDIR}/include
             PATH_SUFFIXES python${VER} python${VER}m
             NO_DEFAULT_PATH)
  endif(PYTHON_ROOTDIR)
  find_path(PYTHON_INCLUDE_PATH
    NAMES Python.h
    PATHS
      /usr/include /usr/local/include /opt/share/casa/current/release/include
    PATH_SUFFIXES
      python${VER} python${VER}m
  )
  if(PYTHON_ROOTDIR)
   find_library(PYTHON_LIBRARY
                NAMES python${VER} python${VER}m
                PATHS ${PYTHON_ROOTDIR}/lib
                NO_DEFAULT_PATH)
  endif(PYTHON_ROOTDIR)

  find_library(PYTHON_LIBRARY
    NAMES python${VER} python${VER}m
    PATHS
      /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib /opt/share/casa/current/release/lib
    PATH_SUFFIXES
      python${_CURRENT_VERSION}/config
  )

  # numpy include path
  message(STATUS PYTHON_INCLUDE_PATH=${PYTHON_INCLUDE_PATH})
  if(PYTHON_INCLUDE_PATH)
    #message("PYTHON_INCLUDE_PATH is ${PYTHON_INCLUDE_PATH}")
    #message("PYTHON_LIBRARY is ${PYTHON_LIBRARY}")
    execute_process(COMMAND dirname ${PYTHON_LIBRARY}
                    COMMAND tr -d '\n'
                    OUTPUT_VARIABLE PYTHON_LIBRARY_DIR)
    #message("PYTHON_LIBRARY_DIR is ${PYTHON_LIBRARY_DIR}")
    find_path(NUMPY_INCLUDE_DIR
              NAMES numpy/arrayobject.h
              PATHS ${PYTHON_LIBRARY_DIR}/python${VER}/site-packages/numpy/core/include
                    ${PYTHON_LIBRARY_DIR}/python${VER}/dist-packages/numpy/core/include
              NO_DEFAULT_PATH)
    #message("NUMPY_INCLUDE_DIR is ${NUMPY_INCLUDE_DIR}")
    find_path(NUMPY_INCLUDE_DIR
              NAMES numpy/arrayobject.h
              PATH_PREFIXES python${VER})
    #message("NUMPY_INCLUDE_DIR is ${NUMPY_INCLUDE_DIR}")
  endif(PYTHON_INCLUDE_PATH)

  if (PYTHON_INCLUDE_PATH AND NUMPY_INCLUDE_DIR)
    set(PYTHON_VERSION "${VER}" CACHE STRING "Python version")
    break()
  endif()
 endforeach(VER)
endif(FALSE)

message(STATUS "PYTHON_VERSION = ${PYTHON_VERSION}")

if(EIGEN3_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${EIGEN3_INCLUDE_DIR}")
endif(EIGEN3_FOUND)

project(libsakura-python)
set(CMAKE_INSTALL_PREFIX /usr/local)

option(BUILD_DOC "Use Doxygen to create the HTML based API documentation" OFF)

message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
#set(CMAKE_BUILD_TYPE RelWithDebInfo)

set(CMAKE_CXX_FLAGS "-std=c++11 -m64 -Wall -fPIC")
if(APPLE)
    set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ -Wall")
endif(APPLE)
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

#set(ADDITIONAL_CXXFLAGS_Default "-mtune=generic -DARCH_SUFFIX=Default -DARCH_DEFAULT=1")
#set(ADDITIONAL_CXXFLAGS_Default "-march=core2 -mtune=generic -DARCH_SUFFIX=Default -DARCH_DEFAULT=1")
set(ADDITIONAL_CXXFLAGS_Default "-march=corei7 -DARCH_SUFFIX=Default -DARCH_DEFAULT=1")

set(SOURCE python-binding.cc)

configure_file (
  "${PROJECT_SOURCE_DIR}/libsakura/config-python.h.in"
  "${PROJECT_BINARY_DIR}/libsakura/config-python.h"
  @ONLY)

set(INCLUDE_DIR ${EIGEN3_INCLUDE_DIR}
	    ${PYTHON_INCLUDE_PATH}
	    "${PROJECT_SOURCE_DIR}/../src" "${PROJECT_BINARY_DIR}/../bin"
	    "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}")
if(LOG4CXX_FOUND)
    set(INCLUDE_DIR ${LOG4CXX_INCLUDE_DIR} "${INCLUDE_DIR}")
endif(LOG4CXX_FOUND)

if (NUMPY_INCLUDE_DIR)
    set(INCLUDE_DIR "${INCLUDE_DIR}" ${NUMPY_INCLUDE_DIR})
endif(NUMPY_INCLUDE_DIR)

include_directories(${INCLUDE_DIR})

#find_library(pthread pthread)
add_library(sakurapy SHARED ${SOURCE})
set(public_headers
	libsakura/sakura-python.h
	${PROJECT_BINARY_DIR}/libsakura/config-python.h
	)

set(libsakura_VERSION_MAJOR 0)
set(libsakura_VERSION_MINOR 1)

set_target_properties(sakurapy PROPERTIES
	VERSION ${libsakura_VERSION_MAJOR}.${libsakura_VERSION_MINOR}
	SOVERSION ${libsakura_VERSION_MAJOR}
	PUBLIC_HEADER "${public_headers}"
	SUFFIX .so
)

# test
if(WITH_CASA)
  message(STATUS "Replace PYTHON_LIBRARY for framework with ordinary library link")
  set(PYTHON_LIBRARY -lpython${CASA_APP_PYTHON_VERSION})
endif(WITH_CASA)

# no need to link Python library except macOS
if (APPLE)
  if(LOG4CXX_FOUND)
      target_link_libraries(sakurapy sakura log4cxx ${PYTHON_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
  else(LOG4CXX_FOUND)
      target_link_libraries(sakurapy sakura ${PYTHON_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
  endif(LOG4CXX_FOUND)
else()
  message(STATUS "PYTHON Version ${PYTHON_VERSION}: not linking libpython")
  if(LOG4CXX_FOUND)
      target_link_libraries(sakurapy sakura log4cxx ${CMAKE_THREAD_LIBS_INIT})
  else(LOG4CXX_FOUND)
      target_link_libraries(sakurapy sakura ${CMAKE_THREAD_LIBS_INIT})
  endif(LOG4CXX_FOUND)
endif()

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
	add_custom_target(pyapidoc
		mkdir -p "${PROJECT_BINARY_DIR}/docs"
		COMMAND rm -f "${PROJECT_BINARY_DIR}/libsakura/config.h"
		COMMAND ln -s "${PROJECT_BINARY_DIR}/../bin/libsakura/config.h" "${PROJECT_BINARY_DIR}/libsakura"
		COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
		COMMAND rm -f "${PROJECT_BINARY_DIR}/libsakura/config.h"
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
		COMMENT "Generating the API documentation with Doxygen")
endif(BUILD_DOC)

add_custom_target(copyTest
	COMMAND cp -f "${PROJECT_SOURCE_DIR}/test.*" "${PROJECT_BINARY_DIR}"
	DEPENDS "${PROJECT_SOURCE_DIR}/test.*"
)

add_custom_target(pytest
	COMMAND ./test.sh test.py
	WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
	)

add_dependencies(pytest copyTest)

install(TARGETS sakurapy
	LIBRARY DESTINATION lib
	PUBLIC_HEADER DESTINATION include/libsakura)

if(WITH_CASA)
  message(STATUS "edit dependency on Python")
  add_custom_command(TARGET sakurapy
                     POST_BUILD
                     COMMAND ${PROJECT_SOURCE_DIR}/../bin/osx_python_to_casa.sh )
endif(WITH_CASA)
