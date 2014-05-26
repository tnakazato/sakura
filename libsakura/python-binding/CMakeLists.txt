# $Id$
#
cmake_minimum_required(VERSION 2.6)
find_package(Threads)
find_package(OpenMP)
find_package(Eigen3 REQUIRED)
find_package(Log4Cxx)
if(FALSE)
  find_package(PythonLibs REQUIRED)
else(FALSE)
 foreach(VER 3.4 3.3 3.2 3.1 3.0 2.7 2.6)
  find_path(PYTHON_INCLUDE_PATH
    NAMES Python.h
    PATHS
      /usr/include /usr/local/include /opt/share/casa/current/release/include
    PATH_SUFFIXES
      python${VER}
  )
  find_library(PYTHON_LIBRARY
    NAMES python${VER}
    PATHS
      /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib /opt/share/casa/current/release/lib64
    PATH_SUFFIXES
      python${_CURRENT_VERSION}/config
  )
 endforeach(VER)
endif(FALSE)

if(EIGEN3_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${EIGEN3_INCLUDE_DIR}")
endif(EIGEN3_FOUND)

project(libsakura-python)
set(CMAKE_INSTALL_PREFIX /usr/local)

option(BUILD_DOC "Use Doxygen to create the HTML based API documentation" OFF)

message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
#set(CMAKE_BUILD_TYPE RelWithDebInfo)

set(CMAKE_CXX_FLAGS "-std=c++11 -m64 -Wall")
message("default CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
message("default CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "-funroll-all-loops -O3 -DNDEBUG ${VECT_VERB} ${PROF_MODE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -funroll-all-loops -O3 -DNDEBUG ${VECT_VERB} ${PROF_MODE}")
message("CMAKE_CXX_FLAGS_RELEASE was set to ${CMAKE_CXX_FLAGS_RELEASE}")

#set(ADDITIONAL_CXXFLAGS_Default "-mtune=generic -DARCH_SUFFIX=Default -DARCH_DEFAULT=1")
#set(ADDITIONAL_CXXFLAGS_Default "-march=core2 -mtune=generic -DARCH_SUFFIX=Default -DARCH_DEFAULT=1")
set(ADDITIONAL_CXXFLAGS_Default "-march=corei7 -DARCH_SUFFIX=Default -DARCH_DEFAULT=1")

set(SOURCE python-binding.cc)

configure_file (
  "${PROJECT_SOURCE_DIR}/libsakura/config-python.h.in"
  "${PROJECT_BINARY_DIR}/libsakura/config-python.h"
  @ONLY)

include_directories(${LOG4CXX_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}
	${PYTHON_INCLUDE_PATH}
	"${PROJECT_SOURCE_DIR}/../src" "${PROJECT_BINARY_DIR}/../bin"
	"${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}")

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
)

target_link_libraries(sakurapy sakura log4cxx ${PYTHON_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

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
		COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
		COMMENT "Generating the API documentation with Doxygen")
endif(BUILD_DOC)

add_custom_target(copyTest
	COMMAND cp -f "${PROJECT_SOURCE_DIR}/test.*" "${PROJECT_BINARY_DIR}" 
	DEPENDS "${PROJECT_SOURCE_DIR}/test.*"
)

add_custom_target(pytest
	COMMAND "./test.sh"
	WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
	)

add_dependencies(pytest copyTest)

install(TARGETS sakurapy
	LIBRARY DESTINATION lib
	PUBLIC_HEADER DESTINATION include/libsakura)
