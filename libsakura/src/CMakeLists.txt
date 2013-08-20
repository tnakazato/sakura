# $Id$
#
cmake_minimum_required(VERSION 2.6)
find_package(Threads)
find_package(OpenMP)
find_package(Eigen3 REQUIRED)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}
		${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)
if(EIGEN3_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${EIGEN3_INCLUDE_DIR}")
endif(EIGEN3_FOUND)


project(libsakura)
set(libsakura_VERSION_MAJOR 0)
set(libsakura_VERSION_MINOR 1)
set(libsakura_PREFIX sakura)
set(CMAKE_INSTALL_PREFIX /usr/local)

option(BUILD_DOC "Use Doxygen to create the HTML based API documentation" OFF)
set(DOC_LANG "Japanese" CACHE STRING "Language for the document generation")

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "-std=c++11 -m64")
message("default CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE "-funroll-all-loops -Ofast -DNDEBUG -ftree-vectorizer-verbose=3")
message("CMAKE_CXX_FLAGS_RELEASE was set to ${CMAKE_CXX_FLAGS_RELEASE}")
set(ADDITIONAL_CXXFLAGS_Default "-march=corei7 -DARCH_SUFFIX=Default")
#set(ADDITIONAL_CXXFLAGS_AfterSandyBridge "-march=corei7-avx -DARCH_SUFFIX=AfterSandyBridge")
set(ADDITIONAL_CXXFLAGS_AfterSandyBridge "-march=core2 -DARCH_SUFFIX=AfterSandyBridge")
set(MULTI_SOURCE baseline.cc bit_operation.cc convolution.cc gridding.cc interpolation.cc logical_operation.cc numeric_operation.cc statistics.cc)
set(UNI_SOURCE concurrent.cc optimized_implementation_factory.cc gen_baseline.cc gen_bit_operation.cc gen_convolution.cc gen_gridding.cc gen_interpolation.cc gen_logical_operation.cc gen_numeric_operation.cc gen_statistics.cc gen_util.c)

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

include_directories(${EIGEN3_INCLUDE_DIR} "${PROJECT_BINARY_DIR}" ${PROJECT_SOURCE_DIR})

foreach(_file ${MULTI_SOURCE})
	get_filename_component(_base ${_file} NAME_WE)
	get_filename_component(_ext ${_file} EXT)

	# default
	add_custom_command(
		OUTPUT "${PROJECT_BINARY_DIR}/${_base}Default${_ext}"
		COMMAND ln -sf "${PROJECT_SOURCE_DIR}/${_file}" "${PROJECT_BINARY_DIR}/${_base}Default${_ext}"
	)
	set_property(SOURCE ${_base}Default${_ext} PROPERTY
		COMPILE_FLAGS ${ADDITIONAL_CXXFLAGS_Default}
	)
	set(OPTIMIZED_SOURCES ${OPTIMIZED_SOURCES} ${_base}Default${_ext})

	# after sandybridge
	add_custom_command(
		OUTPUT "${PROJECT_BINARY_DIR}/${_base}AfterSandyBridge${_ext}"
		COMMAND ln -sf "${PROJECT_SOURCE_DIR}/${_file}" "${PROJECT_BINARY_DIR}/${_base}AfterSandyBridge${_ext}"
	)
	set_property(SOURCE ${_base}AfterSandyBridge${_ext} PROPERTY
		COMPILE_FLAGS "${ADDITIONAL_CXXFLAGS_AfterSandyBridge}"
	)
	set(OPTIMIZED_SOURCES ${OPTIMIZED_SOURCES} ${_base}AfterSandyBridge${_ext})
endforeach()
message(${OPTIMIZED_SOURCES})

#find_library(pthread pthread)
add_library(sakura SHARED ${UNI_SOURCE} ${OPTIMIZED_SOURCES})
set(public_headers
	libsakura/concurrent.h
	libsakura/sakura.h
	libsakura/optimized_implementation_factory.h
	libsakura/optimized_implementation_factory_impl.h
	${PROJECT_BINARY_DIR}/libsakura/config.h
	libsakura/interpolation.h
	)

set_target_properties(sakura PROPERTIES
	VERSION ${libsakura_VERSION_MAJOR}.${libsakura_VERSION_MINOR}
	SOVERSION ${libsakura_VERSION_MAJOR}
	PUBLIC_HEADER "${public_headers}"
)

target_link_libraries(sakura ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS sakura
	LIBRARY DESTINATION lib
	PUBLIC_HEADER DESTINATION include/libsakura)
install(DIRECTORY "${PROJECT_BINARY_DIR}/docs/api" DESTINATION share/doc/libsakura)
