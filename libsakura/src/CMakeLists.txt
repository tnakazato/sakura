# $Id$
#
cmake_minimum_required(VERSION 2.6)
find_package(Threads)
find_package(OpenMP)
find_package(Eigen3 REQUIRED)
find_package(Log4Cxx)
find_package(FFTW3)

if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}
		${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)
if(EIGEN3_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${EIGEN3_INCLUDE_DIR}")
endif(EIGEN3_FOUND)
if(FFTW3_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${FFTW3_INCLUDE_DIR}")
endif(FFTW3_FOUND)


project(libsakura)
set(libsakura_VERSION_MAJOR 0)
set(libsakura_VERSION_MINOR 1)
set(libsakura_PREFIX sakura)
set(CMAKE_INSTALL_PREFIX /usr/local)

option(SCALAR "Disable auto-vectorization by compiler for AfterSandyBridge" OFF)
option(BUILD_DOC "Use Doxygen to create the HTML based API documentation" OFF)
set(DOC_LANG "Japanese" CACHE STRING "Language for the document generation")
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(VECT_VERB "" CACHE STRING "vectorization verbose compile option, such as -fopt-info-vec-optimized -fopt-info-vec-missed -ftree-vectorizer-verbose=2")

set(CMAKE_CXX_FLAGS "-std=c++11 -m64 -Wall")
message("default CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
message("default CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
#set(CMAKE_CXX_FLAGS_RELEASE "-funroll-all-loops -Ofast -fno-finite-math-only -DNDEBUG ${VECT_VERB}")
set(CMAKE_CXX_FLAGS_RELEASE "-funroll-all-loops -O3 -DNDEBUG ${VECT_VERB}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -funroll-all-loops -O3 -DNDEBUG ${VECT_VERB}")
message("CMAKE_CXX_FLAGS_RELEASE was set to ${CMAKE_CXX_FLAGS_RELEASE}")

#set(ADDITIONAL_CXXFLAGS_Default "-mtune=generic -DARCH_SUFFIX=Default -DARCH_DEFAULT=1")
#set(ADDITIONAL_CXXFLAGS_Default "-march=core2 -mtune=generic -DARCH_SUFFIX=Default -DARCH_DEFAULT=1")
set(ADDITIONAL_CXXFLAGS_Default "-march=corei7 -DARCH_SUFFIX=Default -DARCH_DEFAULT=1")
if(SCALAR)
	set(ADDITIONAL_CXXFLAGS_AfterSandyBridge "-march=corei7-avx -DARCH_SUFFIX=AfterSandyBridge -DARCH_SCALAR=1 -fno-tree-vectorize")
else(SCALAR)
	set(ADDITIONAL_CXXFLAGS_AfterSandyBridge "-march=corei7-avx -DARCH_SUFFIX=AfterSandyBridge -DARCH_AFTER_SANDY_BRIDGE=1")
endif(SCALAR)

set(MULTI_SOURCE bit_operation.cc  bool_filter_collection.cc convolution.cc
	gridding.cc interpolation.cc logical_operation.cc numeric_operation.cc
	baseline.cc statistics.cc apply_calibration.cc)
set(UNI_SOURCE gen_util.cc concurrent.cc optimized_implementation_factory.cc
    gen_apply_calibration.cc
	gen_bit_operation.cc gen_bool_filter_collection.cc gen_convolution.cc
	gen_gridding.cc gen_interpolation.cc gen_logical_operation.cc
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

include_directories(${LOG4CXX_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR} "${PROJECT_BINARY_DIR}" ${PROJECT_SOURCE_DIR})

foreach(_file ${MULTI_SOURCE})
	get_filename_component(_base ${_file} NAME_WE)
	get_filename_component(_ext ${_file} EXT)

	# default
	add_custom_command(
		OUTPUT "${PROJECT_BINARY_DIR}/${_base}Default${_ext}"
		COMMAND ln -sf "${PROJECT_SOURCE_DIR}/${_file}" "${PROJECT_BINARY_DIR}/${_base}Default${_ext}"
	)
	set_property(SOURCE ${_base}Default${_ext} PROPERTY
		COMPILE_FLAGS "${ADDITIONAL_CXXFLAGS_Default} \${VECT_VERB}"
	)
	set(OPTIMIZED_SOURCES ${OPTIMIZED_SOURCES} ${_base}Default${_ext})

	# after sandybridge
	add_custom_command(
		OUTPUT "${PROJECT_BINARY_DIR}/${_base}AfterSandyBridge${_ext}"
		COMMAND ln -sf "${PROJECT_SOURCE_DIR}/${_file}" "${PROJECT_BINARY_DIR}/${_base}AfterSandyBridge${_ext}"
	)
	set_property(SOURCE ${_base}AfterSandyBridge${_ext} PROPERTY
		COMPILE_FLAGS "${ADDITIONAL_CXXFLAGS_AfterSandyBridge} \${VECT_VERB}"
	)
	set(OPTIMIZED_SOURCES ${OPTIMIZED_SOURCES} ${_base}AfterSandyBridge${_ext})
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

target_link_libraries(sakura log4cxx fftw3f ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS sakura
	LIBRARY DESTINATION lib
	PUBLIC_HEADER DESTINATION include/libsakura)
install(DIRECTORY "${PROJECT_BINARY_DIR}/docs/api" DESTINATION share/doc/libsakura)
