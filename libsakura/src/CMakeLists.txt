# $Id$
#
cmake_minimum_required(VERSION 2.6)
find_package(Threads)
find_package(OpenMP)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}
		${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)


project(libsakura)
set(libsakura_VERSION_MAJOR 0)
set(libsakura_VERSION_MINOR 1)
set(CMAKE_INSTALL_PREFIX /usr/local)

set(CMAKE_BUILD_TYPE Release)
set(ADDITIONAL_CXXFLAGS_Default "-DARCH_SUFFIX=Default")
set(ADDITIONAL_CXXFLAGS_AfterSandyBridge "-march=corei7-avx -DARCH_SUFFIX=AfterSandyBridge")
set(MULTI_SOURCE gridding.cc)
set(UNI_SOURCE concurrent.cc OptimizedImplementationFactory.cc)

configure_file (
  "${PROJECT_SOURCE_DIR}/libsakura/config.h.in"
  "${PROJECT_BINARY_DIR}/libsakura/config.h"
  @ONLY
  )

include_directories("${PROJECT_BINARY_DIR}" ${PROJECT_SOURCE_DIR})

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
	libsakura/OptimizedImplementationFactory.h
	libsakura/OptimizedImplementationFactoryImpl.h
	${PROJECT_BINARY_DIR}/libsakura/config.h
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
