cmake_minimum_required(VERSION 2.6)
find_package(Threads)

project(libsakura)
set(CMAKE_BUILD_TYPE Release)
set(ADDITIONAL_CXXFLAGS_Default "-DARCH_SUFFIX=Default")
set(ADDITIONAL_CXXFLAGS_AfterSandyBridge "-march=corei7-avx -DARCH_SUFFIX=AfterSandyBridge")
set(libsakura_VERSION_MAJOR 0)
set(libsakura_VERSION_MINOR 1)
set(MULTI_SOURCE gridding.cc)
set(UNI_SOURCE concurrent.cc OptimizedImplementationFactory.cc)

configure_file (
  "${PROJECT_SOURCE_DIR}/libsakura/config.h.in"
  "${PROJECT_BINARY_DIR}/libsakura/config.h"
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
set_target_properties(sakura PROPERTIES
	VERSION ${libsakura_VERSION_MAJOR}.${libsakura_VERSION_MINOR}
	SOVERSION ${libsakura_VERSION_MAJOR}
)

target_link_libraries(sakura ${CMAKE_THREAD_LIBS_INIT})
