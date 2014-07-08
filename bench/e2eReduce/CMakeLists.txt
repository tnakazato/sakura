# $Id$
#
cmake_minimum_required(VERSION 2.6)
enable_language(C CXX)
find_package(Threads)
find_package(Log4CXX)
find_package(XDispatch)
find_package(CasaCore)
find_package(LibSakura)

if(CASACORE_FOUND)
	include_directories("${CASACORE_INCLUDE_DIR}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CASACORE_CXX_FLAGS} -I${CASACORE_INCLUDE_DIR}/casacore")
endif(CASACORE_FOUND)
if(LIBSAKURA_FOUND)
	include_directories("${LIBSAKURA_INCLUDE_DIR}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBSAKURA_CXX_FLAGS}")
endif(LIBSAKURA_FOUND)
if(LOG4CXX_FOUND)
	include_directories("${LOG4CXX_INCLUDE_DIR}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LOG4CXX_CXX_FLAGS}")
endif(LOG4CXX_FOUND)
if(XDISPATCH_FOUND)
	include_directories("${XDISPATCH_INCLUDE_DIR}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${XDISPATCH_CXX_FLAGS}")
endif(XDISPATCH_FOUND)

project(e2eReduce CXX)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -m64")
message("default CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -O3")
#set(CMAKE_CXX_FLAGS_DEBUG   "-g -O3 -march=corei7-avx")
message("CMAKE_CXX_FLAGS_RELEASE was set to ${CMAKE_CXX_FLAGS_RELEASE}")

message("PROJECT_SOURCE_DIR is ${PROJECT_SOURCE_DIR}")
message("PROJECT_BINARY_DIR is ${PROJECT_BINARY_DIR}")

add_executable(main src/main.cc)

add_custom_target(logConfig
	COMMAND cp -f "${PROJECT_SOURCE_DIR}/config.log4j" "${PROJECT_BINARY_DIR}"
	COMMAND cp -f "${PROJECT_SOURCE_DIR}/libsakura.log4j" "${PROJECT_BINARY_DIR}"
	DEPENDS "${PROJECT_SOURCE_DIR}/config.log4j" "${PROJECT_SOURCE_DIR}/libsakura.log4j"
)

add_custom_target(e2eConfig
	COMMAND cp -f "${PROJECT_SOURCE_DIR}/e2etest.config" "${PROJECT_SOURCE_DIR}/e2eana.config" "${PROJECT_BINARY_DIR}" 
	DEPENDS "${PROJECT_SOURCE_DIR}/e2etest.config"
)

add_dependencies(main logConfig)
add_dependencies(main e2eConfig)

if(CASACORE_FOUND)
	target_link_libraries (main ${CASACORE_EXE_LINKER_FLAGS})
endif(CASACORE_FOUND)
if(LIBSAKURA_FOUND)
	target_link_libraries (main ${LIBSAKURA_EXE_LINKER_FLAGS})
endif(LIBSAKURA_FOUND)
if(XDISPATCH_FOUND)
	target_link_libraries (main ${XDISPATCH_EXE_LINKER_FLAGS})
endif(XDISPATCH_FOUND)
if(LOG4CXX_FOUND)
	target_link_libraries (main ${LOG4CXX_EXE_LINKER_FLAGS})
endif(LOG4CXX_FOUND)

target_link_libraries (main -L${PROJECT_BINARY_DIR}/bin -Wl,-rpath,${PROJECT_BINARY_DIR}/bin -lpthread)

#add_custom_target(run      COMMAND ./main      DEPENDS main)
add_custom_target(run COMMAND ./main DEPENDS main)
add_custom_target(runana COMMAND ./main e2eana.config DEPENDS main)
