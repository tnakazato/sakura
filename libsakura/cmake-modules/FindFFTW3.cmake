# cf. http://www.cmake.org/cmake/help
#This module will be called when find_package was called in CMakeLists.txt
#FFTW3_FOUND this will be set to TRUE when it was found
#FFTW3_INCLUDE_DIRS where is include file
#FFTW3_LIBRARIES where is library
#FFTW3_DEFINITIONS other options

message(STATUS "FFTW3_ROOT_DIR=${FFTW3_ROOT_DIR}")

#to find header file and library of fftw3
find_path(FFTW3_INCLUDE_DIR fftw3.h
	PATHS ${FFTW3_ROOT_DIR}
	PATH_SUFFIXES include
	NO_DEFAULT_PATH)
find_path(FFTW3_INCLUDE_DIR fftw3.h)
message(STATUS "FFTW3_INCLUDE_DIR=${FFTW3_INCLUDE_DIR}")
# support for multiple fftw3 libraries
# set multiple library names to _components, e.g., set(_components fftw3f fftw3)
set(_components fftw3)
foreach(_comp ${_components})
	find_library(${_comp}_LIBRARY ${_comp}
		PATHS ${FFTW3_ROOT_DIR} 
		PATH_SUFFIXES lib64 lib
		NO_DEFAULT_PATH)
	find_library(${_comp}_LIBRARY ${_comp})
	if (${_comp}_LIBRARY)
		list(APPEND FFTW3_LIBRARIES ${${_comp}_LIBRARY}) 
	endif(${_comp}_LIBRARY)
endforeach(_comp ${_components})

#to use FindPackageHandleStandardArgs function 
include(FindPackageHandleStandardArgs)
# if it was succeeded, FFTW3_FOUND will be set to TRUE.
find_package_handle_standard_args(FFTW3 DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDE_DIR)

#to store cache
mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARIES)
