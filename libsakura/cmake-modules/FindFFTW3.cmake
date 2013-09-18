# cf. http://www.cmake.org/cmake/help
#This module will be called when find_package was called in CMakeLists.txt
#FFTW3_FOUND this will be set to TRUE when it was found
#FFTW3_INCLUDE_DIRS where is include file
#FFTW3_LIBRARIES where is library
#FFTW3_DEFINITIONS other options

#to find header file and library of fftw3
find_path(FFTW3_INCLUDE_DIR fftw3.h)
find_library(FFTW3_LIBRARIES NAMES fftw3)

#to use FindPackageHandleStandardArgs function 
include(FindPackageHandleStandardArgs)
# if it was succeeded, FFTW3_FOUND will be set to TRUE.
find_package_handle_standard_args(FFTW3 DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDE_DIR)

#to store cache
mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARIES)
