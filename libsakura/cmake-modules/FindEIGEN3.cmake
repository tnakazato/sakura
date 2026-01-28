find_path(EIGEN3_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
      PATHS ${CMAKE_INSTALL_PREFIX}/include /usr/include /opt/eigen/include /opt/eigen
      PATH_SUFFIXES eigen3)

if(EIGEN3_INCLUDE_DIR)
	if(EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h")
		file(STRINGS "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h" EIGEN3_HDR REGEX "^#define[ 	]+EIGEN_[A-Z]+_VERSION[ 	]")

    	string(REGEX REPLACE "^.*EIGEN_WORLD_VERSION[ 	]+([0-9]+).*$" "\\1" EIGEN3_WORLD_VERSION "${EIGEN3_HDR}")
		string(REGEX REPLACE "^.*EIGEN_MAJOR_VERSION[ 	]+([0-9]+).*$" "\\1" EIGEN3_MAJOR_VERSION "${EIGEN3_HDR}")
		string(REGEX REPLACE "^.*EIGEN_MINOR_VERSION[ 	]+([0-9]+).*$" "\\1" EIGEN3_MINOR_VERSION "${EIGEN3_HDR}")
		set(EIGEN3_VERSION_STRING "${EIGEN3_WORLD_VERSION}.${EIGEN3_MAJOR_VERSION}.${EIGEN3_MINOR_VERSION}")
	endif(EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h")

	include(FindPackageHandleStandardArgs)
	find_package_handle_standard_args(EIGEN3
		REQUIRED_VARS EIGEN3_INCLUDE_DIR
		VERSION_VAR EIGEN3_VERSION_STRING)
	mark_as_advanced(EIGEN3_INCLUDE_DIR)
endif(EIGEN3_INCLUDE_DIR)
