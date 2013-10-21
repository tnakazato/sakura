set(xdispatch_install_prefix "/opt/xdispatch")
find_path(XDISPATCH_INCLUDE_DIR xdispatch/dispatch PATHS ${xdispatch_install_prefix}/include)
find_library(XDISPATCH_LIBRARY NAMES xdispatch dispatch PATHS ${xdispatch_install_prefix}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(XDISPATCH DEFAULT_MSG XDISPATCH_LIBRARY XDISPATCH_INCLUDE_DIR)

if(XDISPATCH_FOUND)
  set(XDISPATCH_LIBRARIES ${XDISPATCH_LIBRARY})
  get_filename_component(XDISPATCH_LIBRARY_PATH "${XDISPATCH_LIBRARY}" PATH)
  SET(XDISPATCH_EXE_LINKER_FLAGS "-L ${XDISPATCH_LIBRARY_PATH} -Wl,-rpath,${XDISPATCH_LIBRARY_PATH} -lxdispatch")
else(XDISPATCH_FOUND)
  set(XDISPATCH_LIBRARIES)
endif(XDISPATCH_FOUND)

mark_as_advanced(XDISPATCH_INCLUDE_DIR XDISPATCH_LIBRARIES)
