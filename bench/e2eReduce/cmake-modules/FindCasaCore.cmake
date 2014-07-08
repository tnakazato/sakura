# Here it is assumed that casacore is located in /work/sakura_e2e/casa/
SET(SAKURA_E2E_DIR /work/sakura_e2e/casa)
FIND_PATH(CASACORE_INCLUDE_DIR casacore/casa/Arrays.h PATHS ${SAKURA_E2E_DIR}/include /opt/share/casa/current/release/include)
FIND_LIBRARY(CASACORE_CASA_LIBRARY NAMES casa_casa PATHS ${SAKURA_E2E_DIR}/lib /opt/share/casa/current/release/lib64)
FIND_LIBRARY(CASACORE_TABLE_LIBRARY NAMES casa_tables PATHS ${SAKURA_E2E_DIR}/lib /opt/share/casa/current/release/lib64)
SET(CASACORE_LIBRARY ${CASACORE_CASA_LIBRARY} ${CASACORE_TABLE_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CASACORE DEFAULT_MSG CASACORE_LIBRARY CASACORE_INCLUDE_DIR)

IF(CASACORE_FOUND)
  SET(CASACORE_LIBRARIES ${CASACORE_LIBRARY})
  get_filename_component(CASACORE_LIBRARY_PATH "${CASACORE_CASA_LIBRARY}" PATH)
  SET(CASACORE_EXE_LINKER_FLAGS "-L ${CASACORE_LIBRARY_PATH} -Wl,-rpath,${CASACORE_LIBRARY_PATH} -lcasa_casa -lcasa_tables")
ELSE(CASACORE_FOUND)
  SET(CASACORE_LIBRARIES)
ENDIF(CASACORE_FOUND)

MARK_AS_ADVANCED(CASACORE_INCLUDE_DIR CASACORE_LIBRARIES)
