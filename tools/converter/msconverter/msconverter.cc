#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <cassert>
#include <cstring>
#include <string>
#include <memory>
#include <iostream>

#include <casa/aips.h>
#include <casa/Inputs/Input.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/OS/Directory.h>
#include <tables/Tables/Table.h>
#include <ms/MeasurementSets.h>
#include "SQLite.h"

#define unique_ptr auto_ptr
#define enter() do { cout << "Enter: " << __FUNCTION__ << endl; } while (0)

using namespace std;
using namespace casa;
using namespace sqlite;

namespace {
union {
  int i;
  float f;
  double d;
  void *p;
  int64_t i64;
} DUMMY_AREA[1];

template<typename T>
void bindArrayAsBlob(PreparedStatement *stmt, int pos, Array<T> const &v) throw (SQLException) {
  Bool deleteIt = false;
  T const *data = v.getStorage(deleteIt);
  try {
    size_t elements = v.nelements();
    cout << "data: " << data <<", " << elements << endl;;
    if (data == NULL) {
      assert(elements == 0);
      data = (T const *)DUMMY_AREA;
    }
    stmt->setTransientBlob(pos, data, sizeof(T) * elements);
  } catch (...) {
    v.freeStorage(data, deleteIt);
    throw;
  }
  v.freeStorage(data, deleteIt);
}

template<>
void bindArrayAsBlob<Complex>(PreparedStatement *stmt, int pos, Array<Complex> const &v) throw (SQLException) {
  size_t elements = v.nelements();
  float *const blob = new float[elements * 2];
  float *real = blob;
  float *imag = &blob[elements];
  try {
    Array<Complex>::const_iterator end = v.end();
    size_t i = 0;
    for (Array<Complex>::const_iterator it = v.begin(); it != end; ++it, ++i) {
      real[i] = it->real();
      imag[i] = it->imag();
    }
    stmt->setTransientBlob(pos, blob, 2 * sizeof(float) * elements);
  } catch (...) {
    delete [] blob;
    throw;
  }
  delete [] blob;
}

void savePol(Connection *con, MeasurementSet &ms) {
  enter();
  MSPolarization &pol = ms.polarization();
  ROMSPolarizationColumns const cols(pol);

  char const *dbcols[] = {
    "POLARIZATION_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "NUM_CORR", // INTEGER NOT NULL,
    "CORR_TYPE", // BLOB NOT NULL,
    "CORR_PRODUCT", // BLOB NOT NULL,
    "FLAG_ROW", // INTEGER NOT NULL
    NULL
  };

  string sql = "insert into POLARIZATION (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = pol.nrow();
  for (uInt i = 0; i < nrow; i++) {
    int pos = 0;
    stmt->setInt(++pos, i); // POLARIZATION_ID
    stmt->setInt(++pos, cols.numCorr()(i)); // NUM_CORR
    bindArrayAsBlob<Int>(stmt.get(), ++pos, cols.corrType()(i));
    bindArrayAsBlob<Int>(stmt.get(), ++pos, cols.corrProduct()(i));
    stmt->setInt(++pos, cols.flagRow()(i)); // FLAG_ROW
    assert(pos == stmt->getParameterCount());
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
}

void saveSw(Connection *con, MeasurementSet &ms) {
  enter();
  MSSpectralWindow &tab = ms.spectralWindow();
  ROMSSpWindowColumns const cols(tab);

  char const *dbcols[] = {
    "SPECTRAL_WINDOW_ID", // INTEGER NOT NULL,
    "NUM_CHAN", // INTEGER NOT NULL,
    "NAME", // TEXT NOT NULL,
    "REF_FREQUENCY", // REAL NOT NULL,
    "CHAN_FREQ", // BLOB NOT NULL,
    "CHAN_WIDTH", // BLOB NOT NULL,
    "MEAS_FREQ_REF", // INTEGER NOT NULL,
    "EFFECTIVE_BW", // BLOB NOT NULL,
    "RESOLUTION", // BLOB NOT NULL,
    "TOTAL_BANDWIDTH", // REAL NOT NULL,
    "NET_SIDEBAND", // INTEGER NOT NULL,
    "BBC_NO", // INTEGER DEFAULT NULL,
    "BBC_SIDEBAND", // INTEGER DEFAULT NULL,
    "IF_CONV_CHAIN", // INTEGER NOT NULL,
    "RECEIVER_ID", // INTEGER DEFAULT NULL,
    "FREQ_GROUP", // INTEGER NOT NULL,
    "FREQ_GROUP_NAME", // TEXT NOT NULL,
    "DOPPLER_ID", // INTEGER DEFAULT NULL,
    "FLAG_ROW", // INTEGER NOT NULL
    NULL
  };

  string sql = "insert into SPECTRAL_WINDOW (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    int pos = 0;
    stmt->setInt(++pos, i); // SPECTRAL_WINDOW_ID
    stmt->setInt(++pos, cols.numChan()(i)); // NUM_CHAN
    stmt->setTransientString(++pos, cols.name()(i).c_str()); // NAME
    stmt->setDouble(++pos, cols.refFrequency()(i)); // REF_FREQUENCY
    bindArrayAsBlob<Double>(stmt.get(), ++pos, cols.chanFreq()(i));
    bindArrayAsBlob<Double>(stmt.get(), ++pos, cols.chanWidth()(i));
    stmt->setInt(++pos, cols.measFreqRef()(i)); // MEAS_FREQ_REF
    bindArrayAsBlob<Double>(stmt.get(), ++pos, cols.effectiveBW()(i));
    bindArrayAsBlob<Double>(stmt.get(), ++pos, cols.resolution()(i));
    stmt->setDouble(++pos, cols.totalBandwidth()(i)); // TOTAL_BANDWIDTH
    stmt->setInt(++pos, cols.netSideband()(i)); // NET_SIDEBAND

    if (cols.bbcNo().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.bbcNo()(i)); //BBC_NO
    }
    if (cols.bbcSideband().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.bbcSideband()(i)); // BBC_SIDEBAND
    }
    stmt->setInt(++pos, cols.ifConvChain()(i)); // IF_CONV_CHAIN
    if (cols.receiverId().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.receiverId()(i)); // RECEIVER_ID
    }
    stmt->setInt(++pos, cols.freqGroup()(i)); // FREQ_GROUP
    stmt->setTransientString(++pos, cols.freqGroupName()(i).c_str()); // FREQ_GROUP_NAME
    if (cols.dopplerId().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.dopplerId()(i)); // DOPPLER_ID
    }
    stmt->setInt(++pos, cols.flagRow()(i)); // FLAG_ROW
    assert(pos == stmt->getParameterCount());
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
}

void saveDataDesc(Connection *con, MeasurementSet &ms) {
  enter();
  MSDataDescription &tab = ms.dataDescription();
  ROMSDataDescColumns const cols(tab);

  char const *dbcols[] = {
    "DATA_DESC_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "SPECTRAL_WINDOW_ID", // INTEGER NOT NULL,
    "POLARIZATION_ID", //INTEGER NOT NULL,
    "LAG_ID", //INTEGER DEFAULT NULL,
    "FLAG_ROW", //INTEGER NOT NULL,
    NULL
  };

  string sql = "insert into DATA_DESCRIPTION (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    int pos = 0;
    stmt->setInt(++pos, i); // DATA_DESC_ID
    stmt->setInt(++pos, cols.spectralWindowId()(i)); // SPECTRAL_WINDOW_ID
    stmt->setInt(++pos, cols.polarizationId()(i)); // POLARIZATION_ID
    if (cols.lagId().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.lagId()(i)); // LAG_ID
    }
    stmt->setInt(++pos, cols.flagRow()(i)); // FLAG_ROW
    assert(pos == stmt->getParameterCount());
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
}

void saveFreqOffset(Connection *con, MeasurementSet &ms) {
  enter();
  MSFreqOffset &tab = ms.freqOffset();
  ROMSFreqOffsetColumns const cols(tab);
  
  // return if optional FREQ_OFFSET table doesn't exist
  if (tab.isNull()) return;

  char const *dbcols[] = {
    "FREQ_OFFSET_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "ANTENNA1", // INTEGER NOT NULL,
    "ANTENNA2", // INTEGER NOT NULL,
    "FEED_ID", //INTEGER NOT NULL,
    "SPECTRAL_WINDOW_ID", //INTEGER DEFAULT NULL,
    "TIME", //REAL NOT NULL,
    "INTERVAL", //REAL NOT NULL,
    "OFFSET", //REAL NOT NULL,
    NULL
  };

  string sql = "insert into FREQ_OFFSET (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    int pos = 0;
    stmt->setInt(++pos, i); // FREQ_OFFSET_ID
    stmt->setInt(++pos, cols.antenna1()(i)); // ANTENNA1
    stmt->setInt(++pos, cols.antenna2()(i)); // ANTENNA2
    stmt->setInt(++pos, cols.feedId()(i)); // FEED_ID
    stmt->setInt(++pos, cols.spectralWindowId()(i)); // SPECTRAL_WINDOW_ID
    stmt->setDouble(++pos, cols.time()(i)); // TIME
    stmt->setDouble(++pos, cols.interval()(i)); // INTERVAL
    stmt->setDouble(++pos, cols.offset()(i)); // OFFSET
    assert(pos == stmt->getParameterCount());
    int result = stmt->executeUpdate();
  }  
}

void saveState(Connection *con, MeasurementSet &ms) {
  enter();
  MSState &tab = ms.state();
  ROMSStateColumns const cols(tab);

  char const *dbcols[] = {
    "STATE_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "SIG", // INTEGER NOT NULL,
    "REF", // INTEGER NOT NULL,
    "CAL", // REAL NOT NULL,
    "LOAD", // REAL NOT NULL,
    "SUB_SCAN", // INTEGER NOT NULL,
    "OBS_MODE", // TEXT NOT NULL,
    "FLAG_ROW", // INTEGER NOT NULL
    NULL
  };

  enum {
    STATE_ID = 1,
    SIG,
    REF,
    CAL,
    LOAD,
    SUB_SCAN,
    OBS_MODE,
    FLAG_ROW
  };
  string sql = "insert into STATE (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    stmt->setInt(STATE_ID, i);
    stmt->setInt(SIG, cols.sig()(i));
    stmt->setInt(REF, cols.ref()(i));
    stmt->setDouble(CAL, cols.cal()(i));
    stmt->setDouble(LOAD, cols.load()(i));
    stmt->setInt(SUB_SCAN, cols.subScan()(i));
    stmt->setTransientString(OBS_MODE, cols.obsMode()(i).c_str());
    stmt->setInt(FLAG_ROW, cols.flagRow()(i));
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
}

void saveField(Connection *con, MeasurementSet &ms) {
  enter();
  MSField &tab = ms.field();
  ROMSFieldColumns const cols(tab);

  // FIELD table
  {
  char const *dbcols[] = {
    "FIELD_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "NAME", // TEXT NOT NULL,
    "CODE", // TEXT NOT NULL,
    "TIME", // REAL NOT NULL,
    "NUM_POLY", // INTEGER NOT NULL,
    "SOURCE_ID", // INTEGER NOT NULL,
    "EPHEMERIS_ID", // INTEGER DEFAULT NULL,
    "FLAG_ROW", // INTEGER NOT NULL
    NULL
  };

  string sql = "insert into FIELD (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    int pos = 0;
    stmt->setInt(++pos, i); // FIELD_ID
    stmt->setTransientString(++pos, cols.name()(i).c_str()); // NAME
    stmt->setTransientString(++pos, cols.code()(i).c_str()); // CODE
    stmt->setDouble(++pos, cols.time()(i)); // TIME
    stmt->setInt(++pos, cols.numPoly()(i)); // NUM_POLY
    stmt->setInt(++pos, cols.sourceId()(i)); // SOURCE_ID
    if (cols.ephemerisId().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.ephemerisId()(i)); // EPHEMERIS_ID
    }
    stmt->setInt(++pos, cols.flagRow()(i)); // FLAG_ROW
    assert(pos == stmt->getParameterCount());
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
  }

  // FIELD_POLY table
  {
  char const *dbcols[] = {
    "FIELD_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "IDX", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "DELAY_DIRX", // REAL NOT NULL,
    "DELAY_DIRY", // REAL NOT NULL,
    "PHASE_DIRX", // REAL NOT NULL,
    "PHASE_DIRY", // REAL NOT NULL,
    "REFERENCE_DIRX", // REAL NOT NULL,
    "REFERENCE_DIRY", // REAL NOT NULL,
    NULL
  };
  
  string sql = "insert into FIELD_POLY (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  uInt id = 0;
  for (uInt i = 0; i < nrow; i++) {
    uInt npoly = cols.numPoly()(i)+1; // NUM_POLY+1
    Matrix<Double> delayDir = cols.delayDir()(i);
    Matrix<Double> phaseDir = cols.phaseDir()(i);
    Matrix<Double> referenceDir = cols.referenceDir()(i);
    for (uInt j = 0 ; j < npoly ; j++ ) {
      int pos = 0;
      stmt->setInt(++pos, i); // FIELD_ID
      stmt->setInt(++pos, id++); // IDX
      stmt->setDouble(++pos, delayDir(0,j)); // DELAY_DIRX
      stmt->setDouble(++pos, delayDir(1,j)); // DELAY_DIRY
      stmt->setDouble(++pos, phaseDir(0,j)); // PHASE_DIRX
      stmt->setDouble(++pos, phaseDir(1,j)); // PHASE_DIRY
      stmt->setDouble(++pos, referenceDir(0,j)); // REFERENCE_DIRX
      stmt->setDouble(++pos, referenceDir(1,j)); // REFERENCE_DIRY
      assert(pos == stmt->getParameterCount());
      int result = stmt->executeUpdate();
    }
  }
  }
}

void saveFeed(Connection *con, MeasurementSet &ms) {
  enter();
  MSFeed &tab = ms.feed();
  ROMSFeedColumns const cols(tab);

  // FEED table
  {
  char const *dbcols[] = {
    "FEED_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "ANTENNA_ID", // INTEGER NOT NULL,
    "SPECTRAL_WINDOW_ID", // INTEGER NOT NULL,
    "TIME", // REAL NOT NULL,
    "INTERVAL", // REAL NOT NULL,
    "NUM_RECEPTORS", // INTEGER NOT NULL,
    "BEAM_ID", // INTEGER NOT NULL,
    "FOCUS_LENGTH", // REAL DEFAULT NULL,
    "PHASED_FEED_ID", // INTEGER DEFAULT NULL,
    "POSITIONX", // REAL NOT NULL,
    "POSITIONY", // REAL NOT NULL,
    "POSITIONZ", // REAL NOT NULL,
    NULL
  };

  string sql = "insert into FEED (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    int pos = 0;
    stmt->setInt(++pos, i); // FEED_ID
    stmt->setInt(++pos, cols.antennaId()(i)); // ANTENNA_ID
    stmt->setInt(++pos, cols.spectralWindowId()(i)); // SPECTRAL_WINDOW_ID
    stmt->setDouble(++pos, cols.time()(i)); // TIME
    stmt->setDouble(++pos, cols.interval()(i)); // INTERVAL
    stmt->setInt(++pos, cols.numReceptors()(i)); // NUM_RECEPTORS
    stmt->setInt(++pos, cols.beamId()(i)); // BEAM_ID
    if (cols.focusLength().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setDouble(++pos, cols.focusLength()(i)); // FOCUS_LENGTH
    }
    if (cols.phasedFeedId().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.phasedFeedId()(i)); // PHASED_FEED_ID
    }
    {
      Array< Double > const &v = cols.position()(i); // POSITION
      const IPosition &shape = v.shape();
      Double const *data = v.data();
      for (size_t i = 0; i < 3; i++) {
	//"POSITIONX", "POSITIONY", "POSITIONZ"
	stmt->setDouble(++pos, data[i]);
      }
    }
    assert(pos == stmt->getParameterCount());
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
  }

  // FEED_RECEPTOR table
  {
  char const *dbcols[] = {
    "FEED_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "IDX", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "BEAM_OFFSETX", // REAL NOT NULL,
    "BEAM_OFFSETY", // REAL NOT NULL,
    "POLARIZATION_TYPE", // TEXT NOT NULL,
    "POL_RESPONSE", // BLOB NOT NULL,
    "RECEPTOR_ANGLE", // REAL NOT NULL,
    NULL
  };

  string sql = "insert into FEED_RECEPTOR (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  uInt id = 0 ;
  for (uInt i = 0; i < nrow; i++) {
    Int nrec = cols.numReceptors()(i);
    Matrix<Double> beamOffset = cols.beamOffset()(i);
    Vector<String> polType = cols.polarizationType()(i);
    Matrix<Complex> polResponse = cols.polResponse()(i);
    Vector<Double> recAngle = cols.receptorAngle()(i);
    for (Int j = 0 ; j < nrec ; j++ ) {
      int pos = 0;
      stmt->setInt(++pos, i); // FEED_ID
      stmt->setInt(++pos, id++); // IDX
      stmt->setDouble(++pos, beamOffset(0,j)); // BEAMOFFSETX
      stmt->setDouble(++pos, beamOffset(1,j)); // BEAMOFFSETY
      stmt->setTransientString(++pos, polType[j].c_str()); // POLARIZATION_TYPE
      bindArrayAsBlob<Complex>(stmt.get(), ++pos, polResponse.row(j)); // POL_RESPONSE
      stmt->setDouble(++pos, recAngle(j)); // RECEPTOR_ANGLE
      assert(pos == stmt->getParameterCount());
      int result = stmt->executeUpdate();
    }
  }
  }
}

void saveAntenna(Connection *con, MeasurementSet &ms) {
  enter();
  MSAntenna &tab = ms.antenna();
  ROMSAntennaColumns const cols(tab);

  char const *dbcols[] = {
    "ANTENNA_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "NAME", // TEXT NOT NULL,
    "STATION", // TEXT NOT NULL,
    "ANTENNA_TYPE", // TEXT NOT NULL,
    "MOUNT", // TEXT NOT NULL,
    "POSITIONX", // REAL NOT NULL,
    "POSITIONY", // REAL NOT NULL,
    "POSITIONZ", // REAL NOT NULL,
    "OFFSETX", // REAL NOT NULL,
    "OFFSETY", // REAL NOT NULL,
    "OFFSETZ", // REAL NOT NULL,
    "DISH_DIAMETER", // REAL NOT NULL,
    "ORBIT_ID", // INTEGER DEFAULT NULL,
    "MEAN_ORBIT", // BLOB DEFAULT NULL,
    "PHASED_ARRAY_ID", // INTEGER DEFAULT NULL,
    "FLAG_ROW", // INTEGER NOT NULL
    NULL
  };

  string sql = "insert into ANTENNA (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    int pos = 0;
    stmt->setInt(++pos, i); // ANTENNA_ID
    stmt->setTransientString(++pos, cols.name()(i).c_str()); // NAME
    stmt->setTransientString(++pos, cols.station()(i).c_str()); // STATION
    stmt->setTransientString(++pos, cols.type()(i).c_str()); // ANTENNA_TYPE
    stmt->setTransientString(++pos, cols.mount()(i).c_str()); // MOUNT
    {
      Array< Double > const &v = cols.position()(i); // POSITION
      const IPosition &shape = v.shape();
      Double const *data = v.data();
      for (size_t i = 0; i < 3; i++) {
	//"POSITIONX", "POSITIONY", "POSITIONZ"
	stmt->setDouble(++pos, data[i]);
      }
    }
    {
      Array< Double > const &v = cols.position()(i); // OFFSET
      const IPosition &shape = v.shape();
      Double const *data = v.data();
      for (size_t i = 0; i < 3; i++) {
	// "OFFSETX", "OFFSETY", "OFFSETZ"
	stmt->setDouble(++pos, data[i]);
      }
    }
    stmt->setDouble(++pos, cols.dishDiameter()(i)); // DISH_DIAMETER
    if (cols.orbitId().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.orbitId()(i)); // ORBIT_ID
    }
    if (cols.meanOrbit().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Double>(stmt.get(), ++pos, cols.meanOrbit()(i));
    }
    if (cols.phasedArrayId().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.phasedArrayId()(i));
    }
    stmt->setInt(++pos, cols.flagRow()(i)); // FLAG_ROW
    assert(pos == stmt->getParameterCount());
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
}

void saveMain(Connection *con, MeasurementSet &ms) {
  enter();
  ROMSColumns const cols(ms);

  char const *dbcols[] = {
    "MAIN_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "TIME", // REAL NOT NULL,
    "TIME_EXTRA_PREC", // REAL DEFAULT NULL,
    "ANTENNA1", // INTEGER NOT NULL,
    "ANTENNA2", // INTEGER NOT NULL,
    "ANTENNA3", // INTEGER DEFAULT NULL,
    "FEED1", // INTEGER NOT NULL,
    "FEED2", // INTEGER NOT NULL,
    "FEED3", // INTEGER DEFAULT NULL,
    "DATA_DESC_ID", // INTEGER NOT NULL,
    "PROCESSOR_ID", // INTEGER NOT NULL,
    "PHASE_ID", // INTEGER DEFAULT NULL,
    "FIELD_ID", // INTEGER NOT NULL,
    "INTERVAL", // REAL NOT NULL,
    "EXPOSURE", // REAL NOT NULL,
    "TIME_CENTROID", // REAL NOT NULL,
    "PULSAR_BIN", // INTEGER DEFAULT NULL,
    "PULSAR_GATE_ID", // INTEGER DEFAULT NULL,
    "SCAN_NUMBER", // INTEGER NOT NULL,
    "ARRAY_ID", // INTEGER NOT NULL,
    "OBSERVATION_ID", // INTEGER NOT NULL,
    "STATE_ID", // INTEGER NOT NULL,
    "BASELINE_REF", // INTEGER DEFAULT NULL,
    "U", // REAL NOT NULL,
    "V", // REAL NOT NULL,
    "W", // REAL NOT NULL,
    "U2", // REAL DEFAULT NULL,
    "V2", // REAL DEFAULT NULL,
    "W2", // REAL DEFAULT NULL,
    "DATA", // BLOB DEFAULT NULL,
    "FLOAT_DATA", // BLOB DEFAULT NULL,
    "VIDEO_POINT", // BLOB DEFAULT NULL,
    "LAG_DATA", // BLOB DEFAULT NULL,
    "SIGMA", // BLOB NOT NULL,
    "SIGMA_SPECTRUM", // BLOB DEFAULT NULL,
    "WEIGHT", // BLOB NOT NULL,
    "WEIGHT_SPECTRUM", // BLOB DEFAULT NULL,
    "FLAG", // BLOB NOT NULL,
    "FLAG_CATEGORY", // BLOB NOT NULL,
    "FLAG_ROW", // INTEGER NOT NULL,
    NULL
  };

  string sql = "insert into MAIN (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = ms.nrow();
  for (uInt i = 0; i < nrow; i++) {
    int pos = 0;
    stmt->setInt(++pos, i); // MAIN_ID
    stmt->setDouble(++pos, cols.time()(i)); // TIME
    if (cols.timeExtraPrec().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setDouble(++pos, cols.timeExtraPrec()(i)); // TIME_EXTRA_PREC
    }
    stmt->setInt(++pos, cols.antenna1()(i)); // ANTENNA1
    stmt->setInt(++pos, cols.antenna2()(i)); // ANTENNA2
    if (cols.antenna3().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.antenna3()(i)); // ANTENNA3
    }
    stmt->setInt(++pos, cols.feed1()(i)); // FEED1
    stmt->setInt(++pos, cols.feed2()(i)); // FEED2
    if (cols.feed3().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.feed3()(i)); // FEED3
    }
    stmt->setInt(++pos, cols.dataDescId()(i)); // DATA_DESC_ID
    stmt->setInt(++pos, cols.processorId()(i)); // PROCESSOR_ID
    if (cols.phaseId().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.phaseId()(i)); // PHASE_ID
    }
    stmt->setInt(++pos, cols.fieldId()(i)); // FIELD_ID
    stmt->setDouble(++pos, cols.interval()(i)); // INTERVAL
    stmt->setDouble(++pos, cols.exposure()(i)); // EXPOSURE
    stmt->setDouble(++pos, cols.timeCentroid()(i)); // TIME_CENTROID
    if (cols.pulsarBin().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.pulsarBin()(i)); // PULSAR_BIN
    }
    if (cols.pulsarGateId().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.pulsarGateId()(i)); // PULSAR_GATE_ID
    }
    stmt->setInt(++pos, cols.scanNumber()(i)); // SCAN_NUMBER
    stmt->setInt(++pos, cols.arrayId()(i)); // ARRAY_ID
    stmt->setInt(++pos, cols.observationId()(i)); // OBSERVATION_ID
    stmt->setInt(++pos, cols.stateId()(i)); // STATE_ID
    if (cols.baselineRef().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.baselineRef()(i)); // BASELINE_REF
    }
    {
      Array< Double > const &v = cols.uvw()(i);
      const IPosition &shape = v.shape();
      Double const *data = v.data();
      for (size_t i = 0; i < 3; i++) {
	//"U", "V", "W"
	stmt->setDouble(++pos, data[i]);
      }
    }
    if (cols.uvw2().isNull()) {
      for (size_t i = 0; i < 3; i++) {
	stmt->setNull(++pos);
      }
    } else {
      Array< Double > const &v = cols.uvw2()(i);
      const IPosition &shape = v.shape();
      Double const *data = v.data();
      for (size_t i = 0; i < 3; i++) {
	//"U2", "V2", "W2"
	stmt->setDouble(++pos, data[i]);
      }
    }
    if (cols.data().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Complex>(stmt.get(), ++pos, cols.data()(i));
    }
    if (cols.floatData().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.floatData()(i));
    }
    if (cols.videoPoint().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Complex>(stmt.get(), ++pos, cols.videoPoint()(i));
    }
    if (cols.lagData().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Complex>(stmt.get(), ++pos, cols.lagData()(i));
    }
    bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.sigma()(i));
    if (cols.sigmaSpectrum().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.sigmaSpectrum()(i));
    }
    bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.weight()(i));
    if (cols.weightSpectrum().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.weightSpectrum()(i));
    }
    bindArrayAsBlob<Bool>(stmt.get(), ++pos, cols.flag()(i));
    bindArrayAsBlob<Bool>(stmt.get(), ++pos, cols.flagCategory()(i));
    stmt->setInt(++pos, cols.flagRow()(i)); // FLAG_ROW
    assert(pos == stmt->getParameterCount());
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
}

void saveFlagCmd(Connection *con, MeasurementSet &ms) {
  enter();
  MSFlagCmd &tab = ms.flagCmd();
  ROMSFlagCmdColumns const cols(tab);

  char const *dbcols[] = {
    "FLAG_CMD_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "TIME", // NOT NULL,
    "INTERVAL", // REAL NOT NULL,
    "FLAG_TYPE", // TEXT NOT NULL,
    "REASON", // TEXT NOT NULL,
    "LEVEL", // INTEGER NOT NULL,
    "SEVERITY", // INTEGER NOT NULL,
    "APPLIED", // INTEGER NOT NULL, bool
    "COMMAND", // TEXT NOT NULL
    NULL
  };

  enum {
    FLAG_CMD_ID = 1,
    TIME,
    INTERVAL,
    FLAG_TYPE,
    REASON,
    LEVEL,
    SEVERITY,
    APPLIED,
    COMMAND
  };

  string sql = "insert into FLAG_CMD (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));

  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    stmt->setInt(FLAG_CMD_ID, i);
    stmt->setDouble(TIME, cols.time()(i));
    stmt->setDouble(INTERVAL, cols.interval()(i));
    stmt->setTransientString(FLAG_TYPE, cols.type()(i).c_str());
    stmt->setTransientString(REASON, cols.reason()(i).c_str());
    stmt->setInt(LEVEL, cols.level()(i));
    stmt->setInt(SEVERITY, cols.severity()(i));
    stmt->setInt(APPLIED, cols.applied()(i) ? 1 : 0);
    stmt->setTransientString(COMMAND, cols.command()(i).c_str());
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
}

void mssave(Connection *con, char const*filename) {
  enter();
  static void (*funcs[])(Connection *, MeasurementSet &) = {
    saveDataDesc, saveField, saveFeed, saveFreqOffset, 
    saveState, saveFlagCmd,
    NULL
  };
  MeasurementSet ms(filename);
  for (size_t i = 0; funcs[i] != NULL; i++) {
    con->execute("BEGIN");
    funcs[i](con, ms);
    con->execute("COMMIT");
  }
}

double currenttime() {
  struct timeval tv;
  int result = gettimeofday(&tv, NULL);
  return tv.tv_sec + ((double)tv.tv_usec) / 1000000.;
}

char *readFileContent(char const *filename) {
  char *buf = NULL;
  ifstream fs(filename);
  try {
    if (! fs) {
      throw "Could not open file.";
    }
    streampos begPos = fs.tellg();
    fs.seekg(0, fstream::end);
    streampos eofPos = fs.tellg();
    size_t size = eofPos - begPos;

    fs.clear();
    fs.seekg(0, fstream::beg);
 
    buf = new char[size + 1];
    memset(buf, 0, size + 1);
    if (! fs.read(buf, size)) {
      throw "failed to read.";
    }
    if (size != static_cast<size_t>(fs.gcount())) {
      throw "failed to read.";
    }
  } catch (...) {
    delete[] buf;
    fs.close();
    throw;
  }
  fs.close();
  return buf;
}

void conv(char const *prefix, char const *msfile, char const *basename) {
  string msm = basename;
  msm += "m.db";
  cout << msm << endl;
  string mst = basename;
  mst += "t.db";

  char const sql[] = "/sql/";
  // create empty master db from msm
  {
    char const ddl_file[] = "MSM.ddl";
    size_t ddl_path_size = strlen(prefix) + strlen(sql) + strlen(ddl_file) + 1;
    char ddl_path[ddl_path_size];
    snprintf(ddl_path, ddl_path_size, "%s%s%s", prefix, sql, ddl_file);
    cerr << ddl_path << endl;
    char *msm_ddl = readFileContent(ddl_path);
    try {
      string dbfile = "file:";
      dbfile += msm;
      unique_ptr<Connection> con(Connection::open(dbfile.c_str()));
      con->execute(msm_ddl);
    } catch (...) {
      delete[] msm_ddl;
      throw;
    }
    delete[] msm_ddl;
  }
    
  // create empty transaction db from mst
  {
    char const ddl_file[] = "MST.ddl";
    size_t ddl_path_size = strlen(prefix) + strlen(sql) + strlen(ddl_file) + 1;
    char ddl_path[ddl_path_size];
    snprintf(ddl_path, ddl_path_size, "%s%s%s", prefix, sql, ddl_file);
    cerr << ddl_path << endl;
    char *mst_ddl = readFileContent(ddl_path);
    try {
      string dbfile = "file:";
      dbfile += mst;
      unique_ptr<Connection> con(Connection::open(dbfile.c_str()));
      con->execute(mst_ddl);
    } catch (...) {
      delete[] mst_ddl;
      throw;
    }
    delete[] mst_ddl;
  }

  // open empty transaction db attached with master db.
  {
    string dbfile = "file:";
    dbfile += mst;
    unique_ptr<Connection> con(Connection::open(dbfile.c_str()));
    // attach master db
    unique_ptr<PreparedStatement> stmt(con->prepare("ATTACH DATABASE ? AS msm;"));
    stmt->setTransientString(1, msm.c_str());
    stmt->executeUpdate();
    mssave(con.get(), msfile);
  }
}

char const *progName = "";

void usage() {
  cerr << "Usage: " << progName << " [options] MSFile basename\n";
  cerr << "options:: \n";
  cerr << "\t--prefix path\tA path where sakura is installed.\n";
  cerr << "\t-p path\n";
}

}

int main(int argc, char const * const argv[]) { 
  progName = argv[0];
  char const *prefix = "/opt/sakura";

  static struct option const long_options[] = {
    {"prefix", 1, NULL, 'p'},
    {0, 0, NULL, 0}
  };

  for (;;) {
    int option_index = 0;
    int optCh = getopt_long (argc, const_cast<char *const *>(argv), "p:",
			     long_options, &option_index);
    if (optCh == -1) {
      break;
    }
    switch (optCh) {
    case 'p':
      prefix = optarg;
      break;
    case '?':
      usage();
      return 1;
    default:
      assert(false);
      return 1;
    }
  }

  int argStart = optind;
  if (argc - argStart != 2) {
    usage();
    return 1;
  }
  cerr << prefix << endl;
  double start = currenttime();
  conv(prefix, argv[argStart], argv[argStart + 1]);
  double end = currenttime();
  cout << end - start << "sec\n";
  return 0;
}
