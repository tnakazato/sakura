#include <sys/time.h>
#include <unistd.h>
#include <getopt.h>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <cstdlib>
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
char const ENV_VAR_SAKURA_ROOT[] = "SAKURA_ROOT";

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
    T const *blobData = data;
    //cout << "data: " << data <<", " << elements << endl;;
    if (data == NULL) {
      assert(elements == 0);
      blobData = (T const *)DUMMY_AREA;
    }
    stmt->setTransientBlob(pos, blobData, sizeof(T) * elements);
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

  // SPECTRAL_WINDOW table
  {
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

  enum {
    SPECTRAL_WINDOW_ID = 1,
    NUM_CHAN,
    NAME,
    REF_FREQUENCY,
    CHAN_FREQ,
    CHAN_WIDTH,
    MEAS_FREQ_REF,
    EFFECTIVE_BW,
    RESOLUTION,
    TOTAL_BANDWIDTH,
    NET_SIDEBAND,
    BBC_NO,
    BBC_SIDEBAND,
    IF_CONV_CHAIN,
    RECEIVER_ID,
    FREQ_GROUP,
    FREQ_GROUP_NAME,
    DOPPLER_ID,
    FLAG_ROW
  };

  string sql = "insert into SPECTRAL_WINDOW (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    stmt->setInt(SPECTRAL_WINDOW_ID, i);
    stmt->setInt(NUM_CHAN, cols.numChan()(i));
    stmt->setTransientString(NAME, cols.name()(i).c_str());
    stmt->setDouble(REF_FREQUENCY, cols.refFrequency()(i));
    bindArrayAsBlob<Double>(stmt.get(), CHAN_FREQ, cols.chanFreq()(i));
    bindArrayAsBlob<Double>(stmt.get(), CHAN_WIDTH, cols.chanWidth()(i));
    stmt->setInt(MEAS_FREQ_REF, cols.measFreqRef()(i));
    bindArrayAsBlob<Double>(stmt.get(), EFFECTIVE_BW, cols.effectiveBW()(i));
    bindArrayAsBlob<Double>(stmt.get(), RESOLUTION, cols.resolution()(i));
    stmt->setDouble(TOTAL_BANDWIDTH, cols.totalBandwidth()(i));
    stmt->setInt(NET_SIDEBAND, cols.netSideband()(i));

    if (cols.bbcNo().isNull()) {
      stmt->setNull(BBC_NO);
    } else {
      stmt->setInt(BBC_NO, cols.bbcNo()(i));
    }
    if (cols.bbcSideband().isNull()) {
      stmt->setNull(BBC_SIDEBAND);
    } else {
      stmt->setInt(BBC_SIDEBAND, cols.bbcSideband()(i));
    }
    stmt->setInt(IF_CONV_CHAIN, cols.ifConvChain()(i));
    if (cols.receiverId().isNull()) {
      stmt->setNull(RECEIVER_ID);
    } else {
      stmt->setInt(RECEIVER_ID, cols.receiverId()(i));
    }
    stmt->setInt(FREQ_GROUP, cols.freqGroup()(i));
    stmt->setTransientString(FREQ_GROUP_NAME, cols.freqGroupName()(i).c_str());
    if (cols.dopplerId().isNull()) {
      stmt->setNull(DOPPLER_ID);
    } else {
      stmt->setInt(DOPPLER_ID, cols.dopplerId()(i));
    }
    stmt->setInt(FLAG_ROW, cols.flagRow()(i));
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
  }

  //SPECTRAL_WINDOW_ASSOC table
  {
  // return if optional ASSOC_SPW_ID and ASSOC_NATURE columns don't exist
  if (cols.assocSpwId().isNull() && cols.assocNature().isNull()) return;

  char const *dbcols[] = {
    "SPECTRAL_WINDOW_ID", // INTEGER NOT NULL,
    "ASSOC_SPW_ID", // INTEGER NOT NULL,
    "ASSOC_NATURE", // TEXT NOT NULL,
    NULL
  };

  enum {
    SPECTRAL_WINDOW_ID = 1,
    ASSOC_SPW_ID,
    ASSOC_NATURE
  };

  string sql = "insert into SPECTRAL_WINDOW_ASSOC (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  Bool const isId = (!cols.assocSpwId().isNull());
  Bool const isNat = (!cols.assocNature().isNull());
  for (uInt i = 0; i < nrow; i++) {
    Vector< Int > assocSpwId;
    Vector< String > assocNature;
    size_t nelem = 0;
    Bool const isIdRow = isId && cols.assocSpwId().isDefined(i);
    Bool const isNatRow = isNat && cols.assocNature().isDefined(i);
    if (isIdRow) {
      assocSpwId = cols.assocSpwId()(i); // ASSOC_SPW_ID
      nelem = assocSpwId.nelements();
    }
    if (isNatRow) {
      assocNature = cols.assocNature()(i); // ASSOC_NATURE
      if (isIdRow) assert(assocNature.nelements() == nelem);
      else nelem = assocNature.nelements();
    }
    for (size_t j = 0; j < nelem; j++) {
      stmt->setInt(SPECTRAL_WINDOW_ID, i);
      if (isIdRow) {
	stmt->setInt(ASSOC_SPW_ID, assocSpwId[j]);
      } else {
	stmt->setNull(ASSOC_SPW_ID);
      }
      if (isNatRow) {
	stmt->setTransientString(ASSOC_NATURE, assocNature[j].c_str());
      } else {
	stmt->setNull(ASSOC_NATURE);
      }
      int result = stmt->executeUpdate();
      assert(result == 1);
    }
  }
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

void saveSource(Connection *con, MeasurementSet &ms) {
  enter();
  MSSource &tab = ms.source();
  ROMSSourceColumns const cols(tab);

  // return if optional SOURCE table doesn't exist
  if (tab.isNull()) return;

  // SOURCE table
  {
  char const *dbcols[] = {
    "SOURCE_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "TIME", // REAL NOT NULL,
    "INTERVAL", // REAL NOT NULL,
    "SPECTRAL_WINDOW_ID", // INTEGER NOT NULL,
    "NUM_LINES", // INTEGER NOT NULL,
    "NAME", // TEXT NOT NULL,
    "CALIBRATION_GROUP", // INTEGER NOT NULL,
    "CODE", // TEXT NOT NULL,
    "DIRECTIONX", // REAL NOT NULL,
    "DIRECTIONY", // REAL NOT NULL,
    "POSITIONX", // REAL,
    "POSITIONY", // REAL,
    "POSITIONZ", // REAL,
    "PROPER_MOTIONX", // REAL NOT NULL,
    "PROPER_MOTIONY", // REAL NOT NULL,
    "PULSAR_ID", // INTEGER,
    NULL
  };

  enum {
    SOURCE_ID = 1,
    TIME,
    INTERVAL,
    SPECTRAL_WINDOW_ID,
    NUM_LINES,
    NAME,
    CALIBRATION_GROUP,
    CODE,
    DIRECTIONX,
    DIRECTIONY,
    POSITIONX,
    POSITIONY,
    POSITIONZ,
    PROPER_MOTIONX,
    PROPER_MOTIONY,
    PULSAR_ID
  };

  string sql = "insert into SOURCE (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    stmt->setInt(SOURCE_ID, cols.sourceId()(i));
    stmt->setDouble(TIME, cols.time()(i));
    stmt->setDouble(INTERVAL, cols.interval()(i));
    stmt->setInt(SPECTRAL_WINDOW_ID, cols.spectralWindowId()(i));
    stmt->setInt(NUM_LINES, cols.numLines()(i));
    stmt->setTransientString(NAME, cols.name()(i).c_str());
    stmt->setInt(CALIBRATION_GROUP, cols.calibrationGroup()(i));
    stmt->setTransientString(CODE, cols.code()(i).c_str());
    {
      Array< Double > const &v = cols.direction()(i); // DIRECTION
      Double const *data = v.data();
      assert(v.nelements() == 2);
      stmt->setDouble(DIRECTIONX, data[0]);
      stmt->setDouble(DIRECTIONY, data[1]);
    }
    {
      if (cols.position().isNull()) {
	stmt->setNull(POSITIONX);
	stmt->setNull(POSITIONY);
	stmt->setNull(POSITIONZ);
      } else {
	Array< Double > const &v = cols.position()(i); // POSITION
	assert(v.nelements() == 3);
	Double const *data = v.data();
	stmt->setDouble(POSITIONX, data[0]);
	stmt->setDouble(POSITIONY, data[1]);
	stmt->setDouble(POSITIONZ, data[2]);
      }
    }
    {
      Array< Double > const &v = cols.properMotion()(i); // PROPER_MOTION
      Double const *data = v.data();
      assert(v.nelements() == 2);
      stmt->setDouble(PROPER_MOTIONX, data[0]);
      stmt->setDouble(PROPER_MOTIONY, data[1]);
    }
    if (cols.pulsarId().isNull()) {
      stmt->setNull(PULSAR_ID);
    } else {
      stmt->setInt(PULSAR_ID, cols.pulsarId()(i));
    }
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
  }

  // SOURCE_LINE_ATTR table
  {
  // skip if none of optional TRANSITION, REST_FREQUENCY,
  // and SYSVEL columns exists
  if ( !cols.transition().isNull() || !cols.restFrequency().isNull() ||
       !cols.sysvel().isNull() ) {

    char const *dbcols[] = {
      "SOURCE_ID", // INTEGER NOT NULL,
      "IDX",  // INTEGER NOT NULL (Starts with 0),
      "TRANSITION", // TEXT,
      "REST_FREQUENCY", // REAL,
      "SYSVEL", // REAL,
    NULL
    };

    enum {
      SOURCE_ID = 1,
      IDX,
      TRANSITION,
      REST_FREQUENCY,
      SYSVEL
    };

    string sql = "insert into SOURCE_LINE_ATTR (";
    sql += SQL::join(dbcols) + ") values (";
    sql += SQL::bindChars(dbcols) + ")";
    unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
    uInt nrow = tab.nrow();
    Bool const isTr = (!cols.transition().isNull());
    Bool const isRf = (!cols.restFrequency().isNull());
    Bool const isSys = (!cols.sysvel().isNull());
    Vector< String > transition;
    Vector< Double > restFrequency, sysvel;
    for (uInt i = 0; i < nrow; i++) {
      uInt const nLine = cols.numLines()(i);
      if (nLine > 0) {
	if ( isTr ) {
	  transition = cols.transition()(i); // TRANSITION
	  assert(transition.nelements() == nLine);
	}
	if ( isRf ) {
	  restFrequency = cols.restFrequency()(i); // REST_FREQUENCY
	  assert(restFrequency.nelements() == nLine);
	}
	if ( isSys ) {
	  sysvel = cols.sysvel()(i); // SYSVEL
	  assert(sysvel.nelements() == nLine);
	}
	for (uInt j = 0; j < nLine; j++) {
	  stmt->setInt(SOURCE_ID, cols.sourceId()(i));
	  stmt->setInt(IDX, j);
	  if (isTr) {
	    stmt->setTransientString(TRANSITION, transition[j].c_str());
	  } else {
	    stmt->setNull(TRANSITION);
	  }
	  if (isRf) {
	    stmt->setDouble(REST_FREQUENCY, restFrequency[j]);
	  } else {
	    stmt->setNull(REST_FREQUENCY);
	  }
	  if (isSys) {
	    stmt->setDouble(SYSVEL, sysvel[j]);
	  } else {
	    stmt->setNull(SYSVEL);
	  }
	  int result = stmt->executeUpdate();
	  assert(result == 1);
	}
      }
    }
  }
  }

//   // SOURCE_MODEL table
//   {
//   // return if optional SOURCE_MODEL column doesn't exist
//   if (cols.sourceModel().isNull()) return;
//   char const *dbcols[] = {
//     "SOURCE_ID", // INTEGER NOT NULL,
//     "MODEL_KEY", // TEXT NOT NULL,
//     "MODEL", // NONE (Table Record),
//     NULL
//   };

//   enum {
//     SOURCE_ID = 1,
//     MODEL_KEY,
//     MODEL
//   };

//   string sql = "insert into SOURCE_MODEL (";
//   sql += SQL::join(dbcols) + ") values (";
//   sql += SQL::bindChars(dbcols) + ")";
//   unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
//   uInt nrow = tab.nrow();
//   for (uInt i = 0; i < nrow; i++) {
//     TableRecord const &v = cols.sourceModel()(i);
    
//     if (v.nfields() > 0) {
//       for (uInt j = 0; j < v.nfields(); j++) {
// 	stmt->setInt(SOURCE_ID, cols.sourceId()(i));
// 	cout << "FIELD = " << j << ":" << v.name(j).c_str() << endl;
// 	//stmt->setTransientString(MODEL_KEY, v.name(j).c_str());
// 	cout << "data type = " << v.type(j) << endl;

// 	int result = stmt->executeUpdate();
// 	assert(result == 1);
//       }
//     }
//   }
//   }

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
    assert(result == 1);
  }  
}

void saveColumnKeyword(Connection *con, MeasurementSet &ms) {
  //enter();
  //MS &tab = ms.();
  //ROMS const cols(tab);
}

void saveTableKeyword(Connection *con, MeasurementSet &ms) {
  enter();
  //cout << "[ " << ms.keywordSet() << " ]" << flush;
}

void saveObs(Connection *con, MeasurementSet &ms) {
  enter();
  MSObservation &tab = ms.observation();
  ROMSObservationColumns const cols(tab);

  // OBSERVATION table
  {
  char const *dbcols[] = {
    "OBSERVATION_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "TELESCOPE_NAME", // TEXT NOT NULL,
    "TIME_RANGE_START", // REAL NOT NULL,
    "TIME_RANGE_END", // REAL NOT NULL,
    "OBSERVER", // TEXT NOT NULL,
    "SCHEDULE_TYPE", // TEXT NOT NULL,
    "PROJECT", // TEXT NOT NULL,
    "RELEASE_DATE", // REAL NOT NULL,
    "FLAG_ROW", // INTEGER NOT NULL
    NULL
  };
  enum {
    OBSERVATION_ID = 1,
    TELESCOPE_NAME,
    TIME_RANGE_START,
    TIME_RANGE_END,
    OBSERVER,
    SCHEDULE_TYPE,
    PROJECT,
    RELEASE_DATE,
    FLAG_ROW
  };
  string sql = "insert into OBSERVATION (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    stmt->setInt(OBSERVATION_ID, i);
    stmt->setTransientString(TELESCOPE_NAME, cols.telescopeName()(i).c_str());

    Array<Double> const &v = cols.timeRange()(i);
    assert(v.nelements() == 2);
    Double const *data = v.data();
    stmt->setDouble(TIME_RANGE_START, data[0]);
    stmt->setDouble(TIME_RANGE_END, data[1]);
    stmt->setTransientString(OBSERVER, cols.observer()(i).c_str());
    stmt->setTransientString(SCHEDULE_TYPE, cols.scheduleType()(i).c_str());
    stmt->setTransientString(PROJECT, cols.project()(i).c_str());
    stmt->setDouble(RELEASE_DATE, cols.releaseDate()(i));
    stmt->setInt(FLAG_ROW, cols.flagRow()(i));
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
  }

  // OBSERVATION_LOG table
  {
  char const *dbcols[] = {
    "OBSERVATION_ID", // INTEGER NOT NULL,
    "IDX", // INTEGER NOT NULL,
    "LOG", // TEXT NOT NULL,
    NULL
  };
  enum {
    OBSERVATION_ID = 1,
    IDX, 
    LOG
  };
  string sql = "insert into OBSERVATION_LOG (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  try {
    for (uInt i = 0; i < nrow; i++) {
      Array<String> const &v = cols.log()(i);
      String const *data = v.data();
      uInt nlog = v.nelements();
      for (uInt j = 0; j < nlog; j++) {
	stmt->setInt(OBSERVATION_ID, i);
	stmt->setInt(IDX, j);
	stmt->setTransientString(LOG, data[j].c_str());
	int result = stmt->executeUpdate();
	assert(result == 1);
      }
    }
  } catch (AipsError x) {
    // do not create OBSERVATION_LOG table if LOG contains no data.
  }



  }

  // OBSERVATION_SCHEDULE table
  {
  char const *dbcols[] = {
    "OBSERVATION_ID", // INTEGER NOT NULL,
    "IDX", // INTEGER NOT NULL,
    "SCHEDULE", // TEXT NOT NULL,
    NULL
  };
  enum {
    OBSERVATION_ID = 1,
    IDX,
    SCHEDULE
  };
  string sql = "insert into OBSERVATION_SCHEDULE (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  try {
    for (uInt i = 0; i < nrow; i++) {
      Array<String> const &schedule = cols.schedule()(i);
      String const *data = schedule.data();
      for (uInt j = 0; j < schedule.nelements(); j++) {
	stmt->setInt(OBSERVATION_ID, i);
	stmt->setInt(IDX, j);
	stmt->setTransientString(SCHEDULE, data[j].c_str());
	int result = stmt->executeUpdate();
	assert(result == 1);
      }
    }
  } catch (AipsError x) {
    // do not create OBSERVATION_SCHEDULE table if SCHEDULE contains no data.
  }

  }

}

void saveProcessor(Connection *con, MeasurementSet &ms) {
  enter();
  MSProcessor &tab = ms.processor();
  ROMSProcessorColumns const cols(tab);

  char const *dbcols[] = {
    "PROCESSOR_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "PROCESSOR_TYPE", // TEXT NOT NULL,
    "SUB_TYPE", //  TEXT NOT NULL,
    "TYPE_ID", // INTEGER NOT NULL,
    "MODE_ID", // INTEGER NOT NULL,
    "PASS_ID", // INTEGER,
    "FLAG_ROW", // INTEGER NOT NULL
    NULL
  };

  enum {
    PROCESSOR_ID = 1,
    PROCESSOR_TYPE,
    SUB_TYPE,
    TYPE_ID,
    MODE_ID,
    PASS_ID,
    FLAG_ROW
  };
  string sql = "insert into PROCESSOR (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    stmt->setInt(PROCESSOR_ID, i);
    stmt->setTransientString(PROCESSOR_TYPE, cols.type()(i).c_str());
    stmt->setTransientString(SUB_TYPE, cols.subType()(i).c_str());
    stmt->setInt(TYPE_ID, cols.typeId()(i));
    stmt->setInt(MODE_ID, cols.modeId()(i));
    if (cols.passId().isNull()) {
      stmt->setNull(PASS_ID);
    } else {
      stmt->setInt(PASS_ID, cols.passId()(i));
    }
    stmt->setInt(FLAG_ROW, cols.flagRow()(i));
    int result = stmt->executeUpdate();
    assert(result == 1);
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
    Int sourceId = cols.sourceId()(i); // SOURCE_ID
    if ( sourceId < 0 )
      stmt->setNull(++pos);
    else
      stmt->setInt(++pos, sourceId);
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
  for (uInt i = 0; i < nrow; i++) {
    uInt npoly = cols.numPoly()(i)+1; // NUM_POLY+1
    Matrix<Double> delayDir = cols.delayDir()(i);
    Matrix<Double> phaseDir = cols.phaseDir()(i);
    Matrix<Double> referenceDir = cols.referenceDir()(i);
    for (uInt j = 0 ; j < npoly ; j++ ) {
      int pos = 0;
      stmt->setInt(++pos, i); // FIELD_ID
      stmt->setInt(++pos, j); // IDX
      stmt->setDouble(++pos, delayDir(0,j)); // DELAY_DIRX
      stmt->setDouble(++pos, delayDir(1,j)); // DELAY_DIRY
      stmt->setDouble(++pos, phaseDir(0,j)); // PHASE_DIRX
      stmt->setDouble(++pos, phaseDir(1,j)); // PHASE_DIRY
      stmt->setDouble(++pos, referenceDir(0,j)); // REFERENCE_DIRX
      stmt->setDouble(++pos, referenceDir(1,j)); // REFERENCE_DIRY
      assert(pos == stmt->getParameterCount());
      int result = stmt->executeUpdate();
      assert(result == 1);
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
      //const IPosition &shape = v.shape();
      Double const *data = v.data();
      for (size_t j = 0; j < 3; j++) {
	//"POSITIONX", "POSITIONY", "POSITIONZ"
	stmt->setDouble(++pos, data[j]);
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
  for (uInt i = 0; i < nrow; i++) {
    Int nrec = cols.numReceptors()(i);
    Matrix<Double> beamOffset = cols.beamOffset()(i);
    Vector<String> polType = cols.polarizationType()(i);
    Matrix<Complex> polResponse = cols.polResponse()(i);
    Vector<Double> recAngle = cols.receptorAngle()(i);
    for (Int j = 0 ; j < nrec ; j++ ) {
      int pos = 0;
      stmt->setInt(++pos, i); // FEED_ID
      stmt->setInt(++pos, j); // IDX
      stmt->setDouble(++pos, beamOffset(0,j)); // BEAMOFFSETX
      stmt->setDouble(++pos, beamOffset(1,j)); // BEAMOFFSETY
      stmt->setTransientString(++pos, polType[j].c_str()); // POLARIZATION_TYPE
      bindArrayAsBlob<Complex>(stmt.get(), ++pos, polResponse.row(j)); // POL_RESPONSE
      stmt->setDouble(++pos, recAngle(j)); // RECEPTOR_ANGLE
      assert(pos == stmt->getParameterCount());
      int result = stmt->executeUpdate();
      assert(result == 1);
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
      //const IPosition &shape = v.shape();
      Double const *data = v.data();
      for (size_t j = 0; j < 3; j++) {
	//"POSITIONX", "POSITIONY", "POSITIONZ"
	stmt->setDouble(++pos, data[j]);
      }
    }
    {
      Array< Double > const &v = cols.offset()(i); // OFFSET
      //const IPosition &shape = v.shape();
      Double const *data = v.data();
      for (size_t j = 0; j < 3; j++) {
	// "OFFSETX", "OFFSETY", "OFFSETZ"
	stmt->setDouble(++pos, data[j]);
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

void saveWeather(Connection *con, MeasurementSet &ms) {
  enter();
  MSWeather &tab = ms.weather();
  ROMSWeatherColumns const cols(tab);

  char const *dbcols[] = {
    "ANTENNA_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "TIME", // REAL NOT NULL,
    "INTERVAL", // REAL NOT NULL,
    "H2O", // REAL DEFAULT NULL,
    "IONOS_ELECTRON", // REAL DEFAULT NULL,
    "PRESSURE", // REAL DEFAULT NULL,
    "REL_HUMIDITY", // REAL DEFAULT NULL,
    "TEMPERATURE", // REAL DEFAULT NULL,
    "DEW_POINT", // REAL DEFAULT NULL,
    "WIND_DIRECTION", // REAL DEFAULT NULL,
    "WIND_SPEED", // REAL DEFAULT NULL,
    "H2O_FLAG", // INTEGER DEFAULT NULL,
    "IONOS_ELECTRON_FLAG", // INTEGER DEFAULT NULL,
    "PRESSURE_FLAG",  // INTEGER DEFAULT NULL,
    "REL_HUMIDITY_FLAG",  // INTEGER DEFAULT NULL,
    "TEMPERATURE_FLAG",  // INTEGER DEFAULT NULL,
    "DEW_POINT_FLAG",  // INTEGER DEFAULT NULL,
    "WIND_DIRECTION_FLAG",  // INTEGER DEFAULT NULL,
    "WIND_SPEED_FLAG",  // INTEGER DEFAULT NULL,
    NULL
  };

  string sql = "insert into WEATHER (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    int pos = 0;
    stmt->setInt(++pos, i); // ANTENNA_ID
    stmt->setDouble(++pos, cols.time()(i)); // TIME
    stmt->setDouble(++pos, cols.interval()(i)); // INTERVAL

     if (cols.H2O().isNull()) {
       stmt->setNull(++pos);
     } else {
       stmt->setDouble(++pos, cols.H2O()(i)); // H2O
     }

    if (cols.ionosElectron().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setDouble(++pos, cols.ionosElectron()(i)); // IONOS_ELECTRON
    }
    if (cols.pressure().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setDouble(++pos, cols.pressure()(i)); // PRESSURE
    }
    if (cols.relHumidity().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setDouble(++pos, cols.relHumidity()(i)); //REL_HUMIDITY 
    }
    if (cols.temperature().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setDouble(++pos, cols.temperature()(i)); // TEMPERATURE
    }
    if (cols.dewPoint().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setDouble(++pos, cols.dewPoint()(i)); // DEW_POINT
    }
    if (cols.windDirection().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setDouble(++pos, cols.windDirection()(i)); // WIND_DIRECTION
    }
    if (cols.windSpeed().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setDouble(++pos, cols.windSpeed()(i)); // WIND_SPEED
    }
     if (cols.H2OFlag().isNull()) {
       stmt->setNull(++pos);
     } else {
       stmt->setInt(++pos, cols.H2OFlag()(i)); // H2O_FLAG
     }
    if (cols.ionosElectronFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.ionosElectronFlag()(i)); // IONOS_ELECTRON_FLAG
    }
    if (cols.pressureFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.pressureFlag()(i)); // PRESSURE_FLAG
    }
    if (cols.relHumidityFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.relHumidityFlag()(i)); // REL_HUMIDITY_FLAG
    }
    if (cols.temperatureFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.temperatureFlag()(i)); // TEMPERATURE_FLAG
    }
    if (cols.dewPointFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.dewPointFlag()(i)); // DEW_POINT_FLAG
    }
    if (cols.windDirectionFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.windDirectionFlag()(i)); // WIND_DIRECTION_FLAG
    }
    if (cols.windSpeedFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.windSpeedFlag()(i)); // WIND_SPEED_FLAG
    }

    assert(pos == stmt->getParameterCount());
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
}

void savePointing(Connection *con, MeasurementSet &ms) {
  enter();
  MSPointing &tab = ms.pointing();
  ROMSPointingColumns const cols(tab);

  char const *dbcols[] = {
    "ANTENNA_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "TIME", // REAL NOT NULL,
    "INTERVAL", // REAL NOT NULL,
    "NAME", // TEXT NOT NULL,
    "NUM_POLY", // INTEGER NOT NULL,
    "TIME_ORIGIN", // REAL NOT NULL,
    "DIRECTION", // BLOB
    "TARGET", // BLOB
    "POINTING_OFFSET", // BLOB 
    "SOURCE_OFFSET", // BLOB 
    "ENCODERX", // REAL 
    "ENCODERY", // REAL 
    "POINTING_MODEL_ID", // INTEGER DEFAULT NULL,
    "TRACKING", // INTEGER NOT NULL,
    "ON_SOURCE", // INTEGER DEFAULT NULL,
    "OVER_THE_TOP", // INTEGER DEFAULT NULL,
    NULL
  };

  string sql = "insert into POINTING (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    int pos = 0;
    stmt->setInt(++pos, i); // ANTENNA_ID
    stmt->setDouble(++pos, cols.time()(i)); // TIME
    stmt->setDouble(++pos, cols.interval()(i)); // INTERVAL
    stmt->setTransientString(++pos, cols.name()(i).c_str()); // NAME
    stmt->setInt(++pos, cols.numPoly()(i)); // NUM_POLY

    stmt->setDouble(++pos, cols.timeOrigin()(i)); // TIME_ORIGIN
    bindArrayAsBlob<Double>(stmt.get(), ++pos, cols.direction()(i));//DIRECTION
    bindArrayAsBlob<Double>(stmt.get(), ++pos, cols.target()(i));//TARGET

    if (cols.pointingOffset().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Double>(stmt.get(), ++pos, cols.pointingOffset()(i));//POINTING_OFFSET
    }
    if (cols.sourceOffset().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Double>(stmt.get(), ++pos, cols.sourceOffset()(i));//SOURCE_OFFSET
    }
    if (cols.encoder().isNull()){
      stmt->setNull(++pos);
      stmt->setNull(++pos);
    }else {
      Array< Double > const &v = cols.encoder()(i); // ENCODERX ENCODERY
      //const IPosition &shape = v.shape();
      Double const *data = v.data();
      for (size_t j = 0; j < 2; j++) {
	// "ENCODERX", "ENCODERY"
	stmt->setDouble(++pos, data[j]);
      }
    }
    if (cols.pointingModelId().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.pointingModelId()(i)); // POINTING_MODEL_ID
    }

    stmt->setInt(++pos, cols.tracking()(i)); // TRACKING

    if (cols.onSource().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.onSource()(i)); // ON_SOURCE
    }
    if (cols.overTheTop().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.overTheTop()(i)); // OVER_THE_TOP
    }

    assert(pos == stmt->getParameterCount());
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
}

void saveSysCal(Connection *con, MeasurementSet &ms) {
  enter();
  MSSysCal &tab = ms.sysCal();
  ROMSSysCalColumns const cols(tab);

  char const *dbcols[] = {
    "SYSCAL_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "ANTENNA_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "FEED_ID", // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "SPECTRAL_WINDOW_ID", // INTEGER NOT NULL,
    "TIME", // REAL NOT NULL,
    "INTERVAL", // REAL NOT NULL,
    "PHASE_DIFF", // REAL DEFAULT NULL,
    "TCAL", // BLOB
    "TRX", // BLOB
    "TSKY", // BLOB
    "TSYS", // BLOB
    "TANT", // BLOB
    "TANT_TSYS", // BLOB
    "TCAL_SPECTRUM", // BLOB
    "TRX_SPECTRUM", // BLOB
    "TSKY_SPECTRUM", // BLOB
    "TSYS_SPECTRUM", // BLOB
    "TANT_SPECTRUM", // BLOB
    "TANT_TSYS_SPECTRUM", // BLOB
    "PHASE_DIFF_FLAG", // INTEGER DEFAULT NULL,
    "TCAL_FLAG", // INTEGER NOT NULL,
    "TRX_FLAG", // INTEGER NOT NULL,
    "TSKY_FLAG", // INTEGER NOT NULL,
    "TSYS_FLAG", // INTEGER NOT NULL,
    "TANT_FLAG", // INTEGER NOT NULL,
    "TANT_TSYS_FLAG", // INTEGER NOT NULL,
    NULL
  };

  string sql = "insert into SYSCAL (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    int pos = 0;
    stmt->setInt(++pos, i); // SYSCAL_ID
    stmt->setInt(++pos, cols.antennaId()(i)); // ANTENNA_ID
    stmt->setInt(++pos, cols.feedId()(i)); // FEED_ID
    stmt->setInt(++pos, cols.spectralWindowId()(i)); // SPECTRAL_WINDOW_ID

    stmt->setDouble(++pos, cols.time()(i)); // TIME
    stmt->setDouble(++pos, cols.interval()(i)); // INTERVAL

    if (cols.phaseDiff().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setFloat(++pos, cols.phaseDiff()(i)); // PHASE_DIFF
    }

    if (cols.tcal().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.tcal()(i));// TCAL 
    }
    if (cols.trx().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.trx()(i));// TRX
    }
    if (cols.tsky().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.tsky()(i));// TSKY
    }
    if (cols.tsys().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.tsys()(i));// TSYS
    }
    if (cols.tant().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.tant()(i));// TANT
    }
    if (cols.tantTsys().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.tantTsys()(i));// TANT_TSYS
    }
    if (cols.tcalSpectrum().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.tcalSpectrum()(i));// TCAL_SPECTRUM
    }
    if (cols.trxSpectrum().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.trxSpectrum()(i));// TRX_SPECTRUM
    }
    if (cols.tskySpectrum().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.tskySpectrum()(i));// TSKY_SPECTRUM
    }
    if (cols.tsysSpectrum().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.tsysSpectrum()(i));// TSYS_SPECTRUM
    }
    if (cols.tantSpectrum().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.tantSpectrum()(i));// TANT_SPECTRUM
    }
    if (cols.tantTsysSpectrum().isNull()) {
      stmt->setNull(++pos);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), ++pos, cols.tantTsysSpectrum()(i));// TANT_TSYS_SPECTRUM
    }

    if (cols.phaseDiffFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.phaseDiffFlag()(i)); // PHASE_DIFF_FLAG
    }
    if (cols.tcalFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.tcalFlag()(i)); // TCAL_FLAG
    }
    if (cols.trxFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.trxFlag()(i)); // TRX_FLAG
    }
    if (cols.tskyFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.tskyFlag()(i)); // TSKY_FLAG
    }
    if (cols.tsysFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.tsysFlag()(i)); // TSYS_FLAG
    }
    if (cols.tantFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.tantFlag()(i)); // TANT_FLAG
    }
    if (cols.tantTsysFlag().isNull()) {
      stmt->setNull(++pos);
    } else {
      stmt->setInt(++pos, cols.tantTsysFlag()(i)); // TANT_TSYS_FLAG
    }

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
    "BASELINE_REF", // INTEGER DEFAULT NULL, -- bool
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
    "FLAG_ROW", // INTEGER NOT NULL, -- bool
    NULL
  };

  enum {
    MAIN_ID = 1,
    TIME,
    TIME_EXTRA_PREC,
    ANTENNA1,
    ANTENNA2,
    ANTENNA3,
    FEED1,
    FEED2,
    FEED3,
    DATA_DESC_ID,
    PROCESSOR_ID,
    PHASE_ID,
    FIELD_ID,
    INTERVAL,
    EXPOSURE,
    TIME_CENTROID,
    PULSAR_BIN,
    PULSAR_GATE_ID,
    SCAN_NUMBER,
    ARRAY_ID,
    OBSERVATION_ID,
    STATE_ID,
    BASELINE_REF,
    U,
    V,
    W,
    U2,
    V2,
    W2,
    DATA,
    FLOAT_DATA,
    VIDEO_POINT,
    LAG_DATA,
    SIGMA,
    SIGMA_SPECTRUM,
    WEIGHT,
    WEIGHT_SPECTRUM,
    FLAG,
    FLAG_CATEGORY,
    FLAG_ROW
  };

  string sql = "insert into MAIN (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));

  uInt nrow = ms.nrow();
  for (uInt i = 0; i < nrow; i++) {
    stmt->setInt(MAIN_ID, i);
    stmt->setDouble(TIME, cols.time()(i));
    if (cols.timeExtraPrec().isNull()) {
      stmt->setNull(TIME_EXTRA_PREC);
    } else {
      stmt->setDouble(TIME_EXTRA_PREC, cols.timeExtraPrec()(i));
    }
    stmt->setInt(ANTENNA1, cols.antenna1()(i));
    stmt->setInt(ANTENNA2, cols.antenna2()(i));
    if (cols.antenna3().isNull()) {
      stmt->setNull(ANTENNA3);
    } else {
      stmt->setInt(ANTENNA3, cols.antenna3()(i));
    }
    stmt->setInt(FEED1, cols.feed1()(i));
    stmt->setInt(FEED2, cols.feed2()(i));
    if (cols.feed3().isNull()) {
      stmt->setNull(FEED3);
    } else {
      stmt->setInt(FEED3, cols.feed3()(i));
    }
    stmt->setInt(DATA_DESC_ID, cols.dataDescId()(i));
    stmt->setInt(PROCESSOR_ID, cols.processorId()(i));
    if (cols.phaseId().isNull()) {
      stmt->setNull(PHASE_ID);
    } else {
      stmt->setInt(PHASE_ID, cols.phaseId()(i));
    }
    stmt->setInt(FIELD_ID, cols.fieldId()(i));
    stmt->setDouble(INTERVAL, cols.interval()(i));
    stmt->setDouble(EXPOSURE, cols.exposure()(i));
    stmt->setDouble(TIME_CENTROID, cols.timeCentroid()(i));
    if (cols.pulsarBin().isNull()) {
      stmt->setNull(PULSAR_BIN);
    } else {
      stmt->setInt(PULSAR_BIN, cols.pulsarBin()(i));
    }
    if (cols.pulsarGateId().isNull()) {
      stmt->setNull(PULSAR_GATE_ID);
    } else {
      stmt->setInt(PULSAR_GATE_ID, cols.pulsarGateId()(i));
    }
    stmt->setInt(SCAN_NUMBER, cols.scanNumber()(i));
    stmt->setInt(ARRAY_ID, cols.arrayId()(i));
    stmt->setInt(OBSERVATION_ID, cols.observationId()(i));
    stmt->setInt(STATE_ID, cols.stateId()(i));
    if (cols.baselineRef().isNull()) {
      stmt->setNull(BASELINE_REF);
    } else {
      stmt->setInt(BASELINE_REF, cols.baselineRef()(i) ? 1 : 0);
    }
    {
      Array< Double > const &v = cols.uvw()(i);
      // const IPosition &shape = v.shape();
      assert(v.nelements() == 3);
      Double const *data = v.data();
      stmt->setDouble(U, data[0]);
      stmt->setDouble(V, data[1]);
      stmt->setDouble(W, data[2]);
    }
    if (cols.uvw2().isNull()) {
      stmt->setNull(U2);
      stmt->setNull(V2);
      stmt->setNull(W2);
    } else {
      Array< Double > const &v = cols.uvw2()(i);
      // const IPosition &shape = v.shape();
      assert(v.nelements() == 3);
      Double const *data = v.data();
      stmt->setDouble(U2, data[0]);
      stmt->setDouble(V2, data[1]);
      stmt->setDouble(W2, data[2]);
    }
    if (cols.data().isNull()) {
      stmt->setNull(DATA);
    } else {
      bindArrayAsBlob<Complex>(stmt.get(), DATA, cols.data()(i));
    }
    if (cols.floatData().isNull()) {
      stmt->setNull(FLOAT_DATA);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), FLOAT_DATA, cols.floatData()(i));
    }
    if (cols.videoPoint().isNull()) {
      stmt->setNull(VIDEO_POINT);
    } else {
      bindArrayAsBlob<Complex>(stmt.get(), VIDEO_POINT, cols.videoPoint()(i));
    }
    if (cols.lagData().isNull()) {
      stmt->setNull(LAG_DATA);
    } else {
      bindArrayAsBlob<Complex>(stmt.get(), LAG_DATA, cols.lagData()(i));
    }
    bindArrayAsBlob<Float>(stmt.get(), SIGMA, cols.sigma()(i));
    if (cols.sigmaSpectrum().isNull()) {
      stmt->setNull(SIGMA_SPECTRUM);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), SIGMA_SPECTRUM, cols.sigmaSpectrum()(i));
    }
    bindArrayAsBlob<Float>(stmt.get(), WEIGHT, cols.weight()(i));
    if (cols.weightSpectrum().isNull()) {
      stmt->setNull(WEIGHT_SPECTRUM);
    } else {
      bindArrayAsBlob<Float>(stmt.get(), WEIGHT_SPECTRUM, cols.weightSpectrum()(i));
    }
    bindArrayAsBlob<Bool>(stmt.get(), FLAG, cols.flag()(i));
    bindArrayAsBlob<Bool>(stmt.get(), FLAG_CATEGORY, cols.flagCategory()(i));
    stmt->setInt(FLAG_ROW, cols.flagRow()(i) ? 1 : 0);
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
}

void saveMainCoordinate(Connection *con, MeasurementSet &ms) {
  enter();
  // There is no corresponding table in MeasurementSet.
  // Nothing to do.
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

void saveDoppler(Connection *con, MeasurementSet &ms) {
  // Note: DOPPLER table is optional
  enter();
  MSDoppler &tab = ms.doppler();
  ROMSDopplerColumns const cols(tab);

  // return if optional DOPPLER table doesn't exist
  if (tab.isNull()) return;

  char const *dbcols[] = {
    "DOPPLER_ID",  // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "SOURCE_ID",  // INTEGER NOT NULL -> REFERENCES SOURCE (SOURCE_ID)
    "TRANSITION_ID", // INTEGER NOT NULL,
    "VELDEF",  // REAL NOT NULL,
    NULL
  };

  enum {
    DOPPLER_ID = 1,
    SOURCE_ID,
    TRANSITION_ID,
    VELDEF
  };

  string sql = "insert into DOPPLER (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));

  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    stmt->setInt(DOPPLER_ID, i);
    stmt->setInt(SOURCE_ID, cols.sourceId()(i));
    stmt->setInt(TRANSITION_ID, cols.transitionId()(i));
    stmt->setDouble(VELDEF, cols.velDef()(i));

    int result = stmt->executeUpdate();
    assert(result == 1);
  }
}

void saveHistory(Connection *con, MeasurementSet &ms) {
  enter();
  MSHistory &tab = ms.history();
  ROMSHistoryColumns const cols(tab);

  // HISTORY table
  {
  char const *dbcols[] = {
    "HISTORY_ID",  // INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
    "TIME",  // REAL NOT NULL,
    "OBSERVATION_ID",  // INTEGER NOT NULL,
    "MESSAGE",  // TEXT NOT NULL,
    "PRIORITY",  // TEXT NOT NULL,
    "ORIGIN",  // TEXT NOT NULL,
    "OBJECT_ID",  // INTEGER NOT NULL,
    "APPLICATION",  // TEXT NOT NULL
    NULL
  };

  enum {
    HISTORY_ID = 1,
    TIME,
    OBSERVATION_ID,
    MESSAGE,
    PRIORITY,
    ORIGIN,
    OBJECT_ID,
    APPLICATION
  };

  string sql = "insert into HISTORY (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    stmt->setInt(HISTORY_ID, i);
    stmt->setDouble(TIME, cols.time()(i));
    stmt->setInt(OBSERVATION_ID, cols.observationId()(i));
    stmt->setTransientString(MESSAGE, cols.message()(i).c_str());
    stmt->setTransientString(PRIORITY, cols.priority()(i).c_str());
    stmt->setTransientString(ORIGIN, cols.origin()(i).c_str());
    stmt->setInt(OBJECT_ID, cols.objectId()(i));
    stmt->setTransientString(APPLICATION, cols.application()(i).c_str());

    int result = stmt->executeUpdate();
    assert(result == 1);
  }
  }

  // HISTORY_CMD table
  {
  char const *dbcols[] = {
    "HISTORY_ID",  // INTEGER NOT NULL --> REFERENCES HISTORY (HISTORY_ID),
    "IDX",  // INTEGER NOT NULL (Starts with 0),
    "CLI_COMMAND",  // TEXT NOT NULL,
    NULL
  };

  enum {
    HISTORY_ID = 1,
    IDX,
    CLI_COMMAND
  };

  string sql = "insert into HISTORY_COMMAND (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    if ( ! cols.cliCommand().isDefined(i) ) continue;
    Array< String > const &v = cols.cliCommand()(i); // CLI_COMMAND
    String const *data = v.data();
    for (size_t j = 0; j < v.nelements(); j++) {
      stmt->setInt(HISTORY_ID, i);
      stmt->setInt(IDX, j);
      stmt->setTransientString(CLI_COMMAND, data[j].c_str());

      int result = stmt->executeUpdate();
      assert(result == 1);
    }
  }
  }

  // HISTORY_PARAM
  {
  char const *dbcols[] = {
    "HISTORY_ID",  // INTEGER NOT NULL -> REFERENCES HISTORY (HISTORY_ID),
    "IDX",  // INTEGER NOT NULL (Starts with 0),
    "APP_PARAMS",  // TEXT NOT NULL,
    NULL
  };

  enum {
    HISTORY_ID = 1,
    IDX,
    APP_PARAMS
  };

  string sql = "insert into HISTORY_PARAM (";
  sql += SQL::join(dbcols) + ") values (";
  sql += SQL::bindChars(dbcols) + ")";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  uInt nrow = tab.nrow();
  for (uInt i = 0; i < nrow; i++) {
    if ( ! cols.appParams().isDefined(i) ) continue;
    Array< String > const &v = cols.appParams()(i); // APP_PARAMS
    String const *data = v.data();
    for (size_t j = 0; j < v.nelements(); j++) {
      stmt->setInt(HISTORY_ID, i);
      stmt->setInt(IDX, j);
      stmt->setTransientString(APP_PARAMS, data[j].c_str());

      int result = stmt->executeUpdate();
      assert(result == 1);
    }
  }
  }

}

void mssave(Connection *con, char const*filename) {
  enter();
  static void (*funcs[])(Connection *, MeasurementSet &) = {
    saveTableKeyword,
    saveState, saveFlagCmd, saveAntenna, saveSw,
    savePol,
    saveProcessor, 
    saveHistory,
    saveObs,
    //saveSource, // depending on SpectralWindow
    saveWeather, savePointing, // depending on Antenna
    saveDataDesc, // depending on SpectralWindow, Polarization
    saveField, // depending on Source
    saveDoppler, // depending on Source
    saveFeed, // depending on Antenna, SpectralWindow
    saveFreqOffset, saveSysCal, // depending on Feed, Antenna, SpectralWindow
    saveMain, // depending on Processor, State, DataDescription, Field, Feed, Antenna
    saveMainCoordinate, // depending on Main
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
  //int result = gettimeofday(&tv, NULL);
  gettimeofday(&tv, NULL);
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
  msm += ".mdb";
  string mst = basename;
  mst += ".tdb";

  char const sql[] = "/sql/";
  // create empty master db from msm
  {
    char const ddl_file[] = "MSM.ddl";
    size_t ddl_path_size = strlen(prefix) + strlen(sql) + strlen(ddl_file) + 1;
    char ddl_path[ddl_path_size];
    snprintf(ddl_path, ddl_path_size, "%s%s%s", prefix, sql, ddl_file);
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
  char const *prefix = PREFIX;

  char const *root = getenv(ENV_VAR_SAKURA_ROOT);
  if (root != NULL) {
    prefix = root;
  }

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
  cout << "SAKURA_ROOT: " << prefix << endl;
  double start = currenttime();
  conv(prefix, argv[argStart], argv[argStart + 1]);
  double end = currenttime();
  cout << end - start << "sec\n";
  return 0;
}
