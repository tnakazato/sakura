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
#define elementsof(x) (sizeof(x) / sizeof(*(x)))

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

double currenttime() {
  struct timeval tv;
  int result = gettimeofday(&tv, NULL);
  return tv.tv_sec + ((double)tv.tv_usec) / 1000000.;
}

void timing(string msg, double t) {
  cout << "Timing: " << msg << ": " << t << " sec\n";
}

void insertToPointingSQL_(Connection *con,MeasurementSet &ms)
{
  enter();
  // MeasurementSet ms(filename);
  MSPointing &tab = ms.pointing();
  //  ROMSPointingColumns rocols(tab);

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
  //  uInt nrow = tab.nrow();
  uInt nrow =48000;

  for (uInt i = 0; i < nrow; i++) {
//     int pos = 0;
//     stmt->setInt(++pos, i); // ANTENNA_ID
//     stmt->setDouble(++pos, 2012.0); // TIME
//     stmt->setDouble(++pos,3.14 ); // INTERVAL
//     stmt->setTransientString(++pos,"test-kawa" ); // NAME
//     stmt->setInt(++pos, 1); // NUM_POLY

//     stmt->setDouble(++pos, 360000); // TIME_ORIGIN
//     //bindArrayAsBlob<Double>(stmt.get(), ++pos, 3.14);//DIRECTION
//     //bindArrayAsBlob<Double>(stmt.get(), ++pos, 3.141592);//TARGET

//     if (cols.pointingOffset().isNull()) {
//       stmt->setNull(++pos);
//     } else {
//       //bindArrayAsBlob<Double>(stmt.get(), ++pos, 4.154);//POINTING_OFFSET
//     }
//     if (cols.sourceOffset().isNull()) {
//       stmt->setNull(++pos);
//     } else {
//       //bindArrayAsBlob<Double>(stmt.get(), ++pos, 5.134);//SOURCE_OFFSET
//     }
//     if (cols.encoder().isNull()){
//       stmt->setNull(++pos);
//       stmt->setNull(++pos);
//     }else {
//       Array< Double > const &v = cols.encoder()(i); // ENCODERX ENCODERY
//       //const IPosition &shape = v.shape();
//       Double const *data = v.data();
//       for (size_t j = 0; j < 2; j++) {
// 	// "ENCODERX", "ENCODERY"
// 	stmt->setDouble(++pos, 5.14);
//       }
//     }
//     if (cols.pointingModelId().isNull()) {
//       stmt->setNull(++pos);
//     } else {
//       stmt->setInt(++pos, 9); // POINTING_MODEL_ID
//     }

//     stmt->setInt(++pos, 1); // TRACKING

//     if (cols.onSource().isNull()) {
//       stmt->setNull(++pos);
//     } else {
//       stmt->setInt(++pos, 1); // ON_SOURCE
//     }
//     if (cols.overTheTop().isNull()) {
//       stmt->setNull(++pos);
//     } else {
//       stmt->setInt(++pos, 1); // OVER_THE_TOP
//     }

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

  //void insertToPointingSQL(Connection *con, char const*filename) {
  // enter();
  // MeasurementSet ms(filename);
  // insertToPointingSQL_(con->get(),ms);
  //}

void insertToPointingSQL(Connection *con, char const*filename) {
  enter();
  static void (*funcs[])(Connection *, MeasurementSet &) = {
    insertToPointingSQL_, // depending on Antenna
    NULL
  };
  MeasurementSet ms(filename);
  for (size_t i = 0; funcs[i] != NULL; i++) {
    con->execute("BEGIN");
    funcs[i](con, ms);
    con->execute("COMMIT");
  }
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
    insertToPointingSQL(con.get(), msfile);
  }
}

struct Entry {
  char const *option;
  void (*func)(Connection *con, char const*filename);
} entries[] = {
   {"pointing", insertToPointingSQL}
};

char const *progName = "";

void usage() {
  cerr << "Usage: " << progName << " [options] action tdbFile\n";
  cerr << "options:: \n";
  cerr << "\t--prefix path\tA path where sakura is installed.\n";
  cerr << "\t-p path\n";
  cerr << "\taction action is one of followings\n";
  for (size_t i = 0; i < elementsof(entries); i++) {
    cerr << "\t\t" << entries[i].option << endl;
  }
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


// int main(int argc, char const *argv[]) { 
//   progName = argv[0];
//   char const *prefix = PREFIX;

//   char const *root = getenv(ENV_VAR_SAKURA_ROOT);
//   if (root != NULL) {
//     prefix = root;
//   }

//   static struct option const long_options[] = {
//     {"prefix", 1, NULL, 'p'},
//     {0, 0, NULL, 0}
//   };

//   for (;;) {
//     int option_index = 0;
//     int optCh = getopt_long (argc, const_cast<char *const *>(argv), "p:",
// 			     long_options, &option_index);
//     if (optCh == -1) {
//       break;
//     }
//     switch (optCh) {
//     case 'p':
//       prefix = optarg;
//       break;
//     case '?':
//       usage();
//       return 1;
//     default:
//       assert(false);
//       return 1;
//     }
//   }

//   int argStart = optind;
//   if (argc - argStart != 2) {
//     usage();
//     return 1;
//   }
//   cout << "SAKURA_ROOT: " << prefix << endl;
//   void (*func)(Connection *con, char const*filename) = NULL;
//   for (size_t i = 0; i < elementsof(entries); i++) {
//     if (strcmp(entries[i].option, argv[argStart]) == 0) {
//       func = entries[i].func;
//     }
//   }
//   if (func == NULL) {
//     usage();
//     return 1;
//   }
//   {
//     double start = currenttime();
//     string tdb = "file:";
//     tdb += argv[argStart + 1];
//     tdb += "?mode=ro";

//     cout << "Opening: " << tdb << endl;
//     unique_ptr<Connection> con(Connection::open(tdb.c_str()));
//     func(con.get(), argv[argStart + 1]);
//     // conv(prefix, argv[argStart], argv[argStart + 1]);
//     double end = currenttime();
//     timing("Total: ", end - start);
//   }
//   return 0;
// }
