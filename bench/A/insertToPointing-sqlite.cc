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
#include <fstream>

#include "SQLite.h"

#define unique_ptr auto_ptr
#define enter() do { cout << "Enter: " << __FUNCTION__ << endl; } while (0)
#define elementsof(x) (sizeof(x) / sizeof(*(x)))

using namespace std;
using namespace sqlite;

namespace {
char const ENV_VAR_SAKURA_ROOT[] = "SAKURA_ROOT";

double currenttime() {
  struct timeval tv;
  int result = gettimeofday(&tv, NULL);
  return tv.tv_sec + ((double)tv.tv_usec) / 1000000.;
}

void timing(string msg, double t) {
  cout << "Timing: " << msg << ": " << t << " sec\n";
}

void insertToPointingSQL_(Connection *con, int nrow){
  enter();

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

  cout << "row# = " << nrow << endl;

  size_t elements = 2;
  double const data[2]={3.13347,-0.329554};
  double const *blobData = data;
 
  double start_insert = currenttime();
  for (int i = 0; i < nrow; i++) {
    int pos = 0;
    stmt->setInt(++pos, 0); // ANTENNA_ID
    stmt->setDouble(++pos, 9999999999.0+i*0.001); // TIME
    stmt->setDouble(++pos,1.152); // INTERVAL
    stmt->setTransientString(++pos,"Antennae" ); // NAME
    stmt->setInt(++pos, 0); // NUM_POLY
    stmt->setDouble(++pos, 4.8e+09); // TIME_ORIGIN
    stmt->setTransientBlob(++pos, blobData, sizeof(double) * elements);//DIRECTION
    stmt->setTransientBlob(++pos, blobData, sizeof(double) * elements);//TARGET
    stmt->setTransientBlob(++pos, blobData, sizeof(double) * elements);//POINTING_OFFSET
    stmt->setTransientBlob(++pos, blobData, sizeof(double) * elements);//SOURCE_OFFSET

    stmt->setDouble(++pos, 3.14);// ENCODERX
    stmt->setDouble(++pos, 0.4); // ENCODERY

    stmt->setInt(++pos, 0); // POINTING_MODEL_ID
    stmt->setInt(++pos, 1); // TRACKING
    stmt->setInt(++pos, 1); // ON_SOURCE
    stmt->setInt(++pos, 1); // OVER_THE_TOP
  
    assert(pos == stmt->getParameterCount());
    int result = stmt->executeUpdate();
    assert(result == 1);
  }
  double end_insert = currenttime();
  cout << "inserted(SQL):" << end_insert - start_insert <<" sec\n";
}

void insertToPointingSQL(Connection *con,char const *nrow) {
  enter();
  int nrow_ = atoi(nrow);
  con->execute("BEGIN");
  insertToPointingSQL_(con,nrow_);
  con->execute("COMMIT");
}

char *readFileContent(char const *filename) {
  enter();
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

  void conv(char const *prefix,char const *basename,char const *nrow,char const *outfilename) {
  enter();
  string msm = outfilename;
  msm += ".mdb";

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
      insertToPointingSQL(con.get(), nrow);
    } catch (...) {
      delete[] msm_ddl;
      throw;
    }
    delete[] msm_ddl;
  }
}

struct Entry {
  char const *option;
  void (*func)(char const *prefix,char const *basename,char const *nrow, char const *outfilename);
} entries[] = {
   {"pointing", conv}
};

char const *progName = "";

void usage() {
  cerr << "Usage: " << progName << " [options] action num_of_row outfilename\n";
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
  if (argc - argStart != 3) {
    usage();
    return 1;
  }
  cout << "SAKURA_ROOT: " << prefix << endl;
  double start = currenttime();
  conv(prefix,argv[argStart],argv[argStart+1],argv[argStart+2]);
  double end = currenttime();
  cout << "Total: " << end - start << "sec\n";
  return 0;
}
