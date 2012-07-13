#include <cstdio>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <string>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <sys/time.h>
#include <sstream>
#include "SQLite.h"

using namespace std;
using namespace sqlite;

#define unique_ptr auto_ptr
#define enter() do { cout << "Enter: " << __FUNCTION__ << endl; } while (0)
#define elementsof(x) (sizeof(x) / sizeof(*(x)))

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

void fetchStringAndReal( Connection *con, char const *sql )
{
  enter();
  unique_ptr<PreparedStatement> stmt(con->prepare(sql));
  unique_ptr<ResultSet> rs(stmt->executeQuery());
  int n = 0 ;
  int count = 0 ;
  while (rs->next()) {
    const char *c = rs->getTransientString( 1, &n ) ;
    double d = rs->getDouble( 2 ) ;
    //cout << c << endl ;
    //cout << setprecision(14) << d << endl;
    count++ ;
  }
  cout << "fetched " << count << " records" << endl ;
}

void fetchRealAndInt( Connection *con, char const *sql )
{
  enter();
  unique_ptr<PreparedStatement> stmt(con->prepare(sql));
  unique_ptr<ResultSet> rs(stmt->executeQuery());
  int count = 0 ;
  while (rs->next()) {
    double d = rs->getDouble( 1 ) ;
    int i = rs->getInt( 2 ) ;
    //cout << setprecision(14) << d << endl ;
    //cout << i << endl ;
    count++ ;
  }
  cout << "fetched " << count << " records" << endl ;
}

struct Entry {
  char const *option;
  char const *sql;
  char const *key;
  void (*func)(Connection *con, char const *sql);
} entries[] = {
  {"timeAndAntennaByTime", "select time, antenna_id from pointing ?;", "time", fetchRealAndInt},
  {"nameAndIntervalByTime", "select name, interval from pointing ?;", "time", fetchStringAndReal},
  {"nameAndIntervalByTimeOrigin", "select name, interval from pointing ?;", "time_origin", fetchStringAndReal}  
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

int main(int argc, char const *argv[]) { 
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
  void (*func)(Connection *con, char const *sql) = NULL;
  char const *sql = NULL;
  char const *key = NULL;
  for (size_t i = 0; i < elementsof(entries); i++) {
    if (strcmp(entries[i].option, argv[argStart]) == 0) {
      func = entries[i].func;
      sql = entries[i].sql;
      key = entries[i].key;
    }
  }
  if (func == NULL) {
    usage();
    return 1;
  }
  double fromTime = 4819385780.3519993 ; 
  double toTime = 4819386236.5440006 ;
  stringstream oss ;
  oss << " where " << key << " >= " 
      << setprecision(16) << fromTime 
      << " and " << key << " <= " 
      << setprecision(16) << toTime ;
  string sel = oss.str() ;
  {
    double start = currenttime();
    string tdb = "file:";
    tdb += argv[argStart + 1];
    tdb += "?mode=ro";

    string sqlsel = sql ;
    sqlsel.replace( sqlsel.find( "?" ), 1, sel ) ;
    cout << "sqlsel: " << sqlsel << endl ;

    cout << "Opening: " << tdb << endl;
    unique_ptr<Connection> con(Connection::open(tdb.c_str()));
//     func(con.get(), sql);
    func(con.get(), sqlsel.c_str());
    double end = currenttime();
    timing("Total: ", end - start);
  }
  return 0;
}
