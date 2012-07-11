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

void fetchReal(Connection *con, char const *sql) {
  unique_ptr<PreparedStatement> stmt(con->prepare(sql));
  unique_ptr<ResultSet> rs(stmt->executeQuery());
  while (rs->next()) {
    double d = rs->getDouble(1);
    //cout << setprecision(14) << d << endl;
  }
}

void fetchInt(Connection *con, char const *sql) {
  unique_ptr<PreparedStatement> stmt(con->prepare(sql));
  unique_ptr<ResultSet> rs(stmt->executeQuery());
  while (rs->next()) {
    int n = rs->getInt(1);
    //cout << n << endl;
  }
}

struct Entry {
  char const *option;
  char const *sql;
  void (*func)(Connection *con, char const *sql);
} entries[] = {
  {"time", "select time from main;", fetchReal},
  {"timeSort", "select time from main order by time;", fetchReal},
  {"scanNumber", "select scan_number from main;", fetchInt},
  {"scanNumberSort", "select scan_number from main order by scan_number;", fetchInt}
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
  for (size_t i = 0; i < elementsof(entries); i++) {
    if (strcmp(entries[i].option, argv[argStart]) == 0) {
      func = entries[i].func;
      sql = entries[i].sql;
    }
  }
  if (func == NULL) {
    usage();
    return 1;
  }
  {
    double start = currenttime();
    string tdb = "file:";
    tdb += argv[argStart + 1];
    tdb += "?mode=ro";

    cout << "Opening: " << tdb << endl;
    unique_ptr<Connection> con(Connection::open(tdb.c_str()));
    func(con.get(), sql);
    double end = currenttime();
    timing("Total: ", end - start);
  }
  return 0;
}
