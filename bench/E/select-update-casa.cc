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
#include <iomanip>
#include <sstream>

#include <casacore/casa/aips.h>
#include <casacore/casa/Inputs/Input.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Utilities/CountedPtr.h>
#include <casacore/casa/OS/Directory.h>
#include <casacore/tables/Tables/Table.h>
#include <tables/Tables/TableParse.h>
#include <casacore/ms/MeasurementSets.h>

#define unique_ptr auto_ptr
#define enter() do { cout << "Enter: " << __FUNCTION__ << endl; } while (0)
#define elementsof(x) (sizeof(x) / sizeof(*(x)))

using namespace std;
using namespace casa;

namespace {
char const ENV_VAR_SAKURA_ROOT[] = "SAKURA_ROOT";

double currenttime() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + ((double)tv.tv_usec) / 1000000.;
}


Int newAntennaId = 4567;
String newName = "NGC4038/9";

/*
void update_( const char* filename, const String keyCol, double from, double to, const String colName, Int newVal )
{
  String tablename = String(filename) + "/POINTING" ;
  MSPointing msp( tablename ) ;
  stringstream oss ;
  oss << "UPDATE $1 SET " << colName << " = " << newVal 
      << " WHERE " << keyCol << " >= " 
      << setprecision(16) << from
      << " AND " << keyCol << " <= " 
      << setprecision(16) << to ;
  cout << "TaQL: " << oss.str() << endl ;
  const String command( oss.str() ) ;
  tableCommand( command, msp ) ;
}
*/
//void update_( const char* filename, const String keyCol, double from, double to, const String colName, String newVal )
template<class T>
void update_( const char* filename, const String keyCol, double from, double to, const String colName, T newVal, bool isString )
{
  String tablename = String(filename) + "/POINTING" ;
  MSPointing msp( tablename ) ;
  stringstream oss ;
  oss << "UPDATE $1 SET " << colName << " = " 
      << (isString ? "'" : "") << newVal << (isString ? "'" : "")
      << " WHERE " << keyCol << " >= " 
      << setprecision(16) << from
      << " AND " << keyCol << " <= " 
      << setprecision(16) << to ;
  cout << "TaQL: " << oss.str() << endl ;
  const String command( oss.str() ) ;
  tableCommand( command, msp ) ;
}

void updateAntennaByTime( char const*filename,
			  double from,
			  double to )
{
  update_( filename, "TIME", from, to, "ANTENNA_ID", newAntennaId, false ) ;
}
void updateNameByTime ( char const*filename,
			double from,
			double to )
{
  update_( filename, "TIME", from, to, "NAME", newName, true ) ;
}



struct Entry {
  char const *option;
  void (*func)(char const*filename, double, double);
} entries[] = {
  {"antennaByTime", updateAntennaByTime},
  {"nameByTime", updateNameByTime}
};

char const *progName = "";

void usage() {
  cerr << "Usage: " << progName << " [options] action MSFile lowerRange upperRange\n";
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
  if (argc - argStart != 4) {
    usage();
    return 1;
  }
  cout << "SAKURA_ROOT: " << prefix << endl;
  void (*func)(char const*filename, double, double) = NULL;
  for (size_t i = 0; i < elementsof(entries); i++) {
    if (strcmp(entries[i].option, argv[argStart]) == 0) {
      func = entries[i].func;
    }
  }
  if (func == NULL) {
    usage();
    return 1;
  }
  double fromTime = atof(argv[argStart + 2]); 
  double toTime   = atof(argv[argStart + 3]);
  double start = currenttime();
  func(argv[argStart + 1], fromTime, toTime);
  double end = currenttime();
  cout << "Total: " << end - start << "sec\n";
  return 0;
}
