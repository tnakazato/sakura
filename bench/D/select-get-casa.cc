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
  //int result = gettimeofday(&tv, NULL);
  gettimeofday(&tv, NULL);
  return tv.tv_sec + ((double)tv.tv_usec) / 1000000.;
}

template<class T> 
void fetch_( const Table &tab, const String colname )
{
//   double start = currenttime();
  ROScalarColumn<T> col( tab, colname ) ;
  Vector<T> val = col.getColumn() ;
//   uInt len = val.nelements() ;
//   double end = currenttime();
//   cout << colname << ": Fetched " << len << " rows within " 
//        << end - start << "sec\n";
}
 
template<class T, class U>
void fetch2_( const Table &tab, const String name1, const String name2 )
{
//   double start = currenttime();
  ROScalarColumn<T> col1( tab, name1 ) ;
  ROScalarColumn<U> col2( tab, name2 ) ;
  Vector<T> val1 = col1.getColumn() ;
  Vector<U> val2 = col2.getColumn() ;
//   uInt len = val1.nelements() ;
//   double end = currenttime();
//   cout << name1 << "," << name2 << ": Fetched " << len << " rows within " 
//        << end - start << "sec\n";
}

Table select_( char const*filename,
	       double from,
	       double to,
	       const String keycol )
{
  String tablename = String(filename) + "/POINTING" ;
  MSPointing msp( tablename ) ;
  stringstream oss ;
  oss << "SELECT FROM $1 WHERE " << keycol << " >= " 
      << setprecision(16) << from
      << " AND " << keycol << " <= " 
      << setprecision(16) << to ;
  cout << "TaQL: " << oss.str() << endl ;
  const String command( oss.str() ) ;
  Table mspsel = tableCommand( command, msp ) ;
  return mspsel ;
}
		      

void fetchTimeAndAntennaByTime( char const*filename,
				double from,
				double to )
{
  Table tsel = select_( filename, from, to, "TIME" ) ;
//   fetch_<Double>( tsel, "TIME" ) ;
//   fetch_<Int>( tsel, "ANTENNA_ID" ) ;
  fetch2_<Double,Int>( tsel, "TIME", "ANTENNA_ID" ) ;
}

void fetchNameAndIntervalByTime( char const*filename,
				 double from, 
				 double to )
{
  Table tsel = select_( filename, from, to, "TIME" ) ;
//   fetch_<String>( tsel, "NAME" ) ;
//   fetch_<Double>( tsel, "INTERVAL" ) ;
  fetch2_<String,Double>( tsel, "NAME", "INTERVAL" ) ;
}

void fetchNameAndIntervalByTimeOrigin( char const*filename,
				       double from,
				       double to )
{
  MSPointing tsel = select_( filename, from, to, "TIME_ORIGIN" ) ;
//   fetch_<String>( tsel, "NAME" ) ;
//   fetch_<Double>( tsel, "INTERVAL" ) ;
  fetch2_<String,Double>( tsel, "NAME", "INTERVAL" ) ;
}



struct Entry {
  char const *option;
  void (*func)(char const*filename, double, double);
} entries[] = {
  {"timeAndAntennaByTime", fetchTimeAndAntennaByTime},
  {"nameAndIntervalByTime", fetchNameAndIntervalByTime},
  {"nameAndIntervalByTimeOrigin", fetchNameAndIntervalByTimeOrigin}
};

char const *progName = "";

void usage() {
  cerr << "Usage: " << progName << " [options] action MSFile\n";
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
  double fromTime = 4819385780.3519993 ; 
  double toTime = 4819386236.5440006 ;
  double start = currenttime();
  func(argv[argStart + 1], fromTime, toTime);
  double end = currenttime();
  cout << "Total: " << end - start << "sec\n";
  return 0;
}
