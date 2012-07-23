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

//#include <casacore/casa/aips.h>
//#include <casacore/casa/Inputs/Input.h>
#include <casacore/casa/BasicSL/String.h>
//#include <casacore/casa/Utilities/CountedPtr.h>
//#include <casacore/casa/OS/Directory.h>
#include <casacore/tables/Tables/Table.h>
//#include <casacore/ms/MeasurementSets.h>
#include <casacore/tables/Tables/ScaColDesc.h>

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

//////////////////// Actual Test Functions /////////////////////////////

void clearCol_(Table &tab, const String &colname) {
  cout << "---> Trying clear up a column '" << colname << "'." << endl;
  if (tab.canRemoveColumn(colname)){
    tab.removeColumn(colname);
    cout << "---> Column '" << colname << "' is successfully removed." << endl;
  } else {
    cout << "Could not remove column '" << colname << "'." << endl;
    assert(false);
  }
}

void addIntCol(char const*filename) {
  cout << "---> addIntCol: filename = " << filename << endl;
  const String colname( "FILLED_INT" );
  Table tab( filename, Table::Update );

  // Check if colname to add already exists
  TableDesc tdorg = tab.tableDesc();
  if ( tdorg.isColumn(colname) ) {
    cout << "---> Column '" << colname << "' already exists." << endl;
    clearCol_(tab, colname);
  }

  // Now actually add Column
  cout << "---> Add a column '" << colname << "' and fill with 1." << endl;
  double start = currenttime();
  tab.addColumn( ScalarColumnDesc<Int>(colname), False );
  ScalarColumn<Int> newCol( tab, colname );
  newCol.fillColumn(1);
  tab.flush(False,False);
  double end = currenttime();
  cout << "Added: " << end - start << " sec" << endl;

#if 0
  // Value check 
  uInt nRow = tab.nrow();
  cout << "Total rows: " << nRow << endl;
  ROScalarColumn<Int> roCol( tab, colname );
  cout << "---> Checking values in the column ..." << endl;
  for ( uInt i = 0; i < nRow; i++ ){
    assert(roCol.asInt(i)==1);
  }
  cout << "---> ... value check successful!!" << endl;
 
  // Final clear up
  clearCol_(tab, colname);
#endif
}

void addNullCol(char const*filename) {
  cout << "---> addNullCol: filename = " << filename << endl;
  const String colname( "DEFAULT_INT" );
  Table tab( filename, Table::Update );

  // Check if colname to add already exists
  TableDesc tdorg = tab.tableDesc();
  if ( tdorg.isColumn(colname) ) {
    cout << "---> Column '" << colname << "' already exists." << endl;
    clearCol_(tab, colname);
  }

  // Now actually add Column
  int defval = int();
  cout << "---> Add a column '" << colname << "' (initialized by " 
       << defval << ")." << endl;
  double start = currenttime();
  tab.addColumn( ScalarColumnDesc<Int>(colname), False );
  //ScalarColumn<int> newCol( tab, colname );
  //newCol.fillColumn(NULL);
  tab.flush(False,False);
  double end = currenttime();
  cout << "Added: " << end - start << "sec" << endl;

#if 0
  // Value check 
  ROScalarColumn<Int> roCol( tab, colname );
  cout << "---> Checking values in the column ..." << endl;
  uInt nRow = tab.nrow();
  for ( uInt i = 0; i < nRow; i++ ){
    assert(roCol.asInt(i)==defval);
  }
  cout << "---> ... value check successful!!" << endl;
 
  // Final clear up
  clearCol_(tab, colname);
#endif
}




////////////////////////////////////////////////////////////////////////
struct Entry {
  char const *option;
  void (*func)(char const*filename);
} entries[] = {
  {"int", addIntCol},
  {"null", addNullCol}
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

  //int iopt = 0;
  for (;;) {
    int option_index = 0;
    int optCh = getopt_long (argc, const_cast<char *const *>(argv), "p:",
			     long_options, &option_index);
    //cout << "Opt " << iopt << ": " << optCh << endl;
    //iopt++
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
  void (*func)(char const*filename) = NULL;
  for (size_t i = 0; i < elementsof(entries); i++) {
    if (strcmp(entries[i].option, argv[argStart]) == 0) {
      func = entries[i].func;
    }
  }
  if (func == NULL) {
    usage();
    return 1;
  }
  double start = currenttime();
  func(argv[argStart + 1]);
  double end = currenttime();
  cout << "Total: " << end - start << "sec\n";
  return 0;
}
