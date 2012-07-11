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

#include <casacore/casa/aips.h>
#include <casacore/casa/Inputs/Input.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Utilities/CountedPtr.h>
#include <casacore/casa/OS/Directory.h>
#include <casacore/tables/Tables/Table.h>
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

void fetchAllTime_(MeasurementSet &ms) {
  enter();
  ROMSColumns const cols(ms);
  const ROScalarColumn< Double > &timeCol = cols.time();
  uInt nrow = ms.nrow();
  double start = currenttime();
  for (uInt i = 0; i < nrow; i++) {
    Double t = timeCol(i);
    //cout << setprecision(14) << t << endl;
  }
  double end = currenttime();
  cout << "Fetched: " << end - start << "sec\n";
}

void fetchAllTime(char const*filename) {
  enter();
  MeasurementSet ms(filename);
  fetchAllTime_(ms);
}

void fetchAllTimeSort(char const*filename) {
  enter();
  MeasurementSet ms_(filename);
  MeasurementSet ms(ms_.sort("TIME"));
  fetchAllTime_(ms);
}

void fetchAllScanNumber_(MeasurementSet &ms) {
  enter();
  ROMSColumns const cols(ms);
  const ROScalarColumn< Int > &scanNumber = cols.scanNumber();
  uInt nrow = ms.nrow();
  double start = currenttime();
  for (uInt i = 0; i < nrow; i++) {
    Int t = scanNumber(i);
    //cout << t << endl;
  }
  double end = currenttime();
  cout << "Fetched: " << end - start << "sec\n";
}

void fetchAllScanNumber(char const*filename) {
  enter();
  MeasurementSet ms(filename);
  fetchAllScanNumber_(ms);
}

void fetchAllScanNumberSort(char const*filename) {
  enter();
  MeasurementSet ms_(filename);
  MeasurementSet ms(ms_.sort("SCAN_NUMBER"));
  fetchAllScanNumber_(ms);
}

struct Entry {
  char const *option;
  void (*func)(char const*filename);
} entries[] = {
  {"time", fetchAllTime},
  {"timeSort", fetchAllTimeSort},
  {"scanNumber", fetchAllScanNumber},
  {"scanNumberSort", fetchAllScanNumberSort}
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
