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

void fetchAllPointing_(MeasurementSet &ms) {
  enter();
  MSPointing &tab = ms.pointing();
  ROMSPointingColumns const cols(tab);
  uInt nrow = tab.nrow();

  const ROScalarColumn< Double > &timeCol = cols.time();
  const ROScalarColumn< Double > &intervalCol = cols.interval();
  const ROScalarColumn< String > &nameCol = cols.name();
  const ROScalarColumn< Int > &numPolyCol = cols.numPoly(); 
  const ROScalarColumn< Double > &timeOrigCol = cols.timeOrigin(); 
  const ROArrayColumn< Double > &directionCol = cols.direction();
  const ROArrayColumn< Double > &targetCol = cols.target();
  const ROArrayColumn< Double > &pointingOffsetCol = cols.pointingOffset(); 
  const ROArrayColumn< Double > &sourceOffsetCol = cols.sourceOffset(); 
  const ROArrayColumn< Double > &encoderCol = cols.encoder();
  const ROScalarColumn< Int > &pointingModelIdCol = cols.pointingModelId(); 
  const ROScalarColumn< bool > &trackingCol = cols.tracking(); 
  const ROScalarColumn< bool > &onSourceCol = cols.onSource(); 
  const ROScalarColumn< bool > &overTheTopCol = cols.overTheTop(); 

  double start = currenttime();
  for (uInt i = 0; i < nrow; i++) {
    Double dt = timeCol(i);// Time
    Double dinte = intervalCol(i);// INTERVAL
    const String& snam = nameCol(i);// NAME
    Int unpol = numPolyCol(i);// NUM_POLY
    Double dto = timeOrigCol(i);// TIME ORIGIN
    const Array< Double >& t_d = directionCol(i);// DIRECTION
    const Array< Double >& t_t = targetCol(i);// TARGET

    if (cols.pointingOffset().isNull()) {
    }else{
      const Array< Double >& t_poff = pointingOffsetCol(i);//POINTING_OFFSET
    }
    if (cols.sourceOffset().isNull()) {
    }else{
      const Array< Double >& t_soff = sourceOffsetCol(i); //SOURCE_OFFSET
    }
    if (cols.encoder().isNull()) {
    }else{
      const Array< Double >& t_enco = encoderCol(i); // ENCODER
    }
    if (cols.pointingModelId().isNull()) {
    }else{
      Int ipmodelid = pointingModelIdCol(i);  // POINTING_MODEL_ID
    }

    Bool itracki = trackingCol(i); // TRACKING

    if (cols.onSource().isNull()) {
    }else{
      Bool ionso = onSourceCol(i); // ON_SOURCE
    }
    if (cols.overTheTop().isNull()) {
    }else{
      Bool iovtop = overTheTopCol(i); // OVER_THE_TOP
    }
  }
  double end = currenttime();
  cout << "Fetched: " << end - start << "sec\n";
}

void fetchAllPointing(char const*filename) {
  enter();
  MeasurementSet ms(filename);
  fetchAllPointing_(ms);
}

struct Entry {
  char const *option;
  void (*func)(char const*filename);
} entries[] = {
  {"pointing", fetchAllPointing}
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
