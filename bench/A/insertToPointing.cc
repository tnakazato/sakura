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

#include <tables/Tables/ExprNode.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/TableIter.h>
#include <tables/Tables/RefRows.h>
#include <tables/Tables/TableRow.h>


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

void insertToPointing_(int nrow) {
  enter();

  casa::MeasurementSet *mstable_ ;
  TableDesc mspointingDesc = MSPointing::requiredTableDesc() ;

  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::ANTENNA_ID) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::TIME) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::INTERVAL ) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::NAME ) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::NUM_POLY ) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::TIME_ORIGIN ) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::DIRECTION ) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::TARGET ) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::POINTING_OFFSET ) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::SOURCE_OFFSET ) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::ENCODER ) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::POINTING_MODEL_ID ) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::TRACKING ) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::ON_SOURCE ) ;
  MSPointing::addColumnToDesc( mspointingDesc, MSPointingEnums::OVER_THE_TOP ) ;

  SetupNewTable pointingTab( "POINTING", mspointingDesc, Table::New ) ;

  cout << "row# = " << nrow << endl;

  MSPointing mspo(pointingTab);
  mspo.addRow( nrow, True );
  MSPointingColumns cols(mspo);

  Matrix<Double> mDirTarg(2,1);
  mDirTarg(0,0)=3.13347;
  mDirTarg(1,0)=-0.329554;
  Vector<Double> vEnco(2);
  vEnco[0]=3.14;
  vEnco[1]=0.40;

  Double intervalvalue=1.152;

  double start = currenttime();
  for (uInt i = 0; i < nrow; i++) {
    cols.antennaId().put(i,0) ; // ANTENNA_ID
    cols.time().put(i,9999999999.0+i*0.1) ;// TIME
    cols.interval().put(i,intervalvalue) ;// INTERVAL
    cols.name().put(i,"Antennae") ;// NAME
    cols.numPoly().put(i,0);// NUM_POLY
    cols.timeOrigin().put(i,4.8e+09);// TIME_ORIGIN
    cols.direction().put(i,mDirTarg);// DIRECTION
    cols.target().put(i,mDirTarg);// TARGET
    cols.pointingOffset().put(i,mDirTarg);//POINTING_OFFSET
    cols.sourceOffset().put(i,mDirTarg); //SOURCE_OFFSET
    cols.encoder().put(i,vEnco); // ENCODER
    cols.pointingModelId().put(i,0);  // POINTING_MODEL_ID
    cols.tracking().put(i,1); // TRACKING
    cols.onSource().put(i,1); // ON_SOURCE
    cols.overTheTop().put(i,1); // OVER_THE_TOP
    }

  // Finalize
  mspo.closeSubTables() ;
  mspo.unlock() ;

  double end = currenttime();
  cout << "Inserted: " << end - start << "sec\n";
}

void insertToPointing(char const *nrow) {
  int nrow_ = atoi(nrow);
  insertToPointing_(nrow_);
}

struct Entry {
  char const *option;
  void (*func)(char const *nrow);
} entries[] = {
  {"pointing", insertToPointing}
};

char const *progName = "";

void usage() {
  cerr << "Usage: " << progName << " [options] action num_of_row\n";
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
  void (*func)(char const *nrow) = NULL;
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

