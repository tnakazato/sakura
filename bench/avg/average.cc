#include <cstdio>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <string>
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

namespace {

double currenttime() {
  struct timeval tv;
  int result = gettimeofday(&tv, NULL);
  return tv.tv_sec + ((double)tv.tv_usec) / 1000000.;
}

void timing(string msg, double t) {
  cout << "Timing: " << msg << ": " << t << " sec\n";
}

struct AvgRec {
  size_t Nc;
  size_t Nf;
  int arrayId;
  int dataDescId;
  int antenna1;
  int antenna2;
  int stateId;
  int fieldId;
  int feed1;
  int feed2;
  double time;
  double interval;
  float *data;
  size_t rows;
  AvgRec() : data(NULL) {
  }
};

template <typename T>
void fillZero(T *data, size_t n) {
  for (size_t i = 0; i < n; i++) {
    data[i] = 0;
  }
}

void saveResult(Connection *con, vector<AvgRec> avgs) {
  enter();

  char const *cols[] = {
    "ARRAY_ID", // INTEGER NOT NULL,
    "DATA_DESC_ID", // INTEGER NOT NULL,
    "ANTENNA1", // INTEGER NOT NULL,
    "ANTENNA2", // INTEGER NOT NULL,
    "STATE_ID", // INTEGER NOT NULL,
    "FIELD_ID", // INTEGER NOT NULL,
    "FEED1", // INTEGER NOT NULL,
    "FEED2", // INTEGER NOT NULL,
    "TIME", // REAL NOT NULL,
    "INTERVAL", // REAL NOT NULL,
    "DATA", // BLOB DEFAULT NULL,
    NULL
  };
  enum {
    POS_ARRAY_ID = 1,
    POS_DATA_DESC_ID,
    POS_ANTENNA1,
    POS_ANTENNA2,
    POS_STATE_ID,
    POS_FIELD_ID,
    POS_FEED1,
    POS_FEED2,
    POS_TIME,
    POS_INTERVAL,
    POS_DATA
  };
  char const *zeroCols[] = {
    "PROCESSOR_ID", // INTEGER NOT NULL,
    "EXPOSURE", // REAL NOT NULL,
    "TIME_CENTROID", // REAL NOT NULL,
    "SCAN_NUMBER", // INTEGER NOT NULL,
    "OBSERVATION_ID", // INTEGER NOT NULL,
    "U", // REAL NOT NULL,
    "V", // REAL NOT NULL,
    "W", // REAL NOT NULL,
    "FLAG_ROW", // INTEGER NOT NULL,
    NULL
  };
  char const *emptyCols[] = {
    "SIGMA", // BLOB NOT NULL, 1
    "WEIGHT", // BLOB NOT NULL, 1
    "FLAG", // BLOB NOT NULL, false
    "FLAG_CATEGORY", // BLOB NOT NULL,
    NULL
  };
  string sql = "insert into MAIN (";
  sql += SQL::join(cols) + ", ";
  sql += SQL::join(zeroCols) + ", ";
  sql += SQL::join(emptyCols) + ") values (";
  sql += SQL::bindChars(cols) + ", ";
  sql += SQL::repeat(zeroCols, "0") + ", ";
  sql += SQL::repeat(emptyCols, "''") + ")";

  con->execute("BEGIN");
  unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));
  for (vector<AvgRec>::iterator i = avgs.begin(), end = avgs.end();
       i != end; ++i) {
    stmt->setInt(POS_ARRAY_ID, i->arrayId);
    stmt->setInt(POS_DATA_DESC_ID, i->dataDescId);
    stmt->setInt(POS_ANTENNA1, i->antenna1);
    stmt->setInt(POS_ANTENNA2, i->antenna2);
    stmt->setInt(POS_STATE_ID, i->stateId);
    stmt->setInt(POS_FIELD_ID, i->fieldId);
    stmt->setInt(POS_FEED1, i->feed1);
    stmt->setInt(POS_FEED2, i->feed2);
    stmt->setDouble(POS_TIME, i->time);
    stmt->setDouble(POS_INTERVAL, i->interval);
    stmt->setDouble(POS_INTERVAL, i->interval);

    size_t numData = i->Nc * i->Nf;
    stmt->setTransientBlob(POS_DATA, i->data,
			   2 * sizeof(*i->data) * numData);
    int result = stmt->executeUpdate();

    double sum = 0;
    for (size_t j = 0; j < 2 * numData; j++) {
      sum += i->data[j];
    }
    cout << "   ARRAY_ID=" << i->arrayId << endl ;
    cout << "   DATA_DESC_ID=" << i->dataDescId << endl ;
    cout << "   FIELD_ID=" << i->fieldId << endl ;
    cout << "   STATE_ID=" << i->stateId << endl ;
    cout << "   ANTENNA1=" << i->antenna1 << endl ;
    cout << "   ANTENNA2=" << i->antenna2 << endl ;
    cout << "   FEED1=" << i->feed1 << endl ;
    cout << "   FEED2=" << i->feed2 << endl ;
    cout << "Processed " << i->rows << " rows\n";
    cout << "Sum:  " << setprecision(8) << sum << endl ;
  }
  con->execute("COMMIT");
}

struct Group {
  int data_desc_id,
    array_id,
    antenna1,
    antenna2,
    field_id,
    state_id;
  bool operator ==(Group const & other) const {
    return data_desc_id == other.data_desc_id
      && array_id == other.array_id
      && antenna1 == other.antenna1
      && antenna2 == other.antenna2
      && field_id == other.field_id
      && state_id == other.state_id;
  }
  bool operator !=(Group const & other) const {
    return ! (*this == other);
  }
};

struct DataDesc {
  int Nc, Nf;
};

struct AvgContxt {
  float *weightSum;
  size_t numData;
  size_t nrowChunk;
  double timeCen;
  double interval;

  size_t totalRow;
  size_t avgIndex; // index for avgs[]
  bool extra;
  int lastGid;
  int gid;

  bool isFirst;
  Group grp;
  vector<AvgRec> avgs;
  map<int,DataDesc> dataDesc;
};

void addDataDesc(Connection *con, map<int,DataDesc> &dataDesc, int data_desc_id) {
  char const sql[] = "select P.NUM_CORR, SW.NUM_CHAN "
    "from data_description as D, polarization as P, SPECTRAL_WINDOW as SW "
    "where D.data_desc_id = ? and D.polarization_id = P.polarization_id and "
    "D.spectral_window_id = SW.spectral_window_id "
    "limit 1";
  unique_ptr<PreparedStatement> stmt(con->prepare(sql));
  stmt->setInt(1, data_desc_id);
  unique_ptr<ResultSet> rs(stmt->executeQuery());
  if (rs->next()) {
    DataDesc dd;
    dd.Nc = rs->getInt(1);
    dd.Nf = rs->getInt(2);
    dataDesc[data_desc_id] = dd;
  } else {
    cout << data_desc_id << endl;
    throw "Something wrong1";
  }
}

void avg(Connection *con, AvgContxt &ctx, ResultSet *rs) {
  assert(con != NULL);
  enum {
    DATA_DESC_ID = 1,
    ARRAY_ID,
    ANTENNA1,
    ANTENNA2,
    FIELD_ID,
    STATE_ID,
    FEED1,
    FEED2,
    TIME,
    INTERVAL,
    FLAG,
    DATA,
    FLOAT_DATA,
    WEIGHT
  };
  if (ctx.extra) {
    ctx.gid = ctx.lastGid + 1;
  } else {
    Group grp;
    grp.data_desc_id = rs->getInt(DATA_DESC_ID);
    grp.array_id = rs->getInt(ARRAY_ID);
    grp.antenna1 = rs->getInt(ANTENNA1);
    grp.antenna2 = rs->getInt(ANTENNA2);
    grp.field_id = rs->getInt(FIELD_ID);
    grp.state_id = rs->getInt(STATE_ID);
    if (ctx.isFirst) {
      ctx.isFirst = ! ctx.isFirst;
      addDataDesc(con, ctx.dataDesc, grp.data_desc_id);
      ctx.gid = ctx.lastGid + 1;
      ctx.grp = grp;
      ctx.avgs.push_back(AvgRec());
    } else if (ctx.grp != grp) {
      if (ctx.grp.data_desc_id != grp.data_desc_id) {
	addDataDesc(con, ctx.dataDesc, grp.data_desc_id);
      }
      ctx.gid = ctx.lastGid + 1;
      ctx.grp = grp;
      ctx.avgs.push_back(AvgRec());
    }
  }
  if (ctx.lastGid != ctx.gid) {
    if (ctx.lastGid >= 0) {
      // post process
      // final scaling of data
      size_t p = 0;
      for (size_t iFreq = 0; iFreq < ctx.avgs[ctx.avgIndex].Nf; iFreq++) {
	for (size_t ipol = 0; ipol < ctx.avgs[ctx.avgIndex].Nc; ipol++) {
	  double w = ctx.weightSum[ipol];
	  if (w != 0) {
	    // real
	    ctx.avgs[ctx.avgIndex].data[p] /= w;
	    // image
	    ctx.avgs[ctx.avgIndex].data[p + ctx.numData] /= w;
	  }
	  p++;
	}
      }
      ctx.timeCen /= ctx.nrowChunk;
      ctx.interval /= ctx.nrowChunk;
      
      ctx.avgs[ctx.avgIndex].time = ctx.timeCen;
      ctx.avgs[ctx.avgIndex].interval = ctx.interval;
      ctx.avgs[ctx.avgIndex].rows = ctx.nrowChunk;
      ctx.avgIndex++;
      ctx.totalRow += ctx.nrowChunk ;
      delete[] ctx.weightSum;
      ctx.weightSum = 0;
      if (ctx.extra) {
	return;
	//break;
      }
    }
    // pre process
    ctx.avgs[ctx.avgIndex].dataDescId = ctx.grp.data_desc_id;
    DataDesc dd = ctx.dataDesc[ctx.grp.data_desc_id];
    ctx.avgs[ctx.avgIndex].Nc = dd.Nc;
    ctx.avgs[ctx.avgIndex].Nf = dd.Nf;
    ctx.avgs[ctx.avgIndex].arrayId = ctx.grp.array_id;
    ctx.avgs[ctx.avgIndex].antenna1 = ctx.grp.antenna1;
    ctx.avgs[ctx.avgIndex].antenna2 = ctx.grp.antenna2;
    ctx.avgs[ctx.avgIndex].stateId = ctx.grp.state_id;
    ctx.avgs[ctx.avgIndex].fieldId = ctx.grp.field_id;
    ctx.avgs[ctx.avgIndex].feed1 = rs->getInt(FEED1);
    ctx.avgs[ctx.avgIndex].feed2 = rs->getInt(FEED2);

    ctx.numData = ctx.avgs[ctx.avgIndex].Nc * ctx.avgs[ctx.avgIndex].Nf;
    //cout << "numData set: " << ctx.numData << endl;
    ctx.avgs[ctx.avgIndex].data = new float[2 * ctx.numData];
    fillZero<float>(ctx.avgs[ctx.avgIndex].data, 2 * ctx.numData);
    ctx.weightSum = new float[ctx.avgs[ctx.avgIndex].Nc];
    fillZero<float>(ctx.weightSum, ctx.avgs[ctx.avgIndex].Nc);
    ctx.nrowChunk = 0;
    ctx.timeCen = 0.0;
    ctx.interval = 0.0;
  }
  // row process
  //int mainId = rs->getInt(MAIN_ID);
  //cout << "MAIN_ID: " << mainId << endl;

  bool releaseData = false;
  float *data = NULL;
  if (rs->isNull(DATA)) {
    // float -> complex
    int floatSize;
    float const *floatData =
      (float const *)rs->getTransientBlob(FLOAT_DATA, &floatSize);
    //cout << floatSize << ", " << ctx.numData << endl;
    assert(floatSize == sizeof(*floatData) * ctx.numData);
    data = new float[2 * ctx.numData];
    releaseData = true;
    float *const real = data;
    float *const imag = &data[ctx.numData];
    for (size_t i = 0; i < ctx.numData; i++) {
      real[i] = floatData[i];
      imag[i] = 0;
    }
  } else {
    int dataSize;
    data = (float *)rs->getTransientBlob(DATA, &dataSize);
    assert(dataSize == sizeof(*data) * 2 * ctx.numData);
  }

  int flagSize;
  bool *flag = (bool *)rs->getTransientBlob(FLAG, &flagSize);
  assert(flagSize == sizeof(*flag) * ctx.numData);

  int weightSize;
  float *weight = (float *)rs->getTransientBlob(WEIGHT, &weightSize);
  assert(weightSize == sizeof(*weight) * ctx.avgs[ctx.avgIndex].Nc);

  double time = rs->getDouble(TIME);
  double tint = rs->getDouble(INTERVAL);
  ctx.timeCen += time;
  ctx.interval += tint;
      
  size_t p = 0;
  for (size_t iFreq = 0; iFreq < ctx.avgs[ctx.avgIndex].Nf; iFreq++) {
    for (size_t ipol = 0; ipol < ctx.avgs[ctx.avgIndex].Nc; ipol++) {
      if (!flag[p]) {
	// real
	ctx.avgs[ctx.avgIndex].data[p] += data[p] * weight[ipol] * tint;
	// image
	ctx.avgs[ctx.avgIndex].data[p+ctx.numData] += data[p+ctx.numData] * weight[ipol] * tint;
      }
      p++;
    }
  }
  if (releaseData) {
    delete[] data;
    data = NULL;
  }
  for (size_t ipol = 0; ipol < ctx.avgs[ctx.avgIndex].Nc; ipol++) {
    ctx.weightSum[ipol] += weight[ipol] * tint;
  }
  ctx.nrowChunk++;
  ctx.lastGid = ctx.gid;
}

void msaverage(Connection *con) {
  enter();

  vector<vector<int> *> mainIds;
  {
    static char const sql[] =
      "select MAIN_ID from mst.MAIN "
      "order by DATA_DESC_ID, ARRAY_ID, ANTENNA1, ANTENNA2, FIELD_ID, STATE_ID";
    unique_ptr<PreparedStatement> stmt(con->prepare(sql));
    unique_ptr<ResultSet> rs(stmt->executeQuery());
    vector<int> *ids = NULL;
    vector<int>::size_type count = 0;
    while (rs->next()) {
      if (ids == NULL) {
	ids = new vector<int>(1024*1024, -1);
	mainIds.push_back(ids);
	count = 0;
      }
      (*ids)[count] = rs->getInt(1);
      count++;
      if (count >= ids->size()) {
	ids = NULL;
      }
    }
    if (ids) {
      ids->resize(count);
    }
  }

  struct AvgContxt avgContext = {
    NULL,
    0,
    0,
    0,
    0,
      
    0,
    0, // index for avgs[]
    false,
    -1,
    0,

    true,
    { 0, 0, 0, 0, 0, 0 }
  };

  con->execute("BEGIN");
  try {
    static char const * const cols[] = {
      "DATA_DESC_ID",
      "ARRAY_ID",
      "ANTENNA1",
      "ANTENNA2",
      "FIELD_ID",
      "STATE_ID",
      "FEED1",
      "FEED2",
      "TIME",
      "INTERVAL",
      "FLAG",
      "DATA",
      "FLOAT_DATA",
      "WEIGHT",
      NULL
    };

    string sql = "select ";
    sql += SQL::join(cols, ", ") + " from mst.MAIN where MAIN_ID = ? and FLAG_ROW = 0";
    unique_ptr<PreparedStatement> stmt(con->prepare(sql.c_str()));

    for (vector<vector<int> *>::iterator outer = mainIds.begin(),
	   outerEnd = mainIds.end();
	 outer != outerEnd; ++outer) {
      vector<int> *ids = *outer;
      for (vector<int>::iterator inner = ids->begin(),
	     innerEnd = ids->end();
	   inner != innerEnd; ++inner) {
	int mainId = *inner;
	stmt->setInt(1, mainId);
	unique_ptr<ResultSet> rs(stmt->executeQuery());
	if (rs->next()) {
	  avg(con, avgContext, rs.get());
	} else {
	  assert(false);
	  throw "Something wrong2";
	}
      }
      delete ids;
    }
    avgContext.extra = true;
    avg(con, avgContext, NULL);
    cout << "total number of rows processed: " << avgContext.totalRow << endl;
  } catch (...) {
    vector<AvgRec> &avgs =  avgContext.avgs;
    for (vector<AvgRec>::iterator i = avgs.begin(), end = avgs.end();
	 i != end; ++i) {
      delete[] i->data;
    }
    con->execute("ROLLBACK");
    throw;
  }
  con->execute("COMMIT");

  vector<AvgRec> &avgs =  avgContext.avgs;
  try {
    double start = currenttime();
    saveResult(con, avgs);
    double end = currenttime();
    timing("Saving", end - start);
  } catch (...) {
    for (vector<AvgRec>::iterator i = avgs.begin(), end = avgs.end();
	 i != end; ++i) {
      delete[] i->data;
    }
    throw;
  }
  for (vector<AvgRec>::iterator i = avgs.begin(), end = avgs.end();
       i != end; ++i) {
    delete[] i->data;
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

void createTxTables(Connection *con, char const *prefix) {
  char const sql[] = "/sql/";
  // create empty transaction db from mst
  {
    char const ddl_file[] = "MST-nodrop.ddl";
    size_t ddl_path_size = strlen(prefix) + strlen(sql) + strlen(ddl_file) + 1;
    char ddl_path[ddl_path_size];
    snprintf(ddl_path, ddl_path_size, "%s%s%s", prefix, sql, ddl_file);
    char *mst_ddl = readFileContent(ddl_path);
    try {
      con->execute(mst_ddl);
    } catch (...) {
      delete[] mst_ddl;
      throw;
    }
    delete[] mst_ddl;
  }
}

void copy(Connection * con) {
  enter();
  createTxTables(con, "/home/kohji/sakura/dist");

  char const sql[] =
    "insert into SYSCAL select * from mst.SYSCAL;\n"
    "insert into FLAG_CMD select * from mst.FLAG_CMD;\n"
    "insert into HISTORY select * from mst.HISTORY;\n"
    "insert into HISTORY_PARAM select * from mst.HISTORY_PARAM;\n"
    "insert into HISTORY_COMMAND select * from mst.HISTORY_COMMAND;\n";
  con->execute(sql);
  cout << "Copied\n";
}

void average(char const *srcm, char const *srct, char const *dstt) {
  enter();
  string dstTx = "file:";
  dstTx += dstt;
  //dstTx += "?mode=rw";

  cout << "Opening: " << dstTx << endl;
  unique_ptr<Connection> con(Connection::open(dstTx.c_str()));

  {
    string srcMaster = "file:";
    srcMaster += srcm;
    srcMaster += "?mode=ro";

    cout << "Attaching: " << srcMaster << endl;
    unique_ptr<PreparedStatement> stmt(con->prepare("ATTACH DATABASE ? AS msm;"));
    stmt->setTransientString(1, srcMaster.c_str());
    stmt->executeUpdate();
  }

  {
    string srcTx = "file:";
    srcTx += srct;
    srcTx += "?mode=ro";

    cout << "Attaching: " << srcTx << endl;
    unique_ptr<PreparedStatement> stmt(con->prepare("ATTACH DATABASE ? AS mst;"));
    stmt->setTransientString(1, srcTx.c_str());
    stmt->executeUpdate();
  }

  {
    double start = currenttime();
    //con->execute("PRAGMA temp_store = memory;");
    copy(con.get());
    double end = currenttime();
    timing("Copied", end - start);
  }
  msaverage(con.get());
}

}

int main(int argc, char const *argv[]) { 
  char const *const progName = argv[0];

  if (argc != 4) {
    cerr << "Usage: " << progName << " src.mdb src.tdb dst.tdb\n";
    return 1;
  }
  double start = currenttime();
  try {
    average(argv[1], argv[2], argv[3]);
  } catch (char const *ex) {
    cout << ex << endl;
  }
  double end = currenttime();
  timing("Total", end - start);
  return 0;
}
