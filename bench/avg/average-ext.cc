#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <sys/time.h>
#include <sqlite3ext.h>

SQLITE_EXTENSION_INIT1

using namespace std;

namespace {

template <typename T>
void fillZero(T data[], size_t n) {
  for (size_t i = 0; i < n; i++) {
    data[i] = 0;
  }
}

template <typename T>
void deleteArray(void *data) {
  delete[] (T*)data;
}

extern "C" {

// avg_weight(NC, NF, DATA, FLOAT_DATA, WEIGHT, INTERVAL)
static void avgWeightFunc(sqlite3_context *ctx, int argc, sqlite3_value *argv[]) {
  assert(argc >= 6);
  if (sqlite3_value_type(argv[0]) == SQLITE_NULL
      || sqlite3_value_type(argv[1]) == SQLITE_NULL
      || (sqlite3_value_type(argv[2]) == SQLITE_NULL
	  && sqlite3_value_type(argv[3]) == SQLITE_NULL)
      || sqlite3_value_type(argv[4]) == SQLITE_NULL
      || sqlite3_value_type(argv[5]) == SQLITE_NULL) {
    sqlite3_result_null(ctx);
    return;
  }
  // TODO check type
  int nc = sqlite3_value_int(argv[0]);
  int nf = sqlite3_value_int(argv[1]);
  if (nc < 0 || nf < 0) {
    sqlite3_result_error(ctx, "Invalid NC or NF or both.", -1);
    return;
  }
  size_t const numData = nc * nf;

  // TODO check type
  double tint = sqlite3_value_double(argv[5]);

  float const *weight = NULL;
  {
    float *data = (float *)sqlite3_value_blob(argv[4]);
    int sizeData = sqlite3_value_bytes(argv[4]);
    if (sizeData != sizeof(*data) *  nc) {
      sqlite3_result_error(ctx, "Invalid size of WEIGHT.", -1);
      return;
    }
    weight = data;
  }

  float const *real = NULL;
  float const *imag = NULL;
  bool releaseImag = false;

  if (sqlite3_value_type(argv[2]) != SQLITE_NULL) {
    float const *data = (float *)sqlite3_value_blob(argv[2]);
    int sizeData = sqlite3_value_bytes(argv[2]);
    if (sizeData != sizeof(*data) * 2 * numData) {
      sqlite3_result_error(ctx, "Invalid size of DATA.", -1);
      return;
    }
    real = data;
    imag = &data[numData];
  } else {
    float const *data = (float *)sqlite3_value_blob(argv[3]);
    int sizeData = sqlite3_value_bytes(argv[3]);
    if (sizeData != sizeof(*data) *  numData) {
      sqlite3_result_error(ctx, "Invalid size of FLOAT_DATA.", -1);
      return;
    }
    real = data;
    float *tmp = new float[numData];
    imag = tmp;
    fillZero(tmp, numData);
    releaseImag = true;
  }

  float *result = new float[2 * numData];

  size_t p = 0;
  for (size_t iFreq = 0; iFreq < nf; iFreq++) {
    for (size_t ipol = 0; ipol < nc; ipol++) {
      double w = weight[ipol] * tint;
      result[p] = real[p] * w;
      result[p+numData] = imag[p] * w;
      p++;
    }
  }
  if (releaseImag) {
    delete[] imag;
  }

  sqlite3_result_blob(ctx, result, sizeof(*result) * 2 * numData, deleteArray<float>);
}

// float_mul_scalar(DATA, NUM)
static void float_mul_scalarFunc(sqlite3_context *ctx, int argc, sqlite3_value *argv[]) {
  assert(argc == 2);
  if (sqlite3_value_type(argv[0]) == SQLITE_NULL
      || sqlite3_value_type(argv[1]) == SQLITE_NULL) {
    sqlite3_result_null(ctx);
    return;
  }
  // TODO check type
  double num = sqlite3_value_double(argv[1]);

  float const *data = (float *)sqlite3_value_blob(argv[0]);
  int sizeData = sqlite3_value_bytes(argv[0]);
  if (sizeData % sizeof(*data) != 0) {
    sqlite3_result_error(ctx, "Invalid size of DATA.", -1);
    return;
  }
  size_t numData = sizeData / sizeof(*data);

  float *result = new float[numData];

  for (size_t i = 0; i < numData; i++) {
    result[i] = data[i] * num;
  }

  sqlite3_result_blob(ctx, result, sizeof(*result) * numData, deleteArray<float>);
}

// float_sumUp(DATA)
static void float_sumUpFunc(sqlite3_context *ctx, int argc, sqlite3_value *argv[]) {
  assert(argc == 1);
  if (sqlite3_value_type(argv[0]) == SQLITE_NULL) {
    sqlite3_result_null(ctx);
    return;
  }

  float const *data = (float *)sqlite3_value_blob(argv[0]);
  int sizeData = sqlite3_value_bytes(argv[0]);
  if (sizeData % sizeof(*data) != 0) {
    sqlite3_result_error(ctx, "Invalid size of DATA.", -1);
    return;
  }
  size_t numData = sizeData / sizeof(*data);

  float result = 0;

  for (size_t i = 0; i < numData; i++) {
    result += data[i];
  }

  sqlite3_result_double(ctx, result);
}

// complex_div_weight(NC, NF, DATA, WEIGHT)
static void complex_div_weightFunc(sqlite3_context *ctx, int argc, sqlite3_value *argv[]) {
  assert(argc == 4);
  if (sqlite3_value_type(argv[0]) == SQLITE_NULL
      || sqlite3_value_type(argv[1]) == SQLITE_NULL
      || sqlite3_value_type(argv[2]) == SQLITE_NULL
      || sqlite3_value_type(argv[3]) == SQLITE_NULL) {
    sqlite3_result_null(ctx);
    return;
  }
  // TODO check type
  int nc = sqlite3_value_int(argv[0]);
  int nf = sqlite3_value_int(argv[1]);
  if (nc < 0 || nf < 0) {
    sqlite3_result_error(ctx, "Invalid NC or NF or both.", -1);
    return;
  }
  size_t const numData = nc * nf;

  float const *weight = NULL;
  {
    float *data = (float *)sqlite3_value_blob(argv[3]);
    int sizeData = sqlite3_value_bytes(argv[3]);
    if (sizeData != sizeof(*data) *  nc) {
      sqlite3_result_error(ctx, "Invalid size of WEIGHT.", -1);
      return;
    }
    weight = data;
  }

  float const *real = NULL;
  float const *imag = NULL;

  float const *data = (float *)sqlite3_value_blob(argv[2]);
  {
    int sizeData = sqlite3_value_bytes(argv[2]);
    if (sizeData != sizeof(*data) * 2 * numData) {
      sqlite3_result_error(ctx, "Invalid size of DATA.", -1);
      return;
    }
  }
  real = data;
  imag = &data[numData];

  float *result = new float[2 * numData];

  size_t p = 0;
  for (size_t iFreq = 0; iFreq < nf; iFreq++) {
    for (size_t ipol = 0; ipol < nc; ipol++) {
      if (weight[ipol] != 0) {
	result[p] = real[p] / weight[ipol];
	result[p+numData] = imag[p] / weight[ipol];
      }
      p++;
    }
  }

  sqlite3_result_blob(ctx, result, sizeof(*result) * 2 * numData, deleteArray<float>);
}

struct FloatSumCtx {
  size_t elements;
  size_t count;
  float *sum;
};

static void complexSumStep(sqlite3_context *ctx, int argc, sqlite3_value *argv[]) {
  assert(argc >= 1);
  FloatSumCtx *agCtx =
    (FloatSumCtx *)sqlite3_aggregate_context(ctx, sizeof(FloatSumCtx));
  if (agCtx == NULL) {
    sqlite3_result_error_nomem(ctx);
    return;
  }
  if (sqlite3_value_type(argv[0]) == SQLITE_NULL) {
    return;
  }
  if (argc >= 2 && sqlite3_value_type(argv[1]) == SQLITE_NULL) {
    return;
  }
  float const *a = (float const *)sqlite3_value_blob(argv[0]);
  int sizeA = sqlite3_value_bytes(argv[0]);
  if (sizeA < 0) {
    sqlite3_result_error(ctx, "Size of 1st arg should be >= 0.", -1);
    return;
  }
  if (agCtx->sum == NULL) {
    if (sizeA % (2 * sizeof(*agCtx->sum)) != 0) {
      sqlite3_result_error(ctx, "Wrong size of 1st arg", -1);
      return;
    }
    int elements = sizeA / (2 * sizeof(*agCtx->sum));
    agCtx->elements = elements;
    agCtx->sum = new float[2 * agCtx->elements];
    fillZero(agCtx->sum, 2 * agCtx->elements);
  }
  if (sizeA != 2 * sizeof(*agCtx->sum) * agCtx->elements) {
    sqlite3_result_error(ctx, "Inconsistent BLOB size.", -1);
    return;
  }
  bool masterFlag = true; // true means 'Do operation'
  bool const *flags = NULL;
  if (argc >= 2) {
    flags = (bool const *)sqlite3_value_blob(argv[1]);
    int sizeFlag = sqlite3_value_bytes(argv[1]);
    if (sizeFlag != sizeof(*flags) * agCtx->elements) {
      sqlite3_result_error(ctx, "Size of 2nd arg is inconsistent.", -1);
      return;
    }
    masterFlag = false;
  }
  for (size_t i = 0; i < agCtx->elements; i++) {
    if (masterFlag || !flags[i]) {
      agCtx->sum[i] += a[i];
      agCtx->sum[i + agCtx->elements] += a[i + agCtx->elements];
    }
  }
  agCtx->count++;
}

static void complexSumFinal(sqlite3_context *ctx) {
  FloatSumCtx *agCtx =
    (FloatSumCtx *)sqlite3_aggregate_context(ctx, sizeof(FloatSumCtx));
  if (agCtx == NULL) {
    sqlite3_result_error_nomem(ctx);
    return;
  }
  //cout << __FUNCTION__ << " called " << agCtx << " count: " << agCtx->count << "\n";
  if (agCtx->sum == NULL) {
    // step func has never been called.
    sqlite3_result_null(ctx);
  }
  sqlite3_result_blob(ctx, agCtx->sum, 2 * sizeof(*agCtx->sum) * agCtx->elements, deleteArray<float>);
}

static void floatSumStep(sqlite3_context *ctx, int argc, sqlite3_value *argv[]) {
  assert(argc >= 1);
  FloatSumCtx *agCtx =
    (FloatSumCtx *)sqlite3_aggregate_context(ctx, sizeof(FloatSumCtx));
  if (agCtx == NULL) {
    sqlite3_result_error_nomem(ctx);
    return;
  }
  if (sqlite3_value_type(argv[0]) == SQLITE_NULL) {
    return;
  }
  if (argc >= 2 && sqlite3_value_type(argv[1]) == SQLITE_NULL) {
    return;
  }
  float const *a = (float const *)sqlite3_value_blob(argv[0]);
  int sizeA = sqlite3_value_bytes(argv[0]);
  if (sizeA < 0) {
    sqlite3_result_error(ctx, "Size of 1st arg should be >= 0.", -1);
    return;
  }
  if (agCtx->sum == NULL) {
    if (sizeA % sizeof(*agCtx->sum) != 0) {
      sqlite3_result_error(ctx, "Wrong size of 1st arg", -1);
      return;
    }
    int elements = sizeA / sizeof(*agCtx->sum);
    agCtx->elements = elements;
    agCtx->sum = new float[agCtx->elements];
    fillZero(agCtx->sum, agCtx->elements);
  }
  if (sizeA != sizeof(*agCtx->sum) * agCtx->elements) {
    sqlite3_result_error(ctx, "Inconsistent BLOB size.", -1);
    return;
  }
  bool masterFlag = true; // true means 'Do operation'
  bool const *flags = NULL;
  if (argc >= 2) {
    flags = (bool const *)sqlite3_value_blob(argv[1]);
    int sizeFlag = sqlite3_value_bytes(argv[1]);
    if (sizeFlag != sizeof(*flags) * agCtx->elements) {
      sqlite3_result_error(ctx, "Size of 2nd arg is inconsistent.", -1);
      return;
    }
    masterFlag = false;
  }
  for (size_t i = 0; i < agCtx->elements; i++) {
    if (masterFlag || !flags[i]) {
      agCtx->sum[i] += a[i];
    }
  }
  agCtx->count++;
}

static void floatSumFinal(sqlite3_context *ctx) {
  FloatSumCtx *agCtx =
    (FloatSumCtx *)sqlite3_aggregate_context(ctx, sizeof(FloatSumCtx));
  if (agCtx == NULL) {
    sqlite3_result_error_nomem(ctx);
    return;
  }
  //cout << __FUNCTION__ << " called " << agCtx << " count: " << agCtx->count << "\n";
  if (agCtx->sum == NULL) {
    // step func has never been called.
    sqlite3_result_null(ctx);
  }
  sqlite3_result_blob(ctx, agCtx->sum, sizeof(*agCtx->sum) * agCtx->elements, deleteArray<float>);
}

};

}

extern "C" {
int
sqlite3_extension_init(sqlite3 *db,
		       char **pzErrMsg,
		       const sqlite3_api_routines *pApi) {
  SQLITE_EXTENSION_INIT2(pApi)
    ;
  int result = sqlite3_create_function_v2(db, "avg_weight", 6,
		      SQLITE_ANY,
		      NULL,
		      avgWeightFunc,
		      NULL,
		      NULL,
		      NULL);
  if (result != SQLITE_OK) {
    *pzErrMsg = sqlite3_mprintf("Can't create function(s).");
    return 1;
  }

  result = sqlite3_create_function_v2(db, "float_sumUp", 1,
				      SQLITE_ANY, NULL,
				      float_sumUpFunc, NULL, NULL, NULL);
  if (result != SQLITE_OK) {
    *pzErrMsg = sqlite3_mprintf("Can't create function(s).");
    return 1;
  }

  result = sqlite3_create_function_v2(db, "complex_div_weight", 4,
				      SQLITE_ANY, NULL,
				      complex_div_weightFunc, NULL, NULL, NULL);
  if (result != SQLITE_OK) {
    *pzErrMsg = sqlite3_mprintf("Can't create function(s).");
    return 1;
  }

  result = sqlite3_create_function_v2(db, "float_mul_scalar", 2,
		      SQLITE_ANY, NULL,
		      float_mul_scalarFunc, NULL, NULL, NULL);
  if (result != SQLITE_OK) {
    *pzErrMsg = sqlite3_mprintf("Can't create function(s).");
    return 1;
  }

  // FLOAT_SUM(BLOB, FLAG-opt)
  result = sqlite3_create_function_v2(db, "complex_sum", 1,
		      SQLITE_ANY, NULL,
		      NULL, complexSumStep, complexSumFinal, NULL);
  if (result != SQLITE_OK) {
    *pzErrMsg = sqlite3_mprintf("Can't create function(s).");
    return 1;
  }

  result = sqlite3_create_function_v2(db, "complex_sum", 2,
		      SQLITE_ANY, NULL,
		      NULL, complexSumStep, complexSumFinal, NULL);
  if (result != SQLITE_OK) {
    *pzErrMsg = sqlite3_mprintf("Can't create function(s).");
    return 1;
  }

  result = sqlite3_create_function_v2(db, "float_sum", 1,
		      SQLITE_ANY, NULL,
		      NULL, floatSumStep, floatSumFinal, NULL);
  if (result != SQLITE_OK) {
    *pzErrMsg = sqlite3_mprintf("Can't create function(s).");
    return 1;
  }

  result = sqlite3_create_function_v2(db, "float_sum", 2,
		      SQLITE_ANY, NULL,
		      NULL, floatSumStep, floatSumFinal, NULL);
  if (result != SQLITE_OK) {
    *pzErrMsg = sqlite3_mprintf("Can't create function(s).");
    return 1;
  }
  return 0;
}

};
