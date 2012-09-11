#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <sys/types.h>
#include <sys/time.h>
#include <sqlite3ext.h>

SQLITE_EXTENSION_INIT1

using namespace std;

namespace {

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

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

template<typename T>
int checkSize(sqlite3_context *ctx, void const *data, int dataSize,
	       int dims, int const elements[]) throw (char const *) {
  if (dataSize % sizeof(T) != 0) {
    sqlite3_result_error(ctx, "Invalid size of DATA.", -1);
    throw "error";
  }
  int dataElements = dataSize / sizeof(T);
  int n = elements[0];
  if (n == -1) {
    n = 0;
  }
  int n2 = 1;
  for (int i = 1; i < dims; i++) {
    n *= elements[i];
    n2 *= elements[i];
  }
  if (dataElements < n) {
    sqlite3_result_error(ctx, "Too small size of DATA.", -1);
    throw "error";
  }
  if (dataElements % n2 != 0) {
    sqlite3_result_error(ctx, "Invalid size of DATA.", -1);
    throw "error";
  }
  if (elements[0] >= 0) {
    return elements[0];
  }
  return dataElements / n2;
}

int estimateSize(int formattedLen,
		 int dims, int const elements[],
		 int pairLen, int sepLen) {
  if (elements[0] == 0) {
    return pairLen;
  }
  int len = formattedLen;
  for (int i = dims - 1; i >= 0; i--) {
    assert(elements[i] > 0);
    len = pairLen + len * elements[i] + sepLen * (elements[i] - 1);
  }
  return len;
}

template<typename T, int formattedLen, typename DispType>
void dumpData(char *&buf, char *const bufEnd,
	      T const* data, int &dataIdx,
	      char const format[],
	      DispType (*conv)(T value),
	      int dims, int elements[],
	      char const opener[], char const closer[], char const separator[],
	      int openerLen,
	      int closerLen,
	      int separatorLen) {
  assert(dims > 0);
  *buf = '\0';
  strncat(buf, opener, bufEnd - buf);
  buf += openerLen;
  assert(buf <= bufEnd);
  char const *sep = "";
  int sepLen = 0;
  for (int i = 0; i < elements[0]; i++) {
    *buf = '\0';
    strncat(buf, sep, bufEnd - buf);
    buf += sepLen;

    if (dims == 1) {
      int chars = snprintf(buf, bufEnd - buf, format, conv(data[dataIdx++]));
      assert(chars <= formattedLen);
      if (chars <= bufEnd - buf) {
	buf += chars;
      } else {
	buf = bufEnd;
      }
    } else {
      dumpData<T, formattedLen, DispType>(buf, bufEnd,
		data, dataIdx,
		format, conv,
		dims - 1, &elements[1],
		opener, closer, separator,
		openerLen, closerLen, separatorLen);
    }

    sep = separator;
    sepLen = separatorLen;
  }
  *buf = '\0';
  strncat(buf, closer, bufEnd - buf);
  buf += closerLen;
  assert(buf <= bufEnd);
}

template<typename T, int formattedLen, typename DispType>
void setResult(sqlite3_context *ctx, void const *data, int dataSize,
	       char const format[],
	       void (*conv)(),
	       int dims, int elements[],
	       char const opener[], char const closer[], char const separator[]) {
  try {
    elements[0] = checkSize<T>(ctx, data, dataSize, dims, elements);

    int openerLen = strlen(opener);
    int closerLen = strlen(closer);
    int pairLen = openerLen + closerLen;
    int sepLen = strlen(separator);

    int bufSize = estimateSize(formattedLen, dims, elements, pairLen, sepLen);
    assert(bufSize > 0);
    T const *values = reinterpret_cast<T const *>(data);
    char * const buf = new char[bufSize + 1]; // +1 for '\0'
    try {
      char *ptr = buf;
      char * const bufEnd = &buf[bufSize];
      int dataIdx = 0;
      dumpData<T, formattedLen>(ptr, bufEnd,
		  values, dataIdx,
		  format, reinterpret_cast<DispType (*)(T)>(conv),
		  dims, elements,
		  opener, closer, separator,
		  openerLen, closerLen, sepLen);
      assert(ptr <= bufEnd);
      *ptr = '\0';
      sqlite3_result_text(ctx, buf, ptr - buf, deleteArray<char>);
    } catch (...) {
      delete[] buf;
      throw;
    }
  } catch (...) {
    return;
  }
}

template<typename T>
T getParam(T (*func)(sqlite3_value *), sqlite3_value *value, T valueIfNull) {
  if (sqlite3_value_type(value) == SQLITE_NULL) {
    return valueIfNull;
  }
  return func(value);
}


typedef void (*DummyConvFunc_t)(); // dummy signature

struct DumpFunc {
  char const *format;
  void (*func)(sqlite3_context *ctx, void const *data, int dataSize,
	       char const format[], DummyConvFunc_t conv,
	       int dims, int elements[],
	       char const opener[], char const closer[], char const separator[]);
  DummyConvFunc_t conv;
};

template<typename T, typename DispType>
DispType convertForDisplay(T value) {
  return static_cast<DispType>(value);
}

template<>
char convertForDisplay<bool, char>(bool value) {
  return value ? 'T' : 'F';
}



extern "C" {

// dump_blob_as_sometype(data, opener, closer, separator, [nElements, ...])
static void dump_blob_as(sqlite3_context *ctx, int argc, sqlite3_value *argv[]) {
  assert(argc >= 4);
  if (sqlite3_value_type(argv[0]) == SQLITE_NULL) {
    sqlite3_result_null(ctx);
    return;
  }
  int const dims = MAX(argc - 4, 1);
  int elements[dims];
  elements[0] = -1; // this means actual data dependent.
  for (int i = 4; i < argc; i++) {
    if (sqlite3_value_type(argv[i]) != SQLITE_INTEGER) {
      sqlite3_result_error(ctx, "Invalid number of elements.", -1);
      return;
    }
    elements[i - 4] = sqlite3_value_int(argv[i]);
    if (! ((i == 4 && elements[i - 4] == -1)
	   || elements[i - 4] > 0)) {
      sqlite3_result_error(ctx, "Negative number of elements.", -1);
      return;
    }
  }

  // TODO check type
  void const *data = sqlite3_value_blob(argv[0]);
  int dataSize = sqlite3_value_bytes(argv[0]);

  unsigned char const defaultOpener[] = "{";
  unsigned char const defaultCloser[] = "}";
  unsigned char const defaultSeparator[] = ",";
  char const *opener =
    reinterpret_cast<char const *>(
	getParam<unsigned char const *>(sqlite3_value_text, argv[1],
					defaultOpener));
  char const *closer =
    reinterpret_cast<char const *>(
	getParam(sqlite3_value_text, argv[2], defaultCloser));
  char const *separator =
    reinterpret_cast<char const *>(
	getParam(sqlite3_value_text, argv[3], defaultSeparator));

  DumpFunc const *dumpFunc =
    reinterpret_cast<DumpFunc const *>(sqlite3_user_data(ctx));
  assert(dumpFunc != NULL);
  dumpFunc->func(ctx, data, dataSize,
		 dumpFunc->format, dumpFunc->conv,
		 dims, elements,
		 opener, closer, separator);
}

int registerDumpFuncs(sqlite3 *db,
		      char **pzErrMsg,
		      char const funcName[],
		      DumpFunc const *dumpFunc) {
  for (int i = 4; i < 10; i++) {
    int result =
      sqlite3_create_function_v2(db, funcName, i,
				 SQLITE_ANY,
				 const_cast<DumpFunc *>(dumpFunc),
				 dump_blob_as, NULL, NULL,
				 NULL);
    if (result != SQLITE_OK) {
      *pzErrMsg = sqlite3_mprintf("Can't create function(s).");
      return 1;
    }
  }
  return 0;
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

  {
    double (*func)(float) = convertForDisplay<float, double>; // to instantiate template
    static DumpFunc dumpFunc = {
      "%g", setResult<float, 13, double>,
      reinterpret_cast<DummyConvFunc_t>(func)
    };
    int result = registerDumpFuncs(db, pzErrMsg,
				   "dump_blob_as_float",
				   &dumpFunc);
    if (result != 0) {
      return result;
    }
  }

  {
    char (*func)(bool) = convertForDisplay<bool, char>; // to instantiate template
    static DumpFunc dumpFunc = {
      "%c", setResult<bool, 1, char>,
      reinterpret_cast<DummyConvFunc_t>(func)
    };
    int result = registerDumpFuncs(db, pzErrMsg,
				   "dump_blob_as_bool",
				   &dumpFunc);
    if (result != 0) {
      return result;
    }
  }

  {
    double (*func)(double) = convertForDisplay<double, double>; // to instantiate template
    static DumpFunc dumpFunc = {
      "%g", setResult<double, 13, double>,
      reinterpret_cast<DummyConvFunc_t>(func)
    };
    int result = registerDumpFuncs(db, pzErrMsg,
				   "dump_blob_as_double",
				   &dumpFunc);
    if (result != 0) {
      return result;
    }
  }

  {
    int (*func)(int32_t) = convertForDisplay<int32_t, int>; // to instantiate template
    static DumpFunc dumpFunc = {
      "%d", setResult<int32_t, 11, int>,
      reinterpret_cast<DummyConvFunc_t>(func)
    };
    int result = registerDumpFuncs(db, pzErrMsg,
				   "dump_blob_as_int32",
				   &dumpFunc);
    if (result != 0) {
      return result;
    }
  }

  return 0;
}

}
