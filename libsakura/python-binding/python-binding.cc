/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2016
 * National Astronomical Observatory of Japan
 * 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
 * 
 * This file is part of Sakura.
 * 
 * Sakura is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * Sakura is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
 * @SAKURA_LICENSE_HEADER_END@
 */

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <unistd.h>
#include <memory>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include <libsakura/memory_manager.h>
#include "libsakura/sakura-python.h"

/*
 * Refer to following documents:
 * http://docs.python.jp/2/extending/extending.html
 * http://docs.python.jp/2/c-api/capsule.html
 * http://docs.python.jp/2/c-api/arg.html
 *
 * Refer to following documents for Python 3 support:
 * https://docs.python.org/3.5/howto/cporting.html
 *
 */

#define MODULE_NAME "libsakurapy"

#undef KEEP_GIL

#ifdef KEEP_GIL
# define SAKURA_BEGIN_ALLOW_THREADS {
# define SAKURA_END_ALLOW_THREADS }
#else
# define SAKURA_BEGIN_ALLOW_THREADS Py_BEGIN_ALLOW_THREADS
# define SAKURA_END_ALLOW_THREADS Py_END_ALLOW_THREADS
#endif

namespace {

constexpr size_t kMaxNumberOfDimensions = 5;

typedef PyObject *(*FuncForPython)(PyObject *self, PyObject *args);

constexpr char const kAlignedBufferName[] =
#if 0
		"AlignedBuffer.sakura.nao.ac.jp";
#else
		MODULE_NAME ".AlignedBuffer";
#endif

constexpr char const kConvolve1DContextName[] =
#if 0
		"Convolve1DContext.sakura.nao.ac.jp";
#else
		MODULE_NAME ".Convolve1DContext";
#endif

constexpr char const kBaselineContextName[] =
#if 0
		"BaselineContext.sakura.nao.ac.jp";
#else
		MODULE_NAME ".BaselineContext";
#endif

constexpr char const kAlignedBufferForNPName[] =
#if 0
		"AlignedBufferForNP.sakura.nao.ac.jp";
#else
		MODULE_NAME ".AlignedBufferForNP";
#endif
}

extern "C" {

struct LIBSAKURA_SYMBOL(PyAlignedBuffer) {
	void *original_addr;
	void *aligned_addr;
	void (*destructor)(void *);
	size_t elements[kMaxNumberOfDimensions];
	unsigned dimensions;LIBSAKURA_SYMBOL(PyTypeId) type;
};

}

namespace {

struct module_state {
	PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

#define LONG_FROM_SSIZE_T(s) PyLong_FromSsize_t(s)
#define ASLONG(s) PyLong_AsLong(s)
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;

#define LONG_FROM_SSIZE_T(s) PyInt_FromSsize_t(s)
#define ASLONG(s) PyInt_AsLong(s)
#endif

//PyObject *sakura_error;

void DecrementRef(PyObject *obj) {
	Py_XDECREF(obj);
}

class RefHolder: public std::unique_ptr<PyObject, decltype(&DecrementRef)> {
public:
	explicit RefHolder(PyObject *obj) :
			std::unique_ptr<PyObject, decltype(&DecrementRef)>(obj,
					DecrementRef) {
	}
};

template<typename T, char const *Key, LIBSAKURA_SYMBOL(Status) (*Func)(T*)>
void TCapsuleDestructor(PyObject *obj) {
	if (!PyCapsule_IsValid(obj, Key)) {
		return;
	}
	auto ptr = reinterpret_cast<T *>(PyCapsule_GetPointer(obj, Key));
	Func(ptr);
}

void TCapsuleDestructorForAlignedBuffer(PyObject *obj) {
	if (!PyCapsule_IsValid(obj, kAlignedBufferName)) {
		return;
	}
	auto ptr =
			reinterpret_cast<LIBSAKURA_SYMBOL(PyAlignedBuffer) *>(PyCapsule_GetPointer(
					obj, kAlignedBufferName));
	LIBSAKURA_SYMBOL(PyAlignedBufferDestroy)(ptr);
}

PyObject *Initialize(PyObject *self, PyObject *args) {
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Initialize)(nullptr,
			nullptr);
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		struct module_state *st = GETSTATE(self);
		auto sakura_error = st->error;
		PyErr_SetString(sakura_error, "Failed to initialize libsakura.");
		return nullptr;
	}
	Py_RETURN_NONE; // equivalent to Py_INCREF(Py_None); return Py_None;
}

PyObject *CleanUp(PyObject *self, PyObject *args) {
	LIBSAKURA_SYMBOL(CleanUp)();
	Py_RETURN_NONE;
}

struct AlignedBufferConfiguration {
	LIBSAKURA_SYMBOL(PyTypeId) type;
	std::function<bool(LIBSAKURA_SYMBOL(PyAlignedBuffer) const &)> pred;
};

bool TotalElementsGreaterOrEqual(LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf,
		size_t n) {
	assert(0 < buf.dimensions && buf.dimensions <= kMaxNumberOfDimensions);
	size_t total_elements = 1;
	for (size_t i = 0; i < buf.dimensions; ++i) {
		total_elements *= buf.elements[i];
	}
	return total_elements >= n;
}

void *GetPointer(PyArrayObject *obj) {

	PyObject *base = PyArray_BASE(obj);

	// check if base object is valid capsule object
	if (!PyCapsule_IsValid(base, kAlignedBufferForNPName)) {
		return nullptr;
	}

	// check if pointer hold by capsule is aligned
	void *base_ptr = PyCapsule_GetPointer(base, kAlignedBufferForNPName);
	if (!LIBSAKURA_SYMBOL(IsAligned)(base_ptr)) {
		return nullptr;
	}

	// now we can use a pointer hold by numpy array
	void *ptr = PyArray_DATA(obj);

	return ptr;
}

LIBSAKURA_SYMBOL(PyTypeId) MapFromNumPyType(int nptype) {
	LIBSAKURA_SYMBOL(PyTypeId) type = LIBSAKURA_SYMBOL(PyTypeId_kEnd);
	switch(nptype) {
	case(NPY_BOOL):
		type = LIBSAKURA_SYMBOL(PyTypeId_kBool);
		break;
	case(NPY_INT8): // NPY_BYTE is equivalent
		type = LIBSAKURA_SYMBOL(PyTypeId_kInt8);
		break;
	case(NPY_INT32): // NPY_INT is equivalent
		type = LIBSAKURA_SYMBOL(PyTypeId_kInt32);
		break;
	case(NPY_INT64):
	case(NPY_LONGLONG):
		type = LIBSAKURA_SYMBOL(PyTypeId_kInt64);
		break;
	case(NPY_UINT8): // NPY_BYTE is equivalent
		type = LIBSAKURA_SYMBOL(PyTypeId_kUInt8);
		break;
	case(NPY_UINT32): // NPY_INT is equivalent
		type = LIBSAKURA_SYMBOL(PyTypeId_kUInt32);
		break;
	case(NPY_FLOAT32): // NPY_FLOAT is equivalent
		type = LIBSAKURA_SYMBOL(PyTypeId_kFloat);
		break;
	case(NPY_FLOAT64): // NPY_DOUBLE is equivalent
		type = LIBSAKURA_SYMBOL(PyTypeId_kDouble);
		break;
	case(NPY_LONGDOUBLE):
		type = LIBSAKURA_SYMBOL(PyTypeId_kLongDouble);
		break;
	default:
		break;
	}
	return type;
}

//int MapToNumPyType(LIBSAKURA_SYMBOL(PyTypeId) type) {
//	int nptype = NPY_NOTYPE;
//	switch(type) {
//	case(LIBSAKURA_SYMBOL(PyTypeId_kBool)):
//		nptype = NPY_BOOL;
//		break;
//	case(LIBSAKURA_SYMBOL(PyTypeId_kInt8)):
//		nptype = NPY_INT8;
//		break;
//	case(LIBSAKURA_SYMBOL(PyTypeId_kInt32)):
//		nptype = NPY_INT32;
//		break;
//	case(LIBSAKURA_SYMBOL(PyTypeId_kInt64)):
//		nptype = NPY_INT64;
//		break;
//	case(LIBSAKURA_SYMBOL(PyTypeId_kFloat)):
//		nptype = NPY_FLOAT;
//		break;
//	case(LIBSAKURA_SYMBOL(PyTypeId_kDouble)):
//		nptype = NPY_DOUBLE;
//		break;
//	case(LIBSAKURA_SYMBOL(PyTypeId_kLongDouble)):
//		nptype = NPY_LONGDOUBLE;
//		break;
//	default:
//		break;
//	}
//	return nptype;
//}

bool IsValidAlignedNumPyArray(size_t num, AlignedBufferConfiguration const conf[],
		PyObject *arr[], LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[]) {
	for (size_t i = 0; i < num; ++i) {
		// check if given PyObject is compatible with PyArray_Type
		if (!PyArray_Check(arr[i])) {
			printf("not a numpy array\n");
			return false;
		}

		// check if PyArray has valid type
		PyArrayObject *nparr = (PyArrayObject *)arr[i];
		if (MapFromNumPyType(PyArray_TYPE(nparr)) != conf[i].type) {
			printf("type is not expected (%lu): numpy type %d sakura type %d expected %d\n",
					i, PyArray_TYPE(nparr), MapFromNumPyType(PyArray_TYPE(nparr)), conf[i].type);
			return false;
		}

		// check if given PyObject has valid base object
		// (otherwise, the array might not satisfy sakura requirement)
		void *data_ptr = GetPointer(nparr);
		if (!data_ptr) {
			printf("failed to get pointer\n");
			return false;
		}

		// get number of dimensions
		int ndim = PyArray_NDIM(nparr);
		if (ndim <= 0 || kMaxNumberOfDimensions <= (size_t)ndim ) {
			printf("Invalid number of dimensins\n");
			return false;
		}

		// get length of each dimension
		npy_intp *dims = PyArray_DIMS(nparr);

		// fill in PyAlignedBuffer
		bufs[i].type = conf[i].type;
		bufs[i].aligned_addr = data_ptr;
		bufs[i].original_addr = data_ptr;
		bufs[i].dimensions = ndim;
		for (int j = 0; j < ndim; ++j) {
			bufs[i].elements[j] = dims[j];
		}

		// user defined check function
		if (!conf[i].pred(bufs[i])) {
			printf("(%lu): user-defined verification failed\n", i);
			return false;
		}
	}
	return true;
}

bool IsValidAlignedNumPyArrayAllowNone(size_t num, AlignedBufferConfiguration const conf[],
		bool const allow_none[], PyObject *arr[], LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[]) {
	for (size_t i = 0; i < num; ++i) {
		// if parameter allows to set None, and if its value is actually None,
		// set nullptr instead of decapsulate PyAlignedBuffer
		if (allow_none[i] && arr[i] == Py_None) {
			bufs[i].type = conf[i].type;
			bufs[i].aligned_addr = nullptr;
			bufs[i].original_addr = nullptr;
			bufs[i].dimensions = 0;
			continue;
		}

		// check if given PyObject is compatible with PyArray_Type
		if (!PyArray_Check(arr[i])) {
			printf("not a numpy array\n");
			return false;
		}

		// check if PyArray has valid type
		PyArrayObject *nparr = (PyArrayObject *)arr[i];
		if (MapFromNumPyType(PyArray_TYPE(nparr)) != conf[i].type) {
			printf("type is not expected (%lu): numpy type %d sakura type %d expected %d\n",
					i, PyArray_TYPE(nparr), MapFromNumPyType(PyArray_TYPE(nparr)), conf[i].type);
			return false;
		}

		// check if given PyObject has valid base object
		// (otherwise, the array might not satisfy sakura requirement)
		void *data_ptr = GetPointer(nparr);
		if (!data_ptr) {
			printf("failed to get pointer\n");
			return false;
		}

		// get number of dimensions
		int ndim = PyArray_NDIM(nparr);
		if (ndim <= 0 || kMaxNumberOfDimensions <= (size_t)ndim ) {
			printf("Invalid number of dimensins\n");
			return false;
		}

		// get length of each dimension
		npy_intp *dims = PyArray_DIMS(nparr);

		// fill in PyAlignedBuffer
		bufs[i].type = conf[i].type;
		bufs[i].aligned_addr = data_ptr;
		bufs[i].original_addr = data_ptr;
		bufs[i].dimensions = ndim;
		for (int j = 0; j < ndim; ++j) {
			bufs[i].elements[j] = dims[j];
		}

		// user defined check function
		if (!conf[i].pred(bufs[i])) {
			printf("(%lu): user-defined verification failed\n", i);
			return false;
		}
	}
	return true;
}

//bool isValidAlignedBuffer(size_t num, AlignedBufferConfiguration const conf[],
//		PyObject *capsules[],
//		LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[]) {
//	for (size_t i = 0; i < num; ++i) {
//		LIBSAKURA_SYMBOL(PyAlignedBuffer) *buf = nullptr;
//		if (LIBSAKURA_SYMBOL(PyAlignedBufferDecapsulate)(capsules[i],
//				&buf) != LIBSAKURA_SYMBOL(Status_kOK)) {
//			return false;
//		}
//		assert(buf);
//		if (!(buf->type == conf[i].type && conf[i].pred(*buf))) {
//			return false;
//		}
//		bufs[i] = buf;
//	}
//	return true;
//}
//
//bool isValidAlignedBufferAllowNone(size_t num,
//		AlignedBufferConfiguration const conf[],
//		bool const allow_none[], PyObject *capsules[],
//		LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[]) {
//	for (size_t i = 0; i < num; ++i) {
//		// if parameter allows to set None, and if its value is actually None,
//		// set nullptr instead of decapsulate PyAlignedBuffer
//		if (allow_none[i] && capsules[i] == Py_None) {
//			bufs[i] = nullptr;
//			continue;
//		}
//		LIBSAKURA_SYMBOL(PyAlignedBuffer) *buf = nullptr;
//		if (LIBSAKURA_SYMBOL(PyAlignedBufferDecapsulate)(capsules[i],
//				&buf) != LIBSAKURA_SYMBOL(Status_kOK)) {
//			return false;
//		}
//		assert(buf);
//		if (!(buf->type == conf[i].type && conf[i].pred(*buf))) {
//			return false;
//		}
//		bufs[i] = buf;
//	}
//	return true;
//}

// For numpy array holding aligned pointer
// destructor for PyCapsule
void DestructPyCapsuleForNP(PyObject *obj) {
    if (!PyCapsule_IsValid(obj, kAlignedBufferForNPName)) {
        return;
    }
    void *ptr = PyCapsule_GetPointer(obj, kAlignedBufferForNPName);
    printf("LOG: Deallocate aligned buffer for numpy address is %p\n", ptr);
    free(ptr);
}

PyObject *NewUninitializedAlignedNumPyArray(PyObject *self, PyObject *args) {
	PyObject *shape = nullptr;
	int type;
//	printf("REFCOUNT: 0 shape %ld\n", (long)Py_REFCNT(&args[1]));
	if (!PyArg_ParseTuple(args, "iO", &type, &shape)) {
		return nullptr;
	}

	// shape should be a tuple
//	printf("REFCOUNT: 1 shape %ld\n", (long)Py_REFCNT(shape));
	int is_tuple = PyTuple_Check(shape);
	if (!is_tuple) {
		PyErr_SetString(PyExc_ValueError, "Second argument should be a shape tuple.");
		return nullptr;
	}

	// numpy shape
	Py_ssize_t len = PyTuple_Size(shape);
//	printf("LOG: tuple size %ld\n", (long)len);
	std::unique_ptr<npy_intp[]> dims(new npy_intp[len]);
	for (Py_ssize_t i = 0; i < len; ++i) {
		auto item = PyTuple_GetItem(shape, i);
		dims[i] = (npy_intp)PyLong_AsLong(item);
//		printf("LOG: tuple item (%ld) %ld\n", i, dims[i]);
	}

	// create numpy array from alinged pointer
	// assume type is float at this moment
	ssize_t num_elements = 1;
	for (Py_ssize_t i = 0; i < len; ++i) {
		num_elements *= dims[i];
	}

	// type mapping
	size_t element_size = 0;
	int nptype = -1;
	switch(type) {
	case(LIBSAKURA_SYMBOL(PyTypeId_kBool)): {
		element_size = sizeof(bool);
		nptype = NPY_BOOL;
	}
	break;
	case(LIBSAKURA_SYMBOL(PyTypeId_kInt8)): {
		element_size = sizeof(int8_t);
		nptype = NPY_INT8;
	}
	break;
	case(LIBSAKURA_SYMBOL(PyTypeId_kInt32)): {
		element_size = sizeof(int32_t);
		nptype = NPY_INT32;
	}
	break;
	case(LIBSAKURA_SYMBOL(PyTypeId_kInt64)): {
		element_size = sizeof(int64_t);
		nptype = NPY_INT64;
	}
	break;
	case(LIBSAKURA_SYMBOL(PyTypeId_kUInt8)): {
		element_size = sizeof(uint8_t);
		nptype = NPY_UINT8;
	}
	break;
	case(LIBSAKURA_SYMBOL(PyTypeId_kUInt32)): {
		element_size = sizeof(uint32_t);
		nptype = NPY_UINT32;
	}
	break;
	case(LIBSAKURA_SYMBOL(PyTypeId_kFloat)): {
		element_size = sizeof(float);
		nptype = NPY_FLOAT;
	}
	break;
	case(LIBSAKURA_SYMBOL(PyTypeId_kDouble)): {
		element_size = sizeof(double);
		nptype = NPY_DOUBLE;
	}
	break;
	case(LIBSAKURA_SYMBOL(PyTypeId_kLongDouble)): {
		element_size = sizeof(long double);
		nptype = NPY_LONGDOUBLE;
	}
	break;
	default:
		// unsupported type
		PyErr_SetString(PyExc_ValueError, "Unsupported data type.");
		return nullptr;
	}

	// allocate memory with aligned pointer
	void *p = nullptr;
	int status = posix_memalign(&p, sakura_GetAlignment(), num_elements * element_size);
	printf("LOG: allocate memory address is %p\n", p);
	if (status != 0 || !sakura_IsAligned(p)) {
		PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory.");
		return nullptr;
	}
	// encapsulate pointer
	PyObject *capsule = PyCapsule_New(p, kAlignedBufferForNPName, DestructPyCapsuleForNP);
	RefHolder rh_capsule(capsule);
//	printf("REFCOUNT: capsule %ld\n", (long)Py_REFCNT(capsule));

	// create array
	PyObject *arr = PyArray_SimpleNewFromData(len, dims.get(), nptype, p);
	RefHolder rh_arr(arr);
//	printf("REFCOUNT: arr %ld\n", (long)Py_REFCNT(arr));
	if (!PyArray_Check(arr)) {
		PyErr_SetString(PyExc_RuntimeError, "Failed to create ndarray.");
		return nullptr;
	}

	// set capsule object as a base object of the array
	int s = PyArray_SetBaseObject((PyArrayObject *)arr, capsule);
	if (s != 0) {
		PyErr_SetString(PyExc_RuntimeError, "Failed to set base for ndarray.");
		return nullptr;
	}

	// give ownership of capsule object to arr
	rh_capsule.release();

	// arr is properly created so release the ownership
	rh_arr.release();

//	printf("REFCOUNT: capsule %ld\n", (long)Py_REFCNT(capsule));
//
//	printf("REFCOUNT: 2 shape %ld\n", (long)Py_REFCNT(shape));

	return arr;
}
//

template<LIBSAKURA_SYMBOL(Status) (*Func)(size_t, float const *, bool const *,
LIBSAKURA_SYMBOL(StatisticsResultFloat) *)>
PyObject *ComputeStatisticsTemplateForNumPy(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kData, kIsValid
	};
	PyObject *arr[2];

	if (!PyArg_ParseTuple(args, "nOO", &num_data_py, &arr[kData],
			&arr[kIsValid])) {
		return nullptr;
	}
	printf("parsing arguments done\n");

	auto num_data = static_cast<size_t>(num_data_py);
	auto pred = [num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
	{	return TotalElementsGreaterOrEqual(buf, num_data);};
	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), pred },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), pred },

	};

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[2];
	printf("argument validity check\n");
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}
	printf("run sakura function\n");
	LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = Func(num_data,
				reinterpret_cast<float const*>(bufs[kData].aligned_addr),
				reinterpret_cast<bool const*>(bufs[kIsValid].aligned_addr),
				&result);
		SAKURA_END_ALLOW_THREADS
	printf("check sakura status\n");
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	return Py_BuildValue("{s:n,s:n,s:n,s:f,s:f,s:d,s:d}", "count",
			static_cast<Py_ssize_t>(result.count), "index_of_min",
			result.index_of_min, "index_of_max", result.index_of_max, "min",
			result.min, "max", result.max, "sum", result.sum, "square_sum",
			result.square_sum);

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

constexpr FuncForPython ComputeStatisticsNP = ComputeStatisticsTemplateForNumPy<
		sakura_ComputeStatisticsFloat>;
constexpr FuncForPython ComputeAccurateStatisticsNP = ComputeStatisticsTemplateForNumPy<
		sakura_ComputeAccurateStatisticsFloat>;

PyObject *GridConvolvingNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_spectra_py;
	Py_ssize_t start_spectrum_py;
	Py_ssize_t end_spectrum_py;
	Py_ssize_t support_py;
	Py_ssize_t sampling_py;
	Py_ssize_t num_polarizations_py;
	Py_ssize_t num_channels_py;
	Py_ssize_t num_convolution_table_py;
	Py_ssize_t num_polarizations_for_grid_py;
	Py_ssize_t num_channels_for_grid_py;
	Py_ssize_t width_py;
	Py_ssize_t height_py;
	PyObject *weight_only = nullptr;
	enum {
		kSpectrumMask,
		kX,
		kY,
		kPolarizationMap,
		kChannelMap,
		kMask,
		kValue,
		kWeight,
		kConvolutionTable,
		kWeightSum,
		kWeightOfGrid,
		kGrid,
		kEnd
	};
	PyObject *arr[kEnd];

	if (!PyArg_ParseTuple(args, "nnnOOOnnnOnOOOOOnOnnnnOOO", &num_spectra_py,
			&start_spectrum_py, &end_spectrum_py, &arr[kSpectrumMask],
			&arr[kX], &arr[kY], &support_py, &sampling_py,
			&num_polarizations_py, &arr[kPolarizationMap],
			&num_channels_py, &arr[kChannelMap], &arr[kMask],
			&arr[kValue], &arr[kWeight], &weight_only,
			&num_convolution_table_py, &arr[kConvolutionTable],
			&num_polarizations_for_grid_py, &num_channels_for_grid_py,
			&width_py, &height_py, &arr[kWeightSum],
			&arr[kWeightOfGrid], &arr[kGrid])) {
		return nullptr;
	}
	auto num_spectra = static_cast<size_t>(num_spectra_py);
	auto start_spectrum = static_cast<size_t>(start_spectrum_py);
	auto end_spectrum = static_cast<size_t>(end_spectrum_py);
	auto support = static_cast<size_t>(support_py);
	auto sampling = static_cast<size_t>(sampling_py);
	auto num_polarizations = static_cast<size_t>(num_polarizations_py);
	auto num_channels = static_cast<size_t>(num_channels_py);
	auto num_convolution_table = static_cast<size_t>(num_convolution_table_py);
	auto num_polarizations_for_grid =
			static_cast<size_t>(num_polarizations_for_grid_py);
	auto num_channels_for_grid = static_cast<size_t>(num_channels_for_grid_py);
	auto width = static_cast<size_t>(width_py);
	auto height = static_cast<size_t>(height_py);

	auto pred_num_spectra =
			[=](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_spectra);};
	auto pred_num_polarizations =
			[=](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_polarizations);};
	auto pred_num_channels =
			[=](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_channels);};
	auto pred_sp_pol_ch =
			[=](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_spectra * num_polarizations * num_channels);};
	auto pred_sp_ch = [=](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
	{	return TotalElementsGreaterOrEqual(buf, num_spectra * num_channels);};
	auto pred_num_convolution_table =
			[=](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_convolution_table);};
	auto pred_pol_ch_for_grid =
			[=](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_polarizations_for_grid * num_channels_for_grid);};
	auto pred_h_w_pol_ch_for_grid =
			[=](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, height * width * num_polarizations_for_grid * num_channels_for_grid);};

	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), pred_num_spectra },

	{ LIBSAKURA_SYMBOL(PyTypeId_kDouble), pred_num_spectra },

	{ LIBSAKURA_SYMBOL(PyTypeId_kDouble), pred_num_spectra },

	{ LIBSAKURA_SYMBOL(PyTypeId_kInt32), pred_num_polarizations },

	{ LIBSAKURA_SYMBOL(PyTypeId_kInt32), pred_num_channels },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), pred_sp_pol_ch },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), pred_sp_pol_ch },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), pred_sp_ch },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), pred_num_convolution_table },

	{ LIBSAKURA_SYMBOL(PyTypeId_kDouble), pred_pol_ch_for_grid },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), pred_h_w_pol_ch_for_grid },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), pred_h_w_pol_ch_for_grid },

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status =
				LIBSAKURA_SYMBOL(GridConvolvingFloat)(num_spectra,
						start_spectrum, end_spectrum,
						static_cast<bool const *>(bufs[kSpectrumMask].aligned_addr),
						static_cast<double const *>(bufs[kX].aligned_addr),
						static_cast<double const *>(bufs[kY].aligned_addr),
						support, sampling, num_polarizations,
						static_cast<uint32_t const *>(bufs[kPolarizationMap].aligned_addr),
						num_channels,
						static_cast<uint32_t const *>(bufs[kChannelMap].aligned_addr),
						static_cast<bool const *>(bufs[kMask].aligned_addr),
						static_cast<float const *>(bufs[kValue].aligned_addr),
						static_cast<float const *>(bufs[kWeight].aligned_addr),
						weight_only == Py_True, num_convolution_table,
						static_cast<float const *>(bufs[kConvolutionTable].aligned_addr),
						num_polarizations_for_grid, num_channels_for_grid,
						width, height,
						static_cast<double *>(bufs[kWeightSum].aligned_addr),
						static_cast<float *>(bufs[kWeightOfGrid].aligned_addr),
						static_cast<float *>(bufs[kGrid].aligned_addr));
		//fprintf(stderr, "Grid\n");
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_RETURN_NONE;

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

template<typename FromType, typename ToType,
LIBSAKURA_SYMBOL(PyTypeId) FromTypeId, LIBSAKURA_SYMBOL(PyTypeId) ToTypeId,
LIBSAKURA_SYMBOL(Status) (*Func)(size_t, FromType const *, ToType *)>
PyObject *ConvertArrayForNumPy(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kData, kResult, kEnd
	};
	PyObject *arr[kEnd];

	if (!PyArg_ParseTuple(args, "nOO", &num_data_py, &arr[kData],
			&arr[kResult])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	auto pred = [num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
	{	return TotalElementsGreaterOrEqual(buf, num_data);};
	AlignedBufferConfiguration const conf[] = {

	{ FromTypeId, pred },

	{ ToTypeId, pred },

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = Func(num_data,
				reinterpret_cast<FromType const*>(bufs[kData].aligned_addr),
				reinterpret_cast<ToType *>(bufs[kResult].aligned_addr));
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(arr[kResult]);
	return arr[kResult];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

constexpr FuncForPython Uint8ToBoolNP = ConvertArrayForNumPy<uint8_t, bool,
		LIBSAKURA_SYMBOL(PyTypeId_kInt8), LIBSAKURA_SYMBOL(PyTypeId_kBool),
		LIBSAKURA_SYMBOL(Uint8ToBool)>;

constexpr FuncForPython Uint32ToBoolNP = ConvertArrayForNumPy<uint32_t, bool,
		LIBSAKURA_SYMBOL(PyTypeId_kInt32), LIBSAKURA_SYMBOL(PyTypeId_kBool),
		LIBSAKURA_SYMBOL(Uint32ToBool)>;

constexpr FuncForPython InvertBoolNP = ConvertArrayForNumPy<bool, bool,
		LIBSAKURA_SYMBOL(PyTypeId_kBool), LIBSAKURA_SYMBOL(PyTypeId_kBool),
		LIBSAKURA_SYMBOL(InvertBool)>;

constexpr FuncForPython SetFalseFloatIfNanOrInfNP = ConvertArrayForNumPy<float, bool,
		LIBSAKURA_SYMBOL(PyTypeId_kFloat), LIBSAKURA_SYMBOL(PyTypeId_kBool),
		LIBSAKURA_SYMBOL(SetFalseIfNanOrInfFloat)>;

constexpr FuncForPython ComputeMADNP = ConvertArrayForNumPy<float, float,
		LIBSAKURA_SYMBOL(PyTypeId_kFloat), LIBSAKURA_SYMBOL(PyTypeId_kFloat),
		LIBSAKURA_SYMBOL(ComputeMedianAbsoluteDeviationFloat)>;

PyObject *SortValidDataDenselyFloatNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kMask, kData, kEnd
	};
	PyObject *arr[kEnd];
	if (!PyArg_ParseTuple(args, "nOO", &num_data_py, &arr[kMask],
			&arr[kData])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	auto pred = [num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
	{	return TotalElementsGreaterOrEqual(buf, num_data);};

	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), pred },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), pred },

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}

	LIBSAKURA_SYMBOL(Status) status;
	size_t new_num_data;
	SAKURA_BEGIN_ALLOW_THREADS
		status = LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(num_data,
				reinterpret_cast<bool const*>(bufs[kMask].aligned_addr),
				reinterpret_cast<float *>(bufs[kData].aligned_addr),
				&new_num_data);
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}

	return PyLong_FromSize_t(new_num_data);

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

inline void AlignedBoolAnd(size_t elements, bool const src1[],
bool const src2[], bool dst[]) {
	auto dst_u8 = AssumeAligned(reinterpret_cast<uint8_t *>(dst));
	auto src1_u8 = AssumeAligned(reinterpret_cast<uint8_t const*>(src1));
	auto src2_u8 = AssumeAligned(reinterpret_cast<uint8_t const*>(src2));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);
	STATIC_ASSERT(sizeof(*dst_u8) == sizeof(*dst));
	for (size_t i = 0; i < elements; ++i) {
		dst_u8[i] = src1_u8[i] & src2_u8[i];
	}
}

PyObject *LogicalAndNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kData1, kData2, kResult, kEnd
	};
	PyObject *arr[kEnd];

	if (!PyArg_ParseTuple(args, "nOOO", &num_data_py, &arr[kData1],
			&arr[kData2], &arr[kResult])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	auto pred = [num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
	{	return TotalElementsGreaterOrEqual(buf, num_data);};
	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), pred },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), pred },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), pred },

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}
	AlignedBoolAnd(num_data,
			reinterpret_cast<bool const*>(bufs[kData1].aligned_addr),
			reinterpret_cast<bool const*>(bufs[kData2].aligned_addr),
			reinterpret_cast<bool *>(bufs[kResult].aligned_addr));
	Py_INCREF(arr[kResult]);
	return arr[kResult];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

template<typename Type,
LIBSAKURA_SYMBOL(PyTypeId) TypeId,
LIBSAKURA_SYMBOL(Status) (*Func)(size_t, Type const *, size_t, Type const *,
		Type const *, bool *)>
PyObject *RangeCheckNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	Py_ssize_t num_range_py;
	enum {
		kData, kLower, kUpper, kMask, kEnd
	};
	PyObject *arr[kEnd];

	if (!PyArg_ParseTuple(args, "nOnOOO", &num_data_py, &arr[kData],
			&num_range_py, &arr[kLower], &arr[kUpper],
			&arr[kMask])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	auto num_range = static_cast<size_t>(num_range_py);
	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};
	auto predRange =
			[num_range](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_range);};
	AlignedBufferConfiguration const conf[] = {

	{ TypeId, predData },

	{ TypeId, predRange },

	{ TypeId, predRange },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), predData },

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = Func(num_data,
				reinterpret_cast<Type const*>(bufs[kData].aligned_addr),
				num_range,
				reinterpret_cast<Type const*>(bufs[kLower].aligned_addr),
				reinterpret_cast<Type const*>(bufs[kUpper].aligned_addr),
				reinterpret_cast<bool *>(bufs[kMask].aligned_addr));
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(arr[kMask]);
	return arr[kMask];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

constexpr FuncForPython Int32SetTrueIntInRangesExclusiveNP = RangeCheckNP<int32_t,
		LIBSAKURA_SYMBOL(PyTypeId_kInt32),
		LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveInt)>;

constexpr FuncForPython FloatSetTrueIntInRangesExclusiveNP = RangeCheckNP<float,
		LIBSAKURA_SYMBOL(PyTypeId_kFloat),
		LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveFloat)>;

struct Py2Type {
	static void Cast(PyObject *obj, float *v) {
		*v = static_cast<float>(PyFloat_AsDouble(obj));
	}
	static void Cast(PyObject *obj, int *v) {
		*v = static_cast<int>(PyLong_AsLong(obj));
	}
};

template<typename Type,
LIBSAKURA_SYMBOL(PyTypeId) TypeId,
LIBSAKURA_SYMBOL(Status) (*Func)(size_t, Type const *, Type, bool *)>
PyObject *ThresholdingNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	PyObject *threshold_py;
	enum {
		kData, kMask, kEnd
	};
	PyObject *arr[kEnd];

	if (!PyArg_ParseTuple(args, "nOOO", &num_data_py, &arr[kData],
			&threshold_py, &arr[kMask])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	Type threshold;
	Py2Type::Cast(threshold_py, &threshold);

	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};
	AlignedBufferConfiguration const conf[] = {

	{ TypeId, predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), predData },

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = Func(num_data,
				reinterpret_cast<Type const*>(bufs[kData].aligned_addr),
				threshold, reinterpret_cast<bool *>(bufs[kMask].aligned_addr));
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(arr[kMask]);
	return arr[kMask];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

constexpr FuncForPython FloatSetTrueIfGreaterThanNP = ThresholdingNP<float,
		LIBSAKURA_SYMBOL(PyTypeId_kFloat),
		LIBSAKURA_SYMBOL(SetTrueIfGreaterThanFloat)>;

constexpr FuncForPython IntSetTrueIfGreaterThanNP = ThresholdingNP<int,
		LIBSAKURA_SYMBOL(PyTypeId_kInt32),
		LIBSAKURA_SYMBOL(SetTrueIfGreaterThanInt)>;

constexpr FuncForPython FloatSetTrueIfGreaterThanOrEqualsNP = ThresholdingNP<float,
		LIBSAKURA_SYMBOL(PyTypeId_kFloat),
		LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsFloat)>;

constexpr FuncForPython IntSetTrueIfGreaterThanOrEqualsNP = ThresholdingNP<int,
		LIBSAKURA_SYMBOL(PyTypeId_kInt32),
		LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsInt)>;

constexpr FuncForPython FloatSetTrueIfLessThanNP = ThresholdingNP<float,
		LIBSAKURA_SYMBOL(PyTypeId_kFloat),
		LIBSAKURA_SYMBOL(SetTrueIfLessThanFloat)>;

constexpr FuncForPython IntSetTrueIfLessThanNP = ThresholdingNP<int,
		LIBSAKURA_SYMBOL(PyTypeId_kInt32),
		LIBSAKURA_SYMBOL(SetTrueIfLessThanInt)>;

constexpr FuncForPython FloatSetTrueIfLessThanOrEqualsNP = ThresholdingNP<float,
		LIBSAKURA_SYMBOL(PyTypeId_kFloat),
		LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsFloat)>;

constexpr FuncForPython IntSetTrueIfLessThanOrEqualsNP = ThresholdingNP<int,
		LIBSAKURA_SYMBOL(PyTypeId_kInt32),
		LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsInt)>;

template<typename Type,
LIBSAKURA_SYMBOL(PyTypeId) TypeId,
LIBSAKURA_SYMBOL(Status) (*Func)(Type, size_t, Type const *, bool const *,
		Type *)>
PyObject *BinaryBitOperationNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	unsigned int bitmask_py;

	enum {
		kData, kMask, kResult, kEnd
	};
	PyObject *arr[kEnd];

	if (!PyArg_ParseTuple(args, "InOOO", &bitmask_py, &num_data_py,
			&arr[kData], &arr[kMask], &arr[kResult])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	auto bitmask = static_cast<Type>(bitmask_py);
	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};

	AlignedBufferConfiguration const conf[] = {

	{ TypeId, predData }, { LIBSAKURA_SYMBOL(PyTypeId_kBool), predData }, {
			TypeId, predData }
	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = Func(bitmask, num_data,
				reinterpret_cast<Type const*>(bufs[kData].aligned_addr),
				reinterpret_cast<bool const*>(bufs[kMask].aligned_addr),
				reinterpret_cast<Type *>(bufs[kResult].aligned_addr));
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(arr[kResult]);
	return arr[kResult];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

constexpr FuncForPython Uint8OperateBitsOrNP = BinaryBitOperationNP<uint8_t,
		LIBSAKURA_SYMBOL(PyTypeId_kUInt8),
		LIBSAKURA_SYMBOL(OperateBitwiseOrUint8)>;

constexpr FuncForPython Uint32OperateBitsOrNP = BinaryBitOperationNP<uint32_t,
		LIBSAKURA_SYMBOL(PyTypeId_kUInt32),
		LIBSAKURA_SYMBOL(OperateBitwiseOrUint32)>;

constexpr FuncForPython Uint8OperateBitsAndNP = BinaryBitOperationNP<uint8_t,
		LIBSAKURA_SYMBOL(PyTypeId_kUInt8),
		LIBSAKURA_SYMBOL(OperateBitwiseAndUint8)>;

constexpr FuncForPython Uint32OperateBitsAndNP = BinaryBitOperationNP<uint32_t,
		LIBSAKURA_SYMBOL(PyTypeId_kUInt32),
		LIBSAKURA_SYMBOL(OperateBitwiseAndUint32)>;

constexpr FuncForPython Uint8OperateBitsXorNP = BinaryBitOperationNP<uint8_t,
		LIBSAKURA_SYMBOL(PyTypeId_kUInt8),
		LIBSAKURA_SYMBOL(OperateBitwiseXorUint8)>;

constexpr FuncForPython Uint32OperateBitsXorNP = BinaryBitOperationNP<uint32_t,
		LIBSAKURA_SYMBOL(PyTypeId_kUInt32),
		LIBSAKURA_SYMBOL(OperateBitwiseXorUint32)>;

template<typename Type,
LIBSAKURA_SYMBOL(PyTypeId) TypeId,
LIBSAKURA_SYMBOL(Status) (*Func)(size_t, Type const *, bool const *, Type *)>
PyObject *UnaryBitOperationNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;

	enum {
		kData, kMask, kResult, kEnd
	};
	PyObject *arr[kEnd];

	if (!PyArg_ParseTuple(args, "nOOO", &num_data_py, &arr[kData],
			&arr[kMask], &arr[kResult])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};

	AlignedBufferConfiguration const conf[] = {

	{ TypeId, predData }, { LIBSAKURA_SYMBOL(PyTypeId_kBool), predData }, {
			TypeId, predData }

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = Func(num_data,
				reinterpret_cast<Type const*>(bufs[kData].aligned_addr),
				reinterpret_cast<bool const*>(bufs[kMask].aligned_addr),
				reinterpret_cast<Type *>(bufs[kResult].aligned_addr));
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(arr[kResult]);
	return arr[kResult];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

constexpr FuncForPython Uint8OperateBitsNotNP = UnaryBitOperationNP<uint8_t,
		LIBSAKURA_SYMBOL(PyTypeId_kUInt8),
		LIBSAKURA_SYMBOL(OperateBitwiseNotUint8)>;

constexpr FuncForPython Uint32OperateBitsNotNP = UnaryBitOperationNP<uint32_t,
		LIBSAKURA_SYMBOL(PyTypeId_kUInt32),
		LIBSAKURA_SYMBOL(OperateBitwiseNotUint32)>;

template<typename Type,
LIBSAKURA_SYMBOL(PyTypeId) TypeId,
LIBSAKURA_SYMBOL(Status) (*Func)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		double const base_position[/*num_base*/], size_t num_array,
		Type const base_data[/*num_base*num_array*/],
		bool const base_mask[/*num_base*num_array*/], size_t num_interpolated,
		double const interpolated_position[/*num_interpolated*/],
		float interpolated_data[/*num_interpolated*num_array*/],
		bool interpolated_mask[/*num_interpolated*num_array*/])>
PyObject *InterpolateAxisNP(PyObject *self, PyObject *args) {
	int interpolate_type_py;
	uint8_t polynomial_order;
	Py_ssize_t num_base_py;
	Py_ssize_t num_array_py;
	Py_ssize_t num_interpolated_py;
	enum {
		kFromAxis, kFromData, kMask, kToAxis, kToData, kToMask, kEnd
	};
	PyObject *arr[kEnd];

	if (!PyArg_ParseTuple(args, "ibnOnOOnOOO", &interpolate_type_py,
			&polynomial_order, &num_base_py, &arr[kFromAxis],
			&num_array_py, &arr[kFromData], &arr[kMask],
			&num_interpolated_py, &arr[kToAxis], &arr[kToData],
			&arr[kToMask])) {
		return nullptr;
	}
	if (!(0 <= interpolate_type_py
			&& interpolate_type_py
					< LIBSAKURA_SYMBOL(InterpolationMethod_kNumElements))) {
		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
		return nullptr;
	}
	auto interpolate_type =
			static_cast<LIBSAKURA_SYMBOL(InterpolationMethod)>(interpolate_type_py);
	auto num_base = static_cast<size_t>(num_base_py);
	auto num_array = static_cast<size_t>(num_array_py);
	auto num_interpolated = static_cast<size_t>(num_interpolated_py);
	auto predFromAxis =
			[num_base](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_base);};
	auto predFromData =
			[num_base, num_array](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_base*num_array);};
	auto predToAxis =
			[num_interpolated](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_interpolated);};
	auto predToData =
			[num_interpolated, num_array](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_interpolated*num_array);};
	auto predMaskIn =
			[num_base, num_array](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_base * num_array);};
	auto predMaskOut =
			[num_interpolated, num_array](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_interpolated * num_array);};
	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kDouble), predFromAxis },

	{ TypeId, predFromData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), predMaskIn },

	{ LIBSAKURA_SYMBOL(PyTypeId_kDouble), predToAxis },

	{ TypeId, predToData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), predMaskOut }, };
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = Func(interpolate_type, polynomial_order, num_base,
				reinterpret_cast<double const*>(bufs[kFromAxis].aligned_addr),
				num_array,
				reinterpret_cast<Type const*>(bufs[kFromData].aligned_addr),
				reinterpret_cast<bool const*>(bufs[kMask].aligned_addr),
				num_interpolated,
				reinterpret_cast<double const*>(bufs[kToAxis].aligned_addr),
				reinterpret_cast<Type *>(bufs[kToData].aligned_addr),
				reinterpret_cast<bool *>(bufs[kToMask].aligned_addr));
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(arr[kToData]);
	return arr[kToData];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

constexpr FuncForPython InterpolateFloatYAxisNP = InterpolateAxisNP<float,
		LIBSAKURA_SYMBOL(PyTypeId_kFloat),
		LIBSAKURA_SYMBOL(InterpolateYAxisFloat)>;

constexpr FuncForPython InterpolateFloatXAxisNP = InterpolateAxisNP<float,
		LIBSAKURA_SYMBOL(PyTypeId_kFloat),
		LIBSAKURA_SYMBOL(InterpolateXAxisFloat)>;

PyObject *CalibrateDataWithArrayScalingFloatNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kFactor, kData, kOff, kResult, kEnd
	};
	PyObject *arr[kEnd];
	if (!PyArg_ParseTuple(args, "nOOOO", &num_data_py, &arr[kFactor],
			&arr[kData], &arr[kOff], &arr[kResult])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};
	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData }

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = LIBSAKURA_SYMBOL(CalibrateDataWithArrayScalingFloat)(num_data,
				reinterpret_cast<float const*>(bufs[kFactor].aligned_addr),
				reinterpret_cast<float const*>(bufs[kData].aligned_addr),
				reinterpret_cast<float const*>(bufs[kOff].aligned_addr),
				reinterpret_cast<float *>(bufs[kResult].aligned_addr));
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(arr[kResult]);
	return arr[kResult];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

PyObject *CreateConvolve1DContextFFTNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_kernel;
	enum {
		kData, kEnd
	};
	PyObject *arr[kEnd];

	if (!PyArg_ParseTuple(args, "nO", &num_kernel, &arr[kData])) {
		return nullptr;
	}

	LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) status;
	auto pred =
			[num_kernel](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_kernel);};
	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), pred }

	};
	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
		return nullptr;
	}
	SAKURA_BEGIN_ALLOW_THREADS
		status =
		LIBSAKURA_SYMBOL(CreateConvolve1DContextFFTFloat)(
				static_cast<size_t>(num_kernel),
				reinterpret_cast<float *>(bufs[kData].aligned_addr), &context);
		SAKURA_END_ALLOW_THREADS

	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError,
				"sakura_CreateConvolve1DContext failed.");
		return nullptr;
	}
	PyCapsule_Destructor destructor = TCapsuleDestructor<
	LIBSAKURA_SYMBOL(Convolve1DContextFloat), kConvolve1DContextName,
			LIBSAKURA_SYMBOL(DestroyConvolve1DContextFloat)>;
	auto capsule = PyCapsule_New(context, kConvolve1DContextName, destructor);
	if (capsule == nullptr) {
		PyErr_SetString(PyExc_MemoryError, "No memory.");
		return nullptr;
	}
	return capsule;
}

PyObject *Convolve1DNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py, num_kernel_py;
	enum {
		kKernel, kData, kMask, kResult, kWeight, kEnd
	};
	PyObject *arr[kEnd];
	if (!PyArg_ParseTuple(args, "nOnOOOO", &num_kernel_py, &arr[kKernel],
			&num_data_py, &arr[kData], &arr[kMask],
			&arr[kResult], &arr[kWeight])) {
		return nullptr;
	}

	auto num_kernel = static_cast<size_t>(num_kernel_py);
	auto num_data = static_cast<size_t>(num_data_py);

	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};
	auto predKernel =
			[num_kernel](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_kernel);};

	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predKernel },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}

	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = LIBSAKURA_SYMBOL(Convolve1DFloat)(num_kernel,
				reinterpret_cast<float const *>(bufs[kKernel].aligned_addr),
				num_data,
				reinterpret_cast<float const *>(bufs[kData].aligned_addr),
				reinterpret_cast<bool const *>(bufs[kMask].aligned_addr),
				reinterpret_cast<float *>(bufs[kResult].aligned_addr),
				reinterpret_cast<float *>(bufs[kWeight].aligned_addr));
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(arr[kResult]);
	return arr[kResult];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

PyObject *Convolve1DFFTNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kData, kResult, kEnd
	};
	PyObject *arr[kEnd];
	PyObject *context_capsule;
	if (!PyArg_ParseTuple(args, "OnOO", &context_capsule, &num_data_py,
			&arr[kData], &arr[kResult])) {
		return nullptr;
	}
	auto context =
			reinterpret_cast<LIBSAKURA_SYMBOL(Convolve1DContextFloat) *>(PyCapsule_GetPointer(
					context_capsule, PyCapsule_GetName(context_capsule)));
//	if (context == nullptr) {
//		goto invalid_arg;
//	};
	auto num_data = static_cast<size_t>(num_data_py);
	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};
	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData }

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)
			|| context == nullptr) {
		goto invalid_arg;
	}

	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = LIBSAKURA_SYMBOL(Convolve1DFFTFloat)(context, num_data,
				reinterpret_cast<float const*>(bufs[kData].aligned_addr),
				reinterpret_cast<float *>(bufs[kResult].aligned_addr));
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(arr[kResult]);
	return arr[kResult];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

PyObject *CreateGaussianKernelFloatNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	float peak_location, kernel_width;
	enum {
		kResult, kEnd
	};
	PyObject *arr[kEnd];
	// TODO: create kernel array internally and return it
	if (!PyArg_ParseTuple(args, "ffnO", &peak_location, &kernel_width,
			&num_data_py, &arr[kResult])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};
	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData }

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}

//	printf("peak %f, width %f\n", peak_location, kernel_width);
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = LIBSAKURA_SYMBOL(CreateGaussianKernelFloat)(peak_location,
				kernel_width, num_data,
				reinterpret_cast<float *>(bufs[kResult].aligned_addr));
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(arr[kResult]);
	return arr[kResult];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

// TODO: CreateBaselineContextCubicSpline
// TODO: CreateBaselineContextSinusoid
PyObject *CreateLSQFitContextPolynomial(PyObject *self, PyObject *args) {
	Py_ssize_t num_data;
	unsigned int order;
	int lsqfit_type;
	if (!PyArg_ParseTuple(args, "iIn", &lsqfit_type, &order, &num_data)) {
		return nullptr;
	}
	if (!(0 <= lsqfit_type
			&& lsqfit_type < LIBSAKURA_SYMBOL(LSQFitType_kNumElements))) {
		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
		return nullptr;
	}

	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextPolynomialFloat)(
				static_cast<LIBSAKURA_SYMBOL(LSQFitType)>(lsqfit_type),
				static_cast<uint16_t>(order), static_cast<size_t>(num_data),
				&context);
		SAKURA_END_ALLOW_THREADS

	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError,
				"sakura_CreateBaselineContext failed.");
		return nullptr;
	}
	PyCapsule_Destructor destructor = TCapsuleDestructor<
	LIBSAKURA_SYMBOL(LSQFitContextFloat), kBaselineContextName,
			LIBSAKURA_SYMBOL(DestroyLSQFitContextFloat)>;
	auto capsule = PyCapsule_New(context, kBaselineContextName, destructor);
	if (capsule == nullptr) {
		PyErr_SetString(PyExc_MemoryError, "No memory.");
		return nullptr;
	}
	return capsule;
}

// TODO: LSQFitCubicSpline
// TODO: LSQFitSinusoid
PyObject *LSQFitPolynomialNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py, num_coeff_py;
	float clip_threshold_sigma;
	unsigned int num_fitting_max, order;
//	uint16_t order = 3; // TODO make this python parameter.
	enum {
		kData, kMask, kCoeff, kBestFit, kResidual, kFinalMask, kEnd
	};
	PyObject *arr[kEnd];
	PyObject *context_capsule;
	if (!PyArg_ParseTuple(args, "OInOOfInOOOO", &context_capsule, &order,
			&num_data_py, &arr[kData], &arr[kMask],
			&clip_threshold_sigma, &num_fitting_max, &num_coeff_py,
			&arr[kCoeff], &arr[kBestFit], &arr[kResidual],
			&arr[kFinalMask])) {
		return nullptr;
	}
	auto context =
			reinterpret_cast<LIBSAKURA_SYMBOL(LSQFitContextFloat) *>(PyCapsule_GetPointer(
					context_capsule, PyCapsule_GetName(context_capsule)));
	if (context == nullptr) {
		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	auto num_coeff = static_cast<size_t>(num_coeff_py);
	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};
	auto predCoeff =
			[num_coeff](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_coeff);};
	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kDouble), predCoeff },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), predData }

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
//	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
//		goto invalid_arg;
//	}
	// indicate if parameter allows None (true) or not (false)
	bool const allow_none[] = {
			false, // kData
			false, // kMask
			true,  // kCoeff
			true,  // kBestFit
			true,  // kResidual
			false  // kFinalMask
			};
	STATIC_ASSERT(ELEMENTSOF(allow_none) == kEnd);
	if (!IsValidAlignedNumPyArrayAllowNone(ELEMENTSOF(conf), conf, allow_none,
			arr, bufs)) {
		goto invalid_arg;
	}

	LIBSAKURA_SYMBOL(Status) status;
	LIBSAKURA_SYMBOL(LSQFitStatus) bl_status;
	float rms;

	SAKURA_BEGIN_ALLOW_THREADS
		status =
				LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context,
						static_cast<uint16_t>(order), num_data,
						reinterpret_cast<float const*>(bufs[kData].aligned_addr),
						reinterpret_cast<bool const*>(bufs[kMask].aligned_addr),
						clip_threshold_sigma,
						static_cast<uint16_t>(num_fitting_max), num_coeff,
						reinterpret_cast<double *>(bufs[kCoeff].aligned_addr),
//						((bufs[kCoeff]) ?
//								reinterpret_cast<double *>(bufs[kCoeff].aligned_addr) :
//								nullptr),
						reinterpret_cast<float *>(bufs[kBestFit].aligned_addr),
//						((bufs[kBestFit]) ?
//								reinterpret_cast<float *>(bufs[kBestFit].aligned_addr) :
//								nullptr),
						reinterpret_cast<float *>(bufs[kResidual].aligned_addr),
//						((bufs[kResidual]) ?
//								reinterpret_cast<float *>(bufs[kResidual].aligned_addr) :
//								nullptr),
						reinterpret_cast<bool *>(bufs[kFinalMask].aligned_addr),
						&rms, &bl_status);
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (bl_status != LIBSAKURA_SYMBOL(LSQFitStatus_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}

	// TODO: return value should be a struct containing all output values
	Py_INCREF(arr[kResidual]);
	return arr[kResidual];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

template<typename T>
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComplementMaskedValue)(
		size_t num_data, T const *data, bool const *mask, T *result) {
	for (size_t i = 0; i < num_data; ++i) {
		T value = data[i];
		if (!mask[i]) {
			value = 0;
		}
		result[i] = value;
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

template<typename Type,
LIBSAKURA_SYMBOL(PyTypeId) TypeId,
LIBSAKURA_SYMBOL(Status) (*Func)(size_t, Type const *, bool const *, Type *)>
PyObject *ComplementMaskedValueNP(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kData, kMask, kResult, kEnd
	};
	PyObject *arr[kEnd];

	if (!PyArg_ParseTuple(args, "nOOO", &num_data_py, &arr[kData],
			&arr[kMask], &arr[kResult])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};

	AlignedBufferConfiguration const conf[] = {

	{ TypeId, predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), predData },

	{ TypeId, predData }

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) bufs[kEnd];
	if (!IsValidAlignedNumPyArray(ELEMENTSOF(conf), conf, arr, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	//SAKURA_BEGIN_ALLOW_THREADS
	status = Func(num_data,
			reinterpret_cast<Type const*>(bufs[kData].aligned_addr),
			reinterpret_cast<bool const*>(bufs[kMask].aligned_addr),
			reinterpret_cast<Type *>(bufs[kResult].aligned_addr));
	//SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(arr[kResult]);
	return arr[kResult];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

constexpr FuncForPython ComplementMaskedValueFloatNP = ComplementMaskedValueNP<
		float, LIBSAKURA_SYMBOL(PyTypeId_kFloat),
		LIBSAKURA_SYMBOL(ComplementMaskedValue)<float> >;

//PyObject *GetElementsOfAlignedBuffer(PyObject *self, PyObject *args) {
//	PyObject *capsule = nullptr;
//
//	if (!PyArg_ParseTuple(args, "O", &capsule)) {
//		return nullptr;
//	}
//	LIBSAKURA_SYMBOL(PyAlignedBuffer) *buf = nullptr;
//	if (LIBSAKURA_SYMBOL(PyAlignedBufferDecapsulate)(capsule,
//			&buf) != LIBSAKURA_SYMBOL(Status_kOK)) {
//		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
//		return nullptr;
//	}
//	assert(buf);
//
//	RefHolder dims(PyTuple_New(buf->dimensions));
//	if (dims.get() == nullptr) {
//		goto out_of_memory;
//	}
//	for (size_t i = 0; i < buf->dimensions; ++i) {
//		auto element = LONG_FROM_SSIZE_T(buf->elements[i]);
//		if (element == nullptr) {
//			goto out_of_memory;
//		}
//		PyTuple_SetItem(dims.get(), i, element);
//	}
//	return dims.release();
//
//	out_of_memory:
//
//	PyErr_SetString(PyExc_MemoryError, "No memory.");
//	return nullptr;
//}

//PyObject *NewUninitializedAlignedBuffer(PyObject *self, PyObject *args) {
//	PyObject *elements;
//	int typeInt;
//	if (!PyArg_ParseTuple(args, "iO", &typeInt, &elements)) {
//		return nullptr;
//	}
//
//	size_t dimensions = PySequence_Length(elements);
//	if (!(0 < dimensions && dimensions <= kMaxNumberOfDimensions && 0 <= typeInt
//			&& typeInt < LIBSAKURA_SYMBOL(PyTypeId_kEnd))) {
//		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
//		return nullptr;
//	}
//	LIBSAKURA_SYMBOL(PyTypeId) type =
//			static_cast<LIBSAKURA_SYMBOL(PyTypeId)>(typeInt);
//
//	size_t elems[dimensions];
//	size_t total_elements = 1;
//	for (size_t i = 0; i < dimensions; ++i) {
//		RefHolder item(PySequence_GetItem(elements, i));
//		auto n = PyNumber_AsSsize_t(item.get(), PyExc_OverflowError);
//		if (PyErr_Occurred()) {
//			return nullptr;
//		}
//		elems[i] = n;
//		total_elements *= n;
//	}
//
//	static size_t const sizes[] = { sizeof(bool), sizeof(int8_t),
//			sizeof(int32_t), sizeof(int64_t), sizeof(uint8_t), sizeof(uint32_t),
//			sizeof(float), sizeof(double),
//			sizeof(long double) };
//	STATIC_ASSERT(LIBSAKURA_SYMBOL(PyTypeId_kEnd) == ELEMENTSOF((sizes)));
//
//	std::unique_ptr<LIBSAKURA_SYMBOL(PyAlignedBuffer),
//			decltype(&LIBSAKURA_SYMBOL(PyAlignedBufferDestroy))> bufPtr(nullptr,
//			LIBSAKURA_SYMBOL(PyAlignedBufferDestroy));
//	try {
//		char *aligned = nullptr;
//		auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<char>(
//				sizes[type] * total_elements, &aligned);
//		LIBSAKURA_SYMBOL(PyAlignedBuffer) *buf = nullptr;
//		LIBSAKURA_SYMBOL(Status) status =
//		LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(type, addr, aligned, dimensions,
//				elems, LIBSAKURA_PREFIX::Memory::Free, &buf);
//		if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
//			throw std::bad_alloc();
//		}
//		bufPtr.reset(buf);
//		assert(status == LIBSAKURA_SYMBOL(Status_kOK));
//	} catch (std::bad_alloc const&e) {
//		PyErr_SetString(PyExc_MemoryError, "No memory.");
//		return nullptr;
//	}
//
//	PyObject *capsule = nullptr;
//	LIBSAKURA_SYMBOL(Status) status =
//	LIBSAKURA_SYMBOL(PyAlignedBufferEncapsulate)(bufPtr.get(), &capsule);
//	if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
//		PyErr_SetString(PyExc_MemoryError, "No memory.");
//		return nullptr;
//	}
//	assert(status == LIBSAKURA_SYMBOL(Status_kOK) && capsule != nullptr);
//	bufPtr.release();
//	return capsule;
//}

//PyObject *NewAlignedBuffer(PyObject *self, PyObject *args) {
//	PyObject *dataSeq;
//	int type;
//	if (!PyArg_ParseTuple(args, "iO", &type, &dataSeq)) {
//		return nullptr;
//	}
//
//	size_t elements[1];
//	auto &len = elements[0] = PySequence_Length(dataSeq);
//	LIBSAKURA_SYMBOL(PyAlignedBuffer) *buf = nullptr;
//	try {
//		switch (type) {
//		case LIBSAKURA_SYMBOL(PyTypeId_kBool): {
//			bool *aligned = nullptr;
//			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
//			bool>(sizeof(bool) * len, &aligned);
//			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
//				RefHolder item(PySequence_GetItem(dataSeq, i));
//				aligned[i] = item.get() == Py_True;
//			}
//			LIBSAKURA_SYMBOL(Status) status =
//			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
//					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
//					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
//			if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
//				throw std::bad_alloc();
//			}
//			assert(status == LIBSAKURA_SYMBOL(Status_kOK));
//		}
//			break;
//
//		case LIBSAKURA_SYMBOL(PyTypeId_kInt8): {
//			uint8_t *aligned = nullptr;
//			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
//					uint8_t>(sizeof(uint8_t) * len, &aligned);
//			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
//				RefHolder item(PySequence_GetItem(dataSeq, i));
//				auto val = ASLONG(item.get());
//				if (PyErr_Occurred()) {
//					return nullptr;
//				}
//				aligned[i] = val;
//			}
//			LIBSAKURA_SYMBOL(Status) status =
//			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
//					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
//					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
//			if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
//				throw std::bad_alloc();
//			}
//			assert(status == LIBSAKURA_SYMBOL(Status_kOK));
//		}
//			break;
//
//		case LIBSAKURA_SYMBOL(PyTypeId_kInt32): {
//			int32_t *aligned = nullptr;
//			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
//					int32_t>(sizeof(int32_t) * len, &aligned);
//			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
//				RefHolder item(PySequence_GetItem(dataSeq, i));
//				auto val = ASLONG(item.get());
//				if (PyErr_Occurred()) {
//					return nullptr;
//				}
//				aligned[i] = val;
//			}
//			LIBSAKURA_SYMBOL(Status) status =
//			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
//					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
//					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
//			if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
//				throw std::bad_alloc();
//			}
//			assert(status == LIBSAKURA_SYMBOL(Status_kOK));
//		}
//			break;
//
//		case LIBSAKURA_SYMBOL(PyTypeId_kFloat): {
//			float *aligned = nullptr;
//			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
//					float>(sizeof(float) * len, &aligned);
//			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
//				RefHolder item(PySequence_GetItem(dataSeq, i));
//				auto val = PyFloat_AsDouble(item.get());
//				if (PyErr_Occurred()) {
//					return nullptr;
//				}
//				aligned[i] = val;
//			}
//			LIBSAKURA_SYMBOL(Status) status =
//			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
//					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
//					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
//			if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
//				throw std::bad_alloc();
//			}
//			assert(status == LIBSAKURA_SYMBOL(Status_kOK));
//		}
//			break;
//
//		case LIBSAKURA_SYMBOL(PyTypeId_kDouble): {
//			double *aligned = nullptr;
//			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
//					double>(sizeof(double) * len, &aligned);
//			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
//				RefHolder item(PySequence_GetItem(dataSeq, i));
//				auto val = PyFloat_AsDouble(item.get());
//				if (PyErr_Occurred()) {
//					return nullptr;
//				}
//				aligned[i] = val;
//			}
//			LIBSAKURA_SYMBOL(Status) status =
//			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
//					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
//					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
//			if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
//				throw std::bad_alloc();
//			}
//			assert(status == LIBSAKURA_SYMBOL(Status_kOK));
//		}
//			break;
//		default:
//			assert(false);
//			break;
//		}
//	} catch (std::bad_alloc const&e) {
//		PyErr_SetString(PyExc_MemoryError, "No memory.");
//		return nullptr;
//	}
//
//	PyObject *capsule = nullptr;
//	LIBSAKURA_SYMBOL(Status) status =
//	LIBSAKURA_SYMBOL(PyAlignedBufferEncapsulate)(buf, &capsule);
//	if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
//		PyErr_SetString(PyExc_MemoryError, "No memory.");
//		return nullptr;
//	}
//	assert(status == LIBSAKURA_SYMBOL(Status_kOK));
//	return capsule;
//}

PyMethodDef module_methods[] =
		{

				{ "initialize", Initialize, METH_VARARGS, "Initializes libsakura." },

				{ "clean_up", CleanUp, METH_VARARGS, "Cleans up libsakura." },

//		{ "compute_statistics", ComputeStatistics, METH_VARARGS,
//				"Computes statistics of unmasked elements." },
//
//		{ "compute_accurate_statistics", ComputeAccurateStatistics,
//		METH_VARARGS, "Computes accurate statistics of unmasked elements." },
//
				{ "compute_statistics", ComputeStatisticsNP, METH_VARARGS,
						"Computes statistics of unmasked elements." },

				{ "compute_accurate_statistics", ComputeAccurateStatisticsNP,
						METH_VARARGS, "Computes accurate statistics of unmasked elements." },

//		{ "compute_mad", ComputeMAD, METH_VARARGS,
//				"Computes median absolute deviation for sorted input array." },
//
				{ "compute_mad", ComputeMADNP, METH_VARARGS,
						"Computes median absolute deviation for sorted input array." },

//		{ "sort_data_densely", SortValidDataDenselyFloat, METH_VARARGS,
//				"Sorts valid data and return number of valid data." },
//
				{ "sort_data_densely", SortValidDataDenselyFloatNP, METH_VARARGS,
						"Sorts valid data and return number of valid data." },

				{ "grid_convolving", GridConvolvingNP, METH_VARARGS,
						"Grids spectra on X-Y plane with convolving." },

//				{ "uint8_to_bool", Uint8ToBool, METH_VARARGS,
//						"Converts uint8 to bool." },
//
//				{ "uint32_to_bool", Uint32ToBool, METH_VARARGS,
//						"Converts uint32 to bool." },
//
//				{ "invert_bool", InvertBool, METH_VARARGS, "Inverts bool." },

				{ "uint8_to_bool", Uint8ToBoolNP, METH_VARARGS,
						"Converts uint8 to bool." },

				{ "uint32_to_bool", Uint32ToBoolNP, METH_VARARGS,
						"Converts uint32 to bool." },

				{ "invert_bool", InvertBoolNP, METH_VARARGS, "Inverts bool." },

//				{ "logical_and", LogicalAnd, METH_VARARGS,
//						"Takes logical conjunction of boolean arrays." },
//
				{ "logical_and", LogicalAndNP, METH_VARARGS,
						"Takes logical conjunction of boolean arrays." },

//				{ "operate_bits_uint8_or", Uint8OperateBitsOr, METH_VARARGS,
//						"Bit operation OR between an uint8 value and uint8 array." },
//
//				{ "operate_bits_uint32_or", Uint32OperateBitsOr, METH_VARARGS,
//						"Bit operation OR between an uint32 value and uint32 array." },
//
//				{ "operate_bits_uint8_and", Uint8OperateBitsAnd, METH_VARARGS,
//						"Bit operation AND between an uint8 value and uint8 array." },
//
//				{ "operate_bits_uint32_and", Uint32OperateBitsAnd, METH_VARARGS,
//						"Bit operation AND between an uint32 value and uint32 array." },
//
//				{ "operate_bits_uint8_xor", Uint8OperateBitsXor, METH_VARARGS,
//						"Bit operation XOR between an uint8 value and uint8 array." },
//
//				{ "operate_bits_uint32_xor", Uint32OperateBitsXor, METH_VARARGS,
//						"Bit operation XOR between an uint32 value and uint32 array." },

				{ "operate_bits_uint8_or", Uint8OperateBitsOrNP, METH_VARARGS,
						"Bit operation OR between an uint8 value and uint8 array." },

				{ "operate_bits_uint32_or", Uint32OperateBitsOrNP, METH_VARARGS,
						"Bit operation OR between an uint32 value and uint32 array." },

				{ "operate_bits_uint8_and", Uint8OperateBitsAndNP, METH_VARARGS,
						"Bit operation AND between an uint8 value and uint8 array." },

				{ "operate_bits_uint32_and", Uint32OperateBitsAndNP, METH_VARARGS,
						"Bit operation AND between an uint32 value and uint32 array." },

				{ "operate_bits_uint8_xor", Uint8OperateBitsXorNP, METH_VARARGS,
						"Bit operation XOR between an uint8 value and uint8 array." },

				{ "operate_bits_uint32_xor", Uint32OperateBitsXorNP, METH_VARARGS,
						"Bit operation XOR between an uint32 value and uint32 array." },

//				{ "operate_bits_uint8_xor", Uint8OperateBitsNot, METH_VARARGS,
//						"Bit operation NOT of uint8 array." },
//
//				{ "operate_bits_uint32_xor", Uint32OperateBitsNot, METH_VARARGS,
//						"Bit operation NOT of uint32 array." },
//
				{ "operate_bits_uint8_not", Uint8OperateBitsNotNP, METH_VARARGS,
						"Bit operation NOT of uint8 array." },

				{ "operate_bits_uint32_not", Uint32OperateBitsNotNP, METH_VARARGS,
						"Bit operation NOT of uint32 array." },

//				{ "set_true_float_in_ranges_exclusive",
//						FloatSetTrueIntInRangesExclusive, METH_VARARGS,
//						"Sets True if the element is in at least one of ranges." },
//
//				{ "set_true_int_in_ranges_exclusive",
//						Int32SetTrueIntInRangesExclusive, METH_VARARGS,
//						"Sets True if the element is in at least one of ranges." },
//
				{ "set_true_float_in_ranges_exclusive",
						FloatSetTrueIntInRangesExclusiveNP, METH_VARARGS,
						"Sets True if the element is in at least one of ranges." },

				{ "set_true_int_in_ranges_exclusive",
						Int32SetTrueIntInRangesExclusiveNP, METH_VARARGS,
						"Sets True if the element is in at least one of ranges." },

//				{ "set_false_float_if_nan_or_inf", SetFalseFloatIfNanOrInf,
//				METH_VARARGS, "set false if float value is NaN or Inf." },
//
				{ "set_false_float_if_nan_or_inf", SetFalseFloatIfNanOrInfNP,
				METH_VARARGS, "set false if float value is NaN or Inf." },

//				{ "set_true_float_if_greater_than", FloatSetTrueIfGreaterThan,
//				METH_VARARGS,
//						"set true if float value is greater than threshold." },
//
//				{ "set_true_int_if_greater_than", IntSetTrueIfGreaterThan,
//				METH_VARARGS, "set true if int value is greater than threshold." },
//
//				{ "set_true_float_if_greater_than_or_equal",
//						FloatSetTrueIfGreaterThanOrEquals,
//						METH_VARARGS,
//						"set true if float value is greater than or equal to threshold." },
//
//				{ "set_true_int_if_greater_than_or_equal",
//						IntSetTrueIfGreaterThanOrEquals,
//						METH_VARARGS,
//						"set true if int value is greater than or equal to threshold." },
//
//				{ "set_true_float_if_less_than", FloatSetTrueIfLessThan,
//				METH_VARARGS, "set true if float value is less than threshold." },
//
//				{ "set_true_int_if_less_than", IntSetTrueIfLessThan,
//				METH_VARARGS, "set true if int value is less than threshold." },
//
//				{ "set_true_float_if_less_than_or_equal",
//						FloatSetTrueIfLessThanOrEquals,
//						METH_VARARGS,
//						"set true if float value is less than or equal to threshold." },
//
//				{ "set_true_int_if_less_than_or_equal",
//						IntSetTrueIfLessThanOrEquals,
//						METH_VARARGS,
//						"set true if int value is less than or equal to threshold." },
//

				{ "set_true_float_if_greater_than", FloatSetTrueIfGreaterThanNP,
						METH_VARARGS,
						"set true if float value is greater than threshold." },

				{ "set_true_int_if_greater_than", IntSetTrueIfGreaterThanNP,
						METH_VARARGS, "set true if int value is greater than threshold." },

				{ "set_true_float_if_greater_than_or_equal",
						FloatSetTrueIfGreaterThanOrEqualsNP,
						METH_VARARGS,
						"set true if float value is greater than or equal to threshold." },

				{ "set_true_int_if_greater_than_or_equal",
						IntSetTrueIfGreaterThanOrEqualsNP,
						METH_VARARGS,
						"set true if int value is greater than or equal to threshold." },

				{ "set_true_float_if_less_than", FloatSetTrueIfLessThanNP,
						METH_VARARGS, "set true if float value is less than threshold." },

				{ "set_true_int_if_less_than", IntSetTrueIfLessThanNP,
						METH_VARARGS, "set true if int value is less than threshold." },

				{ "set_true_float_if_less_than_or_equal",
						FloatSetTrueIfLessThanOrEqualsNP,
						METH_VARARGS,
						"set true if float value is less than or equal to threshold." },

				{ "set_true_int_if_less_than_or_equal",
						IntSetTrueIfLessThanOrEqualsNP,
						METH_VARARGS,
						"set true if int value is less than or equal to threshold." },

//				{ "interpolate_float_yaxis", InterpolateFloatYAxis,
//				METH_VARARGS, "perform one-dimensional interpolation." },
//
//				{ "interpolate_float_xaxis", InterpolateFloatXAxis,
//				METH_VARARGS, "perform one-dimensional interpolation." },
//
				{ "interpolate_float_yaxis", InterpolateFloatYAxisNP,
						METH_VARARGS, "perform one-dimensional interpolation." },

				{ "interpolate_float_xaxis", InterpolateFloatXAxisNP,
						METH_VARARGS, "perform one-dimensional interpolation." },

//				{ "apply_position_switch_calibration",
//						CalibrateDataWithArrayScalingFloat,
//						METH_VARARGS, "apply position switch calibration." },
//
				{ "apply_position_switch_calibration",
						CalibrateDataWithArrayScalingFloatNP,
						METH_VARARGS, "apply position switch calibration." },

//				{ "create_convolve1d_fft_context", CreateConvolve1DContextFFT,
//				METH_VARARGS, "Creates a context for convolving 1D." },
//
				{ "create_convolve1d_fft_context", CreateConvolve1DContextFFTNP,
						METH_VARARGS, "Creates a context for convolving 1D." },

//				{ "convolve1d_fft", Convolve1DFFT,
//				METH_VARARGS,
//						"perform FFT-based one-dimensional discrete convolution." },

				{ "convolve1d_fft", Convolve1DFFTNP,
						METH_VARARGS,
						"perform FFT-based one-dimensional discrete convolution." },

//				{ "convolve1d", Convolve1D,
//				METH_VARARGS, "perform one-dimensional discrete convolution." },
//
				{ "convolve1d", Convolve1DNP,
						METH_VARARGS, "perform one-dimensional discrete convolution." },

				{ "create_baseline_context", CreateLSQFitContextPolynomial,
						METH_VARARGS, "Creates a context for baseline subtraction." },

//				{ "lsqfit_polynomial", LSQFitPolynomial,
//				METH_VARARGS, "perform baseline subtraction." },

				{ "lsqfit_polynomial", LSQFitPolynomialNP,
						METH_VARARGS, "perform baseline subtraction." },

//				{ "complement_masked_value_float", ComplementMaskedValueFloat,
//				METH_VARARGS, "complement masked value with, tentatively, 0." },

				{ "complement_masked_value_float", ComplementMaskedValueFloatNP,
						METH_VARARGS, "complement masked value with, tentatively, 0." },

//				{ "create_gaussian_kernel_float", CreateGaussianKernelFloat,
//				METH_VARARGS, "create gaussian kernel" },

				{ "create_gaussian_kernel_float", CreateGaussianKernelFloatNP,
						METH_VARARGS, "create gaussian kernel" },

//				{ "get_elements_of_aligned_buffer", GetElementsOfAlignedBuffer,
//				METH_VARARGS, "gets_elements of the aligned buffer." },
//
//				{ "new_uninitialized_aligned_buffer",
//						NewUninitializedAlignedBuffer,
//						METH_VARARGS,
//						"Creates an uninitialized new aligned buffer with supplied elements." },
//
//				{ "new_aligned_buffer", NewAlignedBuffer, METH_VARARGS,
//						"Creates a new aligned buffer." },
//
				{ "new_uninitialized_aligned_ndarray", NewUninitializedAlignedNumPyArray,
						METH_VARARGS, "Creates a new aligned numpy ndarray."},

				{ NULL, NULL, 0, NULL } /* Sentinel */
		};

PyDoc_STRVAR(module_doc, "Python binding of libsakura library.");

}

extern "C" {

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferEncapsulate)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, PyObject **capsule) {
	if (capsule == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (buffer == nullptr) {
		*capsule = nullptr;
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	*capsule = PyCapsule_New(buffer, kAlignedBufferName,
			TCapsuleDestructorForAlignedBuffer);
	if (*capsule == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferDecapsulate)(
		PyObject *capsule,
		LIBSAKURA_SYMBOL(PyAlignedBuffer) **buffer) {
	if (buffer == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (capsule == nullptr || !PyCapsule_IsValid(capsule, kAlignedBufferName)) {
		*buffer = nullptr;
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	*buffer =
			reinterpret_cast<LIBSAKURA_SYMBOL(PyAlignedBuffer) *>(PyCapsule_GetPointer(
					capsule, kAlignedBufferName));
	assert(*buffer);
	return LIBSAKURA_SYMBOL(Status_kOK);
}

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
LIBSAKURA_SYMBOL(PyTypeId) type, void *original_addr, void *aligned_addr,
		size_t dimensions, size_t elements[], void (*destructor)(void *),
		LIBSAKURA_SYMBOL(PyAlignedBuffer) **buffer) {
	if (original_addr == nullptr || aligned_addr == nullptr
			|| elements == nullptr || buffer == nullptr
			|| dimensions > kMaxNumberOfDimensions) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	LIBSAKURA_SYMBOL(PyAlignedBuffer) *buf =
			reinterpret_cast<LIBSAKURA_SYMBOL(PyAlignedBuffer) *>(malloc(
					sizeof(LIBSAKURA_SYMBOL(PyAlignedBuffer))));
	if (buf == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	}
	//printf("%p: alloc\n", original_addr);
	buf->type = type;
	buf->original_addr = original_addr;
	buf->aligned_addr = aligned_addr;
	buf->destructor = destructor;
	buf->dimensions = dimensions;
	memcpy(buf->elements, elements, sizeof(buf->elements[0]) * dimensions);
	*buffer = buf;
	//printf("%p: alloc\n", buf);
	return LIBSAKURA_SYMBOL(Status_kOK);
}

void LIBSAKURA_SYMBOL(PyAlignedBufferDestroy)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer) {
	if (buffer) {
		if (buffer->destructor) {
			//printf("%p: free\n", buffer->original_addr);
			buffer->destructor(buffer->original_addr);
		}
		//printf("%p: free\n", buffer);
		free(buffer);
	}
}

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferDimensions)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, size_t *dimensions) {
	if (buffer == nullptr || dimensions == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	*dimensions = buffer->dimensions;
	return LIBSAKURA_SYMBOL(Status_kOK);
}

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferElements)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, size_t dimensions,
		size_t elements[]) {
	if (buffer == nullptr || elements == nullptr
			|| dimensions > buffer->dimensions) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	memcpy(elements, buffer->elements, sizeof(elements[0]) * dimensions);
	return LIBSAKURA_SYMBOL(Status_kOK);
}

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferAlignedAddr)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, void **aligned_addr) {
	if (buffer == nullptr || aligned_addr == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	*aligned_addr = buffer->aligned_addr;
	return LIBSAKURA_SYMBOL(Status_kOK);
}

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferType)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, LIBSAKURA_SYMBOL(PyTypeId) *type) {
	if (buffer == nullptr || type == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	*type = buffer->type;
	return LIBSAKURA_SYMBOL(Status_kOK);
}

#if PY_MAJOR_VERSION >= 3

static int module_traverse(PyObject *m, visitproc visit, void *arg) {
	Py_VISIT(GETSTATE(m)->error);
	return 0;
}

static int module_clear(PyObject *m) {
	Py_CLEAR(GETSTATE(m)->error);
	return 0;
}

static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	MODULE_NAME,
	module_doc,
	sizeof(struct module_state),
	module_methods,
	NULL,
	module_traverse,
	module_clear,
	NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_libsakurapy(void)

#else

#define INITERROR return

PyMODINIT_FUNC initlibsakurapy(void)
#endif
		{
#if PY_MAJOR_VERSION >= 3
	PyObject *mod = PyModule_Create(&moduledef);
#else
	PyObject *mod = Py_InitModule3(MODULE_NAME, module_methods, module_doc);
#endif
	if (mod == nullptr) {
		INITERROR;
	}

	// to use NumPy C-API
	import_array();

	struct module_state *st = GETSTATE(mod);

	static char excep_name[] = MODULE_NAME ".error";
	static char excep_doc[] = "error on invoking libsakura functions";
	auto sakura_error = st->error;
	sakura_error = PyErr_NewExceptionWithDoc(excep_name, excep_doc, nullptr,
			nullptr);
//	if (sakura_error != nullptr) {
//		Py_INCREF(sakura_error);
//		PyModule_AddObject(mod, "error", sakura_error);
//	}
	if (sakura_error == NULL) {
		DecrementRef(mod);
		INITERROR;
	}
	PyModule_AddIntConstant(mod, "TYPE_BOOL", LIBSAKURA_SYMBOL(PyTypeId_kBool));
	PyModule_AddIntConstant(mod, "TYPE_INT8", LIBSAKURA_SYMBOL(PyTypeId_kInt8));
	PyModule_AddIntConstant(mod, "TYPE_INT32",
			LIBSAKURA_SYMBOL(PyTypeId_kInt32));
	PyModule_AddIntConstant(mod, "TYPE_UINT8", LIBSAKURA_SYMBOL(PyTypeId_kUInt8));
	PyModule_AddIntConstant(mod, "TYPE_UINT32",
			LIBSAKURA_SYMBOL(PyTypeId_kUInt32));
	PyModule_AddIntConstant(mod, "TYPE_FLOAT",
			LIBSAKURA_SYMBOL(PyTypeId_kFloat));
	PyModule_AddIntConstant(mod, "TYPE_DOUBLE",
			LIBSAKURA_SYMBOL(PyTypeId_kDouble));

	// Enum for interpolation methods
	PyModule_AddIntConstant(mod, "INTERPOLATION_METHOD_NEAREST",
			LIBSAKURA_SYMBOL(InterpolationMethod_kNearest));
	PyModule_AddIntConstant(mod, "INTERPOLATION_METHOD_LINEAR",
			LIBSAKURA_SYMBOL(InterpolationMethod_kLinear));
	PyModule_AddIntConstant(mod, "INTERPOLATION_METHOD_POLYNOMIAL",
			LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial));
	PyModule_AddIntConstant(mod, "INTERPOLATION_METHOD_SPLINE",
			LIBSAKURA_SYMBOL(InterpolationMethod_kSpline));
	// Enum for baseline functions
	PyModule_AddIntConstant(mod, "BASELINE_TYPE_POLYNOMIAL",
			LIBSAKURA_SYMBOL(LSQFitType_kPolynomial));
	PyModule_AddIntConstant(mod, "BASELINE_TYPE_CHEBYSHEV",
			LIBSAKURA_SYMBOL(LSQFitType_kChebyshev));

#if PY_MAJOR_VERSION >= 3
	return mod;
#endif
}

}
