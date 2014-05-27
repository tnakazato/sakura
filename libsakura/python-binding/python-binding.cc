#include <Python.h>

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
 */

#define MODULE_NAME "libsakurapy"

namespace {

constexpr size_t kMaxNumberOfDimensions = 5;

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

PyObject *sakura_error;

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
void TCapsuleDesctructor(PyObject *obj) {
	if (!PyCapsule_IsValid(obj, Key)) {
		return;
	}
	auto ptr = reinterpret_cast<T *>(PyCapsule_GetPointer(obj, Key));
	Func(ptr);
}

void TCapsuleDesctructorForAlignedBuffer(PyObject *obj) {
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
		PyErr_SetString(sakura_error, "Failed to initialize libsakura.");
		return nullptr;
	}
	Py_RETURN_NONE; // equivalent to Py_INCREF(Py_None); return Py_None;
}

PyObject *CleanUp(PyObject *self, PyObject *args) {
	LIBSAKURA_SYMBOL(CleanUp)();
	Py_RETURN_NONE;
}

PyObject *GetCurrentTime(PyObject *self, PyObject *args) {
	auto now = LIBSAKURA_SYMBOL(GetCurrentTime)();
	return Py_BuildValue("d", now); // Py_BuildValue builds a value with reference count 1 by default.
}

struct AlignedBufferConfiguration {
	LIBSAKURA_SYMBOL(PyTypeId) type;
	size_t dimensions;
	size_t elements[kMaxNumberOfDimensions];
};

bool isValidAlignedBuffer(size_t num, AlignedBufferConfiguration const conf[],
		PyObject *capsules[],
		LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[]) {
	for (size_t i = 0; i < num; ++i) {
		LIBSAKURA_SYMBOL(PyAlignedBuffer) *buf = nullptr;
		if (LIBSAKURA_SYMBOL(PyAlignedBufferDecapsulate)(capsules[i],
				&buf) != LIBSAKURA_SYMBOL(Status_kOK)) {
			return false;
		}
		assert(buf);
		if (buf->type != conf[i].type) {
			return false;
		}
		if (buf->dimensions != conf[i].dimensions) {
			return false;
		}
		if (memcmp(buf->elements, conf[i].elements,
				sizeof(buf->elements[0]) * buf->dimensions) != 0) {
			return false;
		}
		bufs[i] = buf;
	}
	return true;
}

PyObject *ComputeStatistics(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kData, kIsValid
	};
	PyObject *capsules[2];

	if (!PyArg_ParseTuple(args, "nOO", &num_data_py, &capsules[kData],
			&capsules[kIsValid])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	static AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(TypeId_kFloat), 1, { num_data } },

	{ LIBSAKURA_SYMBOL(TypeId_kBool), 1, { num_data } },

	};

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[2];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(StatisticsResult) result;
	LIBSAKURA_SYMBOL(Status) status;
	Py_BEGIN_ALLOW_THREADS
		status = LIBSAKURA_SYMBOL(ComputeStatistics)(num_data,
				reinterpret_cast<float const*>(bufs[kData]->aligned_addr),
				reinterpret_cast<bool const*>(bufs[kIsValid]->aligned_addr),
				&result);
		Py_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	return Py_BuildValue("{s:n,s:f,s:f,s:f,s:f,s:f,s:f,s:i,s:i}", "count",
			static_cast<Py_ssize_t>(result.count), "sum", result.sum, "mean",
			result.mean, "rms", result.rms, "stddev", result.stddev, "min",
			result.min, "max", result.max, "index_of_min", result.index_of_min,
			"index_of_max", result.index_of_max);

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

PyObject *Uint8ToBool(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kData, kResult
	};
	PyObject *capsules[2];

	if (!PyArg_ParseTuple(args, "nOO", &num_data_py, &capsules[kData],
			&capsules[kResult])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	static AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(TypeId_kInt8), 1, { num_data } },

	{ LIBSAKURA_SYMBOL(TypeId_kBool), 1, { num_data } },

	};

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[2];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	Py_BEGIN_ALLOW_THREADS
		status = LIBSAKURA_SYMBOL(Uint8ToBool)(num_data,
				reinterpret_cast<uint8_t const*>(bufs[kData]->aligned_addr),
				reinterpret_cast<bool *>(bufs[kResult]->aligned_addr));
		Py_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(capsules[kResult]);
	return capsules[kResult];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

PyObject *CreateConvolve1DContext(PyObject *self, PyObject *args) {
	Py_ssize_t num_data;
	int kernel_type;
	Py_ssize_t kernel_width;
	PyObject *use_fft;
	if (!PyArg_ParseTuple(args, "ninO", &num_data, &kernel_type, &kernel_width,
			&use_fft)) {
		return nullptr;
	}
	if (!((0 <= kernel_type
			&& kernel_type < LIBSAKURA_SYMBOL(Convolve1DKernelType_kNumType))
			&& (use_fft == Py_True || use_fft == Py_False))) {
		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
		return nullptr;
	}

	LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) status;
	Py_BEGIN_ALLOW_THREADS
		status =
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
						static_cast<size_t>(num_data),
						static_cast<LIBSAKURA_SYMBOL(Convolve1DKernelType)>(kernel_type),
						static_cast<size_t>(kernel_width), use_fft == Py_True,
						&context);
		Py_END_ALLOW_THREADS

	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError,
				"sakura_CreateConvolve1DContext failed.");
		return nullptr;
	}
	PyCapsule_Destructor destructor = TCapsuleDesctructor<
	LIBSAKURA_SYMBOL(Convolve1DContext), kConvolve1DContextName,
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)>;
	auto capsule = PyCapsule_New(context, kConvolve1DContextName, destructor);
	if (capsule == nullptr) {
		PyErr_SetString(PyExc_MemoryError, "No memory.");
		return nullptr;
	}
	return capsule;
}

PyObject *GetElementsOfAlignedBuffer(PyObject *self, PyObject *args) {
	PyObject *capsule = nullptr;

	if (!PyArg_ParseTuple(args, "O", &capsule)) {
		return nullptr;
	}
	LIBSAKURA_SYMBOL(PyAlignedBuffer) *buf = nullptr;
	if (LIBSAKURA_SYMBOL(PyAlignedBufferDecapsulate)(capsule,
			&buf) != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
		return nullptr;
	}
	assert(buf);

	RefHolder result(PyTuple_New(buf->dimensions));
	if (result.get() == nullptr) {
		goto out_of_memory;
	}
	for (size_t i = 0; i < buf->dimensions; ++i) {
		auto element = PyInt_FromSsize_t(buf->elements[i]);
		if (element == nullptr) {
			goto out_of_memory;
		}
		PyTuple_SetItem(result.get(), i, element);
	}
	return result.release();

	out_of_memory:

	PyErr_SetString(PyExc_MemoryError, "No memory.");
	return nullptr;
}

PyObject *NewUninitializedAlignedBuffer(PyObject *self, PyObject *args) {
	PyObject *elements;
	int typeInt;
	if (!PyArg_ParseTuple(args, "iO", &typeInt, &elements)) {
		return nullptr;
	}

	size_t dimensions = PySequence_Length(elements);
	if (!(0 < dimensions && dimensions <= kMaxNumberOfDimensions && 0 <= typeInt
			&& typeInt < LIBSAKURA_SYMBOL(TypeId_kEnd))) {
		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
		return nullptr;
	}
	LIBSAKURA_SYMBOL(PyTypeId) type =
			static_cast<LIBSAKURA_SYMBOL(PyTypeId)>(typeInt);

	size_t elems[dimensions];
	size_t total_elements = 1;
	for (size_t i = 0; i < dimensions; ++i) {
		RefHolder item(PySequence_GetItem(elements, i));
		auto n = PyNumber_AsSsize_t(item.get(), PyExc_OverflowError);
		if (PyErr_Occurred()) {
			return nullptr;
		}
		elems[i] = n;
		total_elements *= n;
	}

	static size_t const sizes[] = { sizeof(bool), sizeof(int8_t),
			sizeof(int32_t), sizeof(int64_t), sizeof(float), sizeof(double),
			sizeof(long double) };
	STATIC_ASSERT(LIBSAKURA_SYMBOL(TypeId_kEnd) == ELEMENTSOF((sizes)));

	std::unique_ptr<LIBSAKURA_SYMBOL(PyAlignedBuffer),
			decltype(&LIBSAKURA_SYMBOL(PyAlignedBufferDestroy))> bufPtr(nullptr,
	LIBSAKURA_SYMBOL(PyAlignedBufferDestroy));
	try {
		char *aligned = nullptr;
		auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<char>(
				sizes[type] * total_elements, &aligned);
		LIBSAKURA_SYMBOL(PyAlignedBuffer) *buf = nullptr;
		LIBSAKURA_SYMBOL(Status) status =
				LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(type, addr, aligned,
						dimensions, elems, LIBSAKURA_PREFIX::Memory::Free,
						&buf);
		if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
			throw std::bad_alloc();
		}
		bufPtr.reset(buf);
		assert(status == LIBSAKURA_SYMBOL(Status_kOK));
	} catch (std::bad_alloc const&e) {
		PyErr_SetString(PyExc_MemoryError, "No memory.");
		return nullptr;
	}

	PyObject *capsule = nullptr;
	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(PyAlignedBufferEncapsulate)(bufPtr.get(), &capsule);
	if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
		PyErr_SetString(PyExc_MemoryError, "No memory.");
		return nullptr;
	}
	assert(status == LIBSAKURA_SYMBOL(Status_kOK) && capsule != nullptr);
	bufPtr.release();
	return capsule;
}

PyObject *NewAlignedBuffer(PyObject *self, PyObject *args) {
	PyObject *dataSeq;
	int type;
	if (!PyArg_ParseTuple(args, "iO", &type, &dataSeq)) {
		return nullptr;
	}

	size_t elements[1];
	auto &len = elements[0] = PySequence_Length(dataSeq);
	LIBSAKURA_SYMBOL(PyAlignedBuffer) *buf = nullptr;
	try {
		switch (type) {
		case LIBSAKURA_SYMBOL(TypeId_kBool): {
			bool *aligned = nullptr;
			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
			bool>(sizeof(bool) * len, &aligned);
			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
				RefHolder item(PySequence_GetItem(dataSeq, i));
				RefHolder itemInt(PyNumber_Int(item.get()));
				auto val = PyInt_AsLong(itemInt.get());
				aligned[i] = val != 0;
			}
			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
		}
			break;

		case LIBSAKURA_SYMBOL(TypeId_kInt8): {
			uint8_t *aligned = nullptr;
			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
					uint8_t>(sizeof(uint8_t) * len, &aligned);
			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
				RefHolder item(PySequence_GetItem(dataSeq, i));
				RefHolder itemInt(PyNumber_Int(item.get()));
				auto val = PyInt_AsLong(itemInt.get());
				aligned[i] = val;
			}
			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
		}
			break;

		case LIBSAKURA_SYMBOL(TypeId_kFloat): {
			float *aligned = nullptr;
			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
					float>(sizeof(float) * len, &aligned);
			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
				RefHolder item(PySequence_GetItem(dataSeq, i));
				RefHolder itemFloat(PyNumber_Float(item.get()));
				auto val = PyFloat_AsDouble(itemFloat.get());
				aligned[i] = val;
			}
			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
		}
			break;
		default:
			break;
		}
	} catch (std::bad_alloc const&e) {
		PyErr_SetString(PyExc_MemoryError, "No memory.");
		return nullptr;
	}

	PyObject *capsule = nullptr;
	LIBSAKURA_SYMBOL(PyAlignedBufferEncapsulate)(buf, &capsule);
	return capsule;
}

PyMethodDef module_methods[] = {

{ "initialize", Initialize, METH_VARARGS, "Initializes libsakura." },

{ "clean_up", CleanUp, METH_VARARGS, "Cleans up libsakura." },

{ "get_current_time", GetCurrentTime, METH_VARARGS,
		"Gets current time in seconds of type double." },

{ "compute_statistics", ComputeStatistics, METH_VARARGS,
		"Computes statistics of unmasked elements." },

{ "uint8_to_bool", Uint8ToBool, METH_VARARGS,
		"Converts uint8 to bool." },

{ "create_convolve1D_context", CreateConvolve1DContext, METH_VARARGS,
		"Creates a context for convolving 1D." },

{ "get_elements_of_aligned_buffer", GetElementsOfAlignedBuffer, METH_VARARGS,
		"gets_elements of the aligned buffer." },

{ "new_uninitialized_aligned_buffer", NewUninitializedAlignedBuffer,
		METH_VARARGS, "Creates an uninitialized new aligned buffer with supplied elements." },

{ "new_aligned_buffer", NewAlignedBuffer, METH_VARARGS,
		"Creates a new aligned buffer." },

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
			TCapsuleDesctructorForAlignedBuffer);
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
	buf->type = type;
	buf->original_addr = original_addr;
	buf->aligned_addr = aligned_addr;
	buf->destructor = destructor;
	buf->dimensions = dimensions;
	memcpy(buf->elements, elements, sizeof(buf->elements[0]) * dimensions);
	*buffer = buf;
	return LIBSAKURA_SYMBOL(Status_kOK);
}

void LIBSAKURA_SYMBOL(PyAlignedBufferDestroy)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer) {
	if (buffer) {
		if (buffer->destructor) {
			buffer->destructor(buffer->original_addr);
		}
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

PyMODINIT_FUNC initlibsakurapy(void) {
	PyObject *mod = Py_InitModule3(MODULE_NAME, module_methods, module_doc);
	if (mod == nullptr) {
		return;
	}

	static char excep_name[] = MODULE_NAME ".error";
	static char excep_doc[] = "error on invoking libsakura functions";
	sakura_error = PyErr_NewExceptionWithDoc(excep_name, excep_doc, nullptr,
			nullptr);
	if (sakura_error != nullptr) {
		Py_INCREF(sakura_error);
		PyModule_AddObject(mod, "error", sakura_error);
	}
	PyModule_AddIntConstant(mod, "TYPE_BOOL", LIBSAKURA_SYMBOL(TypeId_kBool));
	PyModule_AddIntConstant(mod, "TYPE_INT8", LIBSAKURA_SYMBOL(TypeId_kInt8));
	PyModule_AddIntConstant(mod, "TYPE_INT32", LIBSAKURA_SYMBOL(TypeId_kInt32));
	PyModule_AddIntConstant(mod, "TYPE_FLOAT", LIBSAKURA_SYMBOL(TypeId_kFloat));
	PyModule_AddIntConstant(mod, "TYPE_DOUBLE",
	LIBSAKURA_SYMBOL(TypeId_kDouble));

	PyModule_AddIntConstant(mod, "CONVOLVE1D_KERNEL_TYPE_GAUSSIAN",
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian));
	PyModule_AddIntConstant(mod, "CONVOLVE1D_KERNEL_TYPE_BOXCAR",
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kBoxcar));
	PyModule_AddIntConstant(mod, "CONVOLVE1D_KERNEL_TYPE_HANNING",
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kHanning));
	PyModule_AddIntConstant(mod, "CONVOLVE1D_KERNEL_TYPE_HAMMING",
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kHamming));
}

}
