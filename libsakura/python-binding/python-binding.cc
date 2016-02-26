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
		if (!(buf->type == conf[i].type && conf[i].pred(*buf))) {
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
	auto pred = [num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
	{	return TotalElementsGreaterOrEqual(buf, num_data);};
	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), pred },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), pred },

	};

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[2];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = LIBSAKURA_SYMBOL(ComputeStatisticsFloat)(num_data,
				reinterpret_cast<float const*>(bufs[kData]->aligned_addr),
				reinterpret_cast<bool const*>(bufs[kIsValid]->aligned_addr),
				&result);
		SAKURA_END_ALLOW_THREADS
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

PyObject *GridConvolving(PyObject *self, PyObject *args) {
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
	PyObject *capsules[kEnd];

	if (!PyArg_ParseTuple(args, "nnnOOOnnnOnOOOOOnOnnnnOOO", &num_spectra_py,
			&start_spectrum_py, &end_spectrum_py, &capsules[kSpectrumMask],
			&capsules[kX], &capsules[kY], &support_py, &sampling_py,
			&num_polarizations_py, &capsules[kPolarizationMap],
			&num_channels_py, &capsules[kChannelMap], &capsules[kMask],
			&capsules[kValue], &capsules[kWeight], &weight_only,
			&num_convolution_table_py, &capsules[kConvolutionTable],
			&num_polarizations_for_grid_py, &num_channels_for_grid_py,
			&width_py, &height_py, &capsules[kWeightSum],
			&capsules[kWeightOfGrid], &capsules[kGrid])) {
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

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[kEnd];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status =
				LIBSAKURA_SYMBOL(GridConvolvingFloat)(num_spectra, start_spectrum,
						end_spectrum,
						static_cast<bool const *>(bufs[kSpectrumMask]->aligned_addr),
						static_cast<double const *>(bufs[kX]->aligned_addr),
						static_cast<double const *>(bufs[kY]->aligned_addr),
						support, sampling, num_polarizations,
						static_cast<uint32_t const *>(bufs[kPolarizationMap]->aligned_addr),
						num_channels,
						static_cast<uint32_t const *>(bufs[kChannelMap]->aligned_addr),
						static_cast<bool const *>(bufs[kMask]->aligned_addr),
						static_cast<float const *>(bufs[kValue]->aligned_addr),
						static_cast<float const *>(bufs[kWeight]->aligned_addr),
						weight_only == Py_True, num_convolution_table,
						static_cast<float const *>(bufs[kConvolutionTable]->aligned_addr),
						num_polarizations_for_grid, num_channels_for_grid,
						width, height,
						static_cast<double *>(bufs[kWeightSum]->aligned_addr),
						static_cast<float *>(bufs[kWeightOfGrid]->aligned_addr),
						static_cast<float *>(bufs[kGrid]->aligned_addr));
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
PyObject *ConvertArray(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kData, kResult, kEnd
	};
	PyObject *capsules[kEnd];

	if (!PyArg_ParseTuple(args, "nOO", &num_data_py, &capsules[kData],
			&capsules[kResult])) {
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

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[kEnd];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = Func(num_data,
				reinterpret_cast<FromType const*>(bufs[kData]->aligned_addr),
				reinterpret_cast<ToType *>(bufs[kResult]->aligned_addr));
		SAKURA_END_ALLOW_THREADS
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

constexpr FuncForPython Uint8ToBool = ConvertArray<uint8_t, bool,
LIBSAKURA_SYMBOL(PyTypeId_kInt8), LIBSAKURA_SYMBOL(PyTypeId_kBool),
LIBSAKURA_SYMBOL(Uint8ToBool)>;

constexpr FuncForPython Uint32ToBool = ConvertArray<uint32_t, bool,
LIBSAKURA_SYMBOL(PyTypeId_kInt32), LIBSAKURA_SYMBOL(PyTypeId_kBool),
LIBSAKURA_SYMBOL(Uint32ToBool)>;

constexpr FuncForPython InvertBool = ConvertArray<bool, bool,
LIBSAKURA_SYMBOL(PyTypeId_kBool), LIBSAKURA_SYMBOL(PyTypeId_kBool),
LIBSAKURA_SYMBOL(InvertBool)>;

constexpr FuncForPython SetFalseFloatIfNanOrInf = ConvertArray<float, bool,
LIBSAKURA_SYMBOL(PyTypeId_kFloat), LIBSAKURA_SYMBOL(PyTypeId_kBool),
LIBSAKURA_SYMBOL(SetFalseIfNanOrInfFloat)>;


inline void AlignedBoolAnd(size_t elements, bool const src1[], bool const src2[], bool dst[]) {
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

PyObject *LogicalAnd(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kData1, kData2, kResult, kEnd
	};
	PyObject *capsules[kEnd];

	if (!PyArg_ParseTuple(args, "nOOO", &num_data_py, &capsules[kData1],
			&capsules[kData2], &capsules[kResult])) {
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

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[kEnd];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
		goto invalid_arg;
	}
	AlignedBoolAnd(num_data, reinterpret_cast<bool const*>(bufs[kData1]->aligned_addr),
			reinterpret_cast<bool const*>(bufs[kData2]->aligned_addr),
			reinterpret_cast<bool *>(bufs[kResult]->aligned_addr));
	Py_INCREF(capsules[kResult]);
	return capsules[kResult];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

template<typename Type,
LIBSAKURA_SYMBOL(PyTypeId) TypeId,
LIBSAKURA_SYMBOL(Status) (*Func)(size_t, Type const *, size_t, Type const *,
		Type const *, bool *)>
PyObject *RangeCheck(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	Py_ssize_t num_range_py;
	enum {
		kData, kLower, kUpper, kMask, kEnd
	};
	PyObject *capsules[kEnd];

	if (!PyArg_ParseTuple(args, "nOnOOO", &num_data_py, &capsules[kData],
			&num_range_py, &capsules[kLower], &capsules[kUpper],
			&capsules[kMask])) {
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

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[kEnd];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = Func(num_data,
				reinterpret_cast<Type const*>(bufs[kData]->aligned_addr),
				num_range,
				reinterpret_cast<Type const*>(bufs[kLower]->aligned_addr),
				reinterpret_cast<Type const*>(bufs[kUpper]->aligned_addr),
				reinterpret_cast<bool *>(bufs[kMask]->aligned_addr));
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(capsules[kMask]);
	return capsules[kMask];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

constexpr FuncForPython Int32SetTrueIntInRangesExclusive = RangeCheck<int32_t,
LIBSAKURA_SYMBOL(PyTypeId_kInt32),
LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveInt)>;

constexpr FuncForPython FloatSetTrueIntInRangesExclusive = RangeCheck<float,
LIBSAKURA_SYMBOL(PyTypeId_kFloat),
LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveFloat)>;

template<typename Type,
LIBSAKURA_SYMBOL(PyTypeId) TypeId,
LIBSAKURA_SYMBOL(Status) (*Func)(Type, size_t, Type const *, bool const *, Type *)>
PyObject *BitOperation(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	unsigned int bitmask_py;

	enum {
		kData, kMask, kResult, kEnd
	};
	PyObject *capsules[kEnd];

	if (!PyArg_ParseTuple(args, "InOOO", &bitmask_py, &num_data_py, &capsules[kData],
			&capsules[kMask], &capsules[kResult])) {
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	auto bitmask = static_cast<Type>(bitmask_py);
	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};

	AlignedBufferConfiguration const conf[] = {

	{ TypeId, predData },
	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), predData },
	{ TypeId, predData }

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[kEnd];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = Func(bitmask, num_data,
				reinterpret_cast<Type const*>(bufs[kData]->aligned_addr),
				reinterpret_cast<bool const*>(bufs[kMask]->aligned_addr),
				reinterpret_cast<Type *>(bufs[kResult]->aligned_addr));
		SAKURA_END_ALLOW_THREADS
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

constexpr FuncForPython Uint8OperateBitsOr = BitOperation<uint8_t,
		LIBSAKURA_SYMBOL(PyTypeId_kInt8),
		LIBSAKURA_SYMBOL(OperateBitwiseOrUint8)>;


template<typename Type,
LIBSAKURA_SYMBOL(PyTypeId) TypeId,
LIBSAKURA_SYMBOL(Status) (*Func)(LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
				uint8_t polynomial_order, size_t num_base,
				double const base_position[/*num_base*/], size_t num_array,
				Type const base_data[/*num_base*num_array*/],
				bool const base_mask[/*num_base*num_array*/], size_t num_interpolated,
				double const interpolated_position[/*num_interpolated*/],
				float interpolated_data[/*num_interpolated*num_array*/],
				bool interpolated_mask[/*num_interpolated*num_array*/])>
PyObject *InterpolateAxis(PyObject *self, PyObject *args) {
	int interpolate_type_py;
	uint8_t polynomial_order;
	Py_ssize_t num_base_py;
	Py_ssize_t num_array_py;
	Py_ssize_t num_interpolated_py;
	enum {
		kFromAxis, kFromData, kMask, kToAxis, kToData, kToMask, kEnd
	};
	PyObject *capsules[kEnd];

	if (!PyArg_ParseTuple(args, "ibnOnOOnOOO", &interpolate_type_py, &polynomial_order,
			&num_base_py, &capsules[kFromAxis], &num_array_py, &capsules[kFromData], &capsules[kMask],
			&num_interpolated_py, &capsules[kToAxis], &capsules[kToData], &capsules[kMask])) {
		return nullptr;
	}
	if (!(0 <= interpolate_type_py
			&& interpolate_type_py < LIBSAKURA_SYMBOL(InterpolationMethod_kNumElements))) {
		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
		return nullptr;
	}
	auto interpolate_type = static_cast<LIBSAKURA_SYMBOL(InterpolationMethod)>(interpolate_type_py);
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

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[kEnd];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = Func(interpolate_type, polynomial_order, num_base,
				reinterpret_cast<double const*>(bufs[kFromAxis]->aligned_addr),
				num_array,
				reinterpret_cast<Type const*>(bufs[kFromData]->aligned_addr),
				reinterpret_cast<bool const*>(bufs[kMask]->aligned_addr),
				num_interpolated,
				reinterpret_cast<double const*>(bufs[kToAxis]->aligned_addr),
				reinterpret_cast<Type *>(bufs[kToData]->aligned_addr),
				reinterpret_cast<bool *>(bufs[kToMask]->aligned_addr));
		SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
		return nullptr;
	}
	Py_INCREF(capsules[kToData]);
	return capsules[kToData];

	invalid_arg:

	PyErr_SetString(PyExc_ValueError, "Invalid argument.");
	return nullptr;
}

constexpr FuncForPython InterpolateFloatYAxis = InterpolateAxis<float,
		LIBSAKURA_SYMBOL(PyTypeId_kFloat),
		LIBSAKURA_SYMBOL(InterpolateYAxisFloat)>;

constexpr FuncForPython InterpolateFloatXAxis = InterpolateAxis<float,
		LIBSAKURA_SYMBOL(PyTypeId_kFloat),
		LIBSAKURA_SYMBOL(InterpolateXAxisFloat)>;

PyObject *ApplyPositionSwitchCalibration(PyObject *self, PyObject *args){
	Py_ssize_t num_scaling_factor_py, num_data_py;
	enum {
		kFactor, kData, kOff, kResult, kEnd
	};
	PyObject *capsules[kEnd];
	if (!PyArg_ParseTuple(args, "nOnOOO", &num_scaling_factor_py, &capsules[kFactor],
			&num_data_py, &capsules[kData], &capsules[kOff], &capsules[kResult])) {
		return nullptr;
	}
	auto num_scaling_factor = static_cast<size_t>(num_scaling_factor_py);
	auto num_data = static_cast<size_t>(num_data_py);
	auto predFactor =
			[num_scaling_factor](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_scaling_factor);};
	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};
	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predFactor },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData }

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[kEnd];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = LIBSAKURA_SYMBOL(NormalizeDataAgainstReferenceFloat)(
				num_scaling_factor,
				reinterpret_cast<float const*>(bufs[kFactor]->aligned_addr),
				num_data,
				reinterpret_cast<float const*>(bufs[kData]->aligned_addr),
				reinterpret_cast<float const*>(bufs[kOff]->aligned_addr),
				reinterpret_cast<float *>(bufs[kResult]->aligned_addr));
		SAKURA_END_ALLOW_THREADS
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
			&& kernel_type < LIBSAKURA_SYMBOL(Convolve1DKernelType_kNumElements))
			&& (use_fft == Py_True || use_fft == Py_False))) {
		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
		return nullptr;
	}

	LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status =
				LIBSAKURA_SYMBOL(CreateConvolve1DContextFloat)(
						static_cast<size_t>(num_data),
						static_cast<LIBSAKURA_SYMBOL(Convolve1DKernelType)>(kernel_type),
						static_cast<size_t>(kernel_width), use_fft == Py_True,
						&context);
		SAKURA_END_ALLOW_THREADS

	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError,
				"sakura_CreateConvolve1DContext failed.");
		return nullptr;
	}
	PyCapsule_Destructor destructor = TCapsuleDesctructor<
	LIBSAKURA_SYMBOL(Convolve1DContextFloat), kConvolve1DContextName,
	LIBSAKURA_SYMBOL(DestroyConvolve1DContextFloat)>;
	auto capsule = PyCapsule_New(context, kConvolve1DContextName, destructor);
	if (capsule == nullptr) {
		PyErr_SetString(PyExc_MemoryError, "No memory.");
		return nullptr;
	}
	return capsule;
}

PyObject *Convolve1D(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kData, kResult, kEnd
	};
	PyObject *capsules[kEnd];
	PyObject *context_capsule;
	if (!PyArg_ParseTuple(args, "OnOO", &context_capsule,
			&num_data_py, &capsules[kData], &capsules[kResult])) {
		return nullptr;
	}
	auto context = reinterpret_cast<LIBSAKURA_SYMBOL(Convolve1DContextFloat) *>(
			PyCapsule_GetPointer(context_capsule,
					PyCapsule_GetName(context_capsule)));
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

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[kEnd];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs) || context == nullptr) {
		goto invalid_arg;
	}

	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = LIBSAKURA_SYMBOL(Convolve1DFloat)(
				context, num_data,
				reinterpret_cast<float const*>(bufs[kData]->aligned_addr),
				reinterpret_cast<float *>(bufs[kResult]->aligned_addr));
	SAKURA_END_ALLOW_THREADS
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

PyObject *CreateBaselineContext(PyObject *self, PyObject *args) {
	Py_ssize_t num_data;
	unsigned int order;
	int baseline_type;
	if (!PyArg_ParseTuple(args, "iIn", &baseline_type, &order, &num_data)) {
		return nullptr;
	}
	if (!(0 <= baseline_type
			&& baseline_type < LIBSAKURA_SYMBOL(BaselineType_kNumElements))) {
		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
		return nullptr;
	}

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) status;
	SAKURA_BEGIN_ALLOW_THREADS
		status =
				LIBSAKURA_SYMBOL(CreateBaselineContext)(
						static_cast<LIBSAKURA_SYMBOL(BaselineType)>(baseline_type),
						static_cast<uint16_t>(order),
						static_cast<size_t>(num_data),
						&context);
	SAKURA_END_ALLOW_THREADS

	if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
		PyErr_SetString(PyExc_ValueError,
				"sakura_CreateBaselineContext failed.");
		return nullptr;
	}
	PyCapsule_Destructor destructor = TCapsuleDesctructor<
	LIBSAKURA_SYMBOL(BaselineContext), kBaselineContextName,
	LIBSAKURA_SYMBOL(DestroyBaselineContext)>;
	auto capsule = PyCapsule_New(context, kBaselineContextName, destructor);
	if (capsule == nullptr) {
		PyErr_SetString(PyExc_MemoryError, "No memory.");
		return nullptr;
	}
	return capsule;
}

PyObject *SubtractBaseline(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	float clip_threshold_sigma;
	unsigned int num_fitting_max;
	uint16_t order = 3; // TODO make this python parameter.
	enum {
		kData, kMask, kFinalMask, kResult, kEnd
	};
	PyObject *capsules[kEnd];
	PyObject *context_capsule, *get_residual_py;
	if (!PyArg_ParseTuple(args, "nOOOfIOOO", &num_data_py, &capsules[kData], &capsules[kMask],
			&context_capsule, &clip_threshold_sigma, &num_fitting_max, &get_residual_py,
			&capsules[kFinalMask], &capsules[kResult])) {
		return nullptr;
	}
	auto context = reinterpret_cast<LIBSAKURA_SYMBOL(BaselineContext) *>(
			PyCapsule_GetPointer(context_capsule,
					PyCapsule_GetName(context_capsule)));
	if ((context == nullptr) || !(get_residual_py == Py_True || get_residual_py == Py_False ) ) {
		PyErr_SetString(PyExc_ValueError, "Invalid argument.");
		return nullptr;
	}
	auto num_data = static_cast<size_t>(num_data_py);
	auto predData =
			[num_data](LIBSAKURA_SYMBOL(PyAlignedBuffer) const &buf) -> bool
			{	return TotalElementsGreaterOrEqual(buf, num_data);};
	AlignedBufferConfiguration const conf[] = {

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kBool), predData },

	{ LIBSAKURA_SYMBOL(PyTypeId_kFloat), predData }

	};
	STATIC_ASSERT(ELEMENTSOF(conf) == kEnd);

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[kEnd];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
		goto invalid_arg;
	}

	LIBSAKURA_SYMBOL(Status) status;
	LIBSAKURA_SYMBOL(BaselineStatus) bl_status;
	SAKURA_BEGIN_ALLOW_THREADS
		status = LIBSAKURA_SYMBOL(SubtractBaselineFloat)(
				context, order, num_data,
				reinterpret_cast<float const*>(bufs[kData]->aligned_addr),
				reinterpret_cast<bool const*>(bufs[kMask]->aligned_addr),
				clip_threshold_sigma,
				static_cast<uint16_t>(num_fitting_max), (get_residual_py == Py_True),
				reinterpret_cast<bool *>(bufs[kFinalMask]->aligned_addr),
				reinterpret_cast<float *>(bufs[kResult]->aligned_addr),
				&bl_status);
	SAKURA_END_ALLOW_THREADS
	if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
		goto invalid_arg;
	}
	if (bl_status != LIBSAKURA_SYMBOL(BaselineStatus_kOK)) {
		PyErr_SetString(PyExc_ValueError, "Unexpected error.");
				return nullptr;
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
PyObject *ComplementMaskedValue(PyObject *self, PyObject *args) {
	Py_ssize_t num_data_py;
	enum {
		kData, kMask, kResult, kEnd
	};
	PyObject *capsules[kEnd];

	if (!PyArg_ParseTuple(args, "nOOO", &num_data_py, &capsules[kData],
			&capsules[kMask], &capsules[kResult])) {
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

	LIBSAKURA_SYMBOL(PyAlignedBuffer) *bufs[kEnd];
	if (!isValidAlignedBuffer(ELEMENTSOF(conf), conf, capsules, bufs)) {
		goto invalid_arg;
	}
	LIBSAKURA_SYMBOL(Status) status;
	//SAKURA_BEGIN_ALLOW_THREADS
		status = Func(num_data,
				reinterpret_cast<Type const*>(bufs[kData]->aligned_addr),
				reinterpret_cast<bool const*>(bufs[kMask]->aligned_addr),
				reinterpret_cast<Type *>(bufs[kResult]->aligned_addr));
	//SAKURA_END_ALLOW_THREADS
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

constexpr FuncForPython ComplementMaskedValueFloat = ComplementMaskedValue<
		float,
		LIBSAKURA_SYMBOL(PyTypeId_kFloat),
		LIBSAKURA_SYMBOL(ComplementMaskedValue)<float> >;

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

	RefHolder dims(PyTuple_New(buf->dimensions));
	if (dims.get() == nullptr) {
		goto out_of_memory;
	}
	for (size_t i = 0; i < buf->dimensions; ++i) {
		auto element = PyInt_FromSsize_t(buf->elements[i]);
		if (element == nullptr) {
			goto out_of_memory;
		}
		PyTuple_SetItem(dims.get(), i, element);
	}
	return dims.release();

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
			&& typeInt < LIBSAKURA_SYMBOL(PyTypeId_kEnd))) {
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
	STATIC_ASSERT(LIBSAKURA_SYMBOL(PyTypeId_kEnd) == ELEMENTSOF((sizes)));

	std::unique_ptr<LIBSAKURA_SYMBOL(PyAlignedBuffer),
			decltype(&LIBSAKURA_SYMBOL(PyAlignedBufferDestroy))> bufPtr(nullptr,
	LIBSAKURA_SYMBOL(PyAlignedBufferDestroy));
	try {
		char *aligned = nullptr;
		auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<char>(
				sizes[type] * total_elements, &aligned);
		LIBSAKURA_SYMBOL(PyAlignedBuffer) *buf = nullptr;
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(type, addr, aligned, dimensions,
				elems, LIBSAKURA_PREFIX::Memory::Free, &buf);
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
		case LIBSAKURA_SYMBOL(PyTypeId_kBool): {
			bool *aligned = nullptr;
			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
			bool>(sizeof(bool) * len, &aligned);
			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
				RefHolder item(PySequence_GetItem(dataSeq, i));
				aligned[i] = item.get() == Py_True;
			}
			LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
			if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
				throw std::bad_alloc();
			}
			assert(status == LIBSAKURA_SYMBOL(Status_kOK));
		}
			break;

		case LIBSAKURA_SYMBOL(PyTypeId_kInt8): {
			uint8_t *aligned = nullptr;
			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
					uint8_t>(sizeof(uint8_t) * len, &aligned);
			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
				RefHolder item(PySequence_GetItem(dataSeq, i));
				auto val = PyInt_AsLong(item.get());
				if (PyErr_Occurred()) {
					return nullptr;
				}
				aligned[i] = val;
			}
			LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
			if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
				throw std::bad_alloc();
			}
			assert(status == LIBSAKURA_SYMBOL(Status_kOK));
		}
			break;

		case LIBSAKURA_SYMBOL(PyTypeId_kInt32): {
			int32_t *aligned = nullptr;
			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
					int32_t>(sizeof(int32_t) * len, &aligned);
			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
				RefHolder item(PySequence_GetItem(dataSeq, i));
				auto val = PyInt_AsLong(item.get());
				if (PyErr_Occurred()) {
					return nullptr;
				}
				aligned[i] = val;
			}
			LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
			if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
				throw std::bad_alloc();
			}
			assert(status == LIBSAKURA_SYMBOL(Status_kOK));
		}
			break;

		case LIBSAKURA_SYMBOL(PyTypeId_kFloat): {
			float *aligned = nullptr;
			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
					float>(sizeof(float) * len, &aligned);
			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
				RefHolder item(PySequence_GetItem(dataSeq, i));
				auto val = PyFloat_AsDouble(item.get());
				if (PyErr_Occurred()) {
					return nullptr;
				}
				aligned[i] = val;
			}
			LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
			if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
				throw std::bad_alloc();
			}
			assert(status == LIBSAKURA_SYMBOL(Status_kOK));
		}
			break;

		case LIBSAKURA_SYMBOL(PyTypeId_kDouble): {
			double *aligned = nullptr;
			auto addr = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException<
					double>(sizeof(double) * len, &aligned);
			for (Py_ssize_t i = 0; (size_t) i < len; ++i) {
				RefHolder item(PySequence_GetItem(dataSeq, i));
				auto val = PyFloat_AsDouble(item.get());
				if (PyErr_Occurred()) {
					return nullptr;
				}
				aligned[i] = val;
			}
			LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
					(LIBSAKURA_SYMBOL(PyTypeId)) type, addr, aligned, 1,
					elements, LIBSAKURA_PREFIX::Memory::Free, &buf);
			if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
				throw std::bad_alloc();
			}
			assert(status == LIBSAKURA_SYMBOL(Status_kOK));
		}
			break;
		default:
			assert(false);
			break;
		}
	} catch (std::bad_alloc const&e) {
		PyErr_SetString(PyExc_MemoryError, "No memory.");
		return nullptr;
	}

	PyObject *capsule = nullptr;
	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(PyAlignedBufferEncapsulate)(buf, &capsule);
	if (status == LIBSAKURA_SYMBOL(Status_kNoMemory)) {
		PyErr_SetString(PyExc_MemoryError, "No memory.");
		return nullptr;
	}
	assert(status == LIBSAKURA_SYMBOL(Status_kOK));
	return capsule;
}

PyMethodDef module_methods[] =
		{

		{ "initialize", Initialize, METH_VARARGS, "Initializes libsakura." },

		{ "clean_up", CleanUp, METH_VARARGS, "Cleans up libsakura." },

		{ "get_current_time", GetCurrentTime, METH_VARARGS,
				"Gets current time in seconds of type double." },

		{ "compute_statistics", ComputeStatistics, METH_VARARGS,
				"Computes statistics of unmasked elements." },

		{ "grid_convolving", GridConvolving, METH_VARARGS,
				"Grids spectra on X-Y plane with convolving." },

		{ "uint8_to_bool", Uint8ToBool, METH_VARARGS,
				"Converts uint8 to bool." },

		{ "uint32_to_bool", Uint32ToBool, METH_VARARGS,
				"Converts uint32 to bool." },

		{ "invert_bool", InvertBool, METH_VARARGS, "Inverts bool." },

		{ "logical_and", LogicalAnd, METH_VARARGS,
				"Takes logical conjunction of boolean arrays." },

		{ "operate_bits_uint8_or", Uint8OperateBitsOr, METH_VARARGS,
				"Bit operation OR between an uint8 value and uint8 array." },

		{ "set_true_float_in_ranges_exclusive",
				FloatSetTrueIntInRangesExclusive, METH_VARARGS,
				"Sets True if the element is in at least one of ranges." },

		{ "set_true_int_in_ranges_exclusive",
				Int32SetTrueIntInRangesExclusive, METH_VARARGS,
				"Sets True if the element is in at least one of ranges." },

		{ "set_false_float_if_nan_or_inf", SetFalseFloatIfNanOrInf,
				METH_VARARGS, "set false if float value is NaN or Inf." },

		{ "interpolate_float_yaxis", InterpolateFloatYAxis,
				METH_VARARGS, "perform one-dimensional interpolation." },

		{ "interpolate_float_xaxis", InterpolateFloatXAxis,
						METH_VARARGS, "perform one-dimensional interpolation." },

		{ "apply_position_switch_calibration", ApplyPositionSwitchCalibration,
						METH_VARARGS, "apply position switch calibration." },

		{ "create_convolve1D_context", CreateConvolve1DContext,
				METH_VARARGS, "Creates a context for convolving 1D." },

		{ "convolve1D", Convolve1D,
				METH_VARARGS, "perform one-dimensional discrete convolution." },

		{ "create_baseline_context", CreateBaselineContext,
				METH_VARARGS, "Creates a context for baseline subtraction." },

		{ "subtract_baseline", SubtractBaseline,
				METH_VARARGS, "perform baseline subtraction." },

		{ "complement_masked_value_float", ComplementMaskedValueFloat,
				METH_VARARGS, "complement masked value with, tentatively, 0." },

		{ "get_elements_of_aligned_buffer", GetElementsOfAlignedBuffer,
				METH_VARARGS, "gets_elements of the aligned buffer." },

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
	PyModule_AddIntConstant(mod, "TYPE_BOOL", LIBSAKURA_SYMBOL(PyTypeId_kBool));
	PyModule_AddIntConstant(mod, "TYPE_INT8", LIBSAKURA_SYMBOL(PyTypeId_kInt8));
	PyModule_AddIntConstant(mod, "TYPE_INT32", LIBSAKURA_SYMBOL(PyTypeId_kInt32));
	PyModule_AddIntConstant(mod, "TYPE_FLOAT", LIBSAKURA_SYMBOL(PyTypeId_kFloat));
	PyModule_AddIntConstant(mod, "TYPE_DOUBLE",
	LIBSAKURA_SYMBOL(PyTypeId_kDouble));

	PyModule_AddIntConstant(mod, "CONVOLVE1D_KERNEL_TYPE_GAUSSIAN",
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian));
	PyModule_AddIntConstant(mod, "CONVOLVE1D_KERNEL_TYPE_BOXCAR",
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kBoxcar));
	PyModule_AddIntConstant(mod, "CONVOLVE1D_KERNEL_TYPE_HANNING",
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kHanning));
	PyModule_AddIntConstant(mod, "CONVOLVE1D_KERNEL_TYPE_HAMMING",
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kHamming));
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
	LIBSAKURA_SYMBOL(BaselineType_kPolynomial));
	PyModule_AddIntConstant(mod, "BASELINE_TYPE_CHEBYSHEV",
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev));
	PyModule_AddIntConstant(mod, "BASELINE_TYPE_CUBIC_SPLINE",
	LIBSAKURA_SYMBOL(BaselineType_kCubicSpline));
	PyModule_AddIntConstant(mod, "BASELINE_TYPE_SINUSOID",
	LIBSAKURA_SYMBOL(BaselineType_kSinusoid));
}

}
