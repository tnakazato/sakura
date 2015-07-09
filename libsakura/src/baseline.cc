/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2015
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
/**
 * baseline.cc
 *
 *  Created on: 2013/11/11
 *      Author: wataru
 */

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstddef>
#include <memory>
#include <sstream>
#include <stdexcept>

#if defined(__AVX__) && !defined(ARCH_SCALAR)
#	include <immintrin.h>
#endif

#include "libsakura/localdef.h"
#include "libsakura/logger.h"
#include "libsakura/memory_manager.h"
namespace {
#include "libsakura/packed_operation.h"
}
#include "libsakura/packed_type.h"
#include "libsakura/sakura.h"
#include "baseline.h"

namespace {

auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("baseline");
constexpr size_t kNumBasesCubicSpline = 4;

inline void AllocateMemoryForBasisData(
LIBSAKURA_SYMBOL(BaselineContext) *context) {
	size_t num_total_basis_data = context->num_bases * context->num_basis_data;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> work_basis_data_storage(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*context->basis_data) * num_total_basis_data,
					&context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	context->basis_data_storage = work_basis_data_storage.release();
}

inline size_t GetNumberOfCubicSplineLsqBases(size_t num_pieces) {
	return kNumBasesCubicSpline - 1 + num_pieces;
}

inline size_t GetNumberOfContextBases(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order) {
	size_t num_bases = 0;
	switch (baseline_type) {
	case LIBSAKURA_SYMBOL(BaselineType_kPolynomial):
	case LIBSAKURA_SYMBOL(BaselineType_kChebyshev):
		num_bases = order + 1;
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kCubicSpline):
		num_bases = kNumBasesCubicSpline;
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kSinusoid):
		num_bases = 2 * order + 1;
		break;
	default:
		assert(false);
		break;
	}
	return num_bases;
}

inline size_t GetNumberOfBasesFromOrder(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order,
		size_t const num_nwave, size_t const *nwave) {
	size_t num_bases = 0;
	switch (baseline_type) {
	case LIBSAKURA_SYMBOL(BaselineType_kPolynomial):
	case LIBSAKURA_SYMBOL(BaselineType_kChebyshev):
		num_bases = order + 1;
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kCubicSpline):
		if (order < 1) {
			throw std::invalid_argument(
					"order (number of pieces) must be a positive value!");
		}
		num_bases = GetNumberOfCubicSplineLsqBases(order);
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kSinusoid):
		num_bases = (nwave[0] == 0) ? (2 * num_nwave - 1) : (2 * num_nwave);
		break;
	default:
		assert(false);
		break;
	}
	return num_bases;
}

inline void DoSetBasisDataPolynomial(size_t num_bases, double const i_d,
		size_t *idx, double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto out = AssumeAligned(out_arg);

	size_t i = *idx;
	double val = 1.0;
	for (size_t j = 0; j < num_bases; ++j) {
		out[i++] = val;
		val *= i_d;
	}
	*idx = i;
}

inline void SetBasisDataPolynomial(LIBSAKURA_SYMBOL(BaselineContext) *context) {
	assert(0 < context->num_basis_data);
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	auto data = AssumeAligned(context->basis_data);

	size_t num_basis_data = context->num_basis_data;
	size_t num_bases = context->num_bases;
	size_t idx = 0;
	for (size_t i = 0; i < num_basis_data; ++i) {
		DoSetBasisDataPolynomial(num_bases, static_cast<double>(i), &idx, data);
	}
}

inline void SetBasisDataChebyshev(LIBSAKURA_SYMBOL(BaselineContext) *context) {
	assert(0 < context->num_basis_data);
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	auto data = AssumeAligned(context->basis_data);

	size_t num_basis_data = context->num_basis_data;
	size_t num_bases = context->num_bases;
	size_t idx = 0;
	if (num_basis_data == 1) {   // order == 0
		data[idx] = 1.0;
	} else if (num_bases == 1) { // order == 0
		for (size_t i = 0; i < num_basis_data; ++i) {
			data[idx++] = 1.0;
		}
	} else {
		double max_data_x = static_cast<double>(num_basis_data - 1);
		for (size_t i = 0; i < num_basis_data; ++i) {
			data[idx++] = 1.0;
			double x = 2.0 * static_cast<double>(i) / max_data_x - 1.0;
			data[idx++] = x;
			for (size_t j = 2; j < num_bases; ++j) {
				data[idx] = 2.0 * x * data[idx - 1] - data[idx - 2];
				++idx;
			}
		}
	}
}

inline void SetBasisDataSinusoid(LIBSAKURA_SYMBOL(BaselineContext) *context) {
	assert(0 < context->num_basis_data);
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	auto data = AssumeAligned(context->basis_data);

	size_t num_basis_data = context->num_basis_data;
	size_t max_nwave = context->baseline_param;
	size_t idx = 0;
	if (num_basis_data == 1) {
		data[idx] = 1.0;
	} else {
		double max_data_x = static_cast<double>(num_basis_data - 1);
		double norm = 2.0 * M_PI / max_data_x;
		for (size_t i = 0; i < num_basis_data; ++i) {
			data[idx++] = 1.0;
			for (size_t nwave = 1; nwave <= max_nwave; ++nwave) {
				double x = static_cast<double>(i) * nwave * norm;
				data[idx++] = sin(x);
				data[idx++] = cos(x);
			}
		}
	}
}

inline void SetBasisData(size_t const order,
LIBSAKURA_SYMBOL(BaselineContext) *context) {
	context->num_bases = GetNumberOfContextBases(context->baseline_type, order);
	size_t min_num_basis_data =
			(context->baseline_type
					== LIBSAKURA_SYMBOL(BaselineType_kCubicSpline)) ?
					GetNumberOfCubicSplineLsqBases(order) : context->num_bases;
	if (context->num_basis_data < min_num_basis_data) {
		throw std::invalid_argument("num_basis_data is too small!");
	}
	AllocateMemoryForBasisData(context);
	switch (context->baseline_type) {
	case LIBSAKURA_SYMBOL(BaselineType_kPolynomial):
		SetBasisDataPolynomial(context);
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kChebyshev):
		SetBasisDataChebyshev(context);
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kCubicSpline):
		assert(context->num_bases == kNumBasesCubicSpline);
		SetBasisDataPolynomial(context);
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kSinusoid):
		SetBasisDataSinusoid(context);
		break;
	default:
		assert(false);
		break;
	}
}

inline void DestroyBaselineContext(LIBSAKURA_SYMBOL(BaselineContext) *context) {
	if (context != nullptr) {
		if (context->basis_data_storage != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->basis_data_storage);
		}
		LIBSAKURA_PREFIX::Memory::Free(context);
	}
}

inline void CreateBaselineContext(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order,
		uint16_t const npiece, uint16_t const nwave,
		size_t const num_basis_data,
		LIBSAKURA_SYMBOL(BaselineContext) **context) {
	try {
		std::unique_ptr<LIBSAKURA_SYMBOL(BaselineContext),
		LIBSAKURA_PREFIX::Memory> work_context(
				static_cast<LIBSAKURA_SYMBOL(BaselineContext) *>(LIBSAKURA_PREFIX::Memory::Allocate(
						sizeof(LIBSAKURA_SYMBOL(BaselineContext)))));
		if (work_context == nullptr) {
			throw std::bad_alloc();
		}
		work_context->basis_data_storage = nullptr;
		work_context->baseline_type = baseline_type;
		work_context->num_basis_data = num_basis_data;
		switch (baseline_type) {
		case LIBSAKURA_SYMBOL(BaselineType_kPolynomial):
		case LIBSAKURA_SYMBOL(BaselineType_kChebyshev):
			work_context->baseline_param = order;
			break;
		case LIBSAKURA_SYMBOL(BaselineType_kCubicSpline):
			work_context->baseline_param = npiece;
			break;
		case LIBSAKURA_SYMBOL(BaselineType_kSinusoid):
			work_context->baseline_param = nwave;
			break;
		default:
			assert(false);
			break;
		}

		SetBasisData(work_context->baseline_param, work_context.get());
		*context = work_context.release();
	} catch (...) {
		DestroyBaselineContext(*context);
		throw;
	}
}

inline bool IsNWaveUniqueAndAscending(size_t num_nwave, size_t const *nwave) {
	//check if elements in nwave are
	// (1) stored in ascending order and
	// (2) not duplicate.
	bool res = true;
	for (size_t i = 0; i < num_nwave - 1; ++i) {
		if (nwave[i] >= nwave[i + 1]) {
			res = false;
			break;
		}
	}
	return res;
}

inline void OperateSubtractionFloat(size_t num_in, float const *in1_arg,
		double const *in2_arg, float *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in1_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const in1 = AssumeAligned(in1_arg);
	auto const in2 = AssumeAligned(in2_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_in; ++i) {
		out[i] = in1[i] - static_cast<float>(in2[i]);
	}
}

inline void AddMulMatrix(size_t num_coeff, double const *coeff,
		size_t const *use_idx, size_t num_out, size_t num_bases,
		double const *basis, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	constexpr size_t kPackElements = sizeof(__m256d) / sizeof(double);
	size_t const end = (num_out / kPackElements) * kPackElements;
	auto const zero = _mm256_set1_pd(0.);
	size_t const offset1 = num_bases * 1;
	size_t const offset2 = num_bases * 2;
	size_t const offset3 = num_bases * 3;
#if defined(__AVX2__) && 0 // <--- will make this part effective
	// once _mm256_i64gather_pd gets faster in future version.
	// cf. #744 (2015/7/9 WK)
	auto vindex = _mm256_set_epi64x(offset3, offset2, offset1, 0);
#endif
	for (i = 0; i < end; i += kPackElements) {
		auto total = zero;
		auto bases_row = &basis[num_bases * i];
		for (size_t j = 0; j < num_coeff; ++j) {
			auto ce = _mm256_set1_pd(coeff[j]);
			auto idx = use_idx[j];
#if defined(__AVX2__) && 0 // <--- will make this part effective
			// once _mm256_i64gather_pd gets faster in future version.
			// cf. #744(2015/7/9 WK)
			auto bs = _mm256_i64gather_pd(bases_row+idx, vindex, sizeof(double));
#else
			auto bs = _mm256_set_pd(bases_row[idx + offset3],
					bases_row[idx + offset2], bases_row[idx + offset1],
					bases_row[idx]);
#endif
			total = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
			LIBSAKURA_SYMBOL(SimdPacketAVX), double>(ce, bs, total);
		}
		_mm256_store_pd(&out[i], total);
	}
#endif
	for (; i < num_out; ++i) {
		double out_double = 0.0;
		for (size_t j = 0; j < num_coeff; ++j) {
			out_double += coeff[j] * basis[num_bases * i + use_idx[j]];
		}
		out[i] = out_double;
	}
}

inline void AddMulMatrixCubicSpline(size_t num_boundary, double const *boundary,
		double const *coeff_full, size_t num_out, double const *basis,
		double *out) {
	for (size_t i = 0; i < num_boundary; ++i) {
		size_t coffset = kNumBasesCubicSpline * i;
		size_t start_idx = static_cast<size_t>(ceil(boundary[i]));
		size_t end_idx =
				(i < num_boundary - 1) ?
						static_cast<size_t>(ceil(boundary[i + 1])) : num_out;
		size_t j = start_idx;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
		constexpr size_t kPackElements = sizeof(__m256d) / sizeof(double);
		size_t const start =
				(start_idx % kPackElements == 0) ?
						start_idx :
						(kPackElements
								+ (start_idx / kPackElements) * kPackElements);
		for (j = start_idx; j < start; ++j) {
			double out_double = 0.0;
			for (size_t k = 0; k < kNumBasesCubicSpline; ++k) {
				out_double += coeff_full[k + coffset]
						* basis[kNumBasesCubicSpline * j + k];
			}
			out[j] = out_double;
		}
		size_t const end = start
				+ ((end_idx - start) / kPackElements) * kPackElements;
		auto const zero = _mm256_set1_pd(0.);
		size_t const offset1 = kNumBasesCubicSpline * 1;
		size_t const offset2 = kNumBasesCubicSpline * 2;
		size_t const offset3 = kNumBasesCubicSpline * 3;
#if defined(__AVX2__) && 0 // <--- will make this part effective
		// once _mm256_i64gather_pd gets faster in future version.
		// cf. #744 (2015/7/9 WK)
		auto vindex = _mm256_set_epi64x(offset3, offset2, offset1, 0);
#endif
		for (j = start; j < end; j += kPackElements) {
			auto total = zero;
			auto bases_row = &basis[kNumBasesCubicSpline * j];
			for (size_t k = 0; k < kNumBasesCubicSpline; ++k) {
				auto ce = _mm256_set1_pd(coeff_full[k + coffset]);
#if defined(__AVX2__) && 0 // <--- will make this part effective
				// once _mm256_i64gather_pd gets faster in future
				// version. cf. #744 (2015/7/9 WK)
				auto bs = _mm256_i64gather_pd(bases_row, vindex, sizeof(double));
#else
				auto bs = _mm256_set_pd(bases_row[k + offset3],
						bases_row[k + offset2], bases_row[k + offset1],
						bases_row[k]);
#endif
				total = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
				LIBSAKURA_SYMBOL(SimdPacketAVX), double>(ce, bs, total);
			}
			_mm256_store_pd(&out[j], total);
		}
#endif
		for (; j < end_idx; ++j) {
			double out_double = 0.0;
			for (size_t k = 0; k < kNumBasesCubicSpline; ++k) {
				out_double += coeff_full[k + coffset]
						* basis[kNumBasesCubicSpline * j + k];
			}
			out[j] = out_double;
		}
	}
}

inline std::string GetNotEnoughDataMessage(
		uint16_t const idx_erroneous_fitting) {
	std::stringstream ss;
	ss << "SubtractBaseline: available data became too few in the ";
	ss << idx_erroneous_fitting;
	ss << " ";

	uint16_t imod10 = idx_erroneous_fitting % 10;
	std::string isuffix;
	if (imod10 == 1) {
		isuffix = "st";
	} else if (imod10 == 2) {
		isuffix = "nd";
	} else if (imod10 == 3) {
		isuffix = "rd";
	} else {
		isuffix = "th";
	}
	ss << isuffix << " fitting.";

	return ss.str();
}

inline void GetBoundariesOfPiecewiseData(size_t num_mask, bool const *mask_arg,
		size_t num_pieces, double *boundary_arg, size_t *start_arg,
		size_t *end_arg) {
	assert(num_pieces > 0);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(start_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(end_arg));
	auto const mask = AssumeAligned(mask_arg);
	auto boundary = AssumeAligned(boundary_arg);
	auto start = AssumeAligned(start_arg);
	auto end = AssumeAligned(end_arg);

	size_t num_unmasked_data = 0;
	for (size_t i = 0; i < num_mask; ++i) {
		if (mask[i])
			++num_unmasked_data;
	}
	boundary[0] = 0.0; // the first value of boundary[] must always point the first element.
	size_t idx = 1;
	size_t count_unmasked_data = 0;
	for (size_t i = 0; i < num_mask; ++i) {
		if (idx == num_pieces)
			break;
		if (mask[i]) {
			if (count_unmasked_data
					>= static_cast<double>(num_unmasked_data * idx)
							/ static_cast<double>(num_pieces)) {
				boundary[idx] = static_cast<double>(i);
				++idx;
			}
			++count_unmasked_data;
		}
	}
	for (size_t i = 0; i < num_pieces; ++i) {
		start[i] = static_cast<size_t>(ceil(boundary[i]));
		end[i] =
				(i < num_pieces - 1) ?
						static_cast<size_t>(boundary[i + 1]) : num_mask;
	}

}

inline void GetFullCubicSplineCoefficients(size_t num_pieces,
		double const *boundary_arg, double const *coeff_raw_arg,
		double *coeff_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_raw_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	auto const boundary = AssumeAligned(boundary_arg);
	auto const coeff_raw = AssumeAligned(coeff_raw_arg);
	auto coeff = AssumeAligned(coeff_arg);

	for (size_t i = 0; i < kNumBasesCubicSpline; ++i) {
		coeff[i] = coeff_raw[i];
	}
	for (size_t i = 1; i < num_pieces; ++i) {
		size_t ioffset = kNumBasesCubicSpline * i;
		size_t ioffset_prev = ioffset - kNumBasesCubicSpline;
		size_t j = GetNumberOfCubicSplineLsqBases(i);
		auto const c = coeff_raw[j] - coeff[ioffset_prev + 3];
		coeff[ioffset] = coeff[ioffset_prev]
				- boundary[i] * boundary[i] * boundary[i] * c;
		coeff[ioffset + 1] = coeff[ioffset_prev + 1]
				+ 3.0 * boundary[i] * boundary[i] * c;
		coeff[ioffset + 2] = coeff[ioffset_prev + 2] - 3.0 * boundary[i] * c;
		coeff[ioffset + 3] = coeff_raw[j];
	}
}

inline void SetStandardCubicBases(double const i_d, size_t *idx,
		double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto out = AssumeAligned(out_arg);
	DoSetBasisDataPolynomial(kNumBasesCubicSpline, i_d, idx, out);
}

inline void SetAuxiliaryCubicBases(size_t const num_boundary,
		double const *boundary_arg, double const boundary_final,
		double const i_d, size_t *idx, double *out_arg) {
	size_t i = *idx;
	assert(1 <= i);
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const boundary = AssumeAligned(boundary_arg);
	auto out = AssumeAligned(out_arg);

	auto cb = [](double v) {return v * v * v;};
	auto p = [](double v) {return std::max(0.0, v);};
	auto pcb = [&](double v) {return p(cb(v));};

	--i;
	auto boundary_current = boundary[1];
	out[i] -= pcb(i_d - boundary_current);
	++i;
	for (size_t j = 1; j < num_boundary - 1; ++j) {
		auto boundary_next = boundary[j + 1];
		out[i] = p(cb(i_d - boundary_current) - pcb(i_d - boundary_next));
		boundary_current = boundary_next;
		++i;
	}
	out[i] = pcb(i_d - boundary_final);
	++i;
	*idx = i;
}

inline void SetFullCubicSplineBasisData(size_t num_data, size_t num_boundary,
		double const *boundary_arg, double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const boundary = AssumeAligned(boundary_arg);
	auto out = AssumeAligned(out_arg);

	size_t idx = 0;
	if (num_boundary < 2) {
		for (size_t i = 0; i < num_data; ++i) {
			double i_d = static_cast<double>(i);
			SetStandardCubicBases(i_d, &idx, out);
		}
	} else {
		double boundary_final = boundary[num_boundary - 1];
		for (size_t i = 0; i < num_data; ++i) {
			double i_d = static_cast<double>(i);
			SetStandardCubicBases(i_d, &idx, out);
			SetAuxiliaryCubicBases(num_boundary, boundary, boundary_final, i_d,
					&idx, out);
		}
	}

}

inline void GetBestFitModelAndResidual(size_t num_data, float const *data,
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_coeff,
		double const *coeff, size_t const *use_idx, double *best_fit_model,
		float *residual_data) {
	AddMulMatrix(num_coeff, coeff, use_idx, num_data, context->num_bases,
			context->basis_data, best_fit_model);
	OperateSubtractionFloat(num_data, data, best_fit_model, residual_data);
}

inline void GetBestFitModelAndResidualCubicSpline(size_t num_data,
		float const *data, LIBSAKURA_SYMBOL(BaselineContext) const *context,
		size_t num_pieces, double const *boundary, double const *coeff_full,
		double *best_fit_model, float *residual_data) {
	assert(
			context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kCubicSpline));
	assert(context->num_bases == kNumBasesCubicSpline);
	AddMulMatrixCubicSpline(num_pieces, boundary, coeff_full, num_data,
			context->basis_data, best_fit_model);
	OperateSubtractionFloat(num_data, data, best_fit_model, residual_data);
}

inline void ClipData(size_t num_piece, size_t const *piece_start_arg,
		size_t const *piece_end_arg, float const *data_arg,
		bool const *in_mask_arg, float const lower_bound,
		float const upper_bound, bool *out_mask_arg,
		size_t *clipped_indices_arg, size_t *num_clipped) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(piece_start_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(piece_end_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	auto const piece_start = AssumeAligned(piece_start_arg);
	auto const piece_end = AssumeAligned(piece_end_arg);
	auto const data = AssumeAligned(data_arg);
	auto const in_mask = AssumeAligned(in_mask_arg);
	auto out_mask = AssumeAligned(out_mask_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);

	*num_clipped = 0;
	for (size_t ipiece = 0; ipiece < num_piece; ++ipiece) {
		for (size_t i = piece_start[ipiece]; i < piece_end[ipiece]; ++i) {
			out_mask[i] = in_mask[i];
			if (in_mask[i]) {
				if ((data[i] - lower_bound) * (upper_bound - data[i]) < 0.0f) {
					out_mask[i] = false;
					clipped_indices[*num_clipped] = i;
					++(*num_clipped);
				}
			}
		}
	}
}

template<typename FUNC>
inline void DoSubtractBaselineEngine(size_t num_data, float const *data_arg,
bool const *mask_arg, size_t num_context_bases, size_t num_coeff,
		size_t const *use_bases_idx_arg, double const *basis_arg,
		size_t num_pieces, size_t *piece_start_indices_arg,
		size_t *piece_end_indices_arg, uint16_t num_fitting_max,
		float clip_threshold_sigma, bool get_residual, double *coeff_arg,
		bool *final_mask_arg, float *rms, float *residual_data_arg,
		double *best_fit_model_arg, float *out_arg, FUNC func,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(use_bases_idx_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(basis_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(piece_start_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(piece_end_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(residual_data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(best_fit_model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto const use_bases_idx = AssumeAligned(use_bases_idx_arg);
	auto const basis = AssumeAligned(basis_arg);
	auto piece_start_indices = AssumeAligned(piece_start_indices_arg);
	auto piece_end_indices = AssumeAligned(piece_end_indices_arg);
	auto coeff = AssumeAligned(coeff_arg);
	auto final_mask = AssumeAligned(final_mask_arg);
	auto residual_data = AssumeAligned(residual_data_arg);
	auto best_fit_model = AssumeAligned(best_fit_model_arg);
	auto out = AssumeAligned(out_arg);

	size_t num_lsq_matrix = num_coeff * num_coeff;
	double *lsq_matrix = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_matrix(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*lsq_matrix) * num_lsq_matrix, &lsq_matrix));
	size_t num_lsq_vector = num_coeff;
	double *lsq_vector = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_vector(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*lsq_vector) * num_lsq_vector, &lsq_vector));
	size_t *clipped_indices = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_clipped_indices(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*clipped_indices) * num_data, &clipped_indices));

	size_t num_unmasked_data = num_data;
	size_t num_clipped = num_data;
	for (size_t i = 0; i < num_data; ++i) {
		final_mask[i] = mask[i];
	}
	LIBSAKURA_SYMBOL(Status) status;
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kOK);

	try {
		double rms_d = 0.0;
		for (uint16_t i = 1; i <= num_fitting_max; ++i) {
			if (num_unmasked_data < num_coeff) {
				*baseline_status = LIBSAKURA_SYMBOL(
						BaselineStatus_kNotEnoughData);
				throw std::runtime_error(GetNotEnoughDataMessage(i));
			}
			status = LIBSAKURA_SYMBOL(GetLSQCoefficientsDouble)(num_data, data,
					final_mask, num_context_bases, basis, num_coeff,
					use_bases_idx, lsq_matrix, lsq_vector);
			if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
				throw std::runtime_error("failed in GetLSQCoefficients.");
			}
			status = LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLUDouble)(
					num_coeff, lsq_matrix, lsq_vector, coeff);
			if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
				throw std::runtime_error(
						"failed in SolveSimultaneousEquationsByLU.");
			}

			func();

			LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
			status = LIBSAKURA_SYMBOL(ComputeAccurateStatisticsFloat)(num_data,
					residual_data, final_mask, &result);
			if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
				throw std::runtime_error(
						"failed in ComputeAccurateStatisticsFloat.");
			}
			assert(0 < result.count);
			double mean = result.sum / result.count;
			rms_d = sqrt(fabs(result.square_sum / result.count - mean * mean));

			if (i < num_fitting_max) {
				float clip_threshold_abs = clip_threshold_sigma * rms_d;
				float clip_threshold_lower = mean - clip_threshold_abs;
				float clip_threshold_upper = mean + clip_threshold_abs;
				ClipData(num_pieces, piece_start_indices, piece_end_indices,
						residual_data, final_mask, clip_threshold_lower,
						clip_threshold_upper, final_mask, clipped_indices,
						&num_clipped);
				if (num_clipped == 0) {
					break;
				}
				num_unmasked_data = result.count - num_clipped;
			}
		}
		if (out != nullptr) {
			for (size_t i = 0; i < num_data; ++i) {
				out[i] = get_residual ? residual_data[i] : best_fit_model[i];
			}
		}
		*rms = rms_d;
	} catch (...) {
		if (*baseline_status == LIBSAKURA_SYMBOL(BaselineStatus_kOK)) {
			*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kNG);
		}
		throw;
	}
}

inline void DoSubtractBaseline(size_t num_data, float const *data_arg,
bool const *mask_arg, LIBSAKURA_SYMBOL(BaselineContext) const *context,
		float clip_threshold_sigma, uint16_t num_fitting_max_arg,
		size_t num_coeff, size_t const *use_bases_idx_arg, double *coeff_arg,
		bool *final_mask_arg, float *rms, bool get_residual, float *out_arg,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(use_bases_idx_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto const use_bases_idx = AssumeAligned(use_bases_idx_arg);
	auto coeff = AssumeAligned(coeff_arg);
	auto final_mask = AssumeAligned(final_mask_arg);
	auto out = AssumeAligned(out_arg);

	SIMD_ALIGN
	size_t piece_start_index = 0;
	SIMD_ALIGN
	size_t piece_end_index = num_data;
	double *best_fit_model = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_best_fit_model(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*best_fit_model) * num_data, &best_fit_model));
	float *residual_data = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_residual_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*residual_data) * num_data, &residual_data));

	DoSubtractBaselineEngine(num_data, data, mask, context->num_bases,
			num_coeff, use_bases_idx, context->basis_data, 1,
			&piece_start_index, &piece_end_index,
			std::max(static_cast<uint16_t>(1), num_fitting_max_arg),
			clip_threshold_sigma, get_residual, coeff, final_mask, rms,
			residual_data, best_fit_model, out,
			[&]() {GetBestFitModelAndResidual(num_data, data, context,
						num_coeff, coeff, use_bases_idx, best_fit_model, residual_data);},
			baseline_status);
}

inline void DoSubtractBaselineCubicSpline(size_t num_data,
		float const *data_arg, bool const *mask_arg, size_t num_pieces,
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		float clip_threshold_sigma, uint16_t num_fitting_max_arg,
		double *coeff_full_arg, double *boundary_arg, bool *final_mask_arg,
		float *rms, bool get_residual, float *out_arg,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(
			context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kCubicSpline));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_full_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	auto const *data = AssumeAligned(data_arg);
	auto const *mask = AssumeAligned(mask_arg);
	auto *boundary = AssumeAligned(boundary_arg);
	auto *coeff_full = AssumeAligned(coeff_full_arg);
	auto *final_mask = AssumeAligned(final_mask_arg);
	auto *out = AssumeAligned(out_arg);

	size_t num_coeff = GetNumberOfCubicSplineLsqBases(num_pieces);

	size_t *use_bases_idx = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_use_bases_idx(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*use_bases_idx) * num_coeff, &use_bases_idx));
	for (size_t i = 0; i < num_coeff; ++i) {
		use_bases_idx[i] = i;
	}

	size_t *piece_start_indices = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_piece_start_indices(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*piece_start_indices) * num_pieces,
					&piece_start_indices));
	size_t *piece_end_indices = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_piece_end_indices(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*piece_end_indices) * num_pieces,
					&piece_end_indices));
	GetBoundariesOfPiecewiseData(num_data, mask, num_pieces, boundary,
			piece_start_indices, piece_end_indices);
	double *cspline_basis = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_aux_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*cspline_basis) * num_coeff * num_data,
					&cspline_basis));
	SetFullCubicSplineBasisData(num_data, num_pieces, boundary, cspline_basis);
	double *coeff = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_coeff(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*coeff) * num_coeff, &coeff));
	double *best_fit_model = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_best_fit_model(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*best_fit_model) * num_data, &best_fit_model));
	float *residual_data = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_residual_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*residual_data) * num_data, &residual_data));

	DoSubtractBaselineEngine(num_data, data, mask, num_coeff, num_coeff,
			use_bases_idx, cspline_basis, num_pieces, piece_start_indices,
			piece_end_indices,
			std::max(static_cast<uint16_t>(1), num_fitting_max_arg),
			clip_threshold_sigma, get_residual, coeff, final_mask, rms,
			residual_data, best_fit_model, out,
			[&]() {
				GetFullCubicSplineCoefficients(num_pieces, boundary, coeff, coeff_full);
				GetBestFitModelAndResidualCubicSpline(num_data, data, context,
						num_pieces, boundary, coeff_full, best_fit_model, residual_data);
			}, baseline_status);
}

inline void SetSinusoidUseBasesIndex(size_t const num_nwave,
		size_t const *nwave, size_t *use_bases_idx) {
	size_t iuse = 0;
	size_t i = 0;
	if (nwave[0] == 0) {
		use_bases_idx[iuse++] = nwave[i++];
	}
	for (; i < num_nwave; ++i) {
		use_bases_idx[iuse++] = 2 * nwave[i] - 1;
		use_bases_idx[iuse++] = 2 * nwave[i];
	}
}

inline void SubtractBaselineFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *context, uint16_t order,
		size_t const num_nwave, size_t const *nwave, size_t num_data,
		float const *data_arg, bool const *mask_arg, float clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual, bool *final_mask_arg,
		float *rms, float *out_arg,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	if ((context->baseline_type != LIBSAKURA_SYMBOL(BaselineType_kPolynomial))
			&& (context->baseline_type
					!= LIBSAKURA_SYMBOL(BaselineType_kChebyshev))
			&& (context->baseline_type
					!= LIBSAKURA_SYMBOL(BaselineType_kSinusoid))) {
		throw std::invalid_argument(
				"bad baseline context: baseline type must be polynomial, chebyshev, or sinusoid.");
	}
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto final_mask = AssumeAligned(final_mask_arg);
	auto out = AssumeAligned(out_arg);

	size_t num_coeff = GetNumberOfBasesFromOrder(context->baseline_type, order,
			num_nwave, nwave);
	size_t *use_bases_idx = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_use_bases_idx(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*use_bases_idx) * num_coeff, &use_bases_idx));
	if (context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kSinusoid)) {
		SetSinusoidUseBasesIndex(num_nwave, nwave, use_bases_idx);
	} else {
		for (size_t i = 0; i < num_coeff; ++i) {
			use_bases_idx[i] = i;
		}
	}
	double *coeff = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_coeff(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*coeff) * num_coeff, &coeff));
	DoSubtractBaseline(num_data, data, mask, context, clip_threshold_sigma,
			num_fitting_max, num_coeff, use_bases_idx, coeff, final_mask, rms,
			get_residual, out, baseline_status);
}

inline void SubtractBaselineCubicSplineFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_pieces,
		size_t num_data, float const *data_arg,
		bool const *mask_arg, float clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual, bool *final_mask_arg,
		float *rms, double *boundary_arg, float *out_arg,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	if (context->baseline_type != LIBSAKURA_SYMBOL(BaselineType_kCubicSpline)) {
		throw std::invalid_argument(
				"bad baseline context: baseline type must be cubic spline.");
	}
	assert(
			context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kCubicSpline));
	assert(context->num_bases == kNumBasesCubicSpline);
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));  //<-- will uncomment
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto final_mask = AssumeAligned(final_mask_arg);
	auto boundary = AssumeAligned(boundary_arg); //<-- will uncomment
	auto out = AssumeAligned(out_arg);

	double *coeff_full = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_coeff_full(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*coeff_full) * num_pieces * context->num_bases,
					&coeff_full));

	DoSubtractBaselineCubicSpline(num_data, data, mask, num_pieces, context,
			clip_threshold_sigma, num_fitting_max, coeff_full, boundary,
			final_mask, rms, get_residual, out, baseline_status);
}

inline void GetBestFitBaselineCoefficientsFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const *data_arg, bool const *mask_arg, float clip_threshold_sigma,
		uint16_t num_fitting_max, size_t num_nwave, size_t const *nwave,
		size_t num_coeff, double *coeff_arg,
		bool *final_mask_arg, float *rms,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto coeff = AssumeAligned(coeff_arg);
	auto final_mask = AssumeAligned(final_mask_arg);

	size_t *use_bases_idx = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_use_bases_idx(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*use_bases_idx) * num_coeff, &use_bases_idx));
	if (context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kSinusoid)) {
		SetSinusoidUseBasesIndex(num_nwave, nwave, use_bases_idx);
	} else {
		for (size_t i = 0; i < num_coeff; ++i) {
			use_bases_idx[i] = i;
		}
	}
	DoSubtractBaseline(num_data, data, mask, context, clip_threshold_sigma,
			num_fitting_max, num_coeff, use_bases_idx, coeff, final_mask, rms,
			true, nullptr, baseline_status);
}

inline void GetBestFitBaselineCoefficientsCubicSplineFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const *data_arg, bool const *mask_arg, float clip_threshold_sigma,
		uint16_t num_fitting_max, size_t num_pieces, double *coeff_arg,
		bool *final_mask_arg, float *rms, double *boundary_arg,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(
			context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kCubicSpline));
	assert(context->num_bases == kNumBasesCubicSpline);
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto coeff = AssumeAligned(coeff_arg);
	auto final_mask = AssumeAligned(final_mask_arg);
	auto boundary = AssumeAligned(boundary_arg);

	DoSubtractBaselineCubicSpline(num_data, data, mask, num_pieces, context,
			clip_threshold_sigma, num_fitting_max, coeff, boundary, final_mask,
			rms, true, nullptr, baseline_status);
}

inline void SubtractBaselineUsingCoefficientsFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const *data_arg, size_t num_coeff, double const *coeff_arg,
		size_t num_nwave, size_t const *nwave, float *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const data = AssumeAligned(data_arg);
	auto const coeff = AssumeAligned(coeff_arg);
	auto out = AssumeAligned(out_arg);

	size_t *use_bases_idx = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_use_bases_idx(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*use_bases_idx) * num_coeff, &use_bases_idx));
	if (context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kSinusoid)) {
		SetSinusoidUseBasesIndex(num_nwave, nwave, use_bases_idx);
	} else {
		for (size_t i = 0; i < num_coeff; ++i) {
			use_bases_idx[i] = i;
		}
	}
	double *best_fit_model = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_best_fit_model(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*best_fit_model) * num_data, &best_fit_model));
	GetBestFitModelAndResidual(num_data, data, context, num_coeff, coeff,
			use_bases_idx, best_fit_model, out);
}

inline void SubtractBaselineCubicSplineUsingCoefficientsFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const *data_arg, size_t num_pieces, double const *coeff_arg,
		double const *boundary_arg, float *out_arg) {
	if (num_pieces < 1)
		return;
	assert(
			context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kCubicSpline));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const data = AssumeAligned(data_arg);
	auto const coeff = AssumeAligned(coeff_arg);
	auto const boundary = AssumeAligned(boundary_arg);
	auto out = AssumeAligned(out_arg);

	double *best_fit_model = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_best_fit_model(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*best_fit_model) * num_data, &best_fit_model));

	try {
		GetBestFitModelAndResidualCubicSpline(num_data, data, context,
				num_pieces, boundary, coeff, best_fit_model, out);
	} catch (...) {
		throw;
	}
}

} /* anonymous namespace */

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateBaselineContext)(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order,
		uint16_t const npiece, uint16_t const nwave, size_t const num_data,
		LIBSAKURA_SYMBOL(BaselineContext) **context) noexcept {
	CHECK_ARGS(baseline_type < LIBSAKURA_SYMBOL(BaselineType_kNumElements));
	CHECK_ARGS(context != nullptr);
	size_t param;
	switch (baseline_type) {
	case LIBSAKURA_SYMBOL(BaselineType_kPolynomial):
	case LIBSAKURA_SYMBOL(BaselineType_kChebyshev):
		param = order;
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kCubicSpline):
		CHECK_ARGS(0 < npiece);
		param = npiece;
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kSinusoid):
		param = nwave;
		break;
	default:
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
		break;
	}
	size_t num_bases = GetNumberOfContextBases(baseline_type, param);
	CHECK_ARGS(num_bases <= num_data);

	try {
		CreateBaselineContext(baseline_type, order, npiece, nwave, num_data,
				context);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::invalid_argument &e) {
		LOG4CXX_ERROR(logger, "Order must be smaller than num_data.");
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyBaselineContext)(
LIBSAKURA_SYMBOL(BaselineContext) *context) noexcept {
	CHECK_ARGS(context != nullptr);

	try {
		DestroyBaselineContext(context);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

#include <iostream>
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetNumberOfCoefficients)(
LIBSAKURA_SYMBOL(BaselineContext) const *context, uint16_t const order,
		size_t *num_coeff) noexcept {
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(
			(context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kPolynomial)) || (context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kChebyshev)) || (context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kCubicSpline)));
	size_t num_out;
	LIBSAKURA_SYMBOL(BaselineType) type = context->baseline_type;
	num_out = GetNumberOfBasesFromOrder(type, order, 0, nullptr);
	if (context->baseline_type != LIBSAKURA_SYMBOL(BaselineType_kCubicSpline)) {
		CHECK_ARGS(num_out <= context->num_bases);
	}
	*num_coeff = num_out;
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineFloat)(
LIBSAKURA_SYMBOL(BaselineContext) const *context, uint16_t const order,
		size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], float clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual,
		bool final_mask[/*num_data*/], float out[/*num_data*/], float *rms,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(
			(context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kPolynomial)) || (context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kChebyshev)));
	size_t num_bases_from_order = GetNumberOfBasesFromOrder(
			context->baseline_type, order, 0, nullptr);
	CHECK_ARGS(num_bases_from_order <= context->num_bases);
	CHECK_ARGS(num_bases_from_order <= num_data);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(clip_threshold_sigma > 0.0f);
	CHECK_ARGS(final_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(final_mask));
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));
	CHECK_ARGS(baseline_status != nullptr);

	try {
		SubtractBaselineFloat(context, order, 0, nullptr, num_data, data, mask,
				clip_threshold_sigma, num_fitting_max, get_residual, final_mask,
				rms, out, baseline_status);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(
		LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_pieces,
		size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], float clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual,
		bool final_mask[/*num_data*/], float out[/*num_data*/], float *rms,
		double boundary[/*num_pieces*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(
			context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kCubicSpline));
	CHECK_ARGS(context->num_bases == kNumBasesCubicSpline);
	CHECK_ARGS(0 < num_pieces);
	size_t num_data_min = GetNumberOfCubicSplineLsqBases(num_pieces);
	CHECK_ARGS(num_data_min <= num_data);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(clip_threshold_sigma > 0.0f);
	CHECK_ARGS(final_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(final_mask));
	CHECK_ARGS(boundary != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(boundary));
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));
	CHECK_ARGS(baseline_status != nullptr);

	try {
		SubtractBaselineCubicSplineFloat(context, num_pieces, num_data, data,
				mask, clip_threshold_sigma, num_fitting_max, get_residual,
				final_mask, rms, boundary, out, baseline_status);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineSinusoidFloat)(
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		size_t const num_nwave, size_t const nwave[/*num_nwave*/],
		size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], float clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual,
		bool final_mask[/*num_data*/], float out[/*num_data*/], float *rms,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(
			context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kSinusoid));
	CHECK_ARGS(context->num_bases <= num_data);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(0.0f < clip_threshold_sigma);
	CHECK_ARGS(0 < num_nwave);
	CHECK_ARGS(nwave != nullptr);
	CHECK_ARGS(IsNWaveUniqueAndAscending(num_nwave, nwave));
	CHECK_ARGS(nwave[num_nwave - 1] <= context->baseline_param);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(final_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(final_mask));
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));
	CHECK_ARGS(baseline_status != nullptr);

	try {
		SubtractBaselineFloat(context, 0, num_nwave, nwave, num_data, data,
				mask, clip_threshold_sigma, num_fitting_max, get_residual,
				final_mask, rms, out, baseline_status);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsFloat)(
		LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		float clip_threshold_sigma, uint16_t num_fitting_max, size_t num_coeff,
		double coeff[/*num_coeff*/], bool final_mask[/*num_data*/], float *rms,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(
			(context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kPolynomial)) || (context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kChebyshev)));
	CHECK_ARGS(context->num_bases <= num_data);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(clip_threshold_sigma > 0.0f);
	CHECK_ARGS(coeff != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	CHECK_ARGS(final_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(final_mask));
	CHECK_ARGS(baseline_status != nullptr);

	try {
		GetBestFitBaselineCoefficientsFloat(context, num_data, data, mask,
				clip_threshold_sigma, num_fitting_max, 0, nullptr, num_coeff,
				coeff, final_mask, rms, baseline_status);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsCubicSplineFloat)(
		LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		float clip_threshold_sigma, uint16_t num_fitting_max, size_t num_pieces,
		double coeff[/*num_pieces*kNumBasesCubicSpline*/],
		bool final_mask[/*num_data*/], float *rms,
		double boundary[/*num_pieces*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(
			context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kCubicSpline));
	CHECK_ARGS(context->num_bases == kNumBasesCubicSpline);
	CHECK_ARGS(0 < num_pieces);
	size_t num_data_min = GetNumberOfCubicSplineLsqBases(num_pieces);
	CHECK_ARGS(num_data_min <= num_data);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(clip_threshold_sigma > 0.0f);
	CHECK_ARGS(coeff != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	CHECK_ARGS(final_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(final_mask));
	CHECK_ARGS(boundary != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(boundary));
	CHECK_ARGS(baseline_status != nullptr);

	try {
		GetBestFitBaselineCoefficientsCubicSplineFloat(context, num_data, data,
				mask, clip_threshold_sigma, num_fitting_max, num_pieces, coeff,
				final_mask, rms, boundary, baseline_status);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsSinusoidFloat)(
		LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		float clip_threshold_sigma, uint16_t num_fitting_max, size_t num_nwave,
		size_t const nwave[/*num_nwave*/], size_t num_coeff,
		double coeff[/*num_coeff*/], bool final_mask[/*num_data*/], float *rms,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(
			context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kSinusoid));
	CHECK_ARGS(context->num_bases <= num_data);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(0.0f < clip_threshold_sigma);
	CHECK_ARGS(0 < num_nwave);
	CHECK_ARGS(nwave != nullptr);
	CHECK_ARGS(IsNWaveUniqueAndAscending(num_nwave, nwave));
	CHECK_ARGS(nwave[num_nwave - 1] <= context->baseline_param);
	CHECK_ARGS(num_coeff <= context->num_bases);
	CHECK_ARGS(coeff != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	CHECK_ARGS(final_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(final_mask));
	CHECK_ARGS(baseline_status != nullptr);

	try {
		GetBestFitBaselineCoefficientsFloat(context, num_data, data, mask,
				clip_threshold_sigma, num_fitting_max, num_nwave, nwave,
				num_coeff, coeff, final_mask, rms, baseline_status);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineUsingCoefficientsFloat)(
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const data[/*num_data*/], size_t num_coeff,
		double const coeff[/*num_coeff*/], float out[/*num_data*/]) noexcept {
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(
			(context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kPolynomial)) || (context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kChebyshev)));
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(num_coeff != 0);
	CHECK_ARGS(num_coeff <= context->num_bases);
	CHECK_ARGS(num_coeff <= num_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(coeff != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));

	try {
		SubtractBaselineUsingCoefficientsFloat(context, num_data, data,
				num_coeff, coeff, 0, nullptr, out);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const data[/*num_data*/], size_t num_pieces,
		double const coeff[/*num_pieces*kNumBasesCubicSpline*/],
		double const boundary[/*num_pieces*/], float out[/*num_data*/])
				noexcept {
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(
			context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kCubicSpline));
	CHECK_ARGS(context->num_bases == kNumBasesCubicSpline);
	CHECK_ARGS(0 < num_pieces);
	size_t num_data_min = GetNumberOfCubicSplineLsqBases(num_pieces);
	CHECK_ARGS(num_data_min <= num_data);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(num_pieces <= INT_MAX);
	CHECK_ARGS(coeff != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	CHECK_ARGS(boundary != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(boundary));
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));

	try {
		SubtractBaselineCubicSplineUsingCoefficientsFloat(context, num_data,
				data, num_pieces, coeff, boundary, out);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineSinusoidUsingCoefficientsFloat)(
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const data[/*num_data*/], size_t num_nwave,
		size_t const nwave[/*num_nwave*/], size_t num_coeff,
		double const coeff[/*num_coeff*/], float out[/*num_data*/]) noexcept {
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(
			context->baseline_type == LIBSAKURA_SYMBOL(BaselineType_kSinusoid));
	CHECK_ARGS(context->num_bases <= num_data);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(0 < num_nwave);
	CHECK_ARGS(nwave != nullptr);
	CHECK_ARGS(IsNWaveUniqueAndAscending(num_nwave, nwave));
	CHECK_ARGS(nwave[num_nwave - 1] <= context->baseline_param);
	CHECK_ARGS(num_coeff <= context->num_bases);
	CHECK_ARGS(coeff != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));

	try {
		SubtractBaselineUsingCoefficientsFloat(context, num_data, data,
				num_coeff, coeff, num_nwave, nwave, out);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
