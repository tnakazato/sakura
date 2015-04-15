/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2014
 * National Astronomical Observatory of JapaGetNumBasesn
 * 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
 * 
 * This file is part of Sakura.
 * 
 * Sakura is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * Sakura is distriburedefinitionted in the hope that it will be useful, but WITHOUT
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
#include <cmath>
#include <cstddef>
#include <climits>
#include <iostream>
#include <iomanip>
#include <memory>
#include <sstream>

#if defined(__AVX__) && !defined(ARCH_SCALAR)
#	include <immintrin.h>
#endif

#include "libsakura/sakura.h"
#include "libsakura/localdef.h"
#include "libsakura/memory_manager.h"
#include "libsakura/logger.h"
#include "libsakura/packed_type.h"
namespace {
#include "libsakura/packed_operation.h"
}
#include "baseline.h"

namespace {

auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("baseline");

inline size_t GetNumberOfBasesFromOrder(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order) {
	size_t num_bases = 0;
	switch (baseline_type) {
	case LIBSAKURA_SYMBOL(BaselineType_kPolynomial):
		num_bases = order + 1;
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kChebyshev):
		num_bases = order + 1;
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kCubicSpline):
		if (order < 1) {
			throw std::invalid_argument(
					"order (number of pieces) must be a positive value!");
		}
		num_bases = 4;
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

inline void SetBasisDataPolynomial(LIBSAKURA_SYMBOL(BaselineContext) *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	auto data = AssumeAligned(context->basis_data);

	size_t idx = 0;
	for (size_t i = 0; i < context->num_basis_data; ++i) {
		double val = 1.0;
		data[idx] = val;
		idx++;
		for (size_t j = 1; j < context->num_bases; ++j) {
			val *= static_cast<double>(i);
			data[idx] = val;
			idx++;
		}
	}
}

inline void SetBasisDataChebyshev(LIBSAKURA_SYMBOL(BaselineContext) *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	auto data = AssumeAligned(context->basis_data);

	size_t idx = 0;
	double xrange = (double) (context->num_basis_data - 1);
	for (size_t i = 0; i < context->num_basis_data; ++i) {
		double x = 2.0 * (double) i / xrange - 1.0;
		for (size_t j = 0; j < context->num_bases; ++j) {
			double val = 0.0;
			if (j == 0) {
				val = 1.0;
			} else if (j == 1) {
				val = x;
			} else {
				val = 2.0 * x * data[idx - 1] - data[idx - 2];
			}
			data[idx] = val;
			idx++;
		}
	}
}

inline void SetBasisDataCubicSpline(
LIBSAKURA_SYMBOL(BaselineContext) *context) {
	SetBasisDataPolynomial(context);
}

inline void SetBasisDataSinusoid(LIBSAKURA_SYMBOL(BaselineContext) *context) {
	// TODO implement here.
	assert(false);
}

inline void SetBasisData(uint16_t const order,
LIBSAKURA_SYMBOL(BaselineContext) *context) {
	context->num_bases = GetNumberOfBasesFromOrder(context->baseline_type,
			order);
	AllocateMemoryForBasisData(context);
	switch (context->baseline_type) {
	case LIBSAKURA_SYMBOL(BaselineType_kPolynomial):
		SetBasisDataPolynomial(context);
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kChebyshev):
		SetBasisDataChebyshev(context);
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kCubicSpline):
		if (order < 1) {
			throw std::invalid_argument(
					"order (number of pieces) must be a positive value!");
		}
		SetBasisDataCubicSpline(context);
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kSinusoid):
		SetBasisDataSinusoid(context);
		break;
	default:
		assert(false);
		break;
	}
	if (context->num_bases > context->num_basis_data) {
		throw std::invalid_argument("num_bases exceeds num_basis_data!");
	}
}

inline void CreateBaselineContext(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order,
		size_t const num_basis_data,
		LIBSAKURA_SYMBOL(BaselineContext) **context) {
	std::unique_ptr<LIBSAKURA_SYMBOL(BaselineContext),
	LIBSAKURA_PREFIX::Memory> work_context(
			static_cast<LIBSAKURA_SYMBOL(BaselineContext) *>(LIBSAKURA_PREFIX::Memory::Allocate(
					sizeof(LIBSAKURA_SYMBOL(BaselineContext)))));
	if (work_context == nullptr) {
		throw std::bad_alloc();
	}
	work_context->baseline_type = baseline_type;
	work_context->num_basis_data = num_basis_data;
	SetBasisData(order, work_context.get());
	*context = (LIBSAKURA_SYMBOL(BaselineContext) *) work_context.release();
}

inline void DestroyBaselineContext(LIBSAKURA_SYMBOL(BaselineContext) *context) {
	if (context->basis_data_storage != nullptr) {
		LIBSAKURA_PREFIX::Memory::Free(context->basis_data_storage);
	}
	LIBSAKURA_PREFIX::Memory::Free(context);
}

inline void OperateFloatSubtraction(size_t num_in, float const *in1_arg,
		float const *in2_arg, float *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in1_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto in1 = AssumeAligned(in1_arg);
	auto in2 = AssumeAligned(in2_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_in; ++i) {
		out[i] = in1[i] - in2[i];
	}
}

inline void AddMulMatrix(size_t num_coeff, double const *coeff_arg,
		size_t num_out, size_t num_bases, double const *basis_arg,
		float *out_arg) {
	auto coeff = AssumeAligned(coeff_arg);
	auto basis = AssumeAligned(basis_arg);
	auto out = AssumeAligned(out_arg);
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (num_out / pack_elements) * pack_elements;
	auto const zero = _mm256_set1_pd(0.);
	size_t const offset1 = num_bases * 1;
	size_t const offset2 = num_bases * 2;
	size_t const offset3 = num_bases * 3;
	for (i = 0; i < end; i += pack_elements) {
		auto total = zero;
		auto bases_row = &basis[num_bases * i];
		for (size_t j = 0; j < num_coeff; ++j) {
			auto ce = _mm256_set1_pd(coeff[j]);
			auto bs = _mm256_set_pd(bases_row[j + offset3],
					bases_row[j + offset2], bases_row[j + offset1],
					bases_row[j]);
			total = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
			LIBSAKURA_SYMBOL(SimdPacketAVX), double>(ce, bs, total);
		}
		_mm_store_ps(&out[i], _mm256_cvtpd_ps(total));
	}
#endif
	for (; i < num_out; ++i) {
		double out_double = 0.0;
		for (size_t j = 0; j < num_coeff; ++j) {
			out_double += coeff[j] * basis[num_bases * i + j];
		}
		out[i] = out_double;
	}
}

inline void AddMulMatrixCubicSpline(size_t num_pieces, double const *boundary,
		size_t num_bases, double const *coeff_arg, size_t num_out,
		double const *basis_arg, float *out_arg) {
	auto coeff = AssumeAligned(coeff_arg);
	auto basis = AssumeAligned(basis_arg);
	auto out = AssumeAligned(out_arg);
	for (size_t i = 0; i < num_pieces; ++i) {
		size_t coffset = num_bases * i;
		size_t start_idx = static_cast<size_t>(ceil(boundary[i]));
		size_t end_idx =
				(i < num_pieces - 1) ?
						static_cast<size_t>(ceil(boundary[i + 1])) : num_out;
		size_t j = start_idx;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
		size_t const pack_elements = sizeof(__m256d) / sizeof(double);
		size_t const start =
				(start_idx % pack_elements == 0) ?
						start_idx :
						(pack_elements
								+ (start_idx / pack_elements) * pack_elements);
		for (j = start_idx; j < start; ++j) {
			double out_double = 0.0;
			for (size_t k = 0; k < num_bases; ++k) {
				out_double += coeff[k + coffset] * basis[num_bases * j + k];
			}
			out[j] = out_double;
		}
		size_t const end = start
				+ ((end_idx - start) / pack_elements) * pack_elements;
		auto const zero = _mm256_set1_pd(0.);
		size_t const offset1 = num_bases * 1;
		size_t const offset2 = num_bases * 2;
		size_t const offset3 = num_bases * 3;
		for (j = start; j < end; j += pack_elements) {
			auto total = zero;
			auto bases_row = &basis[num_bases * j];
			for (size_t k = 0; k < num_bases; ++k) {
				auto ce = _mm256_set1_pd(coeff[k + coffset]);
				auto bs = _mm256_set_pd(bases_row[k + offset3],
						bases_row[k + offset2], bases_row[k + offset1],
						bases_row[k]);
				total = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
				LIBSAKURA_SYMBOL(SimdPacketAVX), double>(ce, bs, total);
			}
			_mm_store_ps(&out[j], _mm256_cvtpd_ps(total));
		}
#endif
		for (; j < end_idx; ++j) {
			double out_double = 0.0;
			for (size_t k = 0; k < num_bases; ++k) {
				out_double += coeff[k + coffset] * basis[num_bases * j + k];
			}
			out[j] = out_double;
		}
	}
}

inline void GetBestFitBaselineFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *context, uint16_t const order,
		size_t num_data, float const *data_arg, bool const *mask_arg,
		float *out_arg, LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto basis = AssumeAligned(context->basis_data);
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto out = AssumeAligned(out_arg);

	size_t num_coeff = GetNumberOfBasesFromOrder(context->baseline_type, order);
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
	//size_t num_coeff = context->num_bases;
	double *coeff = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_coeff(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*coeff) * num_coeff, &coeff));

	LIBSAKURA_SYMBOL(Status) coeff_status =
	LIBSAKURA_SYMBOL(GetLeastSquareFittingCoefficientsDouble)(num_data, data,
			mask, context->num_bases, context->basis_data, num_coeff,
			lsq_matrix, lsq_vector);
	if (coeff_status != LIBSAKURA_SYMBOL(Status_kOK)) {
		throw std::runtime_error(
				"failed in GetLeastSquareFittingCoefficients.");
	}
	LIBSAKURA_SYMBOL(Status) solve_status =
	LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLUDouble)(num_coeff,
			lsq_matrix, lsq_vector, coeff);
	if (solve_status != LIBSAKURA_SYMBOL(Status_kOK)) {
		throw std::runtime_error("failed in SolveSimultaneousEquationsByLU.");
	}

	AddMulMatrix(num_coeff, coeff, num_data, context->num_bases, basis, out);
}

inline void ClipData(size_t num_data, float const *data_arg,
bool const *in_mask_arg, float const lower_bound, float const upper_bound,
bool *out_mask_arg, size_t *clipped_indices_arg, size_t *num_clipped) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	auto data = AssumeAligned(data_arg);
	auto in_mask = AssumeAligned(in_mask_arg);
	auto out_mask = AssumeAligned(out_mask_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);

	*num_clipped = 0;

	for (size_t i = 0; i < num_data; ++i) {
		out_mask[i] = in_mask[i];
		if (in_mask[i]) {
			if ((data[i] - lower_bound) * (upper_bound - data[i]) < 0.0f) {
				out_mask[i] = false;
				clipped_indices[*num_clipped] = i;
				(*num_clipped)++;
			}
		}
	}
}

inline void ClipDataPiecewise(size_t num_data, float const *data_arg,
bool const *in_mask_arg, float const lower_bound, float const upper_bound,
bool *out_mask_arg, size_t *clipped_indices_arg, size_t num_pieces,
		double const *boundary, size_t *num_clipped_arg,
		size_t *num_clipped_sum) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(num_clipped_arg));
	auto data = AssumeAligned(data_arg);
	auto in_mask = AssumeAligned(in_mask_arg);
	auto out_mask = AssumeAligned(out_mask_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto num_clipped = AssumeAligned(num_clipped_arg);

	*num_clipped_sum = 0;
	for (size_t i = 0; i < num_pieces; ++i) {
		num_clipped[i] = 0;
		size_t start_idx = static_cast<size_t>(ceil(boundary[i]));
		size_t end_idx =
				(i < num_pieces - 1) ?
						static_cast<size_t>(boundary[i + 1]) : num_data;
		for (size_t j = start_idx; j < end_idx; ++j) {
			if (in_mask[j]) {
				if ((data[j] - lower_bound) * (upper_bound - data[j]) < 0.0f) {
					out_mask[j] = false;
					++num_clipped[i];
					clipped_indices[*num_clipped_sum] = i;
					++(*num_clipped_sum);
				}
			}
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

inline void GetBestFitModelAndResidual(size_t num_data, float const *data,
LIBSAKURA_SYMBOL(BaselineContext) const *baseline_context, size_t num_coeff,
		double const *coeff, float *best_fit_model, float *residual_data) {
	AddMulMatrix(num_coeff, coeff, num_data, baseline_context->num_bases,
			baseline_context->basis_data, best_fit_model);
	OperateFloatSubtraction(num_data, data, best_fit_model, residual_data);
}

inline void GetBestFitModelAndResidualCubicSpline(size_t num_data,
		float const *data,
		LIBSAKURA_SYMBOL(BaselineContext) const *baseline_context,
		double const *coeff, size_t num_pieces, double const *boundary,
		float *best_fit_model, float *residual_data) {
	AddMulMatrixCubicSpline(num_pieces, boundary, baseline_context->num_bases,
			coeff, num_data, baseline_context->basis_data, best_fit_model);
	OperateFloatSubtraction(num_data, data, best_fit_model, residual_data);
}

inline void GetBoundariesOfPiecewiseData(size_t num_mask,
bool const *mask, size_t num_pieces, double *boundary) {
	assert(num_pieces > 0);

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
}

inline void GetUnmaskedDataNumbers(size_t num_data, float const *data,
bool const *mask, size_t num_pieces, double const *boundary,
		size_t *num_clipped) {
	for (size_t i = 0; i < num_pieces; ++i) {
		num_clipped[i] = 0;
		size_t start_idx = static_cast<size_t>(ceil(boundary[i]));
		size_t end_idx =
				(i < num_pieces - 1) ?
						static_cast<size_t>(ceil(boundary[i + 1])) : num_data;
		for (size_t j = start_idx; j < end_idx; ++j) {
			if (mask[j]) {
				++num_clipped[i];
			}
		}
	}
}

inline void GetFullCubicSplineCoefficients(size_t num_pieces,
		double const *boundary, double const *coeff_raw, double *coeff) {
	size_t num_bases = 4;
	for (size_t i = 0; i < num_bases; ++i) {
		coeff[i] = coeff_raw[i];
	}
	/*
	 for (size_t i = 1; i < num_pieces; ++i) {
	 size_t ioffset = num_bases * i;
	 size_t ioffset_prev = ioffset - num_bases;
	 size_t j = num_bases - 1 + i;
	 coeff[ioffset] = coeff[ioffset_prev]
	 - boundary[i] * boundary[i] * boundary[i] * coeff_raw[j];
	 coeff[ioffset + 1] = coeff[ioffset_prev + 1]
	 + 3.0 * boundary[i] * boundary[i] * coeff_raw[j];
	 coeff[ioffset + 2] = coeff[ioffset_prev + 2]
	 - 3.0 * boundary[i] * coeff_raw[j];
	 coeff[ioffset + 3] = coeff[ioffset_prev + 3] + coeff_raw[j];
	 }
	 */
	for (size_t i = 1; i < num_pieces; ++i) {
		size_t ioffset = num_bases * i;
		size_t ioffset_prev = ioffset - num_bases;
		size_t j = num_bases - 1 + i;
		auto c = coeff_raw[j] - coeff[ioffset_prev + 3];
		coeff[ioffset] = coeff[ioffset_prev]
				- boundary[i] * boundary[i] * boundary[i] * c;
		coeff[ioffset + 1] = coeff[ioffset_prev + 1]
				+ 3.0 * boundary[i] * boundary[i] * c;
		coeff[ioffset + 2] = coeff[ioffset_prev + 2] - 3.0 * boundary[i] * c;
		coeff[ioffset + 3] = coeff_raw[j];
	}
}

inline void DoSubtractBaselineCubicSpline(size_t num_data, float const *data,
bool const *mask, size_t num_pieces, double const *boundary,
LIBSAKURA_SYMBOL(BaselineContext) const *baseline_context,
		float clip_threshold_sigma, uint16_t num_fitting_max_arg,
		size_t num_coeff, double *coeff, bool *final_mask, bool get_residual,
		float *out,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {

	float *best_fit_model = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_best_fit_model(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*best_fit_model) * num_data, &best_fit_model));
	float *residual_data = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_residual_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*residual_data) * num_data, &residual_data));
	size_t *clipped_indices = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_clipped_indices(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*clipped_indices) * num_data, &clipped_indices));
	size_t num_bases = baseline_context->num_bases - 1 + num_pieces;
	size_t num_lsq_matrix = num_bases * num_bases;
	double *lsq_matrix = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_matrix(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*lsq_matrix) * num_lsq_matrix, &lsq_matrix));
	size_t num_lsq_vector = num_bases;
	double *lsq_vector = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_vector(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*lsq_vector) * num_lsq_vector, &lsq_vector));
	size_t num_coeff_raw = num_bases;
	double *coeff_raw = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_coeff_raw(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*coeff_raw) * num_coeff_raw, &coeff_raw));

	uint16_t num_fitting_max = std::max((uint16_t) 1, num_fitting_max_arg);

	//size_t num_unmasked_data = num_data;
	size_t *num_unmasked_data = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_num_unmasked_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*num_unmasked_data) * num_pieces,
					&num_unmasked_data));
	GetUnmaskedDataNumbers(num_data, data, mask, num_pieces, boundary,
			num_unmasked_data);
	//size_t num_clipped = num_data;
	size_t *num_clipped = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_num_clipped(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*num_clipped) * num_pieces, &num_clipped));
	GetUnmaskedDataNumbers(num_data, data, mask, num_pieces, boundary,
			num_clipped);

	for (size_t i = 0; i < num_data; ++i) {
		final_mask[i] = mask[i];
	}

	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kOK);

	try {
		for (uint16_t i = 0; i < num_fitting_max; ++i) {
			/*
			 if (num_unmasked_data < baseline_context->num_bases) {
			 *baseline_status = LIBSAKURA_SYMBOL(
			 BaselineStatus_kNotEnoughData);
			 throw std::runtime_error(GetNotEnoughDataMessage(i + 1));
			 }
			 */
			for (size_t j = 0; j < num_pieces; ++j) {
				if (num_unmasked_data[j] < baseline_context->num_bases) {
					*baseline_status = LIBSAKURA_SYMBOL(
							BaselineStatus_kNotEnoughData);
					throw std::runtime_error(GetNotEnoughDataMessage(i + 1));
				}
			}

			LIBSAKURA_SYMBOL(Status) coeff_status =
					LIBSAKURA_SYMBOL(GetLeastSquareFittingCoefficientsCubicSplineDouble)(
							num_data, data, final_mask, num_pieces, boundary,
							baseline_context->basis_data, lsq_matrix,
							lsq_vector);
			if (coeff_status != LIBSAKURA_SYMBOL(Status_kOK)) {
				throw std::runtime_error(
						"failed in GetLeastSquareFittingCoefficients.");
			}

			/*
			 for (size_t j = 0; j < num_bases; ++j) {
			 std::cout << " | ";
			 for (size_t k = 0; k < num_bases; ++k) {
			 std::cout << std::setw(8) << std::right
			 << lsq_matrix[j * num_bases + k] << "  ";
			 }
			 std::cout << " | " << std::setw(8) << std::right
			 << lsq_vector[j] << std::endl;
			 }
			 */
			/*
			 if (num_unmasked_data <= num_clipped) {
			 LIBSAKURA_SYMBOL(Status) coeff_status =
			 LIBSAKURA_SYMBOL(GetLeastSquareFittingCoefficientsDouble)(
			 num_data, data, final_mask, baseline_context->num_bases,
			 baseline_context->basis_data, lsq_matrix, lsq_vector);
			 if (coeff_status != LIBSAKURA_SYMBOL(Status_kOK)) {
			 throw std::runtime_error(
			 "failed in GetLeastSquareFittingCoefficients.");
			 }
			 } else {
			 LIBSAKURA_SYMBOL(Status) coeff_status =
			 LIBSAKURA_SYMBOL(UpdateLeastSquareFittingCoefficientsDouble)(
			 num_data, data, num_clipped, clipped_indices,
			 baseline_context->num_bases,
			 baseline_context->basis_data, lsq_matrix, lsq_vector);
			 if (coeff_status != LIBSAKURA_SYMBOL(Status_kOK)) {
			 throw std::runtime_error(
			 "failed in UpdateLeastSquareFittingCoefficients.");
			 }
			 }
			 */
			LIBSAKURA_SYMBOL(Status) solve_status =
			LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLUDouble)(num_bases,
					lsq_matrix, lsq_vector, coeff_raw);
			if (solve_status != LIBSAKURA_SYMBOL(Status_kOK)) {
				throw std::runtime_error(
						"failed in SolveSimultaneousEquationsByLU.");
			}

			/*
			 for (size_t j = 0; j < num_bases; ++j) {
			 std::cout << "#### coeff_raw[" << j << "] = " << coeff_raw[j]
			 << std::endl;
			 }
			 */
			GetFullCubicSplineCoefficients(num_pieces, boundary, coeff_raw,
					coeff);
			/*
			 for (size_t j = 0; j < num_pieces; ++j) {
			 std::cout << "   [ ";
			 for (size_t k = 0; k < 4; ++k) {
			 std::cout << coeff[4 * j + k] << ", ";
			 }
			 std::cout << " ]" << std::endl << std::flush;
			 }
			 */
			GetBestFitModelAndResidualCubicSpline(num_data, data,
					baseline_context, coeff, num_pieces, boundary,
					best_fit_model, residual_data);

			LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
			LIBSAKURA_SYMBOL(Status) stat_status =
			LIBSAKURA_SYMBOL(ComputeAccurateStatisticsFloat)(num_data,
					residual_data, final_mask, &result);
			if (stat_status != LIBSAKURA_SYMBOL(Status_kOK)) {
				throw std::runtime_error(
						"failed in ComputeAccurateStatisticsFloat.");
			}
			float clip_threshold_abs = clip_threshold_sigma * result.stddev;
			float clip_threshold_lower = result.mean - clip_threshold_abs;
			float clip_threshold_upper = result.mean + clip_threshold_abs;
			size_t num_clipped_sum;
			ClipDataPiecewise(num_data, residual_data, final_mask,
					clip_threshold_lower, clip_threshold_upper, final_mask,
					clipped_indices, num_pieces, boundary, num_clipped,
					&num_clipped_sum);
			if (num_clipped_sum == 0) {
				break;
			}
			for (size_t j = 0; j < num_pieces; ++j) {
				num_unmasked_data[j] = result.count - num_clipped[j];
			}
		}
		if (out != nullptr) {
			for (size_t i = 0; i < num_data; ++i) {
				out[i] = get_residual ? residual_data[i] : best_fit_model[i];
			}
		}

	} catch (...) {
		if (*baseline_status == LIBSAKURA_SYMBOL(BaselineStatus_kOK)) {
			*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kNG);
		}
		throw;
	}
}

inline void DoSubtractBaseline(size_t num_data, float const *data,
bool const *mask,
LIBSAKURA_SYMBOL(BaselineContext) const *baseline_context,
		float clip_threshold_sigma, uint16_t num_fitting_max_arg,
		size_t num_coeff, double *coeff, bool *final_mask, bool get_residual,
		float *out,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {

	float *best_fit_model = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_best_fit_model(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*best_fit_model) * num_data, &best_fit_model));
	float *residual_data = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_residual_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*residual_data) * num_data, &residual_data));
	size_t *clipped_indices = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_clipped_indices(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*clipped_indices) * num_data, &clipped_indices));
	//size_t num_bases = baseline_context->num_bases;
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

	uint16_t num_fitting_max = std::max((uint16_t) 1, num_fitting_max_arg);
	size_t num_unmasked_data = num_data;
	size_t num_clipped = num_data;

	for (size_t i = 0; i < num_data; ++i) {
		final_mask[i] = mask[i];
	}

	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kOK);

	try {
		for (uint16_t i = 0; i < num_fitting_max; ++i) {
			if (num_unmasked_data < num_coeff) {
				*baseline_status = LIBSAKURA_SYMBOL(
						BaselineStatus_kNotEnoughData);
				throw std::runtime_error(GetNotEnoughDataMessage(i + 1));
			}

			if (num_unmasked_data <= num_clipped) {
				LIBSAKURA_SYMBOL(Status) coeff_status =
				LIBSAKURA_SYMBOL(GetLeastSquareFittingCoefficientsDouble)(
						num_data, data, final_mask, baseline_context->num_bases,
						baseline_context->basis_data, num_coeff, lsq_matrix,
						lsq_vector);
				if (coeff_status != LIBSAKURA_SYMBOL(Status_kOK)) {
					throw std::runtime_error(
							"failed in GetLeastSquareFittingCoefficients.");
				}
			} else {
				LIBSAKURA_SYMBOL(Status) coeff_status =
				LIBSAKURA_SYMBOL(UpdateLeastSquareFittingCoefficientsDouble)(
						num_data, data, num_clipped, clipped_indices,
						baseline_context->num_bases,
						baseline_context->basis_data, num_coeff, lsq_matrix,
						lsq_vector);
				if (coeff_status != LIBSAKURA_SYMBOL(Status_kOK)) {
					throw std::runtime_error(
							"failed in UpdateLeastSquareFittingCoefficients.");
				}
			}

			LIBSAKURA_SYMBOL(Status) solve_status =
			LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLUDouble)(num_coeff,
					lsq_matrix, lsq_vector, coeff);

			if (solve_status != LIBSAKURA_SYMBOL(Status_kOK)) {
				throw std::runtime_error(
						"failed in SolveSimultaneousEquationsByLU.");
			}

			GetBestFitModelAndResidual(num_data, data, baseline_context,
					num_coeff, coeff, best_fit_model, residual_data);

			LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
			LIBSAKURA_SYMBOL(Status) stat_status =
			LIBSAKURA_SYMBOL(ComputeAccurateStatisticsFloat)(num_data,
					residual_data, final_mask, &result);
			if (stat_status != LIBSAKURA_SYMBOL(Status_kOK)) {
				throw std::runtime_error(
						"failed in ComputeAccurateStatisticsFloat.");
			}
			float clip_threshold_abs = clip_threshold_sigma * result.stddev;
			float clip_threshold_lower = result.mean - clip_threshold_abs;
			float clip_threshold_upper = result.mean + clip_threshold_abs;
			ClipData(num_data, residual_data, final_mask, clip_threshold_lower,
					clip_threshold_upper, final_mask, clipped_indices,
					&num_clipped);
			if (num_clipped == 0) {
				break;
			}
			num_unmasked_data = result.count - num_clipped;
		}
		if (out != nullptr) {
			for (size_t i = 0; i < num_data; ++i) {
				out[i] = get_residual ? residual_data[i] : best_fit_model[i];
			}
		}

	} catch (...) {
		if (*baseline_status == LIBSAKURA_SYMBOL(BaselineStatus_kOK)) {
			*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kNG);
		}
		throw;
	}
}

inline void GetBestFitBaselineCoefficientsFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const *data_arg, bool const *mask_arg, float clip_threshold_sigma,
		uint16_t num_fitting_max, size_t num_coeff, double *coeff_arg,
		bool *final_mask_arg,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto coeff = AssumeAligned(coeff_arg);
	auto final_mask = AssumeAligned(final_mask_arg);

	DoSubtractBaseline(num_data, data, mask, context, clip_threshold_sigma,
			num_fitting_max, num_coeff, coeff, final_mask,
			true, nullptr, baseline_status);

}

inline void GetBestFitBaselineCoefficientsCubicSplineFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const *data_arg, bool const *mask_arg, float clip_threshold_sigma,
		uint16_t num_fitting_max, size_t num_pieces, double *coeff_arg,
		bool *final_mask_arg,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto coeff = AssumeAligned(coeff_arg);
	auto final_mask = AssumeAligned(final_mask_arg);

	double *boundary = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_boundary(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*boundary) * num_pieces, &boundary));
	GetBoundariesOfPiecewiseData(num_data, mask, num_pieces, boundary);
	size_t num_coeff = context->num_bases * num_pieces;

	DoSubtractBaselineCubicSpline(num_data, data, mask, num_pieces, boundary,
			context, clip_threshold_sigma, num_fitting_max, num_coeff, coeff,
			final_mask, true, nullptr, baseline_status);
}

inline void SubtractBaselineFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *baseline_context, uint16_t order,
		size_t num_data, float const *data_arg,
		bool const *mask_arg, float clip_threshold_sigma,
		uint16_t num_fitting_max_arg, bool get_residual, bool *final_mask_arg,
		float *out_arg, LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(baseline_context->basis_data));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto final_mask = AssumeAligned(final_mask_arg);
	auto out = AssumeAligned(out_arg);

	size_t num_coeff = GetNumberOfBasesFromOrder(
			baseline_context->baseline_type, order);
	double *coeff = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_coeff(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*coeff) * num_coeff, &coeff));

	DoSubtractBaseline(num_data, data, mask, baseline_context,
			clip_threshold_sigma, num_fitting_max_arg, num_coeff, coeff,
			final_mask, get_residual, out, baseline_status);

}

inline void SubtractBaselineCubicSplineFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *baseline_context, size_t num_pieces,
		size_t num_data, float const *data_arg,
		bool const *mask_arg, float clip_threshold_sigma,
		uint16_t num_fitting_max_arg, bool get_residual, bool *final_mask_arg,
		float *out_arg, LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	if (baseline_context->baseline_type
			!= LIBSAKURA_SYMBOL(BaselineType_kCubicSpline)) {
		throw std::invalid_argument(
				"bad baseline context: baseline type must be cubic spline.");
	}
	assert(LIBSAKURA_SYMBOL(IsAligned)(baseline_context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto final_mask = AssumeAligned(final_mask_arg);
	auto out = AssumeAligned(out_arg);

	double *boundary = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_boundary(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*boundary) * num_pieces, &boundary));
	GetBoundariesOfPiecewiseData(num_data, mask, num_pieces, boundary);

	size_t num_coeff = baseline_context->num_bases * num_pieces;
	double *coeff = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_coeff(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*coeff) * num_coeff, &coeff));

	DoSubtractBaselineCubicSpline(num_data, data, mask, num_pieces, boundary,
			baseline_context, clip_threshold_sigma, num_fitting_max_arg,
			num_coeff, coeff, final_mask, get_residual, out, baseline_status);

}

inline void SubtractBaselinePolynomialFloat(size_t num_data,
		float const *data_arg,
		bool const *mask_arg, uint16_t order, float clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual,
		bool *final_mask_arg, float *out_arg,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto final_mask = AssumeAligned(final_mask_arg);
	auto out = AssumeAligned(out_arg);

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	CreateBaselineContext(LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order,
			num_data, &context);
	ScopeGuard guard_for_baseline_context([&]() {
		DestroyBaselineContext(context);
	});

	SubtractBaselineFloat(context, order, num_data, data, mask,
			clip_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, baseline_status);
}

inline void SubtractBaselineUsingCoefficientsFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_data,
		float const *data_arg, size_t num_coeff, double const *coeff_arg,
		float *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto data = AssumeAligned(data_arg);
	auto coeff = AssumeAligned(coeff_arg);
	auto out = AssumeAligned(out_arg);

	float *best_fit_model = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_best_fit_model(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*best_fit_model) * num_data, &best_fit_model));

	try {
		GetBestFitModelAndResidual(num_data, data, context, num_coeff, coeff,
				best_fit_model, out);
	} catch (...) {
		throw;
	}
}

inline void SubtractBaselineCubicSplineUsingCoefficientsFloat(
LIBSAKURA_SYMBOL(BaselineContext) const *baseline_context, size_t num_data,
		float const *data_arg, size_t num_pieces, double const *coeff_arg,
		double const *boundary_arg, float *out_arg) {
	if (num_pieces < 1)
		return;
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto data = AssumeAligned(data_arg);
	auto coeff = AssumeAligned(coeff_arg);
	auto boundary = AssumeAligned(boundary_arg);
	auto out = AssumeAligned(out_arg);

	float *best_fit_model = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_best_fit_model(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*best_fit_model) * num_data, &best_fit_model));

	try {
		GetBestFitModelAndResidualCubicSpline(num_data, data, baseline_context,
				coeff, num_pieces, boundary, best_fit_model, out);
	} catch (...) {
		throw;
	}
}

} /* anonymous namespace */

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateBaselineContext)(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order,
		size_t const num_data,
		LIBSAKURA_SYMBOL(BaselineContext) **context) noexcept {
	if (baseline_type == LIBSAKURA_SYMBOL(BaselineType_kSinusoid)) {
		//--- fails for sinusoid until sinusoidal fitting is implemented (2015/3/25 WK)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (baseline_type >= LIBSAKURA_SYMBOL(BaselineType_kNumElements)) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (context == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	try {
		CreateBaselineContext(baseline_type, order, num_data, context);
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
	if (context == nullptr) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	try {
		DestroyBaselineContext(context);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitBaselineFloat)(
LIBSAKURA_SYMBOL(BaselineContext) const *context, uint16_t const order,
		size_t num_data, float const data[], bool const mask[], float out[],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	if (context == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	size_t num_bases_from_order(
			GetNumberOfBasesFromOrder(context->baseline_type, order));
	if (num_bases_from_order > context->num_bases)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data < num_bases_from_order)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data != context->num_basis_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (baseline_status == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	try {
		GetBestFitBaselineFloat(context, order, num_data, data, mask, out,
				baseline_status);
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
		float const data[], bool const mask[], float clip_threshold_sigma,
		uint16_t num_fitting_max, size_t num_coeff, double coeff[],
		bool final_mask[],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	if (context == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data != context->num_basis_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data < context->num_bases)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (clip_threshold_sigma <= 0.0f)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (coeff == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(coeff)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (final_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(final_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (baseline_status == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	try {
		GetBestFitBaselineCoefficientsFloat(context, num_data, data, mask,
				clip_threshold_sigma, num_fitting_max, num_coeff, coeff,
				final_mask, baseline_status);
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
		double coeff[/*4*num_piece*/], bool final_mask[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	if (context == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data != context->num_basis_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data < context->num_bases)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (clip_threshold_sigma <= 0.0f)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (coeff == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(coeff)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (final_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(final_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (baseline_status == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	try {
		GetBestFitBaselineCoefficientsCubicSplineFloat(context, num_data, data,
				mask, clip_threshold_sigma, num_fitting_max, num_pieces, coeff,
				final_mask, baseline_status);
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

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetNumberOfCoefficients)(
LIBSAKURA_SYMBOL(BaselineContext) const *context, uint16_t const order,
		size_t *num_coeff) noexcept {
	if (context == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	size_t num_out = GetNumberOfBasesFromOrder(context->baseline_type, order);
	if (num_out > context->num_bases)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	*num_coeff = num_out;
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineFloat)(
LIBSAKURA_SYMBOL(BaselineContext) const *context, uint16_t const order,
		size_t num_data, float const data[], bool const mask[],
		float clip_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual, bool final_mask[], float out[],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	if (context == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	size_t num_bases_from_order(
			GetNumberOfBasesFromOrder(context->baseline_type, order));
	if (num_bases_from_order > context->num_bases)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data < num_bases_from_order)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data != context->num_basis_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (clip_threshold_sigma <= 0.0f)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (final_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(final_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (baseline_status == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	try {
		SubtractBaselineFloat(context, order, num_data, data, mask,
				clip_threshold_sigma, num_fitting_max, get_residual, final_mask,
				out, baseline_status);
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
		size_t num_data, float const data[], bool const mask[],
		float clip_threshold_sigma, uint16_t num_fitting_max, bool get_residual,
		bool final_mask[], float out[],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	if (context == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data < context->num_bases)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data != context->num_basis_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (clip_threshold_sigma <= 0.0f)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (final_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(final_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (baseline_status == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	try {
		SubtractBaselineCubicSplineFloat(context, num_pieces, num_data, data,
				mask, clip_threshold_sigma, num_fitting_max, get_residual,
				final_mask, out, baseline_status);
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

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselinePolynomialFloat)(
		size_t num_data, float const data[], bool const mask[], uint16_t order,
		float clip_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual, bool final_mask[], float out[],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	if (num_data <= order)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (clip_threshold_sigma <= 0.0f)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (final_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(final_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (baseline_status == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	try {
		SubtractBaselinePolynomialFloat(num_data, data, mask, order,
				clip_threshold_sigma, num_fitting_max, get_residual, final_mask,
				out, baseline_status);
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
		float const data[], size_t num_coeff, double const coeff[], float out[])
				noexcept {
	if (context == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data != context->num_basis_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data < num_coeff)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_coeff == 0)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_coeff > context->num_bases)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (coeff == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(coeff)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	try {
		SubtractBaselineUsingCoefficientsFloat(context, num_data, data,
				num_coeff, coeff, out);
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
		float const data[], size_t num_pieces, double const coeff[],
		double const boundary[], float out[]) noexcept {
	if (context == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data != context->num_basis_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_data < context->num_bases)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!(num_pieces <= INT_MAX))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (coeff == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(coeff)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (boundary == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(boundary)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

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

