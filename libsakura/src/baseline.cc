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
#include "baseline.h"

namespace {

auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("baseline");
constexpr size_t kNumBasesCubicSpline = 4;

template<typename T>
inline void AllocateMemoryForBasisData(T *context) {
	size_t num_total_basis_data = context->num_bases * context->num_basis_data;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> work_basis_data_storage(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*context->basis_data) * num_total_basis_data,
					&context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	context->basis_data_storage = work_basis_data_storage.release();
}

/**
 * Check if elements of the input @a array are
 * (1) not duplicate and
 * (2) stored in ascending order.
 *
 * @param[in] num_array The number of elements in @a array
 * @param[in] array The input array. Its length must be @a num_array .
 * @return true if the elements of @a array are not duplicate, and also
 * if they are stored in ascending order, otherwise false.
 */
inline bool IsUniqueAndAscendingOrder(size_t num_array, size_t const *array) {
	bool res = true;
	for (size_t i = 1; i < num_array; ++i) {
		if (array[i - 1] >= array[i]) {
			res = false;
			break;
		}
	}
	return res;
}

/**
 * Returns the number of steady bases, which are computable
 * using overall baseline parameter @a order only and don't
 * depend on other conditions specific to each data, for
 * example, mask information etc.
 *
 * @param[in] baseline_type Type of baseline. It must be one
 * of the elements in sakura_BaselineType.
 * @param[in] order Parameter of baseline. It is the maximum
 * polynomial order if baseline_type is
 * BaselineTypeInternal_kPolynomial or
 * BaselineTypeInternal_kChebyshev, or the maximum wave number
 * if baseline_type is BaselineTypeInternal_kSinusoid.
 * @return The returned number should be one of the followings:
 * ( @a order +1) for Polynomial and Chebyshev Polynomial,
 * ( 4 ) for Cubic Spline and
 * ( 2 * @a order +1) for Sinusoid.
 */
inline size_t GetNumberOfContextBases(BaselineTypeInternal const baseline_type,
		uint16_t const order) {
	size_t num_bases = 0;
	switch (baseline_type) {
	case BaselineTypeInternal_kPolynomial:
	case BaselineTypeInternal_kChebyshev:
		num_bases = order + 1;
		break;
	case BaselineTypeInternal_kCubicSpline:
		num_bases = kNumBasesCubicSpline;
		break;
	case BaselineTypeInternal_kSinusoid:
		num_bases = 2 * order + 1;
		break;
	default:
		assert(false);
		break;
	}
	return num_bases;
}

/**
 * Returns the number of all bases for the given baseline type.
 *
 * @param[in] baseline_type Type of baseline. It must be one
 * of the elements in sakura_BaselineType.
 * @param[in] order Parameter of baseline. It is the maximum
 * polynomial order if baseline_type is
 * BaselineTypeInternal_kPolynomial or
 * BaselineTypeInternal_kChebyshev, or the number of spline
 * pieces if baseline_type is BaselineTypeInternal_kCubicSpline,
 * or the maximum wave number if baseline_type is
 * BaselineTypeInternal_kSinusoid.
 * @return The returned number should be one of the followings:
 * ( @a order +1) for Polynomial and Chebyshev Polynomial,
 * ( @a order +3 ) for Cubic Spline and
 * ( 2 * @a order +1) for Sinusoid.
 */
inline size_t GetNumberOfLsqBases(BaselineTypeInternal const baseline_type,
		size_t const order) {
	size_t num_bases = 0;
	switch (baseline_type) {
	case BaselineTypeInternal_kPolynomial:
	case BaselineTypeInternal_kChebyshev:
	case BaselineTypeInternal_kSinusoid:
		num_bases = GetNumberOfContextBases(baseline_type, order);
		break;
	case BaselineTypeInternal_kCubicSpline:
		num_bases = kNumBasesCubicSpline - 1 + order;
		break;
	default:
		assert(false);
		break;
	}
	return num_bases;
}

/**
 * Returns the number of baseline coefficients.
 *
 * @param[in] baseline_type Type of baseline. It must be one
 * of the elements in sakura_BaselineType.
 * @param[in] order It is used as the maximum polynomial order
 * if baseline_type is BaselineTypeInternal_kPolynomial or
 * BaselineTypeInternal_kChebyshev, or it is used as the number
 * of spline pieces if baseline_type is
 * BaselineTypeInternal_kCubicSpline.
 * @param[in] num_nwave The length of @a nwave .
 * @param[in] nwave The maximum wave number. It is used only
 * if baseline_type is BaselineTypeInternal_kSinusoid.
 * @return The returned number should be one of the followings:
 * ( @a order +1) for Polynomial and Chebyshev Polynomial,
 * ( @a order *4 ) for Cubic Spline,
 * ( 2 * @a num_nwave -1) for Sinusoid with @a nwave containing '0', or
 * ( 2 * @a num_nwave ) for Sinusoid with @a nwave not containing '0'.
 */
inline size_t DoGetNumberOfCoefficients(
		BaselineTypeInternal const baseline_type, uint16_t const order,
		size_t const num_nwave, size_t const *nwave) {
	size_t num_bases = 0;
	switch (baseline_type) {
	case BaselineTypeInternal_kPolynomial:
	case BaselineTypeInternal_kChebyshev:
		num_bases = order + 1;
		break;
	case BaselineTypeInternal_kCubicSpline:
		if (order < 1) {
			throw std::invalid_argument(
					"order (number of pieces) must be a positive value!");
		}
		num_bases = kNumBasesCubicSpline * order;
		break;
	case BaselineTypeInternal_kSinusoid:
		if (!IsUniqueAndAscendingOrder(num_nwave, nwave)) {
			throw std::invalid_argument(
					"nwave elements must be in ascending order and not duplicate.");
		}
		num_bases = (nwave[0] == 0) ? (2 * num_nwave - 1) : (2 * num_nwave);
		break;
	default:
		assert(false);
		break;
	}
	return num_bases;
}

/**
 * The function to compute polynomial basis data from order of zero
 * up to the maximum one (@num_bases -1) at a given x-position @a i_d ,
 * then set them in the specific range (from index @idx to ( @idx +
 * @num_bases -1)) in an array @out_arg , i.e., @a 1.0 is set in
 * @a out_arg[idx] , @a i_d in @a out_arg[idx+1] , @a (i_d*i_d) in
 * @a out_arg[idx_2], ..., and @a (i_d^(num_bases-1)) in
 * @a out_arg[idx+num_bases-1] .
 *
 * @tparam U Type of array elements.
 * @param[in] num_bases The number of polynomial bases. Thus the maximum
 * order of polynomial should be ( @a num_bases -1).
 * @param[in] i_d The x-position for which polynomial basis values
 * are to be computed.
 * @param[in,out] idx The starting index of @a out_arg where the first
 * basis value is set. @a idx gains @a num_bases when this function ends.
 * @param[out] out_arg The 1-dimensional array in which polynomial basis
 * values at @i_d are set.
 */
template<typename U>
inline void DoSetBasisDataPolynomial(size_t num_bases, U const i_d, size_t *idx,
		U *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto out = AssumeAligned(out_arg);

	size_t i = *idx;
	U val = static_cast<U>(1.0);
	for (size_t j = 0; j < num_bases; ++j) {
		out[i++] = val;
		val *= i_d;
	}
	*idx = i;
}

/**
 * The function to compute polynomial basis data and set them
 * in the given baseline context. The maximum polynomial order
 * is taken from the given baseline context.
 *
 * @tparam T Type of baseline context.
 * @tparam U Type of basis data.
 */
template<typename T, typename U>
inline void SetBasisDataPolynomial(T *context) {
	assert(0 < context->num_basis_data);
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	auto data = AssumeAligned(context->basis_data);

	size_t num_basis_data = context->num_basis_data;
	size_t num_bases = context->num_bases;
	size_t idx = 0;
	for (size_t i = 0; i < num_basis_data; ++i) {
		DoSetBasisDataPolynomial<U>(num_bases, static_cast<U>(i), &idx, data);
	}
}

/**
 * The function to compute Chebyshev polynomial basis data
 * and set them in the given baseline context. The maximum
 * order is taken from the given baseline context.
 *
 * @tparam T Type of baseline context.
 * @tparam U Type of basis data.
 */
template<typename T, typename U>
inline void SetBasisDataChebyshev(T *context) {
	assert(0 < context->num_basis_data);
	assert(0 < context->num_bases);
	assert(context->num_bases <= context->num_basis_data);
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	auto data = AssumeAligned(context->basis_data);

	size_t num_basis_data = context->num_basis_data;
	size_t num_bases = context->num_bases;
	size_t idx = 0;
	if (num_bases == 1) { // (order == 0) or (num_basis_data == 1)
		for (size_t i = 0; i < num_basis_data; ++i) {
			data[idx++] = 1.0;
		}
	} else {
		U max_data_x = static_cast<U>(num_basis_data - 1);
		assert(0.0 < max_data_x);
		for (size_t i = 0; i < num_basis_data; ++i) {
			data[idx++] = 1.0;
			U x = 2.0 * static_cast<U>(i) / max_data_x - 1.0;
			data[idx++] = x;
			for (size_t j = 2; j < num_bases; ++j) {
				data[idx] = 2.0 * x * data[idx - 1] - data[idx - 2];
				++idx;
			}
		}
	}
}

/**
 * The function to compute sinusoidal basis data and set them
 * in the given baseline context. The maximum wave number is
 * taken from the given baseline context. The value of the
 * @a j -th basis at the @a i -th position is stored at
 * @a context->basis_data[ @a i * ( @a context->num_basis_data
 * ) + @a j ]. As for the order of basis functions, constant
 * (corresponding to zero wave number) comes first, then
 * followed by sine and cosine for wave number of 1, sine and
 * cosine for wave number of 2, and so on.
 *
 * @tparam T Type of baseline context.
 * @tparam U Type of basis data.
 */
template<typename T, typename U>
inline void SetBasisDataSinusoid(T *context) {
	assert(0 < context->num_basis_data);
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	auto data = AssumeAligned(context->basis_data);

	size_t num_basis_data = context->num_basis_data;
	size_t max_nwave = context->baseline_param;
	size_t idx = 0;
	if (num_basis_data == 1) {
		data[idx] = 1.0;
	} else {
		U max_data_x = static_cast<U>(num_basis_data - 1);
		U norm = 2.0 * M_PI / max_data_x;
		for (size_t i = 0; i < num_basis_data; ++i) {
			data[idx++] = 1.0;
			for (size_t nwave = 1; nwave <= max_nwave; ++nwave) {
				U x = static_cast<U>(i) * nwave * norm;
				data[idx++] = sin(x);
				data[idx++] = cos(x);
			}
		}
	}
}

/**
 * It allocates memory space to be used as temporary working
 * area for baseline-related calculation and set them the
 * members of the given baseline context.
 * @param[in,out] context The baseline context.
 */
template<typename T>
inline void AllocateWorkSpaces(T *context) {
	auto const type = context->baseline_type;
	context->num_lsq_bases_max = GetNumberOfLsqBases(type,
			context->baseline_param);
	size_t num_lsq_matrix = context->num_lsq_bases_max
			* context->num_lsq_bases_max;
	context->lsq_matrix = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_matrix(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*context->lsq_matrix) * num_lsq_matrix,
					&context->lsq_matrix));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->lsq_matrix));
	context->lsq_matrix_storage = storage_for_lsq_matrix.release();
	size_t num_lsq_vector = context->num_lsq_bases_max;
	context->lsq_vector = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_vector(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*context->lsq_vector) * num_lsq_vector,
					&context->lsq_vector));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->lsq_vector));
	context->lsq_vector_storage = storage_for_lsq_vector.release();
	context->clipped_indices = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_clipped_indices(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*context->clipped_indices) * context->num_basis_data,
					&context->clipped_indices));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->clipped_indices));
	context->clipped_indices_storage = storage_for_clipped_indices.release();

	context->best_fit_model = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_best_fit_model(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*context->best_fit_model) * context->num_basis_data,
					&context->best_fit_model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->best_fit_model));
	context->best_fit_model_storage = storage_for_best_fit_model.release();
	context->residual_data = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_residual_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*context->residual_data) * context->num_basis_data,
					&context->residual_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->residual_data));
	context->residual_data_storage = storage_for_residual_data.release();

	context->use_bases_idx = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_use_bases_idx(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*context->use_bases_idx)
							* context->num_lsq_bases_max,
					&context->use_bases_idx));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->use_bases_idx));
	if (type != BaselineTypeInternal_kSinusoid) {
		for (size_t i = 0; i < context->num_lsq_bases_max; ++i) {
			context->use_bases_idx[i] = i;
		}
	}
	context->use_bases_idx_storage = storage_for_use_bases_idx.release();

	size_t num_coeff_full_max = context->num_lsq_bases_max;
	if (type == BaselineTypeInternal_kCubicSpline) {
		num_coeff_full_max = DoGetNumberOfCoefficients(type,
				context->baseline_param, 0, nullptr);
	}
	context->coeff_full = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_coeff_full(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*context->coeff_full) * num_coeff_full_max,
					&context->coeff_full));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->coeff_full));
	context->coeff_full_storage = storage_for_coeff_full.release();

	//CubicSpline-specific ones
	context->cspline_basis = nullptr;
	context->cspline_basis_storage = nullptr;
	context->cspline_lsq_coeff = nullptr;
	context->cspline_lsq_coeff_storage = nullptr;
	if (type == BaselineTypeInternal_kCubicSpline) {
		std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_cspline_basis(
				LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
						sizeof(*context->cspline_basis)
								* context->num_lsq_bases_max
								* context->num_basis_data,
						&context->cspline_basis));
		assert(LIBSAKURA_SYMBOL(IsAligned)(context->cspline_basis));
		context->cspline_basis_storage = storage_for_cspline_basis.release();
		std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_cspline_lsq_coeff(
				LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
						sizeof(*context->cspline_lsq_coeff)
								* context->num_lsq_bases_max,
						&context->cspline_lsq_coeff));
		assert(LIBSAKURA_SYMBOL(IsAligned)(context->cspline_lsq_coeff));
		context->cspline_lsq_coeff_storage =
				storage_for_cspline_lsq_coeff.release();
	}
}

/**
 * Set basis data in the given baseline context.
 * @param[in] order The maximum polynomial order for polynomial
 * or Chebyshev polynomial, or the number of spline pieces
 * for cubic spline, or the maximum wave number for sinusoid.
 * @param[in,out] context Pointer to the baseline context.
 */
template<typename T>
inline void SetBasisData(size_t const order, T *context) {
	auto const type = context->baseline_type;
	context->num_bases = GetNumberOfContextBases(type, order);
	size_t min_num_basis_data = GetNumberOfLsqBases(type, order);

	if (context->num_basis_data < min_num_basis_data) {
		throw std::invalid_argument("num_basis_data is too small!");
	}
	AllocateMemoryForBasisData<T>(context);
	switch (type) {
	case BaselineTypeInternal_kPolynomial:
		SetBasisDataPolynomial<T, double>(context);
		break;
	case BaselineTypeInternal_kChebyshev:
		SetBasisDataChebyshev<T, double>(context);
		break;
	case BaselineTypeInternal_kCubicSpline:
		assert(context->num_bases == kNumBasesCubicSpline);
		SetBasisDataPolynomial<T, double>(context);
		break;
	case BaselineTypeInternal_kSinusoid:
		SetBasisDataSinusoid<T, double>(context);
		break;
	default:
		assert(false);
		break;
	}
}

/**
 * Destroy baseline context object.
 * @param[in] context Pointer to baseline context.
 */
template<typename T>
inline void DestroyBaselineContext(T *context) {
	if (context != nullptr) {
		if (context->basis_data_storage != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->basis_data_storage);
		}
		if (context->lsq_matrix_storage != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->lsq_matrix_storage);
		}
		if (context->lsq_vector_storage != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->lsq_vector_storage);
		}
		if (context->clipped_indices != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->clipped_indices_storage);
		}
		if (context->best_fit_model != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->best_fit_model_storage);
		}
		if (context->residual_data != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->residual_data_storage);
		}
		if (context->use_bases_idx != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->use_bases_idx_storage);
		}
		if (context->coeff_full != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->coeff_full_storage);
		}
		if (context->cspline_basis != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->cspline_basis_storage);
		}
		if (context->cspline_lsq_coeff != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(context->cspline_lsq_coeff_storage);
		}
		LIBSAKURA_PREFIX::Memory::Free(context);
	}
}

/**
 * The main function to create baseline context object.
 *
 * @tparam T Type of baseline context.
 * @param[in] baseline_type Type of baseline function. It must be
 * one of those defined in sakura_BaselineType .
 * @param[in] order It is used only if @a baseline_type is
 * BaselineTypeInternal_kPolynomial or BaselineTypeInternal_kChebyshev .
 * @param[in] npiece It is used only if @a baseline_type is
 * BaselineTypeInternal_kCubicSpline .
 * @param[in] nwave It is used only if @a baseline_type is
 * BaselineTypeInternal_kSinusoid .
 * @param[in] num_basis_data Length of the arrays for basis data,
 * input data, etc.
 * @parama[out] Address of pointer to baseline context object.
 */
template<typename T>
inline void CreateBaselineContext(BaselineTypeInternal const baseline_type,
		uint16_t const order, uint16_t const npiece, uint16_t const nwave,
		size_t const num_basis_data, T **context) {
	try {
		std::unique_ptr<T, LIBSAKURA_PREFIX::Memory> work_context(
				static_cast<T *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(T))),
				LIBSAKURA_PREFIX::Memory());
		if (work_context == nullptr) {
			throw std::bad_alloc();
		}
		work_context->basis_data_storage = nullptr;
		work_context->baseline_type = baseline_type;
		work_context->num_basis_data = num_basis_data;
		switch (baseline_type) {
		case BaselineTypeInternal_kPolynomial:
		case BaselineTypeInternal_kChebyshev:
			work_context->baseline_param = order;
			break;
		case BaselineTypeInternal_kCubicSpline:
			work_context->baseline_param = npiece;
			break;
		case BaselineTypeInternal_kSinusoid:
			work_context->baseline_param = nwave;
			break;
		default:
			assert(false);
			break;
		}

		SetBasisData<T>(work_context->baseline_param, work_context.get());
		AllocateWorkSpaces<T>(work_context.get());
		*context = work_context.release();
	} catch (...) {
		DestroyBaselineContext<T>(*context);
		throw;
	}
}

/**
 * Compute difference between two input arrays ( @a in1_arg - @a in2_arg)
 * and set the result in @a out_arg .
 * @tparam T Type of the elements of input/output arrays.
 * @param num_in Length of the input and output arrays.
 */
template<typename T>
inline void OperateSubtraction(size_t num_in, T const *in1_arg,
		T const *in2_arg, T *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in1_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const in1 = AssumeAligned(in1_arg);
	auto const in2 = AssumeAligned(in2_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_in; ++i) {
		out[i] = in1[i] - in2[i];
	}
}

/**
 *
 * @param num_coeff
 * @param coeff
 * @param use_idx
 * @param num_out
 * @param num_bases
 * @param basis
 * @param out
 */
template<typename T, typename U>
inline void AddMulMatrix(size_t num_coeff, U const *coeff,
		size_t const *use_idx, size_t num_out, size_t num_bases, U const *basis,
		T *out) {
	assert(((void)"Not yet implemented", false));
}

template<>
inline void AddMulMatrix<float, double>(size_t num_coeff,
		double const *coeff_arg, size_t const *use_idx_arg, size_t num_out,
		size_t num_bases, double const *basis_arg, float *out_arg) {
	auto coeff = AssumeAligned(coeff_arg);
	auto use_idx = AssumeAligned(use_idx_arg);
	auto basis = AssumeAligned(basis_arg);
	auto out = AssumeAligned(out_arg);

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
			assert(num_bases * i + idx + offset3 < num_bases * num_out);
			auto bs = _mm256_set_pd(bases_row[idx + offset3],
					bases_row[idx + offset2], bases_row[idx + offset1],
					bases_row[idx]);
#endif
			total = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
			LIBSAKURA_SYMBOL(SimdPacketAVX), double>(ce, bs, total);
		}
		//_mm256_store_pd(&out[i], total);
		_mm_store_ps(&out[i], _mm256_cvtpd_ps(total));
	}
#endif
	for (; i < num_out; ++i) {
		double out_double = 0.0;
		for (size_t j = 0; j < num_coeff; ++j) {
			size_t k = num_bases * i + use_idx[j];
			assert(k < num_bases * num_out);
			out_double += coeff[j] * basis[k];
		}
		out[i] = out_double;
	}
}

template<typename T, typename U>
inline void AddMulMatrixCubicSpline(size_t num_boundary, size_t const *boundary,
		U const (*coeff_full)[kNumBasesCubicSpline], size_t num_out,
		U const *basis, T *out) {
	assert(((void)"Not yet implemented", false));
}

template<>
inline void AddMulMatrixCubicSpline<float, double>(size_t num_boundary,
		size_t const *boundary_arg,
		double const (*coeff_full_arg)[kNumBasesCubicSpline], size_t num_out,
		double const *basis_arg, float *out_arg) {
	auto boundary = AssumeAligned(boundary_arg);
	auto coeff_full = AssumeAligned(coeff_full_arg);
	auto basis = AssumeAligned(basis_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_boundary - 1; ++i) {
		size_t start_idx = boundary[i];
		size_t end_idx = boundary[i + 1];
		size_t j = start_idx;
		for (j = start_idx; j < end_idx; ++j) {
#if defined(__AVX__) && !defined(ARCH_SCALAR)
			auto idx = kNumBasesCubicSpline * j;
			auto coeff_packed = _mm256_load_pd(coeff_full[i]);
			auto basis_packed = _mm256_load_pd(&basis[idx]);
			auto out_double_packed = coeff_packed * basis_packed;
			double *out_double = reinterpret_cast<double *>(&out_double_packed);
			out[j] = out_double[0] + out_double[1] + out_double[2]
					+ out_double[3];
#else
			double out_double = 0.0;
			for (size_t k = 0; k < kNumBasesCubicSpline; ++k) {
				size_t l = kNumBasesCubicSpline * j + k;
				assert(l < kNumBasesCubicSpline * num_out);
				out_double += coeff_full[i][k] * basis[l];
			}
			out[j] = out_double;
#endif
		}
	}
}

inline std::string GetNotEnoughDataMessage(
		uint16_t const idx_erroneous_fitting) {
	std::string s;
	s = "SubtractBaseline: available data became too few in the ";
	s += idx_erroneous_fitting;
	s += " ";

	std::string si;
	uint16_t mod100 = idx_erroneous_fitting % 100;
	uint16_t ones_digit = mod100 % 10;
	uint16_t tens_digit = mod100 / 10;
	if (tens_digit == 1) {
		si = "th";
	} else {
		if (ones_digit == 1) {
			si = "st";
		} else if (ones_digit == 2) {
			si = "nd";
		} else if (ones_digit == 3) {
			si = "rd";
		} else {
			si = "th";
		}
	}

	s += si;
	s += " fitting.";

	return s;
}

/**
 * Count the number of data where @a mask_org is true (not masked),
 * then compute the boundary positions @a boundary_arg that divide
 * the unmasked data as evenly as possible into @a num_boundary
 * groups.
 *
 * @tparam U Type of the elements of @a boundary_arg
 * @param[in] num_mask Length of @a mask_arg
 * @param[in] mask_arg Array of mask. If true, the corresponding data
 * with the same index are used for baseline fitting.
 * @param[in] num_boundary Length of @a boundary_arg . It must be
 * the (number of spline pieces + 1).
 * @param[out] boundary_arg Positions of spline piece boundaries. The
 * first element must always point 0, the first x-position of data,
 * and the last element must always point @a num_mask , next of the
 * final x-position of data.
 */
inline void GetBoundariesOfPiecewiseData(size_t num_mask, bool const *mask_arg,
		size_t num_boundary, size_t *boundary_arg) {
	assert(2 <= num_boundary);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	auto const mask = AssumeAligned(mask_arg);
	auto boundary = AssumeAligned(boundary_arg);

	size_t num_unmasked_data = 0;
	for (size_t i = 0; i < num_mask; ++i) {
		if (mask[i])
			++num_unmasked_data;
	}
	size_t boundary_last_idx = num_boundary - 1;
	assert(boundary_last_idx <= num_unmasked_data);
	// the first element of boundary[] must point zero, the first index of mask/data.
	boundary[0] = 0;
	size_t idx = 1;
	size_t count_unmasked_data = 0;
	for (size_t i = 0; i < num_mask; ++i) {
		if (boundary_last_idx <= idx) {
			break;
		}
		if (mask[i]) {
			if (num_unmasked_data * idx
					<= count_unmasked_data * boundary_last_idx) {
				boundary[idx] = i;
				++idx;
			}
			++count_unmasked_data;
		}
	}
	// the last element of boundary[] must point the next of the last index of mask/data.
	boundary[boundary_last_idx] = num_mask;
}

/**
 * Compute the complete set of coefficients of cubic spline
 * curves for all spline pieces.
 *
 * @tparam U Type of coefficient values.
 * @param[in] num_boundary Number of elements of @a boundary_arg.
 * It is (number of spline pieces + 1).
 * @param[in] boundary_arg Boundary positions of spline pieces.
 * The first element is the left end of the first (the left-most)
 * piece, and the last element is the right end of the last (the
 * right-most) piece.
 * @param[in] coeff_raw_arg The result of Least-Square fitting.
 * Its length is (number of spline pieces +3). The first 4
 * elements are the true coefficients of zeroth, first, second and
 * third orders for the left-most piece, and the coefficients for
 * the other pieces can be computed using @a coeff_raw_arg and
 * @a boundary_arg .
 * @param[out] coeff_arg A 2-dimensional array to store
 * the whole cubic spline coefficients. The @a i -th order
 * coefficient for the @a j -th piece from the left side
 * is stored at @a coeff_arg[ @a i ][ @a j ].
 */
template<typename U>
inline void GetFullCubicSplineCoefficients(size_t num_boundary,
		size_t const *boundary_arg, U const *coeff_raw_arg,
		U (*coeff_arg)[kNumBasesCubicSpline]) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_raw_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	auto const boundary = AssumeAligned(boundary_arg);
	auto const coeff_raw = AssumeAligned(coeff_raw_arg);
	auto coeff = AssumeAligned(coeff_arg);

	for (size_t i = 0; i < kNumBasesCubicSpline; ++i) {
		coeff[0][i] = coeff_raw[i];
	}
	U const three = static_cast<U>(3.0);
	for (size_t i = 1; i < num_boundary - 1; ++i) {
		size_t j = GetNumberOfLsqBases(BaselineTypeInternal_kCubicSpline, i);
		auto const c = coeff_raw[j] - coeff[i - 1][3];
		auto const b = static_cast<U>(boundary[i]);
		coeff[i][0] = coeff[i - 1][0] - b * b * b * c;
		coeff[i][1] = coeff[i - 1][1] + three * b * b * c;
		coeff[i][2] = coeff[i - 1][2] - three * b * c;
		coeff[i][3] = coeff_raw[j];
	}
}

/**
 * Compute cubic spline basis data which depend on mask
 * information at a given x-position @a i_d .
 * @tparam U Type of array elements.
 * @param[in] num_boundary Length of @a boundary_arg
 * @param[in] boundary_arg x-positions of boundaries of
 * spline pieces.
 * @param[in] i_d The x-position for which basis values
 * are to be computed.
 * @param[in,out] idx The starting index of @a out_arg
 * where the first basis value is set. @a idx gains
 * ( @a num_boundary -2) or (number of spline pieces -1)
 * when this function ends.
 * @param[out] out_arg The 1-dimensional array in which
 * basis values at @i_d are set.
 */
template<typename U>
inline void SetAuxiliaryCubicBases(size_t const num_boundary,
		size_t const *boundary_arg, U const i_d, size_t *idx, U *out_arg) {
	size_t i = *idx;
	assert(1 <= i);
	assert(i_d <= boundary_arg[num_boundary-1]);
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const boundary = AssumeAligned(boundary_arg);
	auto out = AssumeAligned(out_arg);

	auto cube = [](U v) {return static_cast<U>(v * v * v);};
	auto force_non_negative = [](U v) {return static_cast<U>(std::max(0.0, v));};

	auto cube_prev = cube(i_d - static_cast<U>(boundary[1]));
	out[i - 1] -= force_non_negative(cube_prev);

	for (size_t j = 2; j < num_boundary; ++j) {
		auto cube_current = cube(i_d - static_cast<U>(boundary[j]));
		out[i] = force_non_negative(
				cube_prev - force_non_negative(cube_current));
		cube_prev = cube_current;
		++i;
	}
	*idx = i;
}

/**
 * Compute the whole basis data for cubic spline fitting.
 *
 * @tparam U Type of boundary data.
 * @param[in] num_data Length of a basis data.
 * @param[in] num_boundary Length of @a boundary_arg .
 * @param[in] boundary_arg Boundary positions.
 * @param[out] out_arg A 1D array to store basis data. Its
 * length should be ( @a num_data * ( @a num_boundary + 2)).
 */
template<typename U>
inline void SetFullCubicSplineBasisData(size_t num_data, size_t num_boundary,
		size_t const *boundary_arg, U *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const boundary = AssumeAligned(boundary_arg);
	auto out = AssumeAligned(out_arg);

	size_t idx = 0;
	if (num_boundary <= 2) {
		for (size_t i = 0; i < num_data; ++i) {
			U i_d = static_cast<U>(i);
			DoSetBasisDataPolynomial<U>(kNumBasesCubicSpline, i_d, &idx, out);
		}
	} else {
		for (size_t i = 0; i < num_data; ++i) {
			U i_d = static_cast<U>(i);
			DoSetBasisDataPolynomial<U>(kNumBasesCubicSpline, i_d, &idx, out);
			SetAuxiliaryCubicBases<U>(num_boundary, boundary, i_d, &idx, out);
		}
	}
}

/**
 * Compute the best-fit baseline model using basis data and
 * coefficients, then subtract it from input data.
 *
 * @tparam T Type of input data, best-fit baseline and
 * residual data.
 * @tparam U Type of coefficients.
 * @tparam V Type of baseline context.
 * @param[in] num_data Length of @a data , @a best_fit_model
 * and @a residual_data .
 * @param[in] data Input data.
 * @param[in] context Baseline context.
 * @param[in] num_coeff Length of @a coeff .
 * @param[in] coeff Coefficients.
 * @param[out] best_fit_model The best-fit baseline data.
 * @param[out] residual_data Residual data.
 */
template<typename T, typename U, typename V>
inline void GetBestFitModelAndResidual(size_t num_data, T const *data,
		V const *context, size_t num_coeff, U const *coeff, T *best_fit_model,
		T *residual_data) {
	AddMulMatrix<T, U>(num_coeff, coeff, context->use_bases_idx, num_data,
			context->num_bases, context->basis_data, best_fit_model);
	OperateSubtraction<T>(num_data, data, best_fit_model, residual_data);
}

/**
 * Compute the best-fit cubic spline baseline using basis
 * data, coefficients and boundary positions, then subtract
 * it from input data.
 *
 * @tparam T Type of input data, best-fit baseline and
 * residual data.
 * @tparam U Type of coefficients.
 * @tparam V Type of baseline context.
 * @param[in] num_data Length of @a data , @a best_fit_model
 * and @a residual_data .
 * @param[in] data Input data.
 * @param[in] context Baseline context.
 * @param[in] num_boundary Number of elements of @a boundary .
 * @param[in] boundary Positions of Spline boundaries.
 * @param[in] coeff_full Coefficients. It is a 2D array
 * and the coefficients for @a i -th piece must be stored at
 * @a coeff_full[i][] .
 * @param best_fit_model[out] The best-fit baseline data.
 * @param residual_data[out] Residual data.
 */
template<typename T, typename U, typename V>
inline void GetBestFitModelAndResidualCubicSpline(size_t num_data,
		T const *data, V const *context, size_t num_boundary,
		size_t const *boundary,
		double const (*coeff_full)[kNumBasesCubicSpline], T *best_fit_model,
		T *residual_data) {
	assert(context->baseline_type == BaselineTypeInternal_kCubicSpline);
	assert(context->num_bases == kNumBasesCubicSpline);
	AddMulMatrixCubicSpline<T, U>(num_boundary, boundary, coeff_full, num_data,
			context->basis_data, best_fit_model);
	OperateSubtraction<T>(num_data, data, best_fit_model, residual_data);
}

/**
 * Mask data where the value of @a data_arg is outside the range between
 * @a lower_bound and @a upper_bound .
 *
 * @tparam T Type of data and lower/upper limits.
 * @param[in] num_boundary Number of elements of @a boundary_arg .
 * It must be 2 for baseline types other than Cubic Spline.
 * @param[in] boundary_arg Boundary positions. Its length must be
 * (number of spline pieces +1).
 * @param[in] data_arg Input data.
 * @param[in] in_mask_arg Input mask.
 * @param[in] lower_bound Lower limit of allowed area of @a data_arg values.
 * @param[in] upper_bound Upper limit of allowed area of @a data_arg values.
 * @param[out] out_mask_arg Output mask with positions where @a data_arg is
 * outside the allowed range set false (masked).
 * @param[out] clipped_indices_arg An array to store positions where
 * @a data_arg is outside the allowed range. For safety, its length should
 * not be less than @a boundary_arg[num_boundary-1] .
 * @param[out] num_clipped Number of positions where @a data_arg is outside
 * the allowed range.
 */
template<typename T>
inline void ClipData(size_t num_boundary, size_t *boundary_arg,
		T const *data_arg, bool const *in_mask_arg, T const lower_bound,
		T const upper_bound, bool *out_mask_arg, size_t *clipped_indices_arg,
		size_t *num_clipped) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	auto const boundary = AssumeAligned(boundary_arg);
	auto const data = AssumeAligned(data_arg);
	auto const in_mask = AssumeAligned(in_mask_arg);
	auto out_mask = AssumeAligned(out_mask_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);

	size_t num_clipped_tmp = 0;
	size_t piece_start = 0;
	for (size_t ibound = 1; ibound < num_boundary; ++ibound) {
		size_t piece_end = boundary[ibound] - 1;
		for (size_t i = piece_start; i < piece_end; ++i) {
			bool in_mask_i = in_mask[i];
			bool out_mask_i = in_mask_i;
			if (in_mask_i) {
				auto data_i = data[i];
				if ((data_i - lower_bound) * (upper_bound - data_i) < 0.0) {
					out_mask_i = false;
					clipped_indices[num_clipped_tmp++] = i;
				}
			}
			out_mask[i] = out_mask_i;
		}
		piece_start = piece_end + 1;
	}
	*num_clipped = num_clipped_tmp;
}

/**
 * This function contains the algorithm of baseline fitting itself,
 * and is called from its wrapper functions SubtractBaseline(),
 * SubtractBaselineCubicSpline(), GetBestFitBaselineCoefficients()
 * or GetBestFitBaselineCoefficientsCubicSpline().
 *
 * @tparam T Type of input/output data, best-fit baseline data,
 * residual data, rms of residual data, and threshold level of
 * recursive clipping.
 * @tparam U Type of basis data of baseline model, best-fit
 * coefficients, and boundary data used for cubic spline fitting.
 * @tparam V Type of baseline context.
 * @param[out] out_arg If @a get_residual is true, it is residual
 * (baseline-subtracted) data, otherwise it is the best-fit
 * baseline.
 */
template<typename T, typename U, typename V, typename Func>
inline void DoFitBaseline(V const *context, size_t num_data, T const *data_arg,
bool const *mask_arg, size_t num_context_bases, size_t num_coeff,
		U const *basis_arg, size_t num_boundary, size_t *boundary_arg,
		uint16_t num_fitting_max, T clip_threshold_sigma, bool get_residual,
		U *coeff_arg, T *out_arg, bool *final_mask_arg, T *rms,
		T *residual_data_arg, T *best_fit_model_arg, Func func,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(basis_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(residual_data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(best_fit_model_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto const basis = AssumeAligned(basis_arg);
	auto const boundary = AssumeAligned(boundary_arg);
	auto coeff = AssumeAligned(coeff_arg);
	auto out = AssumeAligned(out_arg);
	auto final_mask = AssumeAligned(final_mask_arg);
	auto residual_data = AssumeAligned(residual_data_arg);
	auto best_fit_model = AssumeAligned(best_fit_model_arg);

	size_t num_unmasked_data = num_data;
	size_t num_clipped = num_data;
	for (size_t i = 0; i < num_data; ++i) {
		final_mask[i] = mask[i];
	}
	LIBSAKURA_SYMBOL(Status) status;

	U rms_d = 0.0;
	assert(0.0 < clip_threshold_sigma);
	for (uint16_t i = 1; i <= num_fitting_max; ++i) {
		if (num_unmasked_data < num_coeff) {
			*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kNotEnoughData);
			throw std::runtime_error(GetNotEnoughDataMessage(i));
		}
		if (num_unmasked_data <= num_clipped) {
			status = LIBSAKURA_SYMBOL(GetLSQCoefficientsDouble)(num_data, data,
					final_mask, num_context_bases, basis, num_coeff,
					context->use_bases_idx, context->lsq_matrix,
					context->lsq_vector);
			if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
				throw std::runtime_error("failed in GetLSQCoefficients.");
			}
		} else {
			status = LIBSAKURA_SYMBOL(UpdateLSQCoefficientsDouble)(num_data,
					data, mask, num_clipped, context->clipped_indices,
					num_context_bases, basis, num_coeff, context->use_bases_idx,
					context->lsq_matrix, context->lsq_vector);
			if (status != LIBSAKURA_SYMBOL(Status_kOK)) {
				throw std::runtime_error("failed in UpdateLSQCoefficients.");
			}
		}
		status = LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLUDouble)(
				num_coeff, context->lsq_matrix, context->lsq_vector, coeff);
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
		U mean = result.sum / result.count;
		rms_d = std::sqrt(
				std::abs(result.square_sum / result.count - mean * mean));

		if (i < num_fitting_max) {
			float clip_threshold_abs = clip_threshold_sigma * rms_d;
			float clip_threshold_lower = mean - clip_threshold_abs;
			float clip_threshold_upper = mean + clip_threshold_abs;
			ClipData<T>(num_boundary, boundary, residual_data, final_mask,
					clip_threshold_lower, clip_threshold_upper, final_mask,
					context->clipped_indices, &num_clipped);
			if (num_clipped == 0) {
				break;
			}
			num_unmasked_data = result.count - num_clipped;
		}
	}
	if (out != nullptr) {
		auto src = AssumeAligned(
				(num_fitting_max == 0) ?
						data : (get_residual ? residual_data : best_fit_model));
		std::copy(src, src + num_data, out);
	}
	*rms = rms_d;
}

/**
 * Convert user-given wave number set used for sinusoidal
 * fitting to the array of indices to access basis data
 * stored in baseline context, then set the array into
 * the baseline context itself.
 *
 * @tparam T Type of baseline context.
 * @param[in] num_nwave Length of nwave.
 * @param[in] nwave User-given set of wave numbers that
 * are used for least-square fitting actually.
 * @param[out] context Baseline context in which the
 * array of indices to access basis data are stored.
 */
template<typename T>
inline void SetSinusoidUseBasesIndex(size_t const num_nwave,
		size_t const *nwave, T *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->use_bases_idx));
	auto use_bases_idx = AssumeAligned(context->use_bases_idx);

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

/**
 * A wrapper function of DoSubtractBaseline() to fit/subtract
 * baseline for polynomial, Chebyshev polynomial and sinusoids.
 * It is called from its C interface sakura_SubtractBaselineFloat().
 *
 * @tparam T Type of input/output data, best-fit baseline data,
 * residual data, rms of residual data, and threshold level of
 * recursive clipping.
 * @tparam U Type of basis data of baseline model, best-fit
 * coefficients, and boundary data used for cubic spline fitting.
 * @tparam V Type of baseline context.
 * @param[out] out_arg If @a get_residual is true, it is residual
 * (baseline-subtracted) data, otherwise it is the best-fit
 * baseline.
 */
template<typename T, typename U, typename V>
inline void SubtractBaseline(V const *context, uint16_t order,
		size_t const num_nwave, size_t const *nwave, size_t num_data,
		T const *data_arg, bool const *mask_arg, T clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual, T *out_arg,
		bool *final_mask_arg, T *rms,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	auto const type = context->baseline_type;
	assert(
			(type == BaselineTypeInternal_kPolynomial)||(type == BaselineTypeInternal_kChebyshev)||(type == BaselineTypeInternal_kSinusoid));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto out = AssumeAligned(out_arg);
	auto final_mask = AssumeAligned(final_mask_arg);

	if (type == BaselineTypeInternal_kSinusoid) {
		SetSinusoidUseBasesIndex<V>(num_nwave, nwave, const_cast<V *>(context));
	}
	size_t const num_coeff = DoGetNumberOfCoefficients(type, order, num_nwave,
			nwave);
	size_t const num_boundary = 2;
	SIMD_ALIGN
	size_t boundary[num_boundary] = { 0, num_data };

	DoFitBaseline<T, U, V>(context, num_data, data, mask, context->num_bases,
			num_coeff, context->basis_data, num_boundary, boundary,
			num_fitting_max, clip_threshold_sigma, get_residual,
			context->coeff_full, out, final_mask, rms, context->residual_data,
			context->best_fit_model,
			[&]() {GetBestFitModelAndResidual<T, U, V>(num_data, data, context,
						num_coeff, context->coeff_full, context->best_fit_model,
						context->residual_data);}, baseline_status);
}

/**
 * A wrapper function of DoSubtractBaseline() to fit/subtract
 * baseline for cubic spline.
 * It is called from its C interface sakura_SubtractBaselineCubicSplineFloat().
 *
 * @tparam T Type of input/output data, best-fit baseline data,
 * residual data, rms of residual data, and threshold level of
 * recursive clipping.
 * @tparam U Type of basis data of baseline model, best-fit
 * coefficients, and boundary data used for cubic spline fitting.
 * @tparam V Type of baseline context.
 * @param[out] out_arg If @a get_residual is true, it is residual
 * (baseline-subtracted) data, otherwise it is the best-fit
 * baseline.
 */
template<typename T, typename U, typename V>
inline void SubtractBaselineCubicSpline(V const *context, size_t num_boundary,
		size_t num_data, T const *data_arg, bool const *mask_arg,
		T clip_threshold_sigma, uint16_t num_fitting_max, bool get_residual,
		T *out_arg, bool *final_mask_arg, float *rms, size_t *boundary_arg,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(context->baseline_type == BaselineTypeInternal_kCubicSpline);
	assert(context->num_bases == kNumBasesCubicSpline);
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto out = AssumeAligned(out_arg);
	auto final_mask = AssumeAligned(final_mask_arg);
	auto boundary = AssumeAligned(boundary_arg);

	GetBoundariesOfPiecewiseData(num_data, mask, num_boundary, boundary);
	SetFullCubicSplineBasisData<U>(num_data, num_boundary, boundary,
			context->cspline_basis);
	size_t num_coeff = GetNumberOfLsqBases(BaselineTypeInternal_kCubicSpline,
			num_boundary - 1);
	auto coeff_full =
			reinterpret_cast<U (*)[kNumBasesCubicSpline]>(context->coeff_full);

	DoFitBaseline<T, U, V>(context, num_data, data, mask, num_coeff, num_coeff,
			context->cspline_basis, num_boundary, boundary, num_fitting_max,
			clip_threshold_sigma, get_residual, context->cspline_lsq_coeff, out,
			final_mask, rms, context->residual_data, context->best_fit_model,
			[&]() {
				GetFullCubicSplineCoefficients(num_boundary, boundary, context->cspline_lsq_coeff, coeff_full);
				GetBestFitModelAndResidualCubicSpline<T, U, V>(num_data, data, context,
						num_boundary, boundary, coeff_full, context->best_fit_model, context->residual_data);
			}, baseline_status);
}

/**
 * A wrapper function of DoSubtractBaseline() to obtain best-fit
 * baseline coefficients for polynomial, Chebyshev polynomial
 * and sinusoids.
 * It is called from its C interface sakura_GetBestFitBaselineCoefficientsFloat().
 *
 * @tparam T Type of input/output data, best-fit baseline data,
 * residual data, rms of residual data, and threshold level of
 * recursive clipping.
 * @tparam U Type of basis data of baseline model, best-fit
 * coefficients, and boundary data used for cubic spline fitting.
 * @tparam V Type of baseline context.
 */
template<typename T, typename U, typename V>
inline void GetBestFitBaselineCoefficients(V const *context, size_t num_data,
		T const *data_arg, bool const *mask_arg, T clip_threshold_sigma,
		uint16_t num_fitting_max, size_t num_nwave, size_t const *nwave,
		size_t num_coeff, U *coeff_arg,
		bool *final_mask_arg, T *rms,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	auto const type = context->baseline_type;
	assert(
			(type == BaselineTypeInternal_kPolynomial)||(type == BaselineTypeInternal_kChebyshev)||(type == BaselineTypeInternal_kSinusoid));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto coeff = AssumeAligned(coeff_arg);
	auto final_mask = AssumeAligned(final_mask_arg);

	if (type == BaselineTypeInternal_kSinusoid) {
		SetSinusoidUseBasesIndex<V>(num_nwave, nwave, const_cast<V *>(context));
	}
	size_t const num_boundary = 2;
	SIMD_ALIGN
	size_t boundary[num_boundary] = { 0, num_data };

	DoFitBaseline<T, U, V>(context, num_data, data, mask, context->num_bases,
			num_coeff, context->basis_data, num_boundary, boundary,
			num_fitting_max, clip_threshold_sigma, true, coeff, nullptr,
			final_mask, rms, context->residual_data, context->best_fit_model,
			[&]() {GetBestFitModelAndResidual<T, U, V>(num_data, data, context,
						num_coeff, coeff, context->best_fit_model,
						context->residual_data);}, baseline_status);
}

/**
 * A wrapper function of DoSubtractBaseline() to obtain best-fit
 * baseline coefficients for cubic spline.
 * It is called from its C interface
 * sakura_GetBestFitBaselineCoefficientsCubicSplineFloat().
 *
 * @tparam T Type of input/output data, best-fit baseline data,
 * residual data, rms of residual data, and threshold level of
 * recursive clipping.
 * @tparam U Type of basis data of baseline model, best-fit
 * coefficients, and boundary data used for cubic spline fitting.
 * @tparam V Type of baseline context.
 */
template<typename T, typename U, typename V>
inline void GetBestFitBaselineCoefficientsCubicSpline(V const *context,
		size_t num_data, T const *data_arg, bool const *mask_arg,
		T clip_threshold_sigma, uint16_t num_fitting_max, size_t num_boundary,
		U (*coeff_arg)[kNumBasesCubicSpline], bool *final_mask_arg, T *rms,
		size_t *boundary_arg,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(context->baseline_type == BaselineTypeInternal_kCubicSpline);
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

	GetBoundariesOfPiecewiseData(num_data, mask, num_boundary, boundary);
	SetFullCubicSplineBasisData<U>(num_data, num_boundary, boundary,
			context->cspline_basis);
	size_t num_coeff = GetNumberOfLsqBases(BaselineTypeInternal_kCubicSpline,
			num_boundary - 1);

	DoFitBaseline<T, U, V>(context, num_data, data, mask, num_coeff, num_coeff,
			context->cspline_basis, num_boundary, boundary, num_fitting_max,
			clip_threshold_sigma, true, context->cspline_lsq_coeff, nullptr,
			final_mask, rms, context->residual_data, context->best_fit_model,
			[&]() {
				GetFullCubicSplineCoefficients(num_boundary, boundary, context->cspline_lsq_coeff, coeff);
				GetBestFitModelAndResidualCubicSpline<T, U, V>(num_data, data, context,
						num_boundary, boundary, coeff, context->best_fit_model, context->residual_data);
			}, baseline_status);
}

/**
 * Subtract baseline using basis data in baseline context and
 * user-given coefficients. This function is for polynomial,
 * Chebyshev polynomial and sinusoid.
 * It is called from C interface
 * sakura_SubtractBaselineUsingCoefficientsFloat().
 *
 * @tparam T Type of input/output data.
 * @tparam U Type of coefficients.
 * @tparam V Type of baseline context.
 * @param[out] out_arg Baseline-subtracted data.
 */
template<typename T, typename U, typename V>
inline void SubtractBaselineUsingCoefficients(V const *context, size_t num_data,
		T const *data_arg, size_t num_coeff, U const *coeff_arg,
		size_t num_nwave, size_t const *nwave, T *out_arg) {
	auto const type = context->baseline_type;
	assert(
			(type == BaselineTypeInternal_kPolynomial)||(type == BaselineTypeInternal_kChebyshev)||(type == BaselineTypeInternal_kSinusoid));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const data = AssumeAligned(data_arg);
	auto const coeff = AssumeAligned(coeff_arg);
	auto out = AssumeAligned(out_arg);

	if (type == BaselineTypeInternal_kSinusoid) {
		SetSinusoidUseBasesIndex<V>(num_nwave, nwave, const_cast<V *>(context));
	}
	GetBestFitModelAndResidual<T, U, V>(num_data, data, context, num_coeff,
			coeff, context->best_fit_model, out);
}

/**
 * Subtract baseline using basis data in baseline context and
 * user-given coefficients. This function is for cubic spline.
 * It is called from C interface
 * sakura_SubtractBaselineCubicSplineUsingCoefficientsFloat().
 *
 * @tparam T Type of input/output data.
 * @tparam U Type of coefficients.
 * @tparam V Type of baseline context.
 * @param[out] out_arg Baseline-subtracted data.
 */
template<typename T, typename U, typename V>
inline void SubtractBaselineCubicSplineUsingCoefficients(V const *context,
		size_t num_data, T const *data_arg, size_t num_boundary,
		U const (*coeff_arg)[kNumBasesCubicSpline], size_t const *boundary_arg,
		T *out_arg) {
	assert(context->baseline_type == BaselineTypeInternal_kCubicSpline);
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(boundary_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const data = AssumeAligned(data_arg);
	auto const coeff = AssumeAligned(coeff_arg);
	auto const boundary = AssumeAligned(boundary_arg);
	auto out = AssumeAligned(out_arg);

	GetBestFitModelAndResidualCubicSpline<T, U, V>(num_data, data, context,
			num_boundary, boundary, coeff, context->best_fit_model, out);
}

} /* anonymous namespace */

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

/**
 * Create a baseline context object for float data to fit.
 * It is only for polynomial and Chebyshev polynomial baseline.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateBaselineContextFloat)(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t order,
		size_t num_data, LIBSAKURA_SYMBOL(BaselineContextFloat) **context)
				noexcept {
	CHECK_ARGS(
			(baseline_type == LIBSAKURA_SYMBOL(BaselineType_kPolynomial)) ||(baseline_type == LIBSAKURA_SYMBOL(BaselineType_kChebyshev)));
	CHECK_ARGS(context != nullptr);
	BaselineTypeInternal baseline_type_internal =
			BaselineTypeInternal_kPolynomial;
	if (baseline_type == LIBSAKURA_SYMBOL(BaselineType_kChebyshev)) {
		baseline_type_internal = BaselineTypeInternal_kChebyshev;
	}
	size_t num_lsq_bases = GetNumberOfLsqBases(baseline_type_internal, order);
	CHECK_ARGS(num_lsq_bases <= num_data);

	try {
		CreateBaselineContext<LIBSAKURA_SYMBOL(BaselineContextFloat)>(
				baseline_type_internal, order, 1, 0, num_data, context);
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

/**
 * Create a baseline context object for float data to fit.
 * It is only for cubic spline baseline.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateBaselineContextCubicSplineFloat)(
		uint16_t npiece, size_t num_data,
		LIBSAKURA_SYMBOL(BaselineContextFloat) **context) noexcept {
	BaselineTypeInternal const baseline_type = BaselineTypeInternal_kCubicSpline;
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(0 < npiece);
	size_t num_lsq_bases = GetNumberOfLsqBases(baseline_type, npiece);
	CHECK_ARGS(num_lsq_bases <= num_data);

	try {
		CreateBaselineContext<LIBSAKURA_SYMBOL(BaselineContextFloat)>(
				baseline_type, 0, npiece, 0, num_data, context);
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

/**
 * Create a baseline context object for float data to fit.
 * It is only for sinusoidal baseline.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateBaselineContextSinusoidFloat)(uint16_t nwave, size_t num_data, LIBSAKURA_SYMBOL(BaselineContextFloat) **context) noexcept {
	BaselineTypeInternal const baseline_type = BaselineTypeInternal_kSinusoid;
	CHECK_ARGS(context != nullptr);
	size_t num_lsq_bases = GetNumberOfLsqBases(baseline_type, nwave);
	CHECK_ARGS(num_lsq_bases <= num_data);

	try {
		CreateBaselineContext<LIBSAKURA_SYMBOL(BaselineContextFloat)>(
				baseline_type, 0, 1, nwave, num_data, context);
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

/**
 * Destroy baseline context object.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyBaselineContextFloat)(
LIBSAKURA_SYMBOL(BaselineContextFloat) *context) noexcept {
	CHECK_ARGS(context != nullptr);

	try {
		DestroyBaselineContext< LIBSAKURA_SYMBOL(BaselineContextFloat)
		>(context);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

/**
 * Computes the number of coefficients for given baseline type and
 * parameters. It can be used for polynomial, Chebyshev, and cubic
 * spline baseline.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetNumberOfCoefficientsFloat)(LIBSAKURA_SYMBOL(BaselineContextFloat) const *context, uint16_t order, size_t *num_coeff) noexcept {
	CHECK_ARGS(context != nullptr);
	auto const type = context->baseline_type;
	CHECK_ARGS(
			(type == BaselineTypeInternal_kPolynomial)
					|| (type == BaselineTypeInternal_kChebyshev)
					|| (type == BaselineTypeInternal_kCubicSpline));
	if (type == BaselineTypeInternal_kCubicSpline) {
		CHECK_ARGS(0 < order);
	}
	CHECK_ARGS(order <= context->baseline_param);

	try {
		*num_coeff = DoGetNumberOfCoefficients(type, order, 0, nullptr);
	} catch (const std::invalid_argument &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

/**
 * Fit baseline using polynomial or Chebyshev polynomial then subtract it.
 *
 * @param[out] out If @a get_residual is true, it is
 * baseline-subtracted data. Otherwise, it is the
 * best-fit baseline data.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineFloat)(
LIBSAKURA_SYMBOL(BaselineContextFloat) const *context, uint16_t order,
		size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], float clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual, float out[/*num_data*/],
		bool final_mask[/*num_data*/], float *rms,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	CHECK_ARGS(baseline_status != nullptr);
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kNG);
	CHECK_ARGS(context != nullptr);
	auto const type = context->baseline_type;
	CHECK_ARGS(
			(type == BaselineTypeInternal_kPolynomial)
					|| (type == BaselineTypeInternal_kChebyshev));
	CHECK_ARGS(order <= context->baseline_param);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(0.0f < clip_threshold_sigma);
	CHECK_ARGS(final_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(final_mask));
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));

	try {
		SubtractBaseline<float, double, LIBSAKURA_SYMBOL(BaselineContextFloat)>(
				context, order, 0, nullptr, num_data, data, mask,
				clip_threshold_sigma, num_fitting_max, get_residual, out,
				final_mask, rms, baseline_status);
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
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kOK);
	return LIBSAKURA_SYMBOL(Status_kOK);
}

/**
 * Fit baseline using cubic spline curve then subtract it.
 *
 * @param[out] out If @a get_residual is true, it is
 * baseline-subtracted data. Otherwise, it is the
 * best-fit baseline data.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(
LIBSAKURA_SYMBOL(BaselineContextFloat) const *context, size_t num_pieces,
		size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], float clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual, float out[/*num_data*/],
		bool final_mask[/*num_data*/], float *rms,
		size_t boundary[/*num_pieces+1*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	CHECK_ARGS(baseline_status != nullptr);
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kNG);
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(context->baseline_type == BaselineTypeInternal_kCubicSpline);
	CHECK_ARGS(0 < num_pieces);
	CHECK_ARGS(num_pieces <= context->baseline_param);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(0.0f < clip_threshold_sigma);
	CHECK_ARGS(final_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(final_mask));
	CHECK_ARGS(boundary != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(boundary));
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));

	try {
		SubtractBaselineCubicSpline<float, double,
		LIBSAKURA_SYMBOL(BaselineContextFloat)>(context, num_pieces + 1,
				num_data, data, mask, clip_threshold_sigma, num_fitting_max,
				get_residual, out, final_mask, rms, boundary, baseline_status);
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
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kOK);
	return LIBSAKURA_SYMBOL(Status_kOK);
}

/**
 * Fit baseline using sinusoids then subtract it.
 *
 * @param[out] out If @a get_residual is true, it is
 * baseline-subtracted data. Otherwise, it is the
 * best-fit baseline data.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineSinusoidFloat)(
LIBSAKURA_SYMBOL(BaselineContextFloat) const *context, size_t num_nwave,
		size_t const nwave[/*num_nwave*/], size_t num_data,
		float const data[/*num_data*/],
		bool const mask[/*num_data*/], float clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual, float out[/*num_data*/],
		bool final_mask[/*num_data*/], float *rms,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	CHECK_ARGS(baseline_status != nullptr);
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kNG);
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(context->baseline_type == BaselineTypeInternal_kSinusoid);
	CHECK_ARGS(0 < num_nwave);
	CHECK_ARGS(nwave != nullptr);
	CHECK_ARGS(IsUniqueAndAscendingOrder(num_nwave, nwave));
	CHECK_ARGS(nwave[num_nwave - 1] <= context->baseline_param);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(0.0f < clip_threshold_sigma);
	CHECK_ARGS(final_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(final_mask));
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));

	try {
		SubtractBaseline<float, double, LIBSAKURA_SYMBOL(BaselineContextFloat)>(
				context, 0, num_nwave, nwave, num_data, data, mask,
				clip_threshold_sigma, num_fitting_max, get_residual, out,
				final_mask, rms, baseline_status);
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
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kOK);
	return LIBSAKURA_SYMBOL(Status_kOK);
}

/**
 * Fit baseline using polynomial or Chebyshev polynomial and
 * compute the best-fit coefficients.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsFloat)(
LIBSAKURA_SYMBOL(BaselineContextFloat) const *context, size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		float clip_threshold_sigma, uint16_t num_fitting_max, size_t num_coeff,
		double coeff[/*num_coeff*/], bool final_mask[/*num_data*/], float *rms,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	CHECK_ARGS(baseline_status != nullptr);
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kNG);
	CHECK_ARGS(context != nullptr);
	auto const type = context->baseline_type;
	CHECK_ARGS(
			(type == BaselineTypeInternal_kPolynomial)
					|| (type == BaselineTypeInternal_kChebyshev));
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(0.0f < clip_threshold_sigma);
	CHECK_ARGS(0 < num_coeff);
	CHECK_ARGS(num_coeff <= context->num_bases);
	CHECK_ARGS(coeff != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	CHECK_ARGS(final_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(final_mask));

	try {
		GetBestFitBaselineCoefficients<float, double,
		LIBSAKURA_SYMBOL(BaselineContextFloat)>(context, num_data, data, mask,
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
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kOK);
	return LIBSAKURA_SYMBOL(Status_kOK);
}

/**
 * Fit baseline using cubic spline and
 * compute the best-fit coefficients.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsCubicSplineFloat)(
LIBSAKURA_SYMBOL(BaselineContextFloat) const *context, size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		float clip_threshold_sigma, uint16_t num_fitting_max, size_t num_pieces,
		double coeff[/*num_pieces*/][kNumBasesCubicSpline],
		bool final_mask[/*num_data*/], float *rms,
		size_t boundary[/*num_pieces+1*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	CHECK_ARGS(baseline_status != nullptr);
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kNG);
	CHECK_ARGS(context != nullptr);
	BaselineTypeInternal const baseline_type = BaselineTypeInternal_kCubicSpline;
	CHECK_ARGS(context->baseline_type == baseline_type);
	CHECK_ARGS(0 < num_pieces);
	CHECK_ARGS(num_pieces <= context->baseline_param);
	CHECK_ARGS(GetNumberOfLsqBases(baseline_type, num_pieces) <= num_data);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(0.0f < clip_threshold_sigma);
	CHECK_ARGS(coeff != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	CHECK_ARGS(final_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(final_mask));
	CHECK_ARGS(boundary != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(boundary));

	try {
		GetBestFitBaselineCoefficientsCubicSpline<float, double,
		LIBSAKURA_SYMBOL(BaselineContextFloat)>(context, num_data, data, mask,
				clip_threshold_sigma, num_fitting_max, num_pieces + 1, coeff,
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
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kOK);
	return LIBSAKURA_SYMBOL(Status_kOK);
}

/**
 * Fit baseline using sinusoids and
 * compute the best-fit coefficients.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsSinusoidFloat)(
LIBSAKURA_SYMBOL(BaselineContextFloat) const *context, size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		float clip_threshold_sigma, uint16_t num_fitting_max, size_t num_nwave,
		size_t const nwave[/*num_nwave*/], size_t num_coeff,
		double coeff[/*num_coeff*/], bool final_mask[/*num_data*/], float *rms,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) noexcept {
	CHECK_ARGS(baseline_status != nullptr);
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kNG);
	CHECK_ARGS(context != nullptr);
	BaselineTypeInternal const baseline_type = BaselineTypeInternal_kSinusoid;
	CHECK_ARGS(context->baseline_type == baseline_type);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(0.0f < clip_threshold_sigma);
	CHECK_ARGS(0 < num_nwave);
	CHECK_ARGS(nwave != nullptr);
	CHECK_ARGS(IsUniqueAndAscendingOrder(num_nwave, nwave));
	CHECK_ARGS(nwave[num_nwave - 1] <= context->baseline_param);
	size_t num_lsq_bases = DoGetNumberOfCoefficients(baseline_type, 0,
			num_nwave, nwave);
	CHECK_ARGS(num_lsq_bases <= num_coeff);
	CHECK_ARGS(num_coeff <= context->num_bases);
	CHECK_ARGS(coeff != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	CHECK_ARGS(final_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(final_mask));

	try {
		GetBestFitBaselineCoefficients<float, double,
		LIBSAKURA_SYMBOL(BaselineContextFloat)>(context, num_data, data, mask,
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
	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kOK);
	return LIBSAKURA_SYMBOL(Status_kOK);
}

/**
 * Subtract baseline using basis data in baseline context and
 * user-given coefficients. This function is for baseline type
 * of polynomial and Chebyshev polynomial.
 *
 * @param[out] out_arg Baseline-subtracted data.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineUsingCoefficientsFloat)(
LIBSAKURA_SYMBOL(BaselineContextFloat) const *context, size_t num_data,
		float const data[/*num_data*/], size_t num_coeff,
		double const coeff[/*num_coeff*/], float out[/*num_data*/]) noexcept {
	CHECK_ARGS(context != nullptr);
	auto const type = context->baseline_type;
	CHECK_ARGS(
			(type == BaselineTypeInternal_kPolynomial)
					|| (type == BaselineTypeInternal_kChebyshev));
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(0 < num_coeff);
	CHECK_ARGS(num_coeff <= context->num_bases);
	CHECK_ARGS(coeff != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));

	try {
		SubtractBaselineUsingCoefficients<float, double,
		LIBSAKURA_SYMBOL(BaselineContextFloat)>(context, num_data, data,
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

/**
 * Subtract baseline using basis data in baseline context and
 * user-given coefficients. This function is for baseline type
 * of cubic spline.
 *
 * @param[out] out_arg Baseline-subtracted data.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
LIBSAKURA_SYMBOL(BaselineContextFloat) const *context, size_t num_data,
		float const data[/*num_data*/], size_t num_pieces,
		double const coeff[/*num_pieces*/][kNumBasesCubicSpline],
		size_t const boundary[/*num_pieces+1*/], float out[/*num_data*/])
				noexcept {
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(context->baseline_type == BaselineTypeInternal_kCubicSpline);
	CHECK_ARGS(0 < num_pieces);
	CHECK_ARGS(num_pieces <= context->baseline_param);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(coeff != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	CHECK_ARGS(boundary != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(boundary));
	CHECK_ARGS(boundary[0] == 0);
	CHECK_ARGS(boundary[num_pieces] == num_data);
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));

	try {
		SubtractBaselineCubicSplineUsingCoefficients<float, double,
		LIBSAKURA_SYMBOL(BaselineContextFloat)>(context, num_data, data,
				num_pieces + 1, coeff, boundary, out);
	} catch (const std::bad_alloc &e) {
		LOG4CXX_ERROR(logger, "Memory allocation failed.");
		return LIBSAKURA_SYMBOL(Status_kNoMemory);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

/**
 * Subtract baseline using basis data in baseline context and
 * user-given coefficients. This function is for baseline type
 * of sinusoids.
 *
 * @param[out] out_arg Baseline-subtracted data.
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineSinusoidUsingCoefficientsFloat)(
LIBSAKURA_SYMBOL(BaselineContextFloat) const *context, size_t num_data,
		float const data[/*num_data*/], size_t num_nwave,
		size_t const nwave[/*num_nwave*/], size_t num_coeff,
		double const coeff[/*num_coeff*/], float out[/*num_data*/]) noexcept {
	CHECK_ARGS(context != nullptr);
	CHECK_ARGS(context->baseline_type == BaselineTypeInternal_kSinusoid);
	CHECK_ARGS(num_data == context->num_basis_data);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(0 < num_nwave);
	CHECK_ARGS(nwave != nullptr);
	CHECK_ARGS(IsUniqueAndAscendingOrder(num_nwave, nwave));
	CHECK_ARGS(nwave[num_nwave - 1] <= context->baseline_param);
	size_t num_lsq_bases = DoGetNumberOfCoefficients(context->baseline_type, 0,
			num_nwave, nwave);
	CHECK_ARGS(num_lsq_bases <= num_coeff);
	CHECK_ARGS(num_coeff <= context->num_bases);
	CHECK_ARGS(coeff != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));

	try {
		SubtractBaselineUsingCoefficients<float, double,
		LIBSAKURA_SYMBOL(BaselineContextFloat)>(context, num_data, data,
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
