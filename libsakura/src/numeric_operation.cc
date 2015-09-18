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
#include <iostream>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>

#if defined(__AVX__) && !defined(ARCH_SCALAR)
#	include <immintrin.h>
#endif

#include <Eigen/Core>
#include <Eigen/LU>
//LM start----------------------------------
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
//LM end----------------------------------

#include "libsakura/sakura.h"
#include "libsakura/localdef.h"
#include "libsakura/logger.h"
#include "libsakura/memory_manager.h"
#include "libsakura/packed_type.h"
namespace {
#include "libsakura/packed_operation.h"
}

using ::Eigen::Map;
using ::Eigen::MatrixXd;
using ::Eigen::VectorXd;
using ::Eigen::Aligned;
using ::Eigen::Stride;
using ::Eigen::Dynamic;
//LM start----------------------------------
using ::Eigen::Matrix;
using ::Eigen::NumericalDiff;
using ::Eigen::LevenbergMarquardt;
//LM end---------------------------------

namespace {

auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("numeric_operation");

#if defined(__AVX__)

template<typename T>
static T NegMultiplyAdd(T const &a, T const &b, T const &c) {
	assert(false); // not defined for this type.
	return T();
}

template<>
__m256d NegMultiplyAdd(__m256d                     const &a, __m256d                     const &b, __m256d                     const &c) {
#if defined(__AVX2__)
	return _mm256_fnmadd_pd(a, b, c);
#else
	return _mm256_sub_pd(c, _mm256_mul_pd(a, b));
#endif
}

#endif

template<size_t kNumBases>
void AddMulVectorTemplate(size_t const *use_idx, double k, double const *vec,
		double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	constexpr size_t kPackElements = sizeof(__m256d) / sizeof(double);
	constexpr size_t kEnd = (kNumBases / kPackElements) * kPackElements;
	auto coeff = _mm256_set1_pd(k);
	for (i = 0; i < kEnd; i += kPackElements) {
		auto v = _mm256_set_pd(vec[use_idx[i + 3]], vec[use_idx[i + 2]],
				vec[use_idx[i + 1]], vec[use_idx[i]]);
		_mm256_storeu_pd(&out[i],
				LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
				LIBSAKURA_SYMBOL(SimdPacketAVX), double>(coeff, v,
						_mm256_loadu_pd(&out[i])));
	}
#endif
	for (; i < kNumBases; ++i) {
		out[i] += k * vec[use_idx[i]];
	}
}

template<size_t kNumBases>
void SubMulVectorTemplate(size_t const *use_idx, double k, double const *vec,
		double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	constexpr size_t kPackElements = sizeof(__m256d) / sizeof(double);
	constexpr size_t kEnd = (kNumBases / kPackElements) * kPackElements;
	auto coeff = _mm256_set1_pd(k);
	for (i = 0; i < kEnd; i += kPackElements) {
		auto v = _mm256_set_pd(vec[use_idx[i + 3]], vec[use_idx[i + 2]],
				vec[use_idx[i + 1]], vec[use_idx[i]]);
		_mm256_storeu_pd(&out[i],
				NegMultiplyAdd(coeff, v, _mm256_loadu_pd(&out[i])));
	}
#endif
	for (; i < kNumBases; ++i) {
		out[i] -= k * vec[use_idx[i]];
	}
}

template<typename Func1, typename Func2>
inline void OperateAddMulVector(size_t const num_model_bases,
		size_t const *use_idx, double k, double const *vec, Func1 func1,
		Func2 func2, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	constexpr size_t kPackElements = sizeof(__m256d) / sizeof(double);
	size_t const end = (num_model_bases / kPackElements) * kPackElements;
	auto coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += kPackElements) {
		auto v = _mm256_set_pd(vec[use_idx[i + 3]], vec[use_idx[i + 2]],
				vec[use_idx[i + 1]], vec[use_idx[i]]);
		_mm256_storeu_pd(&out[i], func1(coeff, v, i));
	}
#endif
	for (; i < num_model_bases; ++i) {
		func2(i);
	}
}

inline void AddMulVector(size_t const num_model_bases, size_t const *use_idx,
		double k, double const *vec, double *out) {
	OperateAddMulVector(num_model_bases, use_idx, k, vec,
#if defined(__AVX__) && !defined(ARCH_SCALAR)
			[&](__m256d coeff, __m256d v, size_t i) {
				return LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
				LIBSAKURA_SYMBOL(SimdPacketAVX), double>(coeff, v,
						_mm256_loadu_pd(&out[i]));},
#else
			[&]() {},
#endif
			[&](size_t i) {out[i] += k * vec[use_idx[i]];}, out);
}

inline void SubMulVector(size_t const num_model_bases, size_t const *use_idx,
		double k, double const *vec, double *out) {
	OperateAddMulVector(num_model_bases, use_idx, k, vec,
#if defined(__AVX__) && !defined(ARCH_SCALAR)
			[&](__m256d coeff, __m256d v, size_t i) {
				return NegMultiplyAdd(coeff, v, _mm256_loadu_pd(&out[i]));},
#else
			[&]() {},
#endif
			[&](size_t i) {out[i] -= k * vec[use_idx[i]];}, out);
}

template<size_t kNumBases>
inline void GetLSQFittingMatrixTemplate(size_t num_mask, bool const *mask_arg,
		size_t num_model_bases, double const *model_arg,
		size_t const *use_bases_idx_arg, double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(use_bases_idx_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const mask = AssumeAligned(mask_arg);
	auto const model = AssumeAligned(model_arg);
	auto const use_bases_idx = AssumeAligned(use_bases_idx_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < kNumBases * kNumBases; ++i) {
		out[i] = 0;
	}
	size_t num_unmasked_data = 0;
	for (size_t i = 0; i < num_mask; ++i) {
		if (mask[i]) {
			auto const model_i = &model[i * num_model_bases];
			for (size_t j = 0; j < kNumBases; ++j) {
				auto out_matrix_j = &out[j * kNumBases];
				AddMulVectorTemplate<kNumBases>(use_bases_idx,
						model_i[use_bases_idx[j]], model_i, out_matrix_j);
			}
			++num_unmasked_data;
		}
	}

	if (num_unmasked_data < kNumBases) {
		throw std::runtime_error(
				"GetLSQFittingMatrixTemplate: too many masked data.");
	}
}

inline void GetLSQFittingMatrix(size_t num_mask, bool const *mask_arg,
		size_t num_model_bases, double const *model_arg, size_t num_lsq_bases,
		size_t const *use_bases_idx_arg, double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(use_bases_idx_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const mask = AssumeAligned(mask_arg);
	auto const model = AssumeAligned(model_arg);
	auto const use_bases_idx = AssumeAligned(use_bases_idx_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_lsq_bases * num_lsq_bases; ++i) {
		out[i] = 0;
	}
	size_t num_unmasked_data = 0;
	for (size_t i = 0; i < num_mask; ++i) {
		if (mask[i]) {
			auto const model_i = &model[i * num_model_bases];
			for (size_t j = 0; j < num_lsq_bases; ++j) {
				auto out_matrix_j = &out[j * num_lsq_bases];
				AddMulVector(num_lsq_bases, use_bases_idx,
						model_i[use_bases_idx[j]], model_i, out_matrix_j);
			}
			++num_unmasked_data;
		}
	}

	if (num_unmasked_data < num_model_bases) {
		throw std::runtime_error(
				"GetLSQFittingMatrix: too many data are masked.");
	}
}

template<size_t kNumBases>
inline void UpdateLSQFittingMatrixTemplate(bool const *mask_arg,
		size_t num_clipped, size_t const *clipped_indices_arg,
		size_t num_model_bases, double const *in_arg, double const *model_arg,
		size_t const *use_bases_idx_arg, double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(use_bases_idx_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const mask = AssumeAligned(mask_arg);
	auto const clipped_indices = AssumeAligned(clipped_indices_arg);
	auto const in = AssumeAligned(in_arg);
	auto const model = AssumeAligned(model_arg);
	auto const use_bases_idx = AssumeAligned(use_bases_idx_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < kNumBases * kNumBases; ++i) {
		out[i] = in[i];
	}
	for (size_t i = 0; i < num_clipped; ++i) {
		if (!mask[clipped_indices[i]])
			continue;
		auto const model_i = &model[clipped_indices[i] * num_model_bases];
		for (size_t j = 0; j < kNumBases; ++j) {
			auto out_matrix_j = &out[j * kNumBases];
			SubMulVectorTemplate<kNumBases>(use_bases_idx,
					model_i[use_bases_idx[j]], model_i, out_matrix_j);
		}
	}
}

inline void UpdateLSQFittingMatrix(bool const *mask_arg, size_t num_clipped,
		size_t const *clipped_indices_arg, size_t num_lsq_bases,
		double const *in_arg, size_t num_model_bases, double const *model_arg,
		size_t const *use_bases_idx_arg, double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(use_bases_idx_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const mask = AssumeAligned(mask_arg);
	auto const clipped_indices = AssumeAligned(clipped_indices_arg);
	auto const in = AssumeAligned(in_arg);
	auto const model = AssumeAligned(model_arg);
	auto const use_bases_idx = AssumeAligned(use_bases_idx_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_lsq_bases * num_lsq_bases; ++i) {
		out[i] = in[i];
	}
	for (size_t i = 0; i < num_clipped; ++i) {
		if (!mask[clipped_indices[i]])
			continue;
		auto const model_i = &model[clipped_indices[i] * num_model_bases];
		for (size_t j = 0; j < num_lsq_bases; ++j) {
			auto out_matrix_j = &out[j * num_lsq_bases];
			SubMulVector(num_lsq_bases, use_bases_idx,
					model_i[use_bases_idx[j]], model_i, out_matrix_j);
		}
	}
}

template<size_t kNumBases, typename T>
inline void GetLSQFittingVectorTemplate(size_t num_data, T const *data_arg,
bool const *mask_arg, size_t num_model_bases, double const *model_arg,
		size_t const *use_bases_idx_arg, double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(use_bases_idx_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto const model = AssumeAligned(model_arg);
	auto const use_bases_idx = AssumeAligned(use_bases_idx_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < kNumBases; ++i) {
		out[i] = 0;
	}
	for (size_t i = 0; i < num_data; ++i) {
		if (mask[i]) {
			auto const model_i = &model[i * num_model_bases];
			auto data_i = data[i];
			AddMulVectorTemplate<kNumBases>(use_bases_idx, data_i, model_i,
					out);
		}
	}
}

template<typename T>
inline void GetLSQFittingVector(size_t num_data, T const *data_arg,
bool const *mask_arg, size_t num_model_bases, double const *model_arg,
		size_t num_lsq_bases, size_t const *use_bases_idx_arg,
		double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(use_bases_idx_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto const model = AssumeAligned(model_arg);
	auto const use_bases_idx = AssumeAligned(use_bases_idx_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_lsq_bases; ++i) {
		out[i] = 0;
	}
	for (size_t i = 0; i < num_data; ++i) {
		if (mask[i]) {
			auto const model_i = &model[i * num_model_bases];
			auto data_i = data[i];
			AddMulVector(num_lsq_bases, use_bases_idx, data_i, model_i, out);
		}
	}
}

template<size_t kNumBases, typename T> inline void UpdateLSQFittingVectorTemplate(
		T const *data_arg, bool const *mask_arg, size_t num_clipped,
		size_t const *clipped_indices_arg, size_t num_model_bases,
		double const *in_arg, double const *model_arg,
		size_t const *use_bases_idx_arg, double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(use_bases_idx_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto const clipped_indices = AssumeAligned(clipped_indices_arg);
	auto const in = AssumeAligned(in_arg);
	auto const model = AssumeAligned(model_arg);
	auto const use_bases_idx = AssumeAligned(use_bases_idx_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < kNumBases; ++i) {
		out[i] = in[i];
	}
	for (size_t i = 0; i < num_clipped; ++i) {
		if (!mask[clipped_indices[i]])
			continue;
		auto const cii = clipped_indices[i];
		auto const model_i = &model[cii * num_model_bases];
		auto data_i = data[cii];
		SubMulVectorTemplate<kNumBases>(use_bases_idx, data_i, model_i, out);
	}
}

template<typename T>
inline void UpdateLSQFittingVector(T const *data_arg, bool const *mask_arg,
		size_t num_clipped, size_t const *clipped_indices_arg,
		size_t num_lsq_bases, double const *in_arg, size_t num_model_bases,
		double const *model_arg, size_t const *use_bases_idx_arg,
		double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(use_bases_idx_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto const data = AssumeAligned(data_arg);
	auto const mask = AssumeAligned(mask_arg);
	auto const clipped_indices = AssumeAligned(clipped_indices_arg);
	auto const in = AssumeAligned(in_arg);
	auto const model = AssumeAligned(model_arg);
	auto const use_bases_idx = AssumeAligned(use_bases_idx_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_lsq_bases; ++i) {
		out[i] = in[i];
	}
	for (size_t i = 0; i < num_clipped; ++i) {
		if (!mask[clipped_indices[i]])
			continue;
		auto const cii = clipped_indices[i];
		auto const model_i = &model[cii * num_model_bases];
		auto data_i = data[cii];
		SubMulVector(num_lsq_bases, use_bases_idx, data_i, model_i, out);
	}
}

template<size_t kNumBases, typename T>
inline void GetLSQCoefficientsTemplate(size_t num_data, T const *data,
bool const *mask, size_t const num_model_bases, double const *basis,
		size_t const num_lsq_bases, size_t const *use_bases_idx,
		double *lsq_matrix, double *lsq_vector) {
	GetLSQFittingMatrixTemplate<kNumBases>(num_data, mask, num_model_bases,
			basis, use_bases_idx, lsq_matrix);
	GetLSQFittingVectorTemplate<kNumBases, T>(num_data, data, mask,
			num_model_bases, basis, use_bases_idx, lsq_vector);
}

template<typename T>
inline void GetLSQCoefficients(size_t num_data, T const *data, bool const *mask,
		size_t const num_model_bases, double const *basis,
		size_t const num_lsq_bases, size_t const *use_bases_idx,
		double *lsq_matrix, double *lsq_vector) {
	GetLSQFittingMatrix(num_data, mask, num_model_bases, basis, num_lsq_bases,
			use_bases_idx, lsq_matrix);
	GetLSQFittingVector<T>(num_data, data, mask, num_model_bases, basis,
			num_lsq_bases, use_bases_idx, lsq_vector);
}

template<size_t kNumBases, typename T>
inline void UpdateLSQCoefficientsTemplate(size_t const num_data, T const *data,
bool const *mask, size_t const num_clipped, size_t const *clipped_indices,
		size_t const num_model_bases, double const *basis,
		size_t const num_lsq_bases, size_t const *use_bases_idx,
		double *lsq_matrix, double *lsq_vector) {
	UpdateLSQFittingMatrixTemplate<kNumBases>(mask, num_clipped,
			clipped_indices, num_model_bases, lsq_matrix, basis, use_bases_idx,
			lsq_matrix);
	UpdateLSQFittingVectorTemplate<kNumBases, T>(data, mask, num_clipped,
			clipped_indices, num_model_bases, lsq_vector, basis, use_bases_idx,
			lsq_vector);
}

template<typename T>
inline void UpdateLSQCoefficients(size_t const num_data, T const *data,
bool const *mask, size_t const num_clipped, size_t const *clipped_indices,
		size_t const num_model_bases, double const *basis,
		size_t const num_lsq_bases, size_t const *use_bases_idx,
		double *lsq_matrix, double *lsq_vector) {
	UpdateLSQFittingMatrix(mask, num_clipped, clipped_indices, num_lsq_bases,
			lsq_matrix, num_model_bases, basis, use_bases_idx, lsq_matrix);
	UpdateLSQFittingVector<T>(data, mask, num_clipped, clipped_indices,
			num_lsq_bases, lsq_vector, num_model_bases, basis, use_bases_idx,
			lsq_vector);
}

template<typename T, typename MatrixT, typename VectorT>
inline void SolveSimultaneousEquationsByLU(size_t num_equations,
		T const *in_matrix_arg, T const *in_vector_arg, T *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_matrix_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_vector_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map<MatrixT, Aligned, Stride<Dynamic, Dynamic> > in_matrix(
			const_cast<T *>(in_matrix_arg), num_equations, num_equations,
			Stride<Dynamic, Dynamic>(1, num_equations));
	Map < VectorT, Aligned
			> in_vector(const_cast<T *>(in_vector_arg), num_equations);

	Map<VectorT>(out, num_equations) = in_matrix.fullPivLu().solve(in_vector);
}

#define RepeatTen(func, start_idx, type) \
		(func<start_idx + 0, type>), \
		(func<start_idx + 1, type>), \
		(func<start_idx + 2, type>), \
		(func<start_idx + 3, type>), \
		(func<start_idx + 4, type>), \
		(func<start_idx + 5, type>), \
		(func<start_idx + 6, type>), \
		(func<start_idx + 7, type>), \
		(func<start_idx + 8, type>), \
		(func<start_idx + 9, type>)

template<typename T>
void GetLSQCoefficientsEntry(size_t const num_data, T const data[/*num_data*/],
bool const mask[/*num_data*/], size_t const num_model_bases,
		double const basis_data[/*num_model_bases*num_data*/],
		size_t const num_lsq_bases,
		size_t const use_bases_idx[/*num_lsq_bases*/],
		double lsq_matrix[/*num_lsq_bases*num_lsq_bases*/],
		double lsq_vector[/*num_lsq_bases*/]) {
	using GetLSQCoefficientsFunc = void (*)(size_t const num_data, T const *data,
			bool const *mask, size_t const num_model_bases, double const *basis_data,
			size_t const num_lsq_bases, size_t const *use_bases_idx,
			double *lsq_matrix, double *lsq_vector);

	static GetLSQCoefficientsFunc const funcs[] = { RepeatTen(
			GetLSQCoefficientsTemplate, 0, T), RepeatTen(
			GetLSQCoefficientsTemplate, 10, T), RepeatTen(
			GetLSQCoefficientsTemplate, 20, T), RepeatTen(
			GetLSQCoefficientsTemplate, 30, T), RepeatTen(
			GetLSQCoefficientsTemplate, 40, T), RepeatTen(
			GetLSQCoefficientsTemplate, 50, T), RepeatTen(
			GetLSQCoefficientsTemplate, 60, T), RepeatTen(
			GetLSQCoefficientsTemplate, 70, T), RepeatTen(
			GetLSQCoefficientsTemplate, 80, T), RepeatTen(
			GetLSQCoefficientsTemplate, 90, T), GetLSQCoefficientsTemplate<100,
			T> };

	if (num_lsq_bases < ELEMENTSOF(funcs)) {
		funcs[num_lsq_bases](num_data, data, mask, num_model_bases, basis_data,
				num_lsq_bases, use_bases_idx, lsq_matrix, lsq_vector);
	} else {
		GetLSQCoefficients<T>(num_data, data, mask, num_model_bases, basis_data,
				num_lsq_bases, use_bases_idx, lsq_matrix, lsq_vector);
	}
}

template<typename T>
void UpdateLSQCoefficientsEntry(size_t const num_data,
		T const data[/*num_data*/], bool const mask[/*num_data*/],
		size_t const num_clipped, size_t const clipped_indices[/*num_data*/],
		size_t const num_model_bases,
		double const basis_data[/*num_model_bases*num_data*/],
		size_t const num_lsq_bases,
		size_t const use_bases_idx[/*num_lsq_bases*/],
		double lsq_matrix[/*num_lsq_bases*num_lsq_bases*/],
		double lsq_vector[/*num_lsq_bases*/]) {
	using UpdateLSQCoefficientsFunc = void (*)(size_t const num_data,
			T const *data, bool const *mask, size_t const num_clipped,
			size_t const *clipped_indices, size_t const num_model_bases,
			double const *basis_data, size_t const num_lsq_bases,
			size_t const *use_bases_idx, double *lsq_matrix,
			double *lsq_vector);

	static UpdateLSQCoefficientsFunc const funcs[] = { RepeatTen(
			UpdateLSQCoefficientsTemplate, 0, T), RepeatTen(
			UpdateLSQCoefficientsTemplate, 10, T), RepeatTen(
			UpdateLSQCoefficientsTemplate, 20, T), RepeatTen(
			UpdateLSQCoefficientsTemplate, 30, T), RepeatTen(
			UpdateLSQCoefficientsTemplate, 40, T), RepeatTen(
			UpdateLSQCoefficientsTemplate, 50, T), RepeatTen(
			UpdateLSQCoefficientsTemplate, 60, T), RepeatTen(
			UpdateLSQCoefficientsTemplate, 70, T), RepeatTen(
			UpdateLSQCoefficientsTemplate, 80, T), RepeatTen(
			UpdateLSQCoefficientsTemplate, 90, T),
			UpdateLSQCoefficientsTemplate<100, T> };

	if (num_lsq_bases < ELEMENTSOF(funcs)) {
		funcs[num_lsq_bases](num_data, data, mask, num_clipped, clipped_indices,
				num_model_bases, basis_data, num_lsq_bases, use_bases_idx,
				lsq_matrix, lsq_vector);
	} else {
		UpdateLSQCoefficients<T>(num_data, data, mask, num_clipped,
				clipped_indices, num_model_bases, basis_data, num_lsq_bases,
				use_bases_idx, lsq_matrix, lsq_vector);
	}
}

//--LM part--------------------------------------
/*
 template<typename T, int Nx = Dynamic, int Ny = Dynamic>
 struct Functor {
 typedef T Scalar;
 enum {
 InputsAtCompileTime = Nx, ValuesAtCompileTime = Ny
 };
 typedef Matrix<Scalar, Nx, 1> InputType;
 typedef Matrix<Scalar, Ny, 1> ValueType;
 typedef Matrix<Scalar, Ny, Nx> JacobianType;
 };

 struct FunctorDouble: Functor<double> {
 FunctorDouble(int inputs, int values, double *x, double *y) :
 inputs_(inputs), values_(values), x(x), y(y) {
 }
 int operator()(const VectorXd &params, VectorXd &values_to_minimize) const {
 return 0;
 }
 const int inputs_;
 const int values_;
 int inputs() const {
 return inputs_;
 }
 int values() const {
 return values_;
 }
 double *x;
 double *y;
 };

 struct GaussianFunctorDouble: FunctorDouble {
 GaussianFunctorDouble(size_t num_components, int inputs, int values,
 double *x, double *y) :
 FunctorDouble(inputs, values, x, y), num_components_(num_components) {
 }
 ;
 int operator()(const VectorXd &params, VectorXd &values_to_minimize) const {
 //params are {peak_amplitude, peak_position, sigma}.
 for (int i = 0; i < values_; ++i) {
 values_to_minimize[i] = 0.0;
 for (size_t iline = 0; iline < num_components_; ++iline) {
 double v = (x[i] - params[3 * iline + 1])
 / params[3 * iline + 2];
 values_to_minimize[i] += params[3 * iline]
 * exp(-0.5 / M_PI * v * v);
 }
 values_to_minimize[i] -= y[i];
 }
 return 0;
 }
 size_t num_components_;
 };

 void DoFitGaussianByLMDouble(size_t const n, size_t const num_data,
 double const x_data[], double const y_data[], double peak_amplitude[],
 double peak_position[], double sigma[]) {
 size_t const num_parameters = 3 * n;
 VectorXd params(num_parameters);
 for (size_t i = 0; i < n; ++i) {
 params << peak_amplitude[i], peak_position[i], sigma[i];
 }

 std::vector<double> x(&x_data[0], &x_data[num_data]);
 std::vector<double> y(&y_data[0], &y_data[num_data]);
 GaussianFunctorDouble functor(n, num_parameters, x.size(), &x[0], &y[0]);
 NumericalDiff<GaussianFunctorDouble> num_diff(functor);
 LevenbergMarquardt < NumericalDiff<GaussianFunctorDouble> > lm(num_diff);
 int info = lm.minimize(params);
 if (info == 0) {
 throw std::runtime_error(
 "FitGaussianByLMDouble: invalid data or parameter for Gaussian fitting.");
 }
 for (size_t i = 0; i < n; ++i) {
 peak_amplitude[i] = params[3 * i];
 peak_position[i] = params[3 * i + 1];
 sigma[i] = params[3 * i + 2];
 }
 }
 */
//--End LM part----------------------------------
} /* anonymous namespace */

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

//--LM part--------------------------------------
/*
 extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FitGaussianByLMFloat)(
 size_t const num_data, float const data[], bool const mask[],
 size_t const num_lines, double amplitude[//num_lines
 ],
 double position[//num_lines
 ],
 double sigma[//num_lines
 ]
 ) noexcept {
 try {
 double data_x[num_data];
 double data_y[num_data];
 size_t num_data_effective = 0;
 for (size_t i = 0; i < num_data; ++i) {
 if (mask[i]) {
 data_x[num_data_effective] = static_cast<double>(i);
 data_y[num_data_effective] = static_cast<double>(data[i]);
 ++num_data_effective;
 }
 }
 DoFitGaussianByLMDouble(num_lines, num_data_effective, data_x, data_y, amplitude, position, sigma);
 } catch (const std::runtime_error &e) {
 LOG4CXX_ERROR(logger, e.what());
 return LIBSAKURA_SYMBOL(Status_kNG);
 } catch (...) {
 assert(false);
 return LIBSAKURA_SYMBOL(Status_kUnknownError);
 }

 return LIBSAKURA_SYMBOL(Status_kOK);
 }
 */
//--End LM part----------------------------------
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetLSQCoefficientsDouble)(
		size_t const num_data, float const data[], bool const mask[],
		size_t const num_model_bases, double const basis_data[],
		size_t const num_lsq_bases, size_t const use_bases_idx[],
		double lsq_matrix[], double lsq_vector[]) noexcept {
	CHECK_ARGS(num_data != 0);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(num_model_bases != 0);
	CHECK_ARGS(num_model_bases <= num_data);
	CHECK_ARGS(basis_data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(basis_data));
	CHECK_ARGS(num_lsq_bases != 0);
	CHECK_ARGS(num_lsq_bases <= num_model_bases);
	CHECK_ARGS(use_bases_idx != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(use_bases_idx));
	CHECK_ARGS(lsq_matrix != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix));
	CHECK_ARGS(lsq_vector != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector));

	try {
		GetLSQCoefficientsEntry<float>(num_data, data, mask, num_model_bases,
				basis_data, num_lsq_bases, use_bases_idx, lsq_matrix,
				lsq_vector);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UpdateLSQCoefficientsDouble)(
		size_t const num_data, float const data[], bool const mask[],
		size_t const num_exclude_indices, size_t const exclude_indices[],
		size_t const num_model_bases, double const basis_data[],
		size_t const num_lsq_bases, size_t const use_bases_idx[],
		double lsq_matrix[], double lsq_vector[]) noexcept {
	CHECK_ARGS(num_data != 0);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(mask));
	CHECK_ARGS(exclude_indices != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(exclude_indices));
	CHECK_ARGS(num_exclude_indices <= num_data);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		CHECK_ARGS(exclude_indices[i] < num_data);
	}
	CHECK_ARGS(num_model_bases != 0);
	CHECK_ARGS(num_model_bases <= num_data);
	CHECK_ARGS(basis_data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(basis_data));
	CHECK_ARGS(num_lsq_bases != 0);
	CHECK_ARGS(num_lsq_bases <= num_model_bases);
	CHECK_ARGS(use_bases_idx != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(use_bases_idx));
	CHECK_ARGS(lsq_matrix != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix));
	CHECK_ARGS(lsq_vector != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector));

	try {
		UpdateLSQCoefficientsEntry<float>(num_data, data, mask,
				num_exclude_indices, exclude_indices, num_model_bases,
				basis_data, num_lsq_bases, use_bases_idx, lsq_matrix,
				lsq_vector);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLUDouble)(
		size_t num_equations, double const in_matrix[],
		double const in_vector[], double out[]) noexcept {
	CHECK_ARGS(in_matrix != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(in_matrix));
	CHECK_ARGS(in_vector != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(in_vector));
	CHECK_ARGS(in_vector != out);
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));

	try {
		SolveSimultaneousEquationsByLU<double, MatrixXd, VectorXd>(
				num_equations, in_matrix, in_vector, out);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}
