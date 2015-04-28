/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2014
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
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>

#if defined(__AVX__) && !defined(ARCH_SCALAR)
#	include <immintrin.h>
#endif

#include <Eigen/Core>
#include <Eigen/LU>

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

namespace {

auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("numeric_operation");

#if defined(__AVX__)

template<typename T>
static T NegMultiplyAdd(T const &a, T const &b, T const &c) {
	assert(false); // not defined for this type.
	return T();
}

template<>
__m256d NegMultiplyAdd(__m256d        const &a, __m256d        const &b, __m256d        const &c) {
#if defined(__AVX2__)
	return _mm256_fnmadd_pd(a, b, c);
#else
	return _mm256_sub_pd(c, _mm256_mul_pd(a, b));
#endif
}

#endif

template<size_t NUM_BASES>
void AddMulVectorTemplate(double k, double const *vec, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (NUM_BASES / pack_elements) * pack_elements;
	auto coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += pack_elements) {
		auto v = _mm256_loadu_pd(&vec[i]);
		_mm256_storeu_pd(&out[i],
				LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
				LIBSAKURA_SYMBOL(SimdPacketAVX), double>(coeff, v,
						_mm256_loadu_pd(&out[i])));
	}
#endif
	for (; i < NUM_BASES; ++i) {
		out[i] += k * vec[i];
	}
}

template<size_t NUM_BASES>
void SubMulVectorTemplate(double k, double const *vec, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (NUM_BASES / pack_elements) * pack_elements;
	auto coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += pack_elements) {
		auto v = _mm256_loadu_pd(&vec[i]);
		_mm256_storeu_pd(&out[i],
				NegMultiplyAdd(coeff, v, _mm256_loadu_pd(&out[i])));
	}
#endif
	for (; i < NUM_BASES; ++i) {
		out[i] -= k * vec[i];
	}
}

inline void AddMulVector(size_t const num_model_bases, double k,
		double const *vec, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (num_model_bases / pack_elements) * pack_elements;
	auto coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += pack_elements) {
		auto v = _mm256_loadu_pd(&vec[i]);
		_mm256_storeu_pd(&out[i],
				LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
				LIBSAKURA_SYMBOL(SimdPacketAVX), double>(coeff, v,
						_mm256_loadu_pd(&out[i])));
	}
#endif
	for (; i < num_model_bases; ++i) {
		out[i] += k * vec[i];
	}
}

inline void SubMulVector(size_t const num_model_bases, double k,
		double const *vec, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (num_model_bases / pack_elements) * pack_elements;
	auto coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += pack_elements) {
		auto v = _mm256_loadu_pd(&vec[i]);
		_mm256_storeu_pd(&out[i],
				NegMultiplyAdd(coeff, v, _mm256_loadu_pd(&out[i])));
	}
#endif
	for (; i < num_model_bases; ++i) {
		out[i] -= k * vec[i];
	}
}

template<size_t NUM_BASES>
inline void GetLSQFittingMatrixTemplate(size_t num_mask, bool const *mask_arg,
		size_t num_model_bases, double const *model_arg, double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto mask = AssumeAligned(mask_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < NUM_BASES * NUM_BASES; ++i) {
		out[i] = 0;
	}
	size_t num_unmasked_data = 0;
	for (size_t i = 0; i < num_mask; ++i) {
		if (mask[i]) {
			auto model_i = &model[i * num_model_bases];
			for (size_t j = 0; j < NUM_BASES; ++j) {
				auto out_matrix_i = &out[j * NUM_BASES];
				AddMulVectorTemplate<NUM_BASES>(model_i[j], model_i,
						out_matrix_i);
			}
			num_unmasked_data++;
		}
	}

	if (num_unmasked_data < NUM_BASES) {
		throw std::runtime_error(
				"GetLSQFittingMatrixTemplate: too many masked data.");
	}
}

inline void GetLSQFittingMatrix(size_t num_mask, bool const *mask_arg,
		size_t num_model_bases, double const *model_arg, size_t num_lsq_bases,
		double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto mask = AssumeAligned(mask_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_lsq_bases * num_lsq_bases; ++i) {
		out[i] = 0;
	}
	size_t num_unmasked_data = 0;
	for (size_t i = 0; i < num_mask; ++i) {
		if (mask[i]) {
			auto model_i = &model[i * num_model_bases];
			for (size_t j = 0; j < num_lsq_bases; ++j) {
				auto out_matrix_j = &out[j * num_lsq_bases];
				AddMulVector(num_lsq_bases, model_i[j], model_i, out_matrix_j);
			}
			num_unmasked_data++;
		}
	}

	if (num_unmasked_data < num_model_bases) {
		throw std::runtime_error(
				"GetLSQFittingMatrix: too many data are masked.");
	}
}

template<size_t NUM_BASES>
inline void UpdateLSQFittingMatrixTemplate(size_t num_clipped,
		size_t const *clipped_indices_arg, size_t num_model_bases,
		double const *in_arg, double const *model_arg, double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto in = AssumeAligned(in_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < NUM_BASES * NUM_BASES; ++i) {
		out[i] = in[i];
	}
	for (size_t i = 0; i < num_clipped; ++i) {
		auto model_i = &model[clipped_indices[i] * num_model_bases];
		for (size_t j = 0; j < NUM_BASES; ++j) {
			auto out_matrix_j = &out[j * NUM_BASES];
			SubMulVectorTemplate<NUM_BASES>(model_i[j], model_i, out_matrix_j);
		}
	}
}

inline void UpdateLSQFittingMatrix(size_t num_clipped,
		size_t const *clipped_indices_arg, size_t num_lsq_bases,
		double const *in_arg, size_t num_model_bases, double const *model_arg,
		double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto in = AssumeAligned(in_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_lsq_bases * num_lsq_bases; ++i) {
		out[i] = in[i];
	}
	for (size_t i = 0; i < num_clipped; ++i) {
		auto model_i = &model[clipped_indices[i] * num_model_bases];
		for (size_t j = 0; j < num_lsq_bases; ++j) {
			auto out_matrix_j = &out[j * num_lsq_bases];
			SubMulVector(num_lsq_bases, model_i[j], model_i, out_matrix_j);
		}
	}
}

template<size_t NUM_BASES>
inline void GetLSQFittingVectorTemplate(size_t num_data, float const *data_arg,
bool const *mask_arg, size_t num_model_bases, double const *model_arg,
		double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < NUM_BASES; ++i) {
		out[i] = 0;
	}
	for (size_t i = 0; i < num_data; ++i) {
		if (mask[i]) {
			auto model_i = &model[i * num_model_bases];
			auto data_i = data[i];
			AddMulVectorTemplate<NUM_BASES>(data_i, model_i, out);
		}
	}
}

inline void GetLSQFittingVector(size_t num_data, float const *data_arg,
bool const *mask_arg, size_t num_model_bases, double const *model_arg,
		size_t num_lsq_bases, double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_lsq_bases; ++i) {
		out[i] = 0;
	}
	for (size_t i = 0; i < num_data; ++i) {
		if (mask[i]) {
			auto model_i = &model[i * num_model_bases];
			auto data_i = data[i];
			AddMulVector(num_lsq_bases, data_i, model_i, out);
		}
	}
}

template<size_t NUM_BASES> inline void UpdateLSQFittingVectorTemplate(
		float const *data_arg, size_t num_clipped,
		size_t const *clipped_indices_arg, size_t num_model_bases,
		double const *in_arg, double const *model_arg, double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto data = AssumeAligned(data_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto in = AssumeAligned(in_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < NUM_BASES; ++i) {
		out[i] = in[i];
	}
	for (size_t i = 0; i < num_clipped; ++i) {
		auto model_i = &model[clipped_indices[i] * num_model_bases];
		auto data_i = data[clipped_indices[i]];
		SubMulVectorTemplate<NUM_BASES>(data_i, model_i, out);
	}
}

inline void UpdateLSQFittingVector(float const *data_arg, size_t num_clipped,
		size_t const *clipped_indices_arg, size_t num_lsq_bases,
		double const *in_arg, size_t num_model_bases, double const *model_arg,
		double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto data = AssumeAligned(data_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto in = AssumeAligned(in_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_lsq_bases; ++i) {
		out[i] = in[i];
	}
	for (size_t i = 0; i < num_clipped; ++i) {
		auto model_i = &model[clipped_indices[i] * num_model_bases];
		auto data_i = data[clipped_indices[i]];
		SubMulVector(num_lsq_bases, data_i, model_i, out);
	}
}

template<size_t NUM_BASES>
inline void GetLSQCoefficientsTemplate(size_t num_data, float const *data_arg,
		bool const *mask_arg, size_t const num_model_bases,
		double const *basis_data, size_t const num_lsq_bases,
		double *lsq_matrix_arg, double *lsq_vector_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector_arg));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto basis = AssumeAligned(basis_data);
	auto lsq_matrix = AssumeAligned(lsq_matrix_arg);
	auto lsq_vector = AssumeAligned(lsq_vector_arg);

	GetLSQFittingMatrixTemplate<NUM_BASES>(num_data, mask, num_model_bases,
			basis, lsq_matrix);
	GetLSQFittingVectorTemplate<NUM_BASES>(num_data, data, mask,
			num_model_bases, basis, lsq_vector);
}

inline void GetLSQCoefficients(size_t num_data, float const *data_arg,
bool const *mask_arg, size_t const num_model_bases, double const *basis_data,
		size_t const num_lsq_bases, double *lsq_matrix_arg,
		double *lsq_vector_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector_arg));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto basis = AssumeAligned(basis_data);
	auto lsq_matrix = AssumeAligned(lsq_matrix_arg);
	auto lsq_vector = AssumeAligned(lsq_vector_arg);

	GetLSQFittingMatrix(num_data, mask, num_model_bases, basis, num_lsq_bases,
			lsq_matrix);
	GetLSQFittingVector(num_data, data, mask, num_model_bases, basis,
			num_lsq_bases, lsq_vector);
}

template<size_t NUM_BASES>
inline void UpdateLSQCoefficientsTemplate(size_t const num_data,
		float const *data_arg, size_t const num_clipped,
		size_t const *clipped_indices_arg, size_t const num_model_bases,
		double const *basis_data, size_t const num_lsq_bases,
		double *lsq_matrix_arg, double *lsq_vector_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector_arg));
	auto data = AssumeAligned(data_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto basis = AssumeAligned(basis_data);
	auto lsq_matrix = AssumeAligned(lsq_matrix_arg);
	auto lsq_vector = AssumeAligned(lsq_vector_arg);

	UpdateLSQFittingMatrixTemplate<NUM_BASES>(num_clipped, clipped_indices,
			num_model_bases, lsq_matrix, basis, lsq_matrix);
	UpdateLSQFittingVectorTemplate<NUM_BASES>(data, num_clipped,
			clipped_indices, num_model_bases, lsq_vector, basis, lsq_vector);
}

inline void UpdateLSQCoefficients(size_t const num_data, float const *data_arg,
		size_t const num_clipped, size_t const *clipped_indices_arg,
		size_t const num_model_bases, double const *basis_data,
		size_t const num_lsq_bases, double *lsq_matrix_arg,
		double *lsq_vector_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector_arg));
	auto data = AssumeAligned(data_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto basis = AssumeAligned(basis_data);
	auto lsq_matrix = AssumeAligned(lsq_matrix_arg);
	auto lsq_vector = AssumeAligned(lsq_vector_arg);

	UpdateLSQFittingMatrix(num_clipped, clipped_indices, num_lsq_bases,
			lsq_matrix, num_model_bases, basis, lsq_matrix);
	UpdateLSQFittingVector(data, num_clipped, clipped_indices, num_lsq_bases,
			lsq_vector, num_model_bases, basis, lsq_vector);
}

inline void SolveSimultaneousEquationsByLUDouble(size_t num_equations,
		double const *in_matrix_arg, double const *in_vector_arg, double *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_matrix_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_vector_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map<MatrixXd, Aligned, Stride<Dynamic, Dynamic> > in_matrix(
			const_cast<double *>(in_matrix_arg), num_equations, num_equations,
			Stride<Dynamic, Dynamic>(1, num_equations));
	Map < VectorXd, Aligned
			> in_vector(const_cast<double *>(in_vector_arg), num_equations);

	Map < VectorXd > (out, num_equations) = in_matrix.fullPivLu().solve(
			in_vector);
}

#define RepeatTen(func, start_idx) \
		(func<start_idx + 0>), \
		(func<start_idx + 1>), \
		(func<start_idx + 2>), \
		(func<start_idx + 3>), \
		(func<start_idx + 4>), \
		(func<start_idx + 5>), \
		(func<start_idx + 6>), \
		(func<start_idx + 7>), \
		(func<start_idx + 8>), \
		(func<start_idx + 9>)

void GetLSQCoefficientsEntry(size_t const num_data,
		float const data[/*num_data*/],
		bool const mask[/*num_data*/], size_t const num_model_bases,
		double const basis_data[/*num_model_bases*num_data*/],
		size_t const num_lsq_bases,
		double lsq_matrix[/*num_lsq_bases*num_lsq_bases*/],
		double lsq_vector[/*num_lsq_bases*/]) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector));

	typedef void (*GetLSQCoefficientsFunc)(size_t const num_data,
			float const *data, bool const *mask, size_t const num_model_bases,
			double const *basis_data, size_t const num_lsq_bases,
			double *lsq_matrix, double *lsq_vector);

	static GetLSQCoefficientsFunc const funcs[] = { RepeatTen(
			GetLSQCoefficientsTemplate, 0), RepeatTen(
			GetLSQCoefficientsTemplate, 10), RepeatTen(
			GetLSQCoefficientsTemplate, 20), RepeatTen(
			GetLSQCoefficientsTemplate, 30), RepeatTen(
			GetLSQCoefficientsTemplate, 40), RepeatTen(
			GetLSQCoefficientsTemplate, 50), RepeatTen(
			GetLSQCoefficientsTemplate, 60), RepeatTen(
			GetLSQCoefficientsTemplate, 70), RepeatTen(
			GetLSQCoefficientsTemplate, 80), RepeatTen(
			GetLSQCoefficientsTemplate, 90), GetLSQCoefficientsTemplate<100> };

	if (num_lsq_bases < ELEMENTSOF(funcs)) {
		funcs[num_lsq_bases](num_data, data, mask, num_model_bases, basis_data,
				num_lsq_bases, lsq_matrix, lsq_vector);
	} else {
		GetLSQCoefficients(num_data, data, mask, num_model_bases, basis_data,
				num_lsq_bases, lsq_matrix, lsq_vector);
	}
}

void UpdateLSQCoefficientsEntry(size_t const num_data,
		float const data[/*num_data*/], size_t const num_clipped,
		size_t const clipped_indices[/*num_data*/],
		size_t const num_model_bases,
		double const basis_data[/*num_model_bases*num_data*/],
		size_t const num_lsq_bases,
		double lsq_matrix[/*num_lsq_bases*num_lsq_bases*/],
		double lsq_vector[/*num_lsq_bases*/]) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices));
	assert(LIBSAKURA_SYMBOL(IsAligned)(basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector));

	typedef void (*UpdateLSQCoefficientsFunc)(size_t const num_data,
			float const *data, size_t const num_clipped,
			size_t const *clipped_indices, size_t const num_model_bases,
			double const *basis_data, size_t const num_lsq_bases,
			double *lsq_matrix, double *lsq_vector);

	static UpdateLSQCoefficientsFunc const funcs[] = { RepeatTen(
			UpdateLSQCoefficientsTemplate, 0), RepeatTen(
			UpdateLSQCoefficientsTemplate, 10), RepeatTen(
			UpdateLSQCoefficientsTemplate, 20), RepeatTen(
			UpdateLSQCoefficientsTemplate, 30), RepeatTen(
			UpdateLSQCoefficientsTemplate, 40), RepeatTen(
			UpdateLSQCoefficientsTemplate, 50), RepeatTen(
			UpdateLSQCoefficientsTemplate, 60), RepeatTen(
			UpdateLSQCoefficientsTemplate, 70), RepeatTen(
			UpdateLSQCoefficientsTemplate, 80), RepeatTen(
			UpdateLSQCoefficientsTemplate, 90), UpdateLSQCoefficientsTemplate<
			100> };

	if (num_lsq_bases < ELEMENTSOF(funcs)) {
		funcs[num_lsq_bases](num_data, data, num_clipped, clipped_indices,
				num_model_bases, basis_data, num_lsq_bases, lsq_matrix,
				lsq_vector);
	} else {
		UpdateLSQCoefficients(num_data, data, num_clipped, clipped_indices,
				num_model_bases, basis_data, num_lsq_bases, lsq_matrix,
				lsq_vector);
	}
}

} /* anonymous namespace */

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetLSQCoefficientsDouble)(
		size_t const num_data, float const data[], bool const mask[],
		size_t const num_model_bases, double const basis_data[],
		size_t const num_lsq_bases, double lsq_matrix[], double lsq_vector[])
				noexcept {
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
	CHECK_ARGS(lsq_matrix != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix));
	CHECK_ARGS(lsq_vector != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector));

	try {
		GetLSQCoefficientsEntry(num_data, data, mask, num_model_bases,
				basis_data, num_lsq_bases, lsq_matrix, lsq_vector);
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
		size_t const num_data, float const data[],
		size_t const num_exclude_indices, size_t const exclude_indices[],
		size_t const num_model_bases, double const basis_data[],
		size_t const num_lsq_bases, double lsq_matrix[], double lsq_vector[])
				noexcept {
	CHECK_ARGS(num_data != 0);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
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
	CHECK_ARGS(lsq_matrix != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix));
	CHECK_ARGS(lsq_vector != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector));

	try {
		UpdateLSQCoefficientsEntry(num_data, data, num_exclude_indices,
				exclude_indices, num_model_bases, basis_data, num_lsq_bases,
				lsq_matrix, lsq_vector);
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
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));

	try {
		SolveSimultaneousEquationsByLUDouble(num_equations, in_matrix,
				in_vector, out);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}
