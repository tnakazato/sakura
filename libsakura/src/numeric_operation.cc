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

using ::Eigen::Map;
using ::Eigen::MatrixXd;
using ::Eigen::VectorXd;
using ::Eigen::Aligned;

namespace {

auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("numeric_operation");

#if defined(__AVX__)
#if defined(__AVX2__)
#define FMAD(a,b,c)	_mm256_fmadd_pd(a, b, c)
#else
#define FMAD(a,b,c)	_mm256_add_pd(_mm256_mul_pd(a, b), c)
#endif
#define IFMSD(a,b,c)	_mm256_sub_pd(c, _mm256_mul_pd(a, b))
#endif

template<size_t NUM_MODEL_BASES>
void AddMulVectorTemplate(double k, double const *vec, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (NUM_MODEL_BASES / pack_elements) * pack_elements;
	__m256d coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += pack_elements) {
		__m256d v = _mm256_loadu_pd(&vec[i]);
		_mm256_storeu_pd(&out[i], FMAD(coeff, v, _mm256_loadu_pd(&out[i])));
	}
#endif
	for (; i < NUM_MODEL_BASES; ++i) {
		out[i] += k * vec[i];
	}
}

template<size_t NUM_MODEL_BASES>
void SubMulVectorTemplate(double k, double const *vec, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (NUM_MODEL_BASES / pack_elements) * pack_elements;
	__m256d coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += pack_elements) {
		__m256d v = _mm256_loadu_pd(&vec[i]);
		_mm256_storeu_pd(&out[i], IFMSD(coeff, v, _mm256_loadu_pd(&out[i])));
	}
#endif
	for (; i < NUM_MODEL_BASES; ++i) {
		out[i] -= k * vec[i];
	}
}

inline void AddMulVector(size_t const num_model_bases, double k,
		double const *vec, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (num_model_bases / pack_elements) * pack_elements;
	__m256d coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += pack_elements) {
		__m256d v = _mm256_loadu_pd(&vec[i]);
		_mm256_storeu_pd(&out[i], FMAD(coeff, v, _mm256_loadu_pd(&out[i])));
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
	__m256d coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += pack_elements) {
		__m256d v = _mm256_loadu_pd(&vec[i]);
		_mm256_storeu_pd(&out[i], IFMSD(coeff, v, _mm256_loadu_pd(&out[i])));
	}
#endif
	for (; i < num_model_bases; ++i) {
		out[i] -= k * vec[i];
	}
}

template<size_t NUM_MODEL_BASES>
inline void GetMatrixCoefficientsForLeastSquareFittingTemplate(size_t num_mask,
bool const *mask_arg, size_t num_model_bases, double const *model_arg,
		double *out_arg) {
	auto mask = AssumeAligned(mask_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < NUM_MODEL_BASES * NUM_MODEL_BASES; ++i) {
		out[i] = 0;
	}
	size_t num_unmasked_data = 0;
	for (size_t i = 0; i < num_mask; ++i) {
		if (mask[i]) {
			auto model_i = &model[i * NUM_MODEL_BASES];
			for (size_t j = 0; j < NUM_MODEL_BASES; ++j) {
				auto out_matrix_i = &out[j * NUM_MODEL_BASES];
				AddMulVectorTemplate<NUM_MODEL_BASES>(model_i[j], model_i,
						out_matrix_i);
			}
			num_unmasked_data++;
		}
	}

	if (num_unmasked_data < NUM_MODEL_BASES) {
		throw std::runtime_error(
				"GetMatrixCoefficientsForLeastSquareFittingTemplate: too many masked data.");
	}
}

inline void GetMatrixCoefficientsForLeastSquareFitting(size_t num_mask,
bool const *mask_arg, size_t num_model_bases, double const *model_arg,
		double *out_arg) {
	auto mask = AssumeAligned(mask_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_model_bases * num_model_bases; ++i) {
		out[i] = 0;
	}
	size_t num_unmasked_data = 0;
	for (size_t i = 0; i < num_mask; ++i) {
		if (mask[i]) {
			auto model_i = &model[i * num_model_bases];
			for (size_t j = 0; j < num_model_bases; ++j) {
				auto out_matrix_j = &out[j * num_model_bases];
				AddMulVector(num_model_bases, model_i[j], model_i,
						out_matrix_j);
			}
			num_unmasked_data++;
		}
	}

	if (num_unmasked_data < num_model_bases) {
		throw std::runtime_error(
				"GetMatrixCoefficientsForLeastSquareFitting: too many data are masked.");
	}
}

template<size_t NUM_MODEL_BASES>
inline void UpdateMatrixCoefficientsForLeastSquareFittingTemplate(
		size_t num_clipped, size_t const *clipped_indices_arg,
		size_t num_model_bases, double const *in_arg, double const *model_arg,
		double *out_arg) {
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto in = AssumeAligned(in_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < NUM_MODEL_BASES * NUM_MODEL_BASES; ++i) {
		out[i] = in[i];
	}
	for (size_t i = 0; i < num_clipped; ++i) {
		auto model_i = &model[clipped_indices[i] * NUM_MODEL_BASES];
		for (size_t j = 0; j < NUM_MODEL_BASES; ++j) {
			auto out_matrix_j = &out[j * NUM_MODEL_BASES];
			SubMulVectorTemplate<NUM_MODEL_BASES>(model_i[j], model_i,
					out_matrix_j);
		}
	}
}

inline void UpdateMatrixCoefficientsForLeastSquareFitting(size_t num_clipped,
		size_t const *clipped_indices_arg, size_t num_model_bases,
		double const *in_arg, double const *model_arg, double *out_arg) {
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto in = AssumeAligned(in_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_model_bases * num_model_bases; ++i) {
		out[i] = in[i];
	}
	for (size_t i = 0; i < num_clipped; ++i) {
		auto model_i = &model[clipped_indices[i] * num_model_bases];
		for (size_t j = 0; j < num_model_bases; ++j) {
			auto out_matrix_j = &out[j * num_model_bases];
			SubMulVector(num_model_bases, model_i[j], model_i, out_matrix_j);
		}
	}
}

template<size_t NUM_MODEL_BASES>
inline void GetVectorCoefficientsForLeastSquareFittingTemplate(size_t num_data,
		float const *data_arg, bool const *mask_arg, size_t num_model_bases,
		double const *model_arg, double *out_arg) {
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < NUM_MODEL_BASES; ++i) {
		out[i] = 0;
	}
	for (size_t i = 0; i < num_data; ++i) {
		if (mask[i]) {
			auto model_i = &model[i * NUM_MODEL_BASES];
			auto data_i = data[i];
			AddMulVectorTemplate<NUM_MODEL_BASES>(data_i, model_i, out);
		}
	}
}

inline void GetVectorCoefficientsForLeastSquareFitting(size_t num_data,
		float const *data_arg, bool const *mask_arg, size_t num_model_bases,
		double const *model_arg, double *out_arg) {
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_model_bases; ++i) {
		out[i] = 0;
	}
	for (size_t i = 0; i < num_data; ++i) {
		if (mask[i]) {
			auto model_i = &model[i * num_model_bases];
			auto data_i = data[i];
			AddMulVector(num_model_bases, data_i, model_i, out);
		}
	}
}

template<size_t NUM_MODEL_BASES>
inline void UpdateVectorCoefficientsForLeastSquareFittingTemplate(
		float const *data_arg, size_t num_clipped,
		size_t const *clipped_indices_arg, size_t num_model_bases,
		double const *in_arg, double const *model_arg, double *out_arg) {
	auto data = AssumeAligned(data_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto in = AssumeAligned(in_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < NUM_MODEL_BASES; ++i) {
		out[i] = in[i];
	}
	for (size_t i = 0; i < num_clipped; ++i) {
		auto model_i = &model[clipped_indices[i] * NUM_MODEL_BASES];
		auto data_i = data[clipped_indices[i]];
		SubMulVectorTemplate<NUM_MODEL_BASES>(data_i, model_i, out);
	}
}

inline void UpdateVectorCoefficientsForLeastSquareFitting(float const *data_arg,
		size_t num_clipped, size_t const *clipped_indices_arg,
		size_t num_model_bases, double const *in_arg, double const *model_arg,
		double *out_arg) {
	auto data = AssumeAligned(data_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto in = AssumeAligned(in_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_model_bases; ++i) {
		out[i] = in[i];
	}
	for (size_t i = 0; i < num_clipped; ++i) {
		auto model_i = &model[clipped_indices[i] * num_model_bases];
		auto data_i = data[clipped_indices[i]];
		SubMulVector(num_model_bases, data_i, model_i, out);
	}
}

template<size_t NUM_MODEL_BASES>
inline void GetLeastSquareFittingCoefficientsTemplate(size_t num_data,
		float const *data_arg, bool const *mask_arg,
		size_t const num_model_bases, double const *basis_data,
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

	GetMatrixCoefficientsForLeastSquareFittingTemplate<NUM_MODEL_BASES>(
			num_data, mask, num_model_bases, basis, lsq_matrix);
	GetVectorCoefficientsForLeastSquareFittingTemplate<NUM_MODEL_BASES>(
			num_data, data, mask, num_model_bases, basis, lsq_vector);
}

inline void GetLeastSquareFittingCoefficients(size_t num_data,
		float const *data_arg, bool const *mask_arg,
		size_t const num_model_bases, double const *basis_data,
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

	GetMatrixCoefficientsForLeastSquareFitting(num_data, mask, num_model_bases,
			basis, lsq_matrix);
	GetVectorCoefficientsForLeastSquareFitting(num_data, data, mask,
			num_model_bases, basis, lsq_vector);
}

template<size_t NUM_MODEL_BASES>
inline void UpdateLeastSquareFittingCoefficientsTemplate(size_t const num_data,
		float const *data_arg, size_t const num_clipped,
		size_t const *clipped_indices_arg, size_t const num_model_bases,
		double const *basis_data, double *lsq_matrix_arg,
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

	UpdateMatrixCoefficientsForLeastSquareFittingTemplate<NUM_MODEL_BASES>(
			num_clipped, clipped_indices, num_model_bases, lsq_matrix, basis,
			lsq_matrix);
	UpdateVectorCoefficientsForLeastSquareFittingTemplate<NUM_MODEL_BASES>(data,
			num_clipped, clipped_indices, num_model_bases, lsq_vector, basis,
			lsq_vector);
}

inline void UpdateLeastSquareFittingCoefficients(size_t const num_data,
		float const *data_arg, size_t const num_clipped,
		size_t const *clipped_indices_arg, size_t const num_model_bases,
		double const *basis_data, double *lsq_matrix_arg,
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

	UpdateMatrixCoefficientsForLeastSquareFitting(num_clipped, clipped_indices,
			num_model_bases, lsq_matrix, basis, lsq_matrix);
	UpdateVectorCoefficientsForLeastSquareFitting(data, num_clipped,
			clipped_indices, num_model_bases, lsq_vector, basis, lsq_vector);
}

inline void SolveSimultaneousEquationsByLUDouble(size_t num_equations,
		double const *in_matrix_arg, double const *in_vector_arg, double *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_matrix_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_vector_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map < MatrixXd, Aligned
			> in_matrix(const_cast<double *>(in_matrix_arg), num_equations,
					num_equations);
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

void GetLeastSquareFittingCoefficientsEntry(size_t const num_data,
		float const data[/*num_data*/],
		bool const mask[/*num_data*/], size_t const num_model_bases,
		double const basis_data[/*num_model_bases*num_data*/],
		double lsq_matrix[/*num_model_bases*num_model_bases*/],
		double lsq_vector[/*num_model_bases*/]) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector));

	typedef void (*GetLeastSquareFittingCoefficientsFunc)(size_t const num_data,
			float const *data, bool const *mask, size_t const num_model_bases,
			double const *basis_data, double *lsq_matrix, double *lsq_vector);

	static GetLeastSquareFittingCoefficientsFunc const funcs[] = { RepeatTen(
			GetLeastSquareFittingCoefficientsTemplate, 0), RepeatTen(
			GetLeastSquareFittingCoefficientsTemplate, 10), RepeatTen(
			GetLeastSquareFittingCoefficientsTemplate, 20), RepeatTen(
			GetLeastSquareFittingCoefficientsTemplate, 30), RepeatTen(
			GetLeastSquareFittingCoefficientsTemplate, 40), RepeatTen(
			GetLeastSquareFittingCoefficientsTemplate, 50), RepeatTen(
			GetLeastSquareFittingCoefficientsTemplate, 60), RepeatTen(
			GetLeastSquareFittingCoefficientsTemplate, 70), RepeatTen(
			GetLeastSquareFittingCoefficientsTemplate, 80), RepeatTen(
			GetLeastSquareFittingCoefficientsTemplate, 90),
			GetLeastSquareFittingCoefficientsTemplate<100> };

	if (num_model_bases < ELEMENTSOF(funcs)) {
		funcs[num_model_bases](num_data, data, mask, num_model_bases,
				basis_data, lsq_matrix, lsq_vector);
	} else {
		::GetLeastSquareFittingCoefficients(num_data, data, mask,
				num_model_bases, basis_data, lsq_matrix, lsq_vector);
	}
}

void UpdateLeastSquareFittingCoefficientsEntry(size_t const num_data,
		float const data[/*num_data*/], size_t const num_clipped,
		size_t const clipped_indices[/*num_data*/],
		size_t const num_model_bases,
		double const basis_data[/*num_model_bases*num_data*/],
		double lsq_matrix[/*num_model_bases*num_model_bases*/],
		double lsq_vector[/*num_model_bases*/]) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices));
	assert(LIBSAKURA_SYMBOL(IsAligned)(basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector));

	typedef void (*UpdateLeastSquareFittingCoefficientsFunc)(
			size_t const num_data, float const *data, size_t const num_clipped,
			size_t const *clipped_indices, size_t const num_model_bases,
			double const *basis_data, double *lsq_matrix, double *lsq_vector);

	static UpdateLeastSquareFittingCoefficientsFunc const funcs[] = { RepeatTen(
			UpdateLeastSquareFittingCoefficientsTemplate, 0), RepeatTen(
			UpdateLeastSquareFittingCoefficientsTemplate, 10), RepeatTen(
			UpdateLeastSquareFittingCoefficientsTemplate, 20), RepeatTen(
			UpdateLeastSquareFittingCoefficientsTemplate, 30), RepeatTen(
			UpdateLeastSquareFittingCoefficientsTemplate, 40), RepeatTen(
			UpdateLeastSquareFittingCoefficientsTemplate, 50), RepeatTen(
			UpdateLeastSquareFittingCoefficientsTemplate, 60), RepeatTen(
			UpdateLeastSquareFittingCoefficientsTemplate, 70), RepeatTen(
			UpdateLeastSquareFittingCoefficientsTemplate, 80), RepeatTen(
			UpdateLeastSquareFittingCoefficientsTemplate, 90),
			UpdateLeastSquareFittingCoefficientsTemplate<100> };

	if (num_model_bases < ELEMENTSOF(funcs)) {
		funcs[num_model_bases](num_data, data, num_clipped, clipped_indices,
				num_model_bases, basis_data, lsq_matrix, lsq_vector);
	} else {
		::UpdateLeastSquareFittingCoefficients(num_data, data, num_clipped,
				clipped_indices, num_model_bases, basis_data, lsq_matrix,
				lsq_vector);
	}
}

} /* anonymous namespace */

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetLeastSquareFittingCoefficientsDouble)(
		size_t const num_data, float const data[], bool const mask[],
		size_t const num_model_bases, double const basis_data[],
		double lsq_matrix[], double lsq_vector[]) {
	if (num_data == 0)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_model_bases == 0)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_model_bases > num_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (basis_data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(basis_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lsq_matrix == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lsq_vector == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lsq_vector)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	try {
		GetLeastSquareFittingCoefficientsEntry(num_data, data, mask,
				num_model_bases, basis_data, lsq_matrix, lsq_vector);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UpdateLeastSquareFittingCoefficientsDouble)(
		size_t const num_data, float const data[],
		size_t const num_exclude_indices, size_t const exclude_indices[],
		size_t const num_model_bases, double const basis_data[],
		double lsq_matrix[], double lsq_vector[]) {
	if (num_data == 0)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (exclude_indices == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(exclude_indices)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_exclude_indices > num_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		if (exclude_indices[i] >= num_data)
			return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (num_model_bases == 0)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (num_model_bases > num_data)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (basis_data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(basis_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lsq_matrix == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lsq_vector == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lsq_vector)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	try {
		UpdateLeastSquareFittingCoefficientsEntry(num_data, data,
				num_exclude_indices, exclude_indices, num_model_bases,
				basis_data, lsq_matrix, lsq_vector);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLUDouble)(
		size_t num_equations, double const in_matrix[],
		double const in_vector[], double out[]) {
	if (in_matrix == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in_matrix)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (in_vector == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in_vector)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	try {
		SolveSimultaneousEquationsByLUDouble(num_equations, in_matrix,
				in_vector, out);
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}
