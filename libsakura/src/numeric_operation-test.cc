/*
 * @SAKURA_LICENSE_HEADER_START@
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

#define TUNE 1

#include <libsakura/localdef.h>
#include <libsakura/memory_manager.h>
#include <libsakura/optimized_implementation_factory_impl.h>
#include <libsakura/sakura.h>

#include <Eigen/Core>
#include <Eigen/LU>

using ::Eigen::Map;
using ::Eigen::MatrixXd;
using ::Eigen::VectorXd;
using ::Eigen::Aligned;

namespace {

template<size_t NUM_MODEL_BASES>
void AddMulVectorTemplate(double k, double const *vec, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR) && TUNE
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (NUM_MODEL_BASES / pack_elements) * pack_elements;
	__m256d coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += pack_elements) {
		__m256d v = _mm256_loadu_pd(&vec[i]);
		_mm256_storeu_pd(&out[i], _mm256_loadu_pd(&out[i]) + coeff * v);
	}
	if (NUM_MODEL_BASES % pack_elements == 0) {
		return;
	}
#endif
	for (; i < NUM_MODEL_BASES; ++i) {
		out[i] += k * vec[i];
	}
}

template<size_t NUM_MODEL_BASES>
void SubMulVectorTemplate(double k, double const *vec, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR) && TUNE
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (NUM_MODEL_BASES / pack_elements) * pack_elements;
	__m256d coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += pack_elements) {
		__m256d v = _mm256_loadu_pd(&vec[i]);
		_mm256_storeu_pd(&out[i], _mm256_loadu_pd(&out[i]) - coeff * v);
	}
	if (NUM_MODEL_BASES % pack_elements == 0) {
		return;
	}
#endif
	for (; i < NUM_MODEL_BASES; ++i) {
		out[i] -= k * vec[i];
	}
}

inline void AddMulVector(size_t const num_model_bases, double k,
		double const *vec, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR) && TUNE
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (num_model_bases / pack_elements) * pack_elements;
	__m256d coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += pack_elements) {
		__m256d v = _mm256_loadu_pd(&vec[i]);
		_mm256_storeu_pd(&out[i], _mm256_loadu_pd(&out[i]) + coeff * v);
	}
#endif
	for (; i < num_model_bases; ++i) {
		out[i] += k * vec[i];
	}
}

inline void SubMulVector(size_t const num_model_bases, double k,
		double const *vec, double *out) {
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR) && TUNE
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (num_model_bases / pack_elements) * pack_elements;
	__m256d coeff = _mm256_set1_pd(k);
	for (i = 0; i < end; i += pack_elements) {
		__m256d v = _mm256_loadu_pd(&vec[i]);
		_mm256_storeu_pd(&out[i], _mm256_loadu_pd(&out[i]) - coeff * v);
	}
#endif
	for (; i < num_model_bases; ++i) {
		out[i] -= k * vec[i];
	}
}

template<size_t NUM_MODEL_BASES>
inline void GetMatrixCoefficientsForLeastSquareFittingUsingTemplate(
bool const update_on_incremental_clipping, size_t num_mask,
bool const *mask_arg, size_t num_clipped, size_t const *clipped_indices_arg,
		size_t num_model_bases, double const *in_arg, double const *model_arg,
		double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));

	auto mask = AssumeAligned(mask_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto in = AssumeAligned(in_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	if (update_on_incremental_clipping) {
		for (size_t i = 0; i < NUM_MODEL_BASES * NUM_MODEL_BASES; ++i) {
			out[i] = in[i];
		}
		for (size_t i = 0; i < num_clipped; ++i) {
			auto model_i = &model[clipped_indices[i] * NUM_MODEL_BASES];
			for (size_t j = 0; j < NUM_MODEL_BASES; ++j) {
				auto out_matrix_i = &out[j * NUM_MODEL_BASES];
				SubMulVectorTemplate<NUM_MODEL_BASES>(model_i[j], model_i,
						out_matrix_i);
			}
		}
	} else {
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
					"GetMatrixCoefficientsForLeastSquareFittingUsingTemplate: too many masked data.");
		}
	}
}

inline void GetMatrixCoefficientsForLeastSquareFitting(
bool const update_on_incremental_clipping, size_t num_mask,
bool const *mask_arg, size_t num_clipped, size_t const *clipped_indices_arg,
		size_t num_model_bases, double const *in_arg, double const *model_arg,
		double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));

	auto mask = AssumeAligned(mask_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto in = AssumeAligned(in_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	if (update_on_incremental_clipping) {
		for (size_t i = 0; i < num_model_bases * num_model_bases; ++i) {
			out[i] = in[i];
		}

		for (size_t i = 0; i < num_clipped; ++i) {
			auto model_i = &model[clipped_indices[i] * num_model_bases];
			for (size_t j = 0; j < num_model_bases; ++j) {
				auto out_matrix_j = &out[j * num_model_bases];
				SubMulVector(num_model_bases, model_i[j], model_i,
						out_matrix_j);
			}
		}
	} else {
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
}

template<size_t NUM_MODEL_BASES>
inline void GetVectorCoefficientsForLeastSquareFittingUsingTemplate(
bool const update_on_incremental_clipping, size_t num_data,
		float const *data_arg, bool const *mask_arg, size_t num_clipped,
		size_t const *clipped_indices_arg, size_t num_model_bases,
		double const *in_arg, double const *model_arg, double *out_arg) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));

	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto in = AssumeAligned(in_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	if (update_on_incremental_clipping) {
		for (size_t i = 0; i < NUM_MODEL_BASES; ++i) {
			out[i] = in[i];
		}
		for (size_t i = 0; i < num_clipped; ++i) {
			auto model_i = &model[clipped_indices[i] * NUM_MODEL_BASES];
			auto data_i = data[clipped_indices[i]];
			for (size_t j = 0; j < NUM_MODEL_BASES; ++j) {
				out[j] -= model_i[j] * data_i;
			}
		}
	} else {
		for (size_t i = 0; i < NUM_MODEL_BASES; ++i) {
			out[i] = 0;
		}
		for (size_t i = 0; i < num_data; ++i) {
			if (mask[i]) {
				auto model_i = &model[i * NUM_MODEL_BASES];
				auto data_i = data[i];
				for (size_t j = 0; j < NUM_MODEL_BASES; ++j) {
					out[j] += model_i[j] * data_i;
				}
			}
		}
	}
}

inline void GetVectorCoefficientsForLeastSquareFitting(
bool const update_on_incremental_clipping, size_t num_data,
		float const *data_arg, bool const *mask_arg, size_t num_clipped,
		size_t const *clipped_indices_arg, size_t num_model_bases,
		double const *in_arg, double const *model_arg, double *out_arg) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));

	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto in = AssumeAligned(in_arg);
	auto model = AssumeAligned(model_arg);
	auto out = AssumeAligned(out_arg);

	if (update_on_incremental_clipping) {
		for (size_t i = 0; i < num_model_bases; ++i) {
			out[i] = in[i];
		}
		for (size_t i = 0; i < num_clipped; ++i) {
			auto model_i = &model[clipped_indices[i] * num_model_bases];
			auto data_i = data[clipped_indices[i]];
			for (size_t j = 0; j < num_model_bases; ++j) {
				out[j] -= model_i[j] * data_i;
			}
		}
	} else {
		for (size_t i = 0; i < num_model_bases; ++i) {
			out[i] = 0;
		}
		for (size_t i = 0; i < num_data; ++i) {
			if (mask[i]) {
				auto model_i = &model[i * num_model_bases];
				auto data_i = data[i];
				for (size_t j = 0; j < num_model_bases; ++j) {
					out[j] += model_i[j] * data_i;
				}
			}
		}
	}
}

inline void SolveSimultaneousEquationsByLU(size_t num_equations,
		double const *in_matrix_arg, double const *in_vector_arg, double *out) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(in_matrix_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_vector_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map<MatrixXd, Aligned> in_matrix(const_cast<double *>(in_matrix_arg),
			num_equations, num_equations);
	Map<VectorXd, Aligned> in_vector(const_cast<double *>(in_vector_arg),
			num_equations);

	Map<VectorXd>(out, num_equations) = in_matrix.fullPivLu().solve(in_vector);
}

} /* anonymous namespace */

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

namespace LIBSAKURA_PREFIX {

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::GetMatrixCoefficientsForLeastSquareFitting(
bool const update_on_incremental_clipping, size_t num_mask,
bool const mask[/*num_mask*/], size_t num_clipped,
		size_t const clipped_indices[/*num_mask*/], size_t num_model_bases,
		double const in[/*num_model_bases*num_model_bases*/],
		double const model[/*num_model_bases*num_mask*/],
		double out[/*num_model_bases*num_model_bases*/]) const {

	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	assert(num_model_bases > 0);

	typedef void (*GetMatrixCoefficientsForLeastSquareFittingFunc)(
	bool const update_on_incremental_clipping, size_t num_mask,
	bool const *mask, size_t num_clipped, size_t const *clipped_indices,
			size_t num_model_bases, double const *in, double const *model,
			double *out);

	static GetMatrixCoefficientsForLeastSquareFittingFunc const funcs[] =
			{ GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<0>,
					GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<1>,
					GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<2>,
					GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<3>,
					GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<4>,
					GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<5>,
					GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<6>,
					GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<7>,
					GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<8>,
					//::GetMatrixCoefficientsForLeastSquareFitting, //non-template version faster for the case num_bases == 8.
					GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<9>,
					RepeatTen(GetMatrixCoefficientsForLeastSquareFittingUsingTemplate, 10),
					RepeatTen(GetMatrixCoefficientsForLeastSquareFittingUsingTemplate, 20),
					RepeatTen(GetMatrixCoefficientsForLeastSquareFittingUsingTemplate, 30),
					RepeatTen(GetMatrixCoefficientsForLeastSquareFittingUsingTemplate, 40),
					RepeatTen(GetMatrixCoefficientsForLeastSquareFittingUsingTemplate, 50),
					RepeatTen(GetMatrixCoefficientsForLeastSquareFittingUsingTemplate, 60),
					RepeatTen(GetMatrixCoefficientsForLeastSquareFittingUsingTemplate, 70),
					RepeatTen(GetMatrixCoefficientsForLeastSquareFittingUsingTemplate, 80),
					RepeatTen(GetMatrixCoefficientsForLeastSquareFittingUsingTemplate, 90),
					GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<100> };

	if (num_model_bases < ELEMENTSOF(funcs)) {
		funcs[num_model_bases](update_on_incremental_clipping, num_mask, mask,
				num_clipped, clipped_indices, num_model_bases, in, model, out);
	} else {
		::GetMatrixCoefficientsForLeastSquareFitting(
				update_on_incremental_clipping, num_mask, mask, num_clipped,
				clipped_indices, num_model_bases, in, model, out);
	}
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::GetVectorCoefficientsForLeastSquareFitting(
bool const update_on_incremental_clipping, size_t num_data,
		float const data[/*num_data*/],
		bool const mask[/*num_data*/], size_t num_clipped,
		size_t const clipped_indices[/*num_data*/], size_t num_model_bases,
		double const in[/*num_model_bases*/],
		double const model[/*num_model_bases*num_data*/],
		double out[/*num_model_bases*/]) const {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	typedef void (*GetVectorCoefficientsForLeastSquareFittingFunc)(
	bool const update_on_incremental_clipping, size_t num_data,
			float const *data, bool const *mask, size_t num_clipped,
			size_t const *clipped_indices, size_t num_model_bases,
			double const *in, double const *model, double *out);
	static GetVectorCoefficientsForLeastSquareFittingFunc const funcs[] = {
	RepeatTen(GetVectorCoefficientsForLeastSquareFittingUsingTemplate, 0),
	RepeatTen(GetVectorCoefficientsForLeastSquareFittingUsingTemplate, 10),
	RepeatTen(GetVectorCoefficientsForLeastSquareFittingUsingTemplate, 20),
	RepeatTen(GetVectorCoefficientsForLeastSquareFittingUsingTemplate, 30),
	RepeatTen(GetVectorCoefficientsForLeastSquareFittingUsingTemplate, 40),
	RepeatTen(GetVectorCoefficientsForLeastSquareFittingUsingTemplate, 50),
	RepeatTen(GetVectorCoefficientsForLeastSquareFittingUsingTemplate, 60),
	RepeatTen(GetVectorCoefficientsForLeastSquareFittingUsingTemplate, 70),
	RepeatTen(GetVectorCoefficientsForLeastSquareFittingUsingTemplate, 80),
	RepeatTen(GetVectorCoefficientsForLeastSquareFittingUsingTemplate, 90),
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<100> };

	if (num_model_bases < ELEMENTSOF(funcs)) {
		funcs[num_model_bases](update_on_incremental_clipping, num_data, data,
				mask, num_clipped, clipped_indices, num_model_bases, in, model,
				out);
	} else {
		::GetVectorCoefficientsForLeastSquareFitting(
				update_on_incremental_clipping, num_data, data, mask,
				num_clipped, clipped_indices, num_model_bases, in, model, out);
	}
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::SolveSimultaneousEquationsByLU(
		size_t num_equations,
		double const in_matrix[/*num_equations*num_equations*/],
		double const in_vector[/*num_equations*/],
		double out[/*num_equations*/]) const {
	::SolveSimultaneousEquationsByLU(num_equations, in_matrix, in_vector, out);
}

}
