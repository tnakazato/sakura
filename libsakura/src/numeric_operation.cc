#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>

#include <libsakura/localdef.h>
#include <libsakura/memory_manager.h>
#include <libsakura/optimized_implementation_factory_impl.h>
#include <libsakura/sakura.h>

#define EIGEN_DENSEBASE_PLUGIN "eigen_binary_visitor_plugin.h"
#include <Eigen/Core>
#include <Eigen/LU>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

template<size_t NUM_MODEL_BASES>
inline void GetMatrixCoefficientsForLeastSquareFittingScalar(
		size_t num_mask, bool const *mask,
		double const *model, double *out) {

	size_t num_out = NUM_MODEL_BASES * NUM_MODEL_BASES;
	for (size_t i = 0; i < num_out; ++i) {
		out[i] = 0.0;
	}

	for (size_t k = 0; k < num_mask ; ++k){
		if (!mask[k]) continue;

		size_t idx_ij = 0;
		size_t idx_j = k;
		for (size_t i = 0; i < NUM_MODEL_BASES; ++i) {
			idx_ij += i;
			size_t idx_diff = 0;
			size_t idx_i = k + num_mask * i;

			for (size_t j = i; j < NUM_MODEL_BASES; ++j) {

				out[idx_ij] += model[idx_i] * model[idx_j];
				out[idx_ij + idx_diff] = out[idx_ij];

				idx_diff += (NUM_MODEL_BASES - 1);
				idx_ij ++;
				idx_i += num_mask;
			}

			idx_j += num_mask;
		}
	}

}

inline void GetMatrixCoefficientsForLeastSquareFitting(
		size_t num_mask, bool const *mask,
		size_t num_model_bases, double const *model, double *out) {

/*
	for (size_t i = 0; i < num_model_bases; ++i) {
		for (size_t j = 0; j < num_model_bases; ++j) {
			size_t idx = num_model_bases * i + j;
			out[idx] = 0.0;

			for (size_t k = 0; k < num_mask ; ++k){
				if (!mask[k]) continue;

				size_t idx_i = num_mask * i + k;
				size_t idx_j = num_mask * j + k;
				out[idx] += model[idx_i] * model[idx_j];
			}
		}
	}
*/

	size_t num_out = num_model_bases * num_model_bases;
	for (size_t i = 0; i < num_out; ++i) {
		out[i] = 0.0;
	}

	for (size_t k = 0; k < num_mask ; ++k){
		if (!mask[k]) continue;

		size_t idx_ij = 0;
		size_t idx_j = k;
		for (size_t i = 0; i < num_model_bases; ++i) {
			idx_ij += i;
			size_t idx_diff = 0;
			size_t idx_i = k + num_mask * i;

			for (size_t j = i; j < num_model_bases; ++j) {

				out[idx_ij] += model[idx_i] * model[idx_j];
				out[idx_ij + idx_diff] = out[idx_ij];

				idx_diff += (num_model_bases - 1);
				idx_ij ++;
				idx_i += num_mask;
			}

			idx_j += num_mask;
		}
	}

}

template<size_t NUM_MODEL_BASES>
inline void GetVectorCoefficientsForLeastSquareFittingScalar(
		size_t num_data, float const *data, bool const *mask,
		double const *model, double *out) {

	for (size_t i = 0; i < NUM_MODEL_BASES; ++i) {
		out[i] = 0.0;
	}

	for (size_t j = 0; j < num_data ; ++j){
		if (!mask[j]) continue;

		size_t k = j;
		for (size_t i = 0; i < NUM_MODEL_BASES; ++i) {
			out[i] += data[j] * model[k];
			k += num_data;
		}
	}
}

inline void GetVectorCoefficientsForLeastSquareFitting(
		size_t num_data, float const *data, bool const *mask,
		size_t num_model_bases, double const *model, double *out) {

/*
    // original
	for (size_t i = 0; i < num_model_bases; ++i) {
		out[i] = 0.0;
		for (size_t k = 0; k < num_data; ++k) {
			if (!mask[k]) continue;

			size_t idx = num_data * i + k;
			out[i] += model[idx] * data[k];
		}
	}
*/

	for (size_t i = 0; i < num_model_bases; ++i) {
		out[i] = 0.0;
	}

	for (size_t k = 0; k < num_data ; ++k){
		if (!mask[k]) continue;

		size_t idx_j = k;
		for (size_t i = 0; i < num_model_bases; ++i) {
			out[i] += model[idx_j] * data[k];
			idx_j += num_data;
		}
	}
}

inline void SolveSimultaneousEquationsByLU(size_t num_equations,
		double const *lsq_matrix0, double const *lsq_vector0, double *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix0));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector0));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map<Array<double, Dynamic, 1>, Aligned> lsq_matrix0_array(const_cast<double *>(lsq_matrix0),
			num_equations*num_equations);
	Map<Array<double, Dynamic, 1>, Aligned> lsq_vector0_array(const_cast<double *>(lsq_vector0),
			num_equations);
	Map<Array<double, Dynamic, 1>, Aligned> out_array(const_cast<double *>(out),
			num_equations);

	::Eigen::MatrixXd lsq_matrix(num_equations, num_equations);
	::Eigen::VectorXd lsq_vector(num_equations);
	for (size_t i = 0; i < num_equations; ++i) {
		for (size_t j = 0; j < num_equations; ++j) {
			lsq_matrix(i, j) = lsq_matrix0_array[num_equations*i+j];
		}
		lsq_vector(i) = lsq_vector0_array[i];
	}

	::Eigen::VectorXd lu_result = lsq_matrix.fullPivLu().solve(lsq_vector);
	for (size_t i = 0; i < num_equations; ++i) {
		out_array[i] = lu_result(i);
	}
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::GetMatrixCoefficientsForLeastSquareFitting(
		size_t num_mask, bool const mask[/*num_mask*/],
		size_t num_model_bases, double const model[/*num_model_bases*num_mask*/],
		double out[/*num_model_bases*num_model_bases*/]) const {

	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	typedef void (*GetMatrixCoefficientsForLeastSquareFittingFunc)(
			size_t num_mask, bool const *mask,
			double const *model, double *out);

	static GetMatrixCoefficientsForLeastSquareFittingFunc const funcs[] = {
			GetMatrixCoefficientsForLeastSquareFittingScalar<0>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<1>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<2>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<3>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<4>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<5>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<6>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<7>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<8>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<9>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<10>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<11>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<12>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<13>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<14>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<15>,
			GetMatrixCoefficientsForLeastSquareFittingScalar<16> };

	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);

	if (num_model_bases < ELEMENTSOF(funcs)) {
		funcs[num_model_bases](num_mask, mask, model, out);
	} else {
		GetMatrixCoefficientsForLeastSquareFitting(
				num_mask, mask, num_model_bases, model, out);
	}
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::GetVectorCoefficientsForLeastSquareFitting(
		size_t num_data, float const data[/*num_data*/], bool const mask[/*num_data*/],
		size_t num_model_bases, double const model[/*num_model_bases*num_data*/],
		double out[/*num_model_bases*/]) const {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	typedef void (*GetVectorCoefficientsForLeastSquareFittingFunc)(
			size_t num_data, float const *data, bool const *mask,
			double const *model, double *out);
	static GetVectorCoefficientsForLeastSquareFittingFunc const funcs[] = {
			GetVectorCoefficientsForLeastSquareFittingScalar<0>,
			GetVectorCoefficientsForLeastSquareFittingScalar<1>,
			GetVectorCoefficientsForLeastSquareFittingScalar<2>,
			GetVectorCoefficientsForLeastSquareFittingScalar<3>,
			GetVectorCoefficientsForLeastSquareFittingScalar<4>,
			GetVectorCoefficientsForLeastSquareFittingScalar<5>,
			GetVectorCoefficientsForLeastSquareFittingScalar<6>,
			GetVectorCoefficientsForLeastSquareFittingScalar<7>,
			GetVectorCoefficientsForLeastSquareFittingScalar<8>,
			GetVectorCoefficientsForLeastSquareFittingScalar<9>,
			GetVectorCoefficientsForLeastSquareFittingScalar<10>,
			GetVectorCoefficientsForLeastSquareFittingScalar<11>,
			GetVectorCoefficientsForLeastSquareFittingScalar<12>,
			GetVectorCoefficientsForLeastSquareFittingScalar<13>,
			GetVectorCoefficientsForLeastSquareFittingScalar<14>,
			GetVectorCoefficientsForLeastSquareFittingScalar<15>,
			GetVectorCoefficientsForLeastSquareFittingScalar<16> };

	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);

	if (num_model_bases < ELEMENTSOF(funcs)) {
		funcs[num_model_bases](num_data, data, mask, model, out);
	} else {
		GetVectorCoefficientsForLeastSquareFitting(
				num_data, data, mask, num_model_bases, model, out);
	}
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::SolveSimultaneousEquationsByLU(
		size_t num_equations,
		double const lsq_matrix0[/*num_equations*num_equations*/],
		double const lsq_vector0[/*num_equations*/],
		double out[/*num_equations*/]) const {
	::SolveSimultaneousEquationsByLU(num_equations, lsq_matrix0, lsq_vector0, out);
}

}
