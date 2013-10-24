#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/localdef.h"
#include "libsakura/memory_manager.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/sakura.h"

#define EIGEN_DENSEBASE_PLUGIN "eigen_binary_visitor_plugin.h"
#include <Eigen/Core>
#include <Eigen/LU>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {
inline void OperateFloatSubtraction(size_t num_in, float const *in1,
		float const *in2, float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	for (size_t i = 0; i < num_in; ++i) {
		out[i] = in1[i] - in2[i];
	}
}

inline void SolveSimultaneousEquationsByLU(size_t num_eqn,
		double const *lsq_matrix0, double const *lsq_vector0, double *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix0));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector0));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map<Array<double, Dynamic, 1>, Aligned> lsq_matrix0_array(const_cast<double *>(lsq_matrix0),
			num_eqn*num_eqn);
	Map<Array<double, Dynamic, 1>, Aligned> lsq_vector0_array(const_cast<double *>(lsq_vector0),
			num_eqn);
	Map<Array<double, Dynamic, 1>, Aligned> out_array(const_cast<double *>(out),
			num_eqn);

	::Eigen::MatrixXd lsq_matrix(num_eqn, num_eqn);
	::Eigen::VectorXd lsq_vector(num_eqn);
	for (size_t i = 0; i < num_eqn; ++i) {
		for (size_t j = 0; j < num_eqn; ++j) {
			lsq_matrix(i, j) = lsq_matrix0_array[num_eqn*i+j];
		}
		lsq_vector(i) = lsq_vector0_array[i];
	}

	::Eigen::VectorXd lu_result = lsq_matrix.fullPivLu().solve(lsq_vector);
	for (size_t i = 0; i < num_eqn; ++i) {
		out_array[i] = lu_result(i);
	}
}

inline void GetCoefficientsForLeastSquareFitting(size_t num_in,
		float const *in_data, bool const *in_mask,
		size_t num_model, double const *model,
		double *out_matrix, double *out_vector) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(in_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_matrix));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_vector));

	for (size_t i = 0; i < num_model; ++i) {
		for (size_t j = 0; j < num_model; ++j) {
			size_t idx = num_model * i + j;
			out_matrix[idx] = 0.0;

			for (size_t k = 0; k < num_in ; ++k){
				if (!in_mask[k]) continue;

				size_t idx_i = num_in * i + k;
				size_t idx_j = num_in * j + k;
				out_matrix[idx] += model[idx_i] * model[idx_j];
			}
		}

		out_vector[i] = 0.0;
		for (size_t k = 0; k < num_in; ++k) {
			if (!in_mask[k]) continue;

			size_t idx = num_in * i + k;
			out_vector[i] += model[idx] * in_data[k];
		}
	}
}

inline void DoGetBestFitModel(size_t num_data, size_t num_eqn,
		double const *model, double const *coeff, float *out) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	for (size_t i = 0; i < num_data; ++i) {
		out[i] = 0.0f;
		for (size_t j = 0; j < num_eqn; ++j) {
			out[i] += coeff[j] * model[num_data * j + i];
		}
	}
}

inline void GetBestFitModel(size_t num_in,
		float const *in_data, bool const *in_mask,
		size_t num_model, double const *model,
		float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	double *lsq_matrix0 = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_model*num_model));
	//if (lsq_matrix0 == nullptr) {}
	double *lsq_vector0 = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_model));
	double *coeff = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_model));

	GetCoefficientsForLeastSquareFitting(num_in, in_data, in_mask,
			num_model, model, lsq_matrix0, lsq_vector0);

	SolveSimultaneousEquationsByLU(num_model,
			lsq_matrix0, lsq_vector0, coeff);

	DoGetBestFitModel(num_in, num_model, model, coeff, out);

	LIBSAKURA_PREFIX::Memory::Free(lsq_matrix0);
	LIBSAKURA_PREFIX::Memory::Free(lsq_vector0);
	LIBSAKURA_PREFIX::Memory::Free(coeff);
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::OperateFloatSubtraction(
		size_t num_in,
		float const in1[/*num_in*/],
		float const in2[/*num_in*/],
		float out[/*num_in*/]) const {
	::OperateFloatSubtraction(num_in, in1, in2, out);
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::GetCoefficientsForLeastSquareFitting(size_t num_in,
		float const in_data[/*num_in*/], bool const in_mask[/*num_in*/],
		size_t num_model, double const model[/*num_model * num_in*/],
		double out_matrix[/*num_model * num_model*/], double out_vector[/*num_model*/]) const {
	::GetCoefficientsForLeastSquareFitting(num_in, in_data, in_mask, num_model, model, out_matrix, out_vector);
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::SolveSimultaneousEquationsByLU(size_t num_eqn,
		double const lsq_matrix0[/*num_eqn * num_eqn*/],
		double const lsq_vector0[/*num_eqn*/], double out[/*num_eqn*/]) const {
	::SolveSimultaneousEquationsByLU(num_eqn, lsq_matrix0, lsq_vector0, out);
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::GetBestFitModel(size_t num_in,
		float const in_data[/*num_in*/], bool const in_mask[/*num_in*/],
		size_t num_model, double const model[/*num_model * num_in*/],
		float out[/*num_in*/]) const {
	::GetBestFitModel(num_in, in_data, in_mask, num_model, model, out);
}

}
