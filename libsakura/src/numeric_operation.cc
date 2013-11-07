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
inline void OperateFloatSubtraction(size_t num_in, float const *in1,
		float const *in2, float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	STATIC_ASSERT(sizeof(in1) == sizeof(in2));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);

	for (size_t i = 0; i < num_in; ++i) {
		out[i] = in1[i] - in2[i];
	}
}

inline void GetCoefficientsForLeastSquareFitting(size_t num_data,
		float const *data, bool const *mask,
		size_t num_model_bases, double const *model,
		double *out_matrix, double *out_vector) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_matrix));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_vector));

	for (size_t i = 0; i < num_model_bases; ++i) {
		for (size_t j = 0; j < num_model_bases; ++j) {
			size_t idx = num_model_bases * i + j;
			out_matrix[idx] = 0.0;

			for (size_t k = 0; k < num_data ; ++k){
				if (!mask[k]) continue;

				size_t idx_i = num_data * i + k;
				size_t idx_j = num_data * j + k;
				out_matrix[idx] += model[idx_i] * model[idx_j];
			}
		}

		out_vector[i] = 0.0;
		for (size_t k = 0; k < num_data; ++k) {
			if (!mask[k]) continue;

			size_t idx = num_data * i + k;
			out_vector[i] += model[idx] * data[k];
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

inline void DoGetBestFitModel(size_t num_data, size_t num_equations,
		double const *model, double const *coeff, float *out) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	for (size_t i = 0; i < num_data; ++i) {
		out[i] = 0.0f;
		for (size_t j = 0; j < num_equations; ++j) {
			out[i] += coeff[j] * model[num_data * j + i];
		}
	}
}

inline void GetBestFitModel(size_t num_data,
		float const *data, bool const *mask,
		size_t num_model_bases, double const *model,
		float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
	size_t num_lsq_matrix0 = num_model_bases * num_model_bases;
	size_t num_arena = num_lsq_matrix0 + sakura_alignment - 1;
	std::unique_ptr<double[]> storage_for_lsq_matrix0(new double[num_arena]);
	double *lsq_matrix0 = sakura_AlignDouble(num_arena,
			storage_for_lsq_matrix0.get(), num_lsq_matrix0);
	//double *lsq_matrix0 = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_model_bases*num_model_bases));

	size_t num_lsq_vector0 = num_model_bases;
	num_arena = num_lsq_vector0 + sakura_alignment - 1;
	std::unique_ptr<double[]> storage_for_lsq_vector0(new double[num_arena]);
	double *lsq_vector0 = sakura_AlignDouble(num_arena,
			storage_for_lsq_vector0.get(), num_lsq_vector0);
	//double *lsq_vector0 = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_model_bases));

	size_t num_coeff = num_model_bases;
	num_arena = num_coeff + sakura_alignment - 1;
	std::unique_ptr<double[]> storage_for_coeff(new double[num_arena]);
	double *coeff = sakura_AlignDouble(num_arena,
			storage_for_coeff.get(), num_coeff);
	//double *coeff = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_model_bases));

	GetCoefficientsForLeastSquareFitting(num_data, data, mask,
			num_model_bases, model, lsq_matrix0, lsq_vector0);

	SolveSimultaneousEquationsByLU(num_model_bases,
			lsq_matrix0, lsq_vector0, coeff);

	DoGetBestFitModel(num_data, num_model_bases, model, coeff, out);

	//LIBSAKURA_PREFIX::Memory::Free(lsq_matrix0);
	//LIBSAKURA_PREFIX::Memory::Free(lsq_vector0);
	//LIBSAKURA_PREFIX::Memory::Free(coeff);
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

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::GetCoefficientsForLeastSquareFitting(
		size_t num_data, float const data[/*num_data*/], bool const mask[/*num_data*/],
		size_t num_model_bases, double const model[/*num_model_bases*num_data*/],
		double out_matrix[/*num_model_bases*num_model_bases*/],
		double out_vector[/*num_model_bases*/]) const {
	::GetCoefficientsForLeastSquareFitting(num_data, data, mask, num_model_bases, model, out_matrix, out_vector);
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::SolveSimultaneousEquationsByLU(
		size_t num_equations,
		double const lsq_matrix0[/*num_equations*num_equations*/],
		double const lsq_vector0[/*num_equations*/],
		double out[/*num_equations*/]) const {
	::SolveSimultaneousEquationsByLU(num_equations, lsq_matrix0, lsq_vector0, out);
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::GetBestFitModel(
		size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		size_t num_model_bases, double const model[/*num_model_bases*num_data*/],
		float out[/*num_data*/]) const {
	::GetBestFitModel(num_data, data, mask, num_model_bases, model, out);
}

}
