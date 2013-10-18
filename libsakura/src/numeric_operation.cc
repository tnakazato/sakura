#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/localdef.h"
#include "libsakura/memory_manager.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/sakura.h"

#define FORCE_EIGEN 0

#if defined(__AVX__) && (! FORCE_EIGEN)
#include <immintrin.h>
#include <cstdint>

namespace {

void OperateFloatSubtractionSimd(size_t num_in, float const in1[],
		float const in2[], float out[]) {
	std::cout << "OperateFloatSubtractionSimd function is called. This function is not implemented yet." << std::endl;
}
void SolveSimultaneousEquationsByLUSimd(size_t num_eqn,
		double const lsq_matrix0[], double const lsq_vector0[], double out[]) {
	std::cout << "SolveSimultaneousEquationsByLUSimd function is called. This function is not implemented yet." << std::endl;
}
void DoGetBestFitModelSimd(size_t num_chan, size_t num_eqn,
		double const model[], double const coeff[], float out[]) {
	std::cout << "DoGetBestFitModelSimd function is called. This function is not implemented yet." << std::endl;
}

} /* anonymous namespace */

#else /* defined(__AVX__) */

#define EIGEN_DENSEBASE_PLUGIN "eigen_binary_visitor_plugin.h"
#include <Eigen/Core>
#include <Eigen/LU>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

inline void OperateFloatSubtractionEigen(size_t num_in, float const *in1,
		float const *in2, float *out) {
	//std::cout << "OperateFloatSubtractionEigen function is called" << std::endl;

	assert(LIBSAKURA_SYMBOL(IsAligned)(in1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map<Array<float, Dynamic, 1>, Aligned> in1_(const_cast<float *>(in1),
			num_in);
	Map<Array<float, Dynamic, 1>, Aligned> in2_(const_cast<float *>(in2),
			num_in);
	Map<Array<float, Dynamic, 1>, Aligned> out_(const_cast<float *>(out),
			num_in);

	out_ = in1_ - in2_;
}

inline void SolveSimultaneousEquationsByLUEigen(size_t num_eqn,
		double const *lsq_matrix0, double const *lsq_vector0, double *out) {
	//std::cout << "SolveSimultaneousEquationsByLUEigen function is called" << std::endl;

	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix0));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector0));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map<Array<double, Dynamic, 1>, Aligned> lsq_matrix0_(const_cast<double *>(lsq_matrix0),
			num_eqn*num_eqn);
	Map<Array<double, Dynamic, 1>, Aligned> lsq_vector0_(const_cast<double *>(lsq_vector0),
			num_eqn);
	Map<Array<double, Dynamic, 1>, Aligned> out_(const_cast<double *>(out),
			num_eqn);

	::Eigen::MatrixXd lsq_matrix(num_eqn, num_eqn);
	::Eigen::VectorXd lsq_vector(num_eqn);
	for (size_t i = 0; i < num_eqn; ++i) {
		for (size_t j = 0; j < num_eqn; ++j) {
			lsq_matrix(i, j) = lsq_matrix0_[num_eqn*i+j];
		}
		lsq_vector(i) = lsq_vector0_[i];
	}

	::Eigen::VectorXd lu_result = lsq_matrix.fullPivLu().solve(lsq_vector);
	for (size_t i = 0; i < num_eqn; ++i) {
		out_[i] = lu_result(i);
	}
}

inline void DoGetBestFitModelEigen(size_t num_chan, size_t num_eqn,
		double const *model, double const *coeff, float *out) {
	//std::cout << "DoGetBestFitModelEigen function is called" << std::endl;

	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map<Array<double, Dynamic, 1>, Aligned> model_(const_cast<double *>(model),
			num_eqn*num_chan);
	Map<Array<double, Dynamic, 1>, Aligned> coeff_(const_cast<double *>(coeff),
			num_eqn);
	Map<Array<float, Dynamic, 1>, Aligned> out_(const_cast<float *>(out),
			num_chan);

	for (size_t i = 0; i < num_chan; ++i) {
		out_[i] = 0.0f;
		for (size_t j = 0; j < num_eqn; ++j) {
			out_[i] += coeff_[j] * model_[num_chan * j + i];
		}
	}
}

} /* anonymous namespace */

#endif /* defined(__AVX__) */

namespace {
inline void GetLeastSquareMatrix(size_t num_in, float const *in_data,
		bool const *in_mask, size_t num_model, double const *model,
		double *out, double *out_vector) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(in_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_vector));

	for (size_t i = 0; i < num_model; ++i) {
		for (size_t j = 0; j < num_model; ++j) {
			size_t idx = num_model * i + j;
			out[idx] = 0.0;

			for (size_t k = 0; k < num_in ; ++k){
				if (!in_mask[k]) continue;

				size_t idx_i = num_in * i + k;
				size_t idx_j = num_in * j + k;
				out[idx] += model[idx_i] * model[idx_j];
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

inline void GetBestFitModel(size_t num_in, float const *in_data,
		bool const *in_mask, size_t num_model, double const *model,
		float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	double *lsq_matrix0 = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_model*num_model));
	//if (lsq_matrix0 == nullptr) {}
	double *lsq_vector0 = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_model));
	double *coeff = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_model));

	LIBSAKURA_SYMBOL(GetLeastSquareMatrix)(num_in, in_data, in_mask,
			num_model, model, lsq_matrix0, lsq_vector0);

	LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLU)(num_model,
			lsq_matrix0, lsq_vector0, coeff);

	LIBSAKURA_SYMBOL(DoGetBestFitModel)(num_in, num_model, model, coeff, out);

	LIBSAKURA_PREFIX::Memory::Free(lsq_matrix0);
	LIBSAKURA_PREFIX::Memory::Free(lsq_vector0);
	LIBSAKURA_PREFIX::Memory::Free(coeff);
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::OperateFloatSubtraction(size_t num_in,
		float const in1[/*num_in*/], float const in2[/*num_in*/],
		float out[/*num_in*/]) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	OperateFloatSubtractionSimd(num_in, in1, in2, out);
#else
	OperateFloatSubtractionEigen(num_in, in1, in2, out);
#endif
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::GetLeastSquareMatrix(size_t num_in,
		float const in_data[/*num_in*/], bool const in_mask[/*num_in*/],
		size_t num_model, double const model[/*num_model * num_in*/],
		double out[/*num_model * num_model*/], double out_vector[/*num_model*/]) const {
	::GetLeastSquareMatrix(num_in, in_data, in_mask, num_model, model, out, out_vector);
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::SolveSimultaneousEquationsByLU(size_t num_eqn,
		double const lsq_matrix0[/*num_eqn * num_eqn*/],
		double const lsq_vector0[/*num_eqn*/], double out[/*num_eqn*/]) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	SolveSimultaneousEquationsByLUSimd(num_eqn, lsq_matrix0, lsq_vector0, out);
#else
	SolveSimultaneousEquationsByLUEigen(num_eqn, lsq_matrix0, lsq_vector0, out);
#endif
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::DoGetBestFitModel(size_t num_chan,
		size_t num_eqn, double const model[/*num_eqn * num_chan*/],
		double const coeff[/*num_eqn*/], float out[/*num_chan*/]) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	DoGetBestFitModelSimd(num_chan, num_eqn, model, coeff, out);
#else
	DoGetBestFitModelEigen(num_chan, num_eqn, model, coeff, out);
#endif
}

void ADDSUFFIX(NumericOperation, ARCH_SUFFIX)::GetBestFitModel(size_t num_in,
		float const in_data[/*num_in*/], bool const in_mask[/*num_in*/],
		size_t num_model, double const model[/*num_model * num_in*/],
		float out[/*num_in*/]) const {
	::GetBestFitModel(num_in, in_data, in_mask, num_model, model, out);
}

}
