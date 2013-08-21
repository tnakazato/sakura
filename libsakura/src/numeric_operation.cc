#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

#define FORCE_EIGEN 0

#if defined(__AVX__) && (! FORCE_EIGEN)
#include <immintrin.h>
#include <cstdint>

namespace {

void OperateFloatSubtractionSimd(size_t num_in, float const in1[],
		float const in2[], float out[]) {
	std::cout << "OperateFloatSubtractionSimd function is called. This function is not implemented yet." << std::endl;
}
void GetLeastSquareMatrixSimd(size_t num_in, float const in_data[],
		bool const in_mask[], size_t num_model, double const model[],
		double out[], double out_vector[]) {
	std::cout << "GetLeastSquareMatrixSimd function is called. This function is not implemented yet." << std::endl;
}
void SolveSimultaneousEquationsByLUSimd(size_t num_eqn,
		double const lsq_matrix0[], double const lsq_vector0[], double out[]) {
	std::cout << "SolveSimultaneousEquationsByLUSimd function is called. This function is not implemented yet." << std::endl;
}
void DoGetBestFitModelSimd(size_t num_chan, size_t num_eqn,
		double const model[], double const coeff[], float out[]) {
	std::cout << "DoGetBestFitModelSimd function is called. This function is not implemented yet." << std::endl;
}
void GetBestFitModelSimd(size_t num_in, float const in_data[], bool const in_mask[],
		size_t num_model, double const model[], float out[]) {
	std::cout << "GetBestFitModelSimd function is called. This function is not implemented yet." << std::endl;
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

	out = in1_ - in2_;
}

inline void GetLeastSquareMatrixEigen(size_t num_in, float const *in_data,
		bool const *in_mask, size_t num_model, double const *model,
		double *out, double *out_vector) {
	//std::cout << "GetLeastSquareMatrixEigen function is called" << std::endl;

	assert(LIBSAKURA_SYMBOL(IsAligned)(in_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_vector));
	Map<Array<float, Dynamic, 1>, Aligned> in_data_(const_cast<float *>(in_data),
			num_in);
	Map<Array<bool, Dynamic, 1>, Aligned> in_mask_(const_cast<bool *>(in_mask),
			num_in);
	Map<Array<double, Dynamic, 1>, Aligned> model_(const_cast<double *>(model),
			num_in);
	Map<Array<double, Dynamic, 1>, Aligned> out_(const_cast<double *>(out),
			num_in);
	Map<Array<double, Dynamic, 1>, Aligned> out_vector_(const_cast<double *>(out_vector),
			num_in);

	for (size_t i = 0; i < num_model; i++) {
		for (size_t j = 0; j < num_model; j++) {
			size_t idx = num_model * i + j;
			out_[idx] = 0.0;

			for (size_t k = 0; k < num_in ; k++){
				if (!in_mask_[k]) continue;

				size_t idx_i = num_in * i + k;
				size_t idx_j = num_in * j + k;
				out_[idx] += model_[idx_i] * model_[idx_j];
			}
		}

		out_vector[i] = 0.0;
		for (size_t k = 0; k < num_in; k++) {
			if (!in_mask_[k]) continue;

			size_t idx = num_in * i + k;
			out_vector_[i] += model_[idx] * in_data_[k];
		}
	}
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
	for (size_t i = 0; i < num_eqn; i++) {
		for (size_t j = 0; j < num_eqn; j++) {
			lsq_matrix(i, j) = lsq_matrix0_[num_eqn*i+j];
		}
		lsq_vector(i) = lsq_vector0_[i];
	}

	::Eigen::FullPivLU< ::Eigen::MatrixXd > lu(lsq_matrix);
	::Eigen::VectorXd lu_result = lu.solve(lsq_vector);
	for (size_t i = 0; i < num_eqn; i++) {
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

	for (size_t i = 0; i < num_chan; i++) {
		out_[i] = 0.0f;
		for (size_t j = 0; j < num_eqn; j++) {
			out_[i] += coeff_[j] * model_[num_chan * j + i];
		}
	}
}

inline void GetBestFitModelEigen(size_t num_in, float const *in_data,
		bool const *in_mask, size_t num_model, double const *model,
		float *out) {
	//std::cout << "GetBestFitModelEigen function is called" << std::endl;

	assert(LIBSAKURA_SYMBOL(IsAligned)(in_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map<Array<float, Dynamic, 1>, Aligned> in_data_(const_cast<float *>(in_data),
			num_in);
	Map<Array<bool, Dynamic, 1>, Aligned> in_mask_(const_cast<bool *>(in_mask),
			num_in);
	Map<Array<double, Dynamic, 1>, Aligned> model_(const_cast<double *>(model),
			num_model*num_in);
	Map<Array<float, Dynamic, 1>, Aligned> out_(const_cast<float *>(out),
			num_in);

	double *lsq_matrix0;
	double *lsq_vector0;
	double *coeff;
	GetLeastSquareMatrixEigen(num_in, in_data, in_mask,
			num_model, model, lsq_matrix0, lsq_vector0);
	SolveSimultaneousEquationsByLUEigen(num_model,
			lsq_matrix0, lsq_vector0, coeff);
	DoGetBestFitModelEigen(num_in, num_model, model, coeff, out);
}

} /* anonymous namespace */

#endif /* defined(__AVX__) */

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
#if defined( __AVX__) && (! FORCE_EIGEN)
	GetLeastSquareMatrixSimd(num_in, in_data, in_mask, num_model, model, out, out_vector);
#else
	GetLeastSquareMatrixEigen(num_in, in_data, in_mask, num_model, model, out, out_vector);
#endif
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
#if defined( __AVX__) && (! FORCE_EIGEN)
	GetBestFitModelSimd(num_in, in_data, in_mask, num_model, model, out);
#else
	GetBestFitModelEigen(num_in, in_data, in_mask, num_model, model, out);
#endif
}

}
