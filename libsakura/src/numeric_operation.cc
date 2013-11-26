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
inline void GetMatrixCoefficientsForLeastSquareFittingUsingTemplate(
		size_t num_mask, bool const *mask_arg,
		size_t num_model_bases, double const *model_arg,
		double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));

	auto mask  = AssumeAligned(mask_arg);
	auto model = AssumeAligned(model_arg);
	auto out   = AssumeAligned(out_arg);

	for (size_t i = 0; i < NUM_MODEL_BASES * NUM_MODEL_BASES; ++i) {
		out[i] = 0;
	}
	for (size_t i = 0; i < num_mask; ++i) {
		if (mask[i]) {
			auto model_h = &model[i * NUM_MODEL_BASES];
			for (size_t j = 0; j < NUM_MODEL_BASES; ++j) {
				auto out_matrix_i = &out[j * NUM_MODEL_BASES];
				auto model_j = model_h[j];
				for (size_t k = 0; k < NUM_MODEL_BASES; ++k) {
					out_matrix_i[k] += model_j * model_h[k];
				}
			}
		}
	}
}

inline void GetMatrixCoefficientsForLeastSquareFitting(
		size_t num_mask, bool const *mask_arg,
		size_t num_model_bases, double const *model_arg,
		double *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));

	auto mask  = AssumeAligned(mask_arg);
	auto model = AssumeAligned(model_arg);
	auto out   = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_model_bases * num_model_bases; ++i) {
		out[i] = 0;
	}
	for (size_t i = 0; i < num_mask; ++i) {
		if (mask[i]) {
			auto model_i = &model[i * num_model_bases];
			for (size_t j = 0; j < num_model_bases; ++j) {
				auto out_matrix_j = &out[j * num_model_bases];
				auto model_j = model_i[j];
				for (size_t k = 0; k < num_model_bases; ++k) {
					out_matrix_j[k] += model_j * model_i[k];
				}
			}
		}
	}
}

template<size_t NUM_MODEL_BASES>
inline void GetVectorCoefficientsForLeastSquareFittingUsingTemplate(
		size_t num_data, float const *data_arg, bool const *mask_arg,
		size_t num_model_bases, double const *model_arg,
		double *out_arg) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));

	auto data  = AssumeAligned(data_arg);
	auto mask  = AssumeAligned(mask_arg);
	auto model = AssumeAligned(model_arg);
	auto out   = AssumeAligned(out_arg);

	for (size_t i = 0; i < NUM_MODEL_BASES; ++i) {
		out[i] = 0;
	}
	for (size_t i = 0; i < num_data; ++i){
		if (mask[i]) {
			auto model_i = &model[i * NUM_MODEL_BASES];
			auto data_i = data[i];
			for (size_t j = 0; j < NUM_MODEL_BASES; ++j) {
				out[j] += model_i[j] * data_i;
			}
		}
	}
}

inline void GetVectorCoefficientsForLeastSquareFitting(
		size_t num_data, float const *data_arg, bool const *mask_arg,
		size_t num_model_bases, double const *model_arg,
		double *out_arg) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));

	auto data  = AssumeAligned(data_arg);
	auto mask  = AssumeAligned(mask_arg);
	auto model = AssumeAligned(model_arg);
	auto out   = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_model_bases; ++i) {
		out[i] = 0;
	}
	for (size_t i = 0; i < num_data; ++i){
		if (mask[i]) {
			auto model_i = &model[i * num_model_bases];
			auto data_i = data[i];
			for (size_t j = 0; j < num_model_bases; ++j) {
				out[j] += model_i[j] * data_i;
			}
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
	assert(num_model_bases > 0);

	typedef void (*GetMatrixCoefficientsForLeastSquareFittingFunc)(
			size_t num_mask, bool const *mask,
			size_t num_model_bases, double const *model,
			double *out);

	static GetMatrixCoefficientsForLeastSquareFittingFunc const funcs[] = {
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<0>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<1>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<2>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<3>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<4>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<5>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<6>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<7>,
			::GetMatrixCoefficientsForLeastSquareFitting,		//non-template version faster for the case num_bases == 8.
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<9>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<10>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<11>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<12>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<13>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<14>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<15>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<16>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<17>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<18>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<19>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<20>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<21>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<22>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<23>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<24>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<25>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<26>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<27>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<28>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<29>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<30>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<31>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<32>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<33>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<34>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<35>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<36>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<37>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<38>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<39>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<40>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<41>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<42>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<43>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<44>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<45>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<46>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<47>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<48>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<49>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<50>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<51>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<52>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<53>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<54>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<55>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<56>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<57>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<58>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<59>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<60>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<61>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<62>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<63>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<64>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<65>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<66>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<67>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<68>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<69>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<70>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<71>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<72>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<73>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<74>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<75>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<76>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<77>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<78>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<79>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<80>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<81>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<82>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<83>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<84>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<85>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<86>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<87>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<88>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<89>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<90>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<91>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<92>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<93>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<94>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<95>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<96>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<97>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<98>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<99>,
			GetMatrixCoefficientsForLeastSquareFittingUsingTemplate<100>
	};

	if (num_model_bases < ELEMENTSOF(funcs)) {
		funcs[num_model_bases](
				num_mask, mask, num_model_bases, model, out);
	} else {
		::GetMatrixCoefficientsForLeastSquareFitting(
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
			size_t num_model_bases, double const *model, double *out);
	static GetVectorCoefficientsForLeastSquareFittingFunc const funcs[] = {
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<0>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<1>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<2>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<3>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<4>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<5>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<6>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<7>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<8>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<9>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<10>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<11>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<12>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<13>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<14>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<15>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<16>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<17>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<18>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<19>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<20>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<21>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<22>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<23>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<24>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<25>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<26>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<27>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<28>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<29>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<30>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<31>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<32>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<33>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<34>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<35>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<36>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<37>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<38>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<39>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<40>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<41>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<42>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<43>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<44>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<45>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<46>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<47>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<48>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<49>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<50>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<51>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<52>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<53>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<54>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<55>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<56>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<57>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<58>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<59>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<60>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<61>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<62>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<63>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<64>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<65>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<66>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<67>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<68>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<69>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<70>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<71>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<72>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<73>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<74>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<75>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<76>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<77>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<78>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<79>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<80>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<81>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<82>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<83>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<84>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<85>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<86>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<87>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<88>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<89>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<90>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<91>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<92>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<93>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<94>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<95>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<96>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<97>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<98>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<99>,
			GetVectorCoefficientsForLeastSquareFittingUsingTemplate<100>
	};

	if (num_model_bases < ELEMENTSOF(funcs)) {
		funcs[num_model_bases](
				num_data, data, mask, num_model_bases, model, out);
	} else {
		::GetVectorCoefficientsForLeastSquareFitting(
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
