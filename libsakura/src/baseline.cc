#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/memory_manager.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

#define FORCE_EIGEN 0

#if defined(__AVX__) && (! FORCE_EIGEN)
#include <immintrin.h>
#include <cstdint>

namespace {
void DoSubtractBaselineSimd(size_t num_chan, float const in_data[], bool const in_mask[], size_t num_model,
		double model_data[], float clipping_threshold_sigma, int clipping_max_iteration, bool get_residual, float out[]) {
	std::cout << "DoSubtractBaselineSimd function is called. This function is not implemented yet." << std::endl;
}
} /* anonymous namespace */

#else /* defined(__AVX__) */

#define EIGEN_DENSEBASE_PLUGIN "eigen_binary_visitor_plugin.h"
#include <Eigen/Core>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

inline void GetBaselineModel(
		size_t num_chan, int order, double *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	size_t num_model = order + 1;

	for (size_t i = 0; i < num_model; i++) {
		for (size_t j = 0; j < num_chan; j++) {
			out[num_chan*i+j] = pow(static_cast<double>(j), static_cast<double>(i));
		}
	}
}

inline void DoSubtractBaselineEigen(
		size_t num_chan, float const *in_data, bool const *in_mask, size_t num_model,
		double *model_data, float clipping_threshold_sigma, int clipping_max_iteration,
		bool get_residual, float *out) {
	//std::cout << "DoSubtractBaselineEigen function is called" << std::endl;

	assert(LIBSAKURA_SYMBOL(IsAligned)(in_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	Map<Array<float, Dynamic, 1>, Aligned> in_data_(const_cast<float *>(in_data),
			num_chan);
	Map<Array<bool, Dynamic, 1>, Aligned> in_mask_(const_cast<bool *>(in_mask),
				num_chan);
	Map<Array<double, Dynamic, 1>, Aligned> model_data_(const_cast<double *>(model_data),
			num_model*num_chan);
	Map<Array<float, Dynamic, 1>, Aligned> out_(const_cast<float *>(out),
			num_chan);

	bool clip_mask[num_chan];
	for (size_t i = 0; i < num_chan; i++) {
		clip_mask[i] = true;
	}
	Map<Array<bool, Dynamic, 1>, Aligned> clip_mask_(const_cast<bool *>(clip_mask),
				num_chan);

	bool composite_mask[num_chan];
	float best_fit_model[num_chan];
	float residual_data[num_chan];
	//bool new_clip_mask[num_chan];

	int num_clipping = 0;
	bool end_condition = false;

	while (!end_condition) {
		LIBSAKURA_SYMBOL(OperateLogicalAnd)(num_chan, in_mask, clip_mask, composite_mask);
		LIBSAKURA_SYMBOL(GetBestFitModel)(num_chan, in_data, composite_mask, num_model, model_data, best_fit_model);
		LIBSAKURA_SYMBOL(OperateFloatSubtraction)(num_chan, in_data, best_fit_model, residual_data);

		//GetRms(); // calculate rms
		//LIBSAKURA_SYMBOL(SetTrueFloatInRangesInclusive)(); //new_clip_mask generated based on rms and residual_data

		/*
		bool no_change_after_clipping = true;
		for (size_t i = 0; i < num_chan; i++) {
			if (new_clip_mask[i] != clip_mask[i]) {
				no_change_after_clipping = false;
				break;
			}
		}
		if (no_change_after_clipping) end_condition = true;
		*/
		if (clipping_max_iteration <= num_clipping) end_condition = true;
		if (end_condition) break;

		/*
		for (size_t i = 0; i < num_chan; i++) {
			clip_mask[i] = new_clip_mask[i];
		}
		*/
		num_clipping ++;
	}

	if (get_residual) {
		for (size_t i = 0; i < num_chan; i++) {
			out[i] = residual_data[i];
		}
	} else {
		for (size_t i = 0; i < num_chan; i++) {
			out[i] = best_fit_model[i];
		}
	}
}

inline void SubtractBaselinePolynomial(
		size_t num_chan, float const *in_data, bool const *in_mask,
		int order, float clipping_threshold_sigma,
		int clipping_max_iteration, bool get_residual,
		float *out) {
	size_t num_model = order + 1;
	double *model_data = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_chan*num_model));

	GetBaselineModel(num_chan, order, model_data);
	DoSubtractBaselineEigen(num_chan, in_data, in_mask, num_model, model_data,
			clipping_threshold_sigma, clipping_max_iteration, get_residual, out);

	LIBSAKURA_PREFIX::Memory::Free(model_data);
}

} /* anonymous namespace */

#endif /* defined(__AVX__) */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Baseline, ARCH_SUFFIX)::SubtractBaselinePolynomial(
		size_t num_chan, float const in_data[/*num_chan*/],
		bool const in_mask[/*num_chan*/], int order,
		float clipping_threshold_sigma, int clipping_max_iteration,
		bool get_residual, float out[/*num_chan*/]) const {
	SubtractBaselinePolynomial(num_chan, in_data, in_mask, order,
			clipping_threshold_sigma, clipping_max_iteration, get_residual, out);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::GetBaselineModel(
		size_t num_chan, int order, double out[/*(order+1)*num_chan*/]) const {
	GetBaselineModel(num_chan, order, out);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::DoSubtractBaseline(
		size_t num_chan, float const in_data[/*num_chan*/],
		bool const in_mask[/*num_chan*/], size_t num_model,
		double model_data[/*num_model * num_chan*/],
		float clipping_threshold_sigma, int clipping_max_iteration,
		bool get_residual, float out[/*num_chan*/]) const {
#if defined( __AVX__) && (! FORCE_EIGEN)
	DoSubtractBaselineSimd(num_chan, in_data, in_mask, num_model, model_data,
		clipping_threshold_sigma, clipping_max_iteration, get_residual, out);
#else
	DoSubtractBaselineEigen(num_chan, in_data, in_mask, num_model, model_data,
		clipping_threshold_sigma, clipping_max_iteration, get_residual, out);
#endif
}

}
