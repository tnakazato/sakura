#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/memory_manager.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

#define EIGEN_DENSEBASE_PLUGIN "eigen_binary_visitor_plugin.h"
#include <Eigen/Core>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

inline void GetBaselineModel(
		size_t num_data, size_t order, double *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	size_t num_model = order + 1;

	for (size_t i = 0; i < num_model; ++i) {
		for (size_t j = 0; j < num_data; ++j) {
			out[num_data*i+j] = pow(static_cast<double>(j), static_cast<double>(i));
		}
	}
}

inline void DoSubtractBaseline(
		size_t num_data, float const *in_data, bool const *in_mask, size_t num_model,
		double const *model_data, float clipping_threshold_sigma, size_t num_fitting_max,
		bool get_residual, float *out) {
	//std::cout << "DoSubtractBaselineEigen function is called" << std::endl;

	assert(LIBSAKURA_SYMBOL(IsAligned)(in_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	bool *clip_mask = reinterpret_cast<bool *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(bool)*num_data));
	for (size_t i = 0; i < num_data; ++i) {
		clip_mask[i] = true;
	}

	bool *composite_mask = reinterpret_cast<bool *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(bool)*num_data));
	float *best_fit_model = reinterpret_cast<float *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(float)*num_data));
	float *residual_data = reinterpret_cast<float *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(float)*num_data));
	//bool *new_clip_mask = reinterpret_cast<bool *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(bool)*num_data));

	for (size_t i = 0; i < num_fitting_max; ++i) {
		LIBSAKURA_SYMBOL(OperateLogicalAnd)(num_data, in_mask, clip_mask, composite_mask);
		LIBSAKURA_SYMBOL(GetBestFitModel)(num_data, in_data, composite_mask, num_model, model_data, best_fit_model);
		LIBSAKURA_SYMBOL(OperateFloatSubtraction)(num_data, in_data, best_fit_model, residual_data);

		//GetRms(); // calculate rms
		//LIBSAKURA_SYMBOL(SetTrueFloatInRangesInclusive)(); //new_clip_mask generated based on rms and residual_data

		/*
		bool mask_changed_after_clipping = false;
		for (size_t j = 0; j < num_data; ++j) {
			if (new_clip_mask[j] != clip_mask[j]) {
				mask_changed_after_clipping = true;
				break;
			}
		}

		if (mask_changed_after_clipping) {
			for (size_t j = 0; j < num_data; ++j) {
				clip_mask[j] = new_clip_mask[j];
			}
		} else {
			break;
		}
		*/
	}

	if (get_residual) {
		for (size_t i = 0; i < num_data; ++i) {
			out[i] = residual_data[i];
		}
	} else {
		for (size_t i = 0; i < num_data; ++i) {
			out[i] = best_fit_model[i];
		}
	}

	LIBSAKURA_PREFIX::Memory::Free(clip_mask);
	LIBSAKURA_PREFIX::Memory::Free(composite_mask);
	LIBSAKURA_PREFIX::Memory::Free(best_fit_model);
	LIBSAKURA_PREFIX::Memory::Free(residual_data);
	//LIBSAKURA_PREFIX::Memory::Free(new_clip_mask);
}

inline void SubtractBaselinePolynomial(
		size_t num_data, float const *in_data, bool const *in_mask,
		size_t order, float clipping_threshold_sigma,
		size_t num_fitting_max, bool get_residual,
		float *out) {
	size_t num_model = order + 1;
	double *model_data = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_data*num_model));

	GetBaselineModel(num_data, order, model_data);
	DoSubtractBaseline(num_data, in_data, in_mask, num_model, model_data,
			clipping_threshold_sigma, num_fitting_max, get_residual, out);

	LIBSAKURA_PREFIX::Memory::Free(model_data);
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Baseline, ARCH_SUFFIX)::SubtractBaselinePolynomial(
		size_t num_data, float const in_data[/*num_data*/],
		bool const in_mask[/*num_data*/], size_t order,
		float clipping_threshold_sigma, size_t num_fitting_max,
		bool get_residual, float out[/*num_data*/]) const {
	SubtractBaselinePolynomial(num_data, in_data, in_mask, order,
			clipping_threshold_sigma, num_fitting_max, get_residual, out);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::GetBaselineModel(
		size_t num_data, size_t order, double out[/*(order+1)*num_data*/]) const {
	GetBaselineModel(num_data, order, out);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::DoSubtractBaseline(
		size_t num_data, float const in_data[/*num_data*/],
		bool const in_mask[/*num_data*/], size_t num_model,
		double const model_data[/*num_model * num_data*/],
		float clipping_threshold_sigma, size_t num_fitting_max,
		bool get_residual, float out[/*num_data*/]) const {
	DoSubtractBaseline(num_data, in_data, in_mask, num_model, model_data,
		clipping_threshold_sigma, num_fitting_max, get_residual, out);
}

}
