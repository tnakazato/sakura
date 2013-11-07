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

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

inline void GetBaselineModelPolynomial(
		size_t num_each_basis, uint16_t order, double *model) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	size_t num_model_bases = order + 1;

	for (size_t j = 0; j < num_each_basis; ++j) {
		double val = 1.0;
		model[j] = val;
		for (size_t i = 1; i < num_model_bases; ++i) {
			val *= static_cast<double>(j);
			model[num_each_basis*i+j] = val;
		}
	}
}

inline void DoSubtractBaseline(
		size_t num_data, float const *data, bool const *mask,
		size_t num_model_bases, double const *model,
		float clipping_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual, float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
	size_t num_arena = num_data + sakura_alignment - 1;

	std::unique_ptr<bool[]> storage_for_clip_mask(new bool[num_data]);
	bool *clip_mask = static_cast<bool *>(LIBSAKURA_SYMBOL(AlignAny)
			(num_arena*sizeof(bool), storage_for_clip_mask.get(),
					num_data*sizeof(bool)));
	//bool *clip_mask = reinterpret_cast<bool *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(bool)*num_data));
	for (size_t i = 0; i < num_data; ++i) {
		clip_mask[i] = true;
	}

	std::unique_ptr<bool[]> storage_for_composite_mask(new bool[num_data]);
	bool *composite_mask = static_cast<bool *>(LIBSAKURA_SYMBOL(AlignAny)
			(num_arena*sizeof(bool), storage_for_composite_mask.get(),
					num_data*sizeof(bool)));
	//bool *composite_mask = reinterpret_cast<bool *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(bool)*num_data));

	std::unique_ptr<float[]> storage_for_best_fit_model(new float[num_data]);
	float *best_fit_model = sakura_AlignFloat(num_arena,
			storage_for_best_fit_model.get(), num_data);
	//float *best_fit_model = reinterpret_cast<float *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(float)*num_data));

	std::unique_ptr<float[]> storage_for_residual_data(new float[num_data]);
	float *residual_data = sakura_AlignFloat(num_arena,
			storage_for_residual_data.get(), num_data);
	//float *residual_data = reinterpret_cast<float *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(float)*num_data));

	/*
	std::unique_ptr<bool[]> storage_for_new_clip_mask(new bool[num_data]);
	bool *new_clip_mask = static_cast<bool *>(LIBSAKURA_SYMBOL(AlignAny)
			(num_arena*sizeof(bool), storage_for_new_clip_mask.get(),
					num_data*sizeof(bool)));
	*/

	for (size_t i = 0; i < num_fitting_max; ++i) {
		LIBSAKURA_SYMBOL(OperateLogicalAnd)(num_data, mask, clip_mask, composite_mask);
		LIBSAKURA_SYMBOL(GetBestFitModel)(num_data, data, composite_mask, num_model_bases, model, best_fit_model);
		LIBSAKURA_SYMBOL(OperateFloatSubtraction)(num_data, data, best_fit_model, residual_data);

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

	//LIBSAKURA_PREFIX::Memory::Free(clip_mask);
	//LIBSAKURA_PREFIX::Memory::Free(composite_mask);
	//LIBSAKURA_PREFIX::Memory::Free(best_fit_model);
	//LIBSAKURA_PREFIX::Memory::Free(residual_data);
}

inline void SubtractBaselinePolynomial(
		size_t num_data, float const *data, bool const *mask,
		uint16_t order, float clipping_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual,
		float *out) {
	size_t num_model_bases = order + 1;

	size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
	size_t num_model = num_data * num_model_bases;
	size_t num_arena = num_model + sakura_alignment - 1;
	std::unique_ptr<double[]> storage_for_model(new double[num_arena]);
	double *model = sakura_AlignDouble(num_arena,
			storage_for_model.get(), num_model);
	//double *model = reinterpret_cast<double *>(LIBSAKURA_PREFIX::Memory::Allocate(sizeof(double)*num_data*num_model_bases));

	GetBaselineModelPolynomial(num_data, order, model);
	DoSubtractBaseline(num_data, data, mask, num_model_bases, model,
			clipping_threshold_sigma, num_fitting_max, get_residual, out);

	//LIBSAKURA_PREFIX::Memory::Free(model);
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Baseline, ARCH_SUFFIX)::GetBaselineModelPolynomial(
		size_t num_model, uint16_t order, double model[/*(order+1)*num_model*/]) const {
	::GetBaselineModelPolynomial(num_model, order, model);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::DoSubtractBaseline(
		size_t num_data, float const data[/*num_data*/], bool const mask[/*num_data*/],
		size_t num_model_bases, double const model[/*num_model_bases * num_data*/],
		float clipping_threshold_sigma, uint16_t num_fitting_max, bool get_residual,
		float out[/*num_data*/]) const {
	::DoSubtractBaseline(num_data, data, mask, num_model_bases, model,
		clipping_threshold_sigma, num_fitting_max, get_residual, out);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::SubtractBaselinePolynomial(
		size_t num_data, float const data[/*num_data*/], bool const mask[/*num_data*/],
		uint16_t order, float clipping_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual, float out[/*num_data*/]) const {
	::SubtractBaselinePolynomial(num_data, data, mask, order,
			clipping_threshold_sigma, num_fitting_max, get_residual, out);
}

}
