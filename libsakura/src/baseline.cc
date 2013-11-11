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

inline void DoGetBestFitBaseline(size_t num_data, size_t num_equations,
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

inline void GetBestFitBaseline(size_t num_data,
		float const *data, bool const *mask,
		size_t num_model_bases, double const *model,
		float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	size_t num_lsq_matrix0 = num_model_bases * num_model_bases;
	double *lsq_matrix0 = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_matrix0(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*lsq_matrix0) * num_lsq_matrix0, &lsq_matrix0));
	size_t num_lsq_vector0 = num_model_bases;
	double *lsq_vector0 = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_vector0(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*lsq_vector0) * num_lsq_vector0, &lsq_vector0));
	size_t num_coeff = num_model_bases;
	double *coeff = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_coeff(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*coeff) * num_coeff, &coeff));

	LIBSAKURA_SYMBOL(GetCoefficientsForLeastSquareFitting)(
			num_data, data, mask,
			num_model_bases, model,
			lsq_matrix0, lsq_vector0);

	LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLU)(
			num_model_bases,
			lsq_matrix0, lsq_vector0, coeff);

	DoGetBestFitBaseline(num_data, num_model_bases, model, coeff, out);
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

	bool *clip_mask = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_clip_mask(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*clip_mask) * num_data, &clip_mask));
	for (size_t i = 0; i < num_data; ++i) {
		clip_mask[i] = true;
	}

	bool *composite_mask = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_composite_mask(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*composite_mask) * num_data, &composite_mask));

	float *best_fit_model = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_best_fit_model(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*best_fit_model) * num_data, &best_fit_model));

	float *residual_data = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_residual_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*residual_data) * num_data, &residual_data));

	for (size_t i = 0; i < num_fitting_max; ++i) {
		LIBSAKURA_SYMBOL(OperateLogicalAnd)(num_data, mask, clip_mask, composite_mask);
		GetBestFitBaseline(num_data, data, composite_mask, num_model_bases, model, best_fit_model);
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
}

inline void SubtractBaselinePolynomial(
		size_t num_data, float const *data, bool const *mask,
		uint16_t order, float clipping_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual,
		float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	size_t num_model_bases = order + 1;
	size_t num_model = num_data * num_model_bases;
	double *model = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_model(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*model) * num_model, &model));

	GetBaselineModelPolynomial(num_data, order, model);
	DoSubtractBaseline(num_data, data, mask, num_model_bases, model,
			clipping_threshold_sigma, num_fitting_max, get_residual, out);
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Baseline, ARCH_SUFFIX)::GetBaselineModelPolynomial(
		size_t num_model, uint16_t order, double model[/*(order+1)*num_model*/]) const {
	::GetBaselineModelPolynomial(num_model, order, model);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::GetBestFitBaseline(
		size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		size_t num_model_bases, double const model[/*num_model_bases*num_data*/],
		float out[/*num_data*/]) const {
	::GetBestFitBaseline(num_data, data, mask, num_model_bases, model, out);
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
