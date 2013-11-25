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
#include "baseline.h"

namespace {

inline void SetBasisDataPolynomial(
		LIBSAKURA_SYMBOL(BaselineContext) *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));

/*
	for (size_t i = 0; i < (context)->num_basis_data; ++i) {
		double val = 1.0;
		context->basis_data[i] = val;
		for (size_t j = 1; j < context->num_bases; ++j) {
			val *= static_cast<double>(i);
			context->basis_data[context->num_basis_data * j + i] = val;
		}
	}
*/
	size_t idx = 0;
	for (size_t i = 0; i < (context)->num_basis_data; ++i) {
		double val = 1.0;
		context->basis_data[idx] = val;
		idx ++;
		for (size_t j = 1; j < context->num_bases; ++j) {
			val *= static_cast<double>(i);
			context->basis_data[idx] = val;
			idx ++;
		}
	}
}

inline void GetChebyshevPolynomial(
		size_t const n, double const x, double *out) {
	double out_ = 0.0;
	if (n < 0) {
		throw "the order must be zero or positive.";
	} else if (n == 0) {
		out_ = 1.0;
	} else if (n == 1) {
		out_ = x;
	} else if ((x < -1.0)||(x > 1.0)) {
		throw "out of definition range (-1 <= x <= 1).";
	} else if (x == 1.0) {
		out_ = 1.0;
	} else if (x == 0.0) {
		out_ = (n%2 == 0) ? ((n%4 == 0) ? 1.0 : -1.0) : 0.0;
	} else if (x == -1.0) {
		out_ = (n%2 == 0) ? 1.0 : -1.0;
	} else {
		double res = 0.0;

		/*
		double prev1 = 0.0;
		GetChebyshevPolynomial(n-1, x, &prev1);
		double prev2 = 0.0;
		GetChebyshevPolynomial(n-2, x, &prev2);
		res = 2.0 * x * prev1 - prev2;
		*/

		for (size_t m = 0; m <= n/2; ++m) {
			double c = 1.0;
			if (m > 0) {
				for (size_t i = 1; i <= m; ++i) {
					c *= (double)(n-2*m+i)/(double)i;
				}
			}
			res += ((m%2 == 0) ? 1.0 : -1.0)*(double)n/(double)(n-m)*pow(2.0*x, (double)(n-2*m))/2.0*c;
		}

		if (fabs(res) > 1.0) {
			throw "bad result: Chebyshev polynomial values must be [-1:1].";
		}

		out_ = res;
	}

	*out = out_;
}

inline void SetBasisDataChebyshev(
		LIBSAKURA_SYMBOL(BaselineContext) *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));

	size_t idx = 0;
	double xrange = (double)(context->num_basis_data - 1);
	for (size_t i = 0; i < context->num_basis_data; ++i) {
		double x = 2.0 * (double)i/xrange - 1.0;
		for (size_t j = 0; j < context->num_bases; ++j) {
			double val = 0.0;
			//GetChebyshevPolynomial(j, x, &val);
			if (j == 0) {
				val = 1.0;
			} else if (j == 1) {
				val = x;
			} else {
				val = 2.0 * x * context->basis_data[idx-1] - context->basis_data[idx-2];
			}
			context->basis_data[idx] = val;
			idx ++;
		}
	}
}

inline void SetBasisDataCubicSpline(
		LIBSAKURA_SYMBOL(BaselineContext) *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	//
}

inline void SetBasisDataSinusoid(
		LIBSAKURA_SYMBOL(BaselineContext) *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	//
}

inline void SetBasisData(
		LIBSAKURA_SYMBOL(BaselineContext) *context) {
	switch (context->baseline_type) {
	case LIBSAKURA_SYMBOL(BaselineType_kPolynomial):
		SetBasisDataPolynomial(context);
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kChebyshev):
		SetBasisDataChebyshev(context);
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kCubicSpline):
		SetBasisDataCubicSpline(context);
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kSinusoid):
		SetBasisDataSinusoid(context);
		break;
	default:
		assert(false);
		break;
	}
}

inline void CreateBaselineContext(
		LIBSAKURA_SYMBOL(BaselineType) const baseline_type,
		uint16_t const order,
		size_t const num_basis_data,
		LIBSAKURA_SYMBOL(BaselineContext) **context) {
	size_t num_bases = 0;
	switch (baseline_type) {
	case LIBSAKURA_SYMBOL(BaselineType_kPolynomial):
		num_bases = order + 1;
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kChebyshev):
		num_bases = order + 1;
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kCubicSpline):
		assert(order > 0);
		num_bases = order + 3;
		break;
	case LIBSAKURA_SYMBOL(BaselineType_kSinusoid):
		num_bases = 2 * order + 1;
		break;
	default:
		assert(false);
		break;
	}

	std::unique_ptr<LIBSAKURA_SYMBOL(BaselineContext),
			LIBSAKURA_PREFIX::Memory> work_context(
				static_cast<LIBSAKURA_SYMBOL(BaselineContext) *>
					(LIBSAKURA_PREFIX::Memory::Allocate(
						sizeof(LIBSAKURA_SYMBOL(BaselineContext)))),
				LIBSAKURA_PREFIX::Memory());
	if (work_context == nullptr) {
		throw std::bad_alloc();
	}

	work_context->baseline_type  = baseline_type;
	work_context->num_bases      = num_bases;
	work_context->num_basis_data = num_basis_data;

	size_t num_total_basis_data = num_bases * num_basis_data;

	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> work_basis_data_storage(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*work_context->basis_data) * num_total_basis_data,
					&work_context->basis_data),
					LIBSAKURA_PREFIX::Memory());
	assert(LIBSAKURA_SYMBOL(IsAligned)(work_context->basis_data));

	SetBasisData(work_context.get());

	work_context->basis_data_storage = work_basis_data_storage.release();
	*context = (LIBSAKURA_SYMBOL(BaselineContext) *)work_context.release();
}

inline void DestroyBaselineContext(
		LIBSAKURA_SYMBOL(BaselineContext) *context) {
	if (context->basis_data_storage != nullptr) {
		LIBSAKURA_PREFIX::Memory::Free(context->basis_data_storage);
	}
	LIBSAKURA_PREFIX::Memory::Free(context);
}

inline void OperateLogicalAnd(
		size_t num_in, bool const *in1, bool const *in2, bool *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	auto src1 = reinterpret_cast<uint8_t const *>(in1);
	auto src2 = reinterpret_cast<uint8_t const *>(in2);
	for (size_t i = 0; i < num_in; ++i) {
		out[i] = static_cast<bool>(src1[i] & src2[i]);
	}
}

inline void OperateFloatSubtraction(
		size_t num_in, float const *in1, float const *in2, float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	for (size_t i = 0; i < num_in; ++i) {
		out[i] = in1[i] - in2[i];
	}
}

inline void DoGetBestFitBaseline(
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		double const *coeff,
		float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	size_t num_data = context->num_basis_data;
	size_t num_bases = context->num_bases;
	for (size_t i = 0; i < num_data; ++i) {
		out[i] = 0.0f;
		for (size_t j = 0; j < num_bases; ++j) {
			out[i] += coeff[j] * context->basis_data[num_bases * i + j];
		}
	}
}

inline void GetBestFitBaseline(
		size_t num_data, float const *data, bool const *mask,
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	assert(num_data == context->num_basis_data);

	size_t num_lsq_matrix0 = context->num_bases * context->num_bases;
	double *lsq_matrix0 = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_matrix0(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*lsq_matrix0) * num_lsq_matrix0, &lsq_matrix0));
	size_t num_lsq_vector0 = context->num_bases;
	double *lsq_vector0 = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_vector0(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*lsq_vector0) * num_lsq_vector0, &lsq_vector0));
	size_t num_coeff = context->num_bases;
	double *coeff = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_coeff(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*coeff) * num_coeff, &coeff));

	LIBSAKURA_SYMBOL(GetMatrixCoefficientsForLeastSquareFitting)(
			num_data, mask, context->num_bases, context->basis_data, lsq_matrix0);
	LIBSAKURA_SYMBOL(GetVectorCoefficientsForLeastSquareFitting)(
			num_data, data, mask, context->num_bases, context->basis_data, lsq_vector0);

	LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLU)(
			context->num_bases,
			lsq_matrix0, lsq_vector0, coeff);

	DoGetBestFitBaseline(context, coeff, out);
}

inline void SubtractBaseline(
		size_t num_data, float const *data, bool const *mask,
		LIBSAKURA_SYMBOL(BaselineContext) const *baseline_context,
		float clipping_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual, bool *final_mask, float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(baseline_context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask));
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
		if (i > 0) {
			OperateLogicalAnd(num_data, mask, clip_mask, composite_mask);
		}
		GetBestFitBaseline(num_data, data, composite_mask, baseline_context, best_fit_model);
		OperateFloatSubtraction(num_data, data, best_fit_model, residual_data);

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

	for (size_t i = 0; i < num_data; ++i) {
		final_mask[i] = composite_mask[i];
		out[i] = get_residual ? residual_data[i] : best_fit_model[i];
	}
}

inline void SubtractBaselinePolynomial(
		size_t num_data, float const *data, bool const *mask,
		uint16_t order, float clipping_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual,
		bool *final_mask, float *out) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kPolynomial),
			order, num_data, &context);
	SubtractBaseline(num_data, data, mask, context,
			clipping_threshold_sigma, num_fitting_max,
			get_residual, final_mask, out);
	DestroyBaselineContext(context);
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::CreateBaselineContext(
		LIBSAKURA_SYMBOL(BaselineType) const baseline_type,
		uint16_t const order, size_t const num_data,
		LIBSAKURA_SYMBOL(BaselineContext) **context) const {
	::CreateBaselineContext(baseline_type, order, num_data, context);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::DestroyBaselineContext(
		LIBSAKURA_SYMBOL(BaselineContext) *context) const {
	::DestroyBaselineContext(context);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::GetBestFitBaseline(
		size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		float out[/*num_data*/]) const {
	::GetBestFitBaseline(num_data, data, mask, context, out);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::SubtractBaseline(
		size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineContext) const *baseline_context,
		float clipping_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual,
		bool final_mask[/*num_data*/], float out[/*num_data*/]) const {
	::SubtractBaseline(num_data, data, mask, baseline_context,
			clipping_threshold_sigma, num_fitting_max,
			get_residual, final_mask, out);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::SubtractBaselinePolynomial(
		size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		uint16_t order, float clipping_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual,
		bool final_mask[/*num_data*/], float out[/*num_data*/]) const {
	::SubtractBaselinePolynomial(num_data, data, mask, order,
			clipping_threshold_sigma, num_fitting_max,
			get_residual, final_mask, out);
}

}
