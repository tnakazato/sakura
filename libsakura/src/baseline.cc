/*
 * baseline.cc
 *
 *  Created on: 2013/11/11
 *      Author: wataru
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <sstream>

#if defined(__AVX__) && !defined(ARCH_SCALAR)
#	include <immintrin.h>
#endif

#include <libsakura/localdef.h>
#include <libsakura/memory_manager.h>
#include <libsakura/optimized_implementation_factory_impl.h>
#include <libsakura/sakura.h>
#include "baseline.h"

namespace {

inline void SetBasisDataPolynomial(LIBSAKURA_SYMBOL(BaselineContext) *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	auto data = AssumeAligned(context->basis_data);

	size_t idx = 0;
	for (size_t i = 0; i < context->num_basis_data; ++i) {
		double val = 1.0;
		data[idx] = val;
		idx++;
		for (size_t j = 1; j < context->num_bases; ++j) {
			val *= static_cast<double>(i);
			data[idx] = val;
			idx++;
		}
	}
}

inline void SetBasisDataChebyshev(LIBSAKURA_SYMBOL(BaselineContext) *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	auto data = AssumeAligned(context->basis_data);

	size_t idx = 0;
	double xrange = (double) (context->num_basis_data - 1);
	for (size_t i = 0; i < context->num_basis_data; ++i) {
		double x = 2.0 * (double) i / xrange - 1.0;
		for (size_t j = 0; j < context->num_bases; ++j) {
			double val = 0.0;
			if (j == 0) {
				val = 1.0;
			} else if (j == 1) {
				val = x;
			} else {
				val = 2.0 * x * data[idx - 1] - data[idx - 2];
			}
			data[idx] = val;
			idx++;
		}
	}
}

inline void SetBasisDataCubicSpline(
LIBSAKURA_SYMBOL(BaselineContext) *context) {
	//
}

inline void SetBasisDataSinusoid(LIBSAKURA_SYMBOL(BaselineContext) *context) {
	//
}

inline void SetBasisData(LIBSAKURA_SYMBOL(BaselineContext) *context) {
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
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order,
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
	if (num_bases > num_basis_data) {
		throw std::invalid_argument("num_bases exceeds num_basis_data!");
	}

	std::unique_ptr<LIBSAKURA_SYMBOL(BaselineContext),
	LIBSAKURA_PREFIX::Memory> work_context(
			static_cast<LIBSAKURA_SYMBOL(BaselineContext) *>(LIBSAKURA_PREFIX::Memory::Allocate(
					sizeof(LIBSAKURA_SYMBOL(BaselineContext)))));
	if (work_context == nullptr) {
		throw std::bad_alloc();
	}

	work_context->baseline_type = baseline_type;
	work_context->num_bases = num_bases;
	work_context->num_basis_data = num_basis_data;

	size_t num_total_basis_data = num_bases * num_basis_data;

	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> work_basis_data_storage(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*work_context->basis_data) * num_total_basis_data,
					&work_context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(work_context->basis_data));

	SetBasisData(work_context.get());

	work_context->basis_data_storage = work_basis_data_storage.release();
	*context = (LIBSAKURA_SYMBOL(BaselineContext) *) work_context.release();
}

inline void DestroyBaselineContext(LIBSAKURA_SYMBOL(BaselineContext) *context) {
	if (context->basis_data_storage != nullptr) {
		LIBSAKURA_PREFIX::Memory::Free(context->basis_data_storage);
	}
	LIBSAKURA_PREFIX::Memory::Free(context);
}

inline void OperateFloatSubtraction(size_t num_in, float const *in1_arg,
		float const *in2_arg, float *out_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(in1_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto in1 = AssumeAligned(in1_arg);
	auto in2 = AssumeAligned(in2_arg);
	auto out = AssumeAligned(out_arg);

	for (size_t i = 0; i < num_in; ++i) {
		out[i] = in1[i] - in2[i];
	}
}

inline void AddMulMatrix(size_t num_bases, double const *coeff_arg,
		size_t num_out, double const *basis_arg, float *out_arg) {
	auto coeff = AssumeAligned(coeff_arg);
	auto basis = AssumeAligned(basis_arg);
	auto out = AssumeAligned(out_arg);
	size_t i = 0;
#if defined(__AVX__) && !defined(ARCH_SCALAR)
	size_t const pack_elements = sizeof(__m256d) / sizeof(double);
	size_t const end = (num_out / pack_elements) * pack_elements;
	__m256d                                                         const zero = _mm256_set1_pd(0.);
	size_t const offset1 = num_bases * 1;
	size_t const offset2 = num_bases * 2;
	size_t const offset3 = num_bases * 3;
	for (i = 0; i < end; i += pack_elements) {
		__m256d total = zero;
		auto bases_row = &basis[num_bases * i];
		for (size_t j = 0; j < num_bases; ++j) {
			__m256d ce = _mm256_set1_pd(coeff[j]);
			__m256d bs = _mm256_set_pd(bases_row[j + offset3],
					bases_row[j + offset2], bases_row[j + offset1],
					bases_row[j]);
			total += ce * bs;
		}
		_mm_store_ps(&out[i], _mm256_cvtpd_ps(total));
	}
#endif
	for (; i < num_out; ++i) {
		double out_double = 0.0;
		for (size_t j = 0; j < num_bases; ++j) {
			out_double += coeff[j] * basis[num_bases * i + j];
		}
		out[i] = out_double;
	}
}

inline void CreateLeastSquareFittingCoefficients(size_t num_data,
		float const *data_arg, bool const *mask_arg,
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		double *lsq_matrix_arg, double *lsq_vector_arg) {
	if (num_data != context->num_basis_data) {
		throw std::runtime_error(
				"num_data and the model data in baseline context does not conform.");
	}
	if (num_data < context->num_bases) {
		throw std::runtime_error("not enough data for baseline fitting.");
	}
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector_arg));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto basis = AssumeAligned(context->basis_data);
	auto lsq_matrix = AssumeAligned(lsq_matrix_arg);
	auto lsq_vector = AssumeAligned(lsq_vector_arg);

	size_t num_bases = context->num_bases;

	LIBSAKURA_SYMBOL(Status) get_matrix_status =
	LIBSAKURA_SYMBOL(GetMatrixCoefficientsForLeastSquareFitting)(num_data, mask,
			num_bases, basis, lsq_matrix);
	if (get_matrix_status != LIBSAKURA_SYMBOL(Status_kOK)) {
		throw std::runtime_error(
				"failed in GetMatrixCoefficientsForLeastSquareFitting.");
	}

	LIBSAKURA_SYMBOL(Status) get_vector_status =
	LIBSAKURA_SYMBOL(GetVectorCoefficientsForLeastSquareFitting)(num_data, data,
			mask, num_bases, basis, lsq_vector);
	if (get_vector_status != LIBSAKURA_SYMBOL(Status_kOK)) {
		throw std::runtime_error(
				"failed in GetVectorCoefficientsForLeastSquareFitting.");
	}
}

inline void UpdateLeastSquareFittingCoefficients(size_t num_data,
		float const *data_arg, size_t num_clipped,
		size_t const *clipped_indices_arg,
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		double *lsq_matrix_arg, double *lsq_vector_arg) {
	if (num_data != context->num_basis_data) {
		throw std::runtime_error(
				"num_data and the model data in baseline context does not conform.");
	}
	if (num_data < context->num_bases) {
		throw std::runtime_error("not enough data for baseline fitting.");
	}
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector_arg));
	auto data = AssumeAligned(data_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);
	auto basis = AssumeAligned(context->basis_data);
	auto lsq_matrix = AssumeAligned(lsq_matrix_arg);
	auto lsq_vector = AssumeAligned(lsq_vector_arg);

	size_t num_bases = context->num_bases;

	LIBSAKURA_SYMBOL(Status) get_matrix_status =
	LIBSAKURA_SYMBOL(UpdateMatrixCoefficientsForLeastSquareFitting)(num_clipped,
			clipped_indices, num_bases, lsq_matrix, basis, lsq_matrix);
	if (get_matrix_status != LIBSAKURA_SYMBOL(Status_kOK)) {
		throw std::runtime_error(
				"failed in UpdateMatrixCoefficientsForLeastSquareFitting.");
	}

	LIBSAKURA_SYMBOL(Status) get_vector_status =
	LIBSAKURA_SYMBOL(UpdateVectorCoefficientsForLeastSquareFitting)(data,
			num_clipped, clipped_indices, num_bases, lsq_vector, basis,
			lsq_vector);
	if (get_vector_status != LIBSAKURA_SYMBOL(Status_kOK)) {
		throw std::runtime_error(
				"failed in UpdateVectorCoefficientsForLeastSquareFitting.");
	}
}

inline void GetBestFitBaseline(size_t num_data, float const *data_arg,
bool const *mask_arg, LIBSAKURA_SYMBOL(BaselineContext) const *context,
		float *out_arg, LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(context->basis_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto basis = AssumeAligned(context->basis_data);
	auto out = AssumeAligned(out_arg);

	size_t num_lsq_matrix = context->num_bases * context->num_bases;
	double *lsq_matrix = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_matrix(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*lsq_matrix) * num_lsq_matrix, &lsq_matrix));
	size_t num_lsq_vector = context->num_bases;
	double *lsq_vector = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_vector(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*lsq_vector) * num_lsq_vector, &lsq_vector));
	size_t num_coeff = context->num_bases;
	double *coeff = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_coeff(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*coeff) * num_coeff, &coeff));

	CreateLeastSquareFittingCoefficients(num_data, data, mask, context,
			lsq_matrix, lsq_vector);
	LIBSAKURA_SYMBOL(Status) solve_status =
	LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLU)(context->num_bases,
			lsq_matrix, lsq_vector, coeff);
	if (solve_status != LIBSAKURA_SYMBOL(Status_kOK)) {
		throw std::runtime_error("failed in SolveSimultaneousEquationsByLU.");
	}

	AddMulMatrix(context->num_bases, coeff, num_data, basis, out);
}

inline size_t ClipData(size_t num_data, float const *data_arg,
bool const *in_mask_arg, float const lower_bound, float const upper_bound,
bool *out_mask_arg, size_t *clipped_indices_arg) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(clipped_indices_arg));
	auto data = AssumeAligned(data_arg);
	auto in_mask = AssumeAligned(in_mask_arg);
	auto out_mask = AssumeAligned(out_mask_arg);
	auto clipped_indices = AssumeAligned(clipped_indices_arg);

	size_t num_clipped = 0;

	for (size_t i = 0; i < num_data; ++i) {
		out_mask[i] = in_mask[i];
		if (in_mask[i]) {
			if ((data[i] - lower_bound) * (upper_bound - data[i]) < 0.0f) {
				out_mask[i] = false;
				clipped_indices[num_clipped] = i;
				num_clipped++;
			}
		}
	}

	return num_clipped;
}

inline std::string GetNotEnoughDataMessage(
		uint16_t const idx_erroneous_fitting) {
	std::stringstream ss;
	ss << "SubtractBaseline: available data became too few in the ";
	ss << idx_erroneous_fitting;
	ss << " ";

	uint16_t imod10 = idx_erroneous_fitting % 10;
	std::string isuffix;
	if (imod10 == 1) {
		isuffix = "st";
	} else if (imod10 == 2) {
		isuffix = "nd";
	} else if (imod10 == 3) {
		isuffix = "rd";
	} else {
		isuffix = "th";
	}
	ss << isuffix << " fitting.";

	return ss.str();
}

inline void SubtractBaseline(size_t num_data, float const *data_arg,
bool const *mask_arg, LIBSAKURA_SYMBOL(BaselineContext) const *baseline_context,
		float clip_threshold_sigma_arg, uint16_t num_fitting_max_arg,
		bool get_residual, bool *final_mask_arg, float *out_arg,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(baseline_context->basis_data));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto final_mask = AssumeAligned(final_mask_arg);
	auto out = AssumeAligned(out_arg);
	auto basis = AssumeAligned(baseline_context->basis_data);

	float *best_fit_model = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_best_fit_model(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*best_fit_model) * num_data, &best_fit_model));

	float *residual_data = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_residual_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*residual_data) * num_data, &residual_data));

	size_t *clipped_indices = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_clipped_indices(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*clipped_indices) * num_data, &clipped_indices));

	size_t num_bases = baseline_context->num_bases;
	size_t num_lsq_matrix = num_bases * num_bases;
	double *lsq_matrix = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_matrix(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*lsq_matrix) * num_lsq_matrix, &lsq_matrix));
	size_t num_lsq_vector = num_bases;
	double *lsq_vector = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_lsq_vector(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*lsq_vector) * num_lsq_vector, &lsq_vector));
	size_t num_coeff = num_bases;
	double *coeff = nullptr;
	std::unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_coeff(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(*coeff) * num_coeff, &coeff));

	float clip_threshold_sigma = fabs((double) clip_threshold_sigma_arg);
	uint16_t num_fitting_max = std::max((uint16_t)1, num_fitting_max_arg);
	size_t num_unmasked_data = num_data;
	size_t num_clipped = num_data;

	for (size_t i = 0; i < num_data; ++i) {
		final_mask[i] = mask[i];
	}

	*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kOK);

	try {
		for (uint16_t i = 0; i < num_fitting_max; ++i) {
			if (num_unmasked_data < baseline_context->num_bases) {
				*baseline_status =
				LIBSAKURA_SYMBOL(BaselineStatus_kNotEnoughData);
				throw std::runtime_error(GetNotEnoughDataMessage(i + 1));
			}

			if (num_unmasked_data <= num_clipped) {
				CreateLeastSquareFittingCoefficients(num_data, data, final_mask,
						baseline_context, lsq_matrix, lsq_vector);
			} else {
				UpdateLeastSquareFittingCoefficients(num_data, data,
						num_clipped, clipped_indices, baseline_context,
						lsq_matrix, lsq_vector);
			}

			LIBSAKURA_SYMBOL(Status) solve_status =
					LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLU)(
							baseline_context->num_bases, lsq_matrix, lsq_vector,
							coeff);
			if (solve_status != LIBSAKURA_SYMBOL(Status_kOK)) {
				throw std::runtime_error(
						"failed in SolveSimultaneousEquationsByLU.");
			}

			AddMulMatrix(baseline_context->num_bases, coeff, num_data, basis,
					best_fit_model);
			OperateFloatSubtraction(num_data, data, best_fit_model,
					residual_data);

			LIBSAKURA_SYMBOL(StatisticsResult) result;
			LIBSAKURA_SYMBOL(Status) stat_status =
			LIBSAKURA_SYMBOL(ComputeStatistics)(num_data, residual_data,
					final_mask, &result);
			if (stat_status != LIBSAKURA_SYMBOL(Status_kOK)) {
				throw std::runtime_error("failed in ComputeStatistics.");
			}
			float clip_threshold_abs = clip_threshold_sigma * result.stddev;
			float clip_threshold_lower = result.mean - clip_threshold_abs;
			float clip_threshold_upper = result.mean + clip_threshold_abs;
			num_clipped = ClipData(num_data, residual_data, final_mask,
					clip_threshold_lower, clip_threshold_upper, final_mask,
					clipped_indices);
			if (num_clipped == 0) {
				break;
			}
			num_unmasked_data = result.count - num_clipped;
		}

		for (size_t i = 0; i < num_data; ++i) {
			out[i] = get_residual ? residual_data[i] : best_fit_model[i];
		}
	} catch (...) {
		if (*baseline_status == LIBSAKURA_SYMBOL(BaselineStatus_kOK)) {
			*baseline_status = LIBSAKURA_SYMBOL(BaselineStatus_kNG);
		}
		throw;
	}
}

inline void SubtractBaselinePolynomial(size_t num_data, float const *data_arg,
bool const *mask_arg, uint16_t order, float clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual,
		bool *final_mask_arg, float *out_arg,
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(final_mask_arg));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_arg));
	auto data = AssumeAligned(data_arg);
	auto mask = AssumeAligned(mask_arg);
	auto final_mask = AssumeAligned(final_mask_arg);
	auto out = AssumeAligned(out_arg);

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order, num_data, &context);
	ScopeGuard guard_for_baseline_context([&]() {
		DestroyBaselineContext(context);
	});

	try {
		SubtractBaseline(num_data, data, mask, context, clip_threshold_sigma,
				num_fitting_max, get_residual, final_mask, out,
				baseline_status);
	} catch (...) {
		throw;
	}
}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::CreateBaselineContext(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order,
		size_t const num_data,
		LIBSAKURA_SYMBOL(BaselineContext) **context) const {
	::CreateBaselineContext(baseline_type, order, num_data, context);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::DestroyBaselineContext(
LIBSAKURA_SYMBOL(BaselineContext) *context) const {
	::DestroyBaselineContext(context);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::GetBestFitBaseline(size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		float out[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) const {
	::GetBestFitBaseline(num_data, data, mask, context, out, baseline_status);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::SubtractBaseline(size_t num_data,
		float const data[/*num_data*/], bool const mask[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineContext) const *baseline_context,
		float clip_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual,
		bool final_mask[/*num_data*/], float out[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) const {
	::SubtractBaseline(num_data, data, mask, baseline_context,
			clip_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, baseline_status);
}

void ADDSUFFIX(Baseline, ARCH_SUFFIX)::SubtractBaselinePolynomial(
		size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], uint16_t order,
		float clip_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual,
		bool final_mask[/*num_data*/], float out[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status) const {
	::SubtractBaselinePolynomial(num_data, data, mask, order,
			clip_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, baseline_status);
}

}
