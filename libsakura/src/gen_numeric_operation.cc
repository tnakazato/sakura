#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/localdef.h"
#include "libsakura/logger.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/sakura.h"

namespace {
auto logger = LIBSAKURA_PREFIX::Logger::GetLogger("numeric_operation");
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetLeastSquareFittingCoefficients)(
		size_t const num_data, float const data[], bool const mask[],
		size_t const num_model_bases, double const basis_data[],
		double lsq_matrix[], double lsq_vector[]) {
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (basis_data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(basis_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lsq_matrix == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lsq_vector == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lsq_vector)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto numop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetNumericOperationImpl();
	try {
		numop->GetLeastSquareFittingCoefficients(num_data, data, mask,
				num_model_bases, basis_data, lsq_matrix, lsq_vector);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UpdateLeastSquareFittingCoefficients)(
		size_t const num_data, float const data[],
		size_t const num_exclude_indices, size_t const exclude_indices[],
		size_t const num_model_bases, double const basis_data[],
		double lsq_matrix[], double lsq_vector[]) {
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (exclude_indices == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(exclude_indices)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		if (exclude_indices[i] >= num_data)
			return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	if (basis_data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(basis_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lsq_matrix == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lsq_vector == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lsq_vector)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto numop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetNumericOperationImpl();
	try {
		numop->UpdateLeastSquareFittingCoefficients(num_data, data,
				num_exclude_indices, exclude_indices, num_model_bases,
				basis_data, lsq_matrix, lsq_vector);
	} catch (const std::runtime_error &e) {
		LOG4CXX_ERROR(logger, e.what());
		return LIBSAKURA_SYMBOL(Status_kNG);
	} catch (...) {
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLU)(
		size_t num_equations, double const in_matrix[],
		double const in_vector[], double out[]) {
	if (in_matrix == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in_matrix)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (in_vector == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in_vector)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto numop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetNumericOperationImpl();
	try {
		numop->SolveSimultaneousEquationsByLU(num_equations, in_matrix,
				in_vector, out);
	} catch (...) {
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}

	return LIBSAKURA_SYMBOL(Status_kOK);
}
