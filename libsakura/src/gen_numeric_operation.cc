#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateFloatSubtraction)(size_t num_in, float const in1[],
		float const in2[], float out[]) {
	if (in1 == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (in2 == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in1)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in2)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto numop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetNumericOperationImpl();
	numop->OperateFloatSubtraction(num_in, in1, in2, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetLeastSquareMatrix)(size_t num_in,
		float const in_data[], bool const in_mask[], size_t num_model, double const model[],
		double out[], double out_vector[]) {
	if (in_data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (in_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (model == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out_vector == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(model)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out_vector)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto numop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetNumericOperationImpl();
	numop->GetLeastSquareMatrix(num_in, in_data, in_mask, num_model, model, out, out_vector);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLU)(size_t num_eqn,
		double const lsq_matrix0[], double const lsq_vector0[], double out[]) {
	if (lsq_matrix0 == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lsq_vector0 == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix0)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lsq_vector0)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto numop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetNumericOperationImpl();
	numop->SolveSimultaneousEquationsByLU(num_eqn, lsq_matrix0, lsq_vector0, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DoGetBestFitModel)(size_t num_chan,
		size_t num_eqn, double const model[], double const coeff[], float out[]) {
	if (model == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (coeff == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(model)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(coeff)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto numop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetNumericOperationImpl();
	numop->DoGetBestFitModel(num_chan, num_eqn, model, coeff, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitModel)(size_t num_in,
		float const in_data[], bool const in_mask[], size_t num_model,
		double const model[], float out[]) {
	if (in_data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (in_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (model == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in_data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(model)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	auto numop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetNumericOperationImpl();
	numop->GetBestFitModel(num_in, in_data, in_mask, num_model, model, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}
