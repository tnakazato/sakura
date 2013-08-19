#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateFloatSubtraction)(size_t num_in, float const in1[],
		float const in2[], float out[]) {
	assert(in1 != nullptr);
	assert(in2 != nullptr);
	assert(out != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(in1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	auto numop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetNumericOperationImpl();
	numop->OperateFloatSubtraction(num_in, in1, in2, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetLeastSquareMatrix)(size_t num_in,
		float const in_data[], bool const in_mask[], size_t num_model, double const model[],
		double out[], double out_vector[]) {
	assert(in_data    != nullptr);
	assert(in_mask    != nullptr);
	assert(model      != nullptr);
	assert(out        != nullptr);
	assert(out_vector != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in_mask));
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out_vector));

	auto getlsmop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetNumericOperationImpl();
	getlsmop->GetLeastSquareMatrix(num_in, in_data, in_mask, num_model, model, out, out_vector);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLU)(size_t num_eqn,
		double const lsq_matrix0[], double const lsq_vector0[], double out[]) {
	assert(lsq_matrix0 != nullptr);
	assert(lsq_vector0 != nullptr);
	assert(out         != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_matrix0));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lsq_vector0));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	auto solveluop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetNumericOperationImpl();
	solveluop->SolveSimultaneousEquationsByLU(num_eqn, lsq_matrix0, lsq_vector0, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DoGetBestFitModel)(size_t num_chan,
		size_t num_eqn, double const model[], double const coeff[], float out[]) {
	assert(model != nullptr);
	assert(coeff != nullptr);
	assert(out   != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(model));
	assert(LIBSAKURA_SYMBOL(IsAligned)(coeff));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	return LIBSAKURA_SYMBOL(Status_kOK);
}
