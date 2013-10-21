#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateLogicalAnd)(size_t num_in, bool const in1[],
		bool const in2[], bool out[]) {
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

	auto logicop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetLogicalOperationImpl();
	logicop->OperateLogicalAnd(num_in, in1, in2, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}
