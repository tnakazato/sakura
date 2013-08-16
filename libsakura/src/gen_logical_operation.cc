#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateLogicalAnd)(size_t num_in, bool const in1[],
		bool const in2[], bool out[]) {
	assert(in1 != nullptr);
	assert(in2 != nullptr);
	assert(out != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(in1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(in2));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	auto logicop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetLogicalOperationImpl();
	logicop->OperateLogicalAnd(num_in, in1, in2, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}
