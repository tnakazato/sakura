#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBoolsAnd)(size_t num_in, bool const in1[],
		bool const in2[], bool out[]) {
	/*
	assert(in != nullptr);
	assert(out != nullptr);
	assert(edit_mask != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	*/
	/* need to include CHECK_ARGS defined in gen_gridding.cc
	CHECK_ARGS(in != nullptr);
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(edit_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(in));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	*/

	auto basicop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBasicOperationImpl();
	basicop->OperateBoolsAnd(num_in, in1, in2, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}
