#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8And)(uint8_t bit_mask,
		size_t num_in, uint8_t const in[], bool const edit_mask[], uint8_t out[]) {
	assert(in != nullptr);
	assert(out != nullptr);
	assert(edit_mask != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	/* need to include CHECK_ARGS defined in gen_gridding.cc
	CHECK_ARGS(in != nullptr);
	CHECK_ARGS(out != nullptr);
	CHECK_ARGS(edit_mask != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(in));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(out));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	*/

	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImpl();
	bitop->OperateBitsAnd(bit_mask, num_in, in, edit_mask, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}






