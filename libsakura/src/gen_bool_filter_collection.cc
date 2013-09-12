#include <cassert>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatInRangesInclusive)(
		size_t num_data, float const data[], size_t num_condition,
		float const lower_bounds[], float const upper_bounds[], bool result[]) {
	assert(data != nullptr);
	assert(result != nullptr);
	assert(lower_bounds != nullptr);
	assert(upper_bounds != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lower_bounds));
	assert(LIBSAKURA_SYMBOL(IsAligned)(upper_bounds));

	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplFloat();
	bfc->SetTrueInRangesInclusive(num_data, data, num_condition, lower_bounds, upper_bounds, result);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntInRangesInclusive)(
		size_t num_data, int const data[], size_t num_condition,
		int const lower_bounds[], int const upper_bounds[], bool result[]) {
	assert(data != nullptr);
	assert(result != nullptr);
	assert(lower_bounds != nullptr);
	assert(upper_bounds != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(lower_bounds));
	assert(LIBSAKURA_SYMBOL(IsAligned)(upper_bounds));

	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplInt();
	bfc->SetTrueInRangesInclusive(num_data, data, num_condition, lower_bounds, upper_bounds, result);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint8ToBool)(
		size_t num_in, uint8_t const in[], bool out[]) {
	assert(in != nullptr);
	assert(out != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplUint8();
	bfc->ToBool(num_in, in, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint32ToBool)(
		size_t num_in, uint32_t const in[], bool out[]) {
	assert(in != nullptr);
	assert(out != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplUint32();
	bfc->ToBool(num_in, in, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InvertBool)(
		size_t num_in, bool const in[], bool out[]) {
	assert(in != nullptr);
	assert(out != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(in));
	assert(LIBSAKURA_SYMBOL(IsAligned)(out));

	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplUint8();
	bfc->InvertBool(num_in, in, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}
