#include <cassert>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatInRangesInclusive)(
		size_t num_data, float const data[], size_t num_condition,
		float const lower_bounds[], float const upper_bounds[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lower_bounds == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (upper_bounds == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lower_bounds)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(upper_bounds)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	// lower_bounds should be smaller or equals to corresponding upper_bounds.
	for (size_t i = 0; i < num_condition; ++i) {
		if (lower_bounds[i] > upper_bounds[i])
			return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	// Now actual operation
	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplFloat();
	bfc->SetTrueInRangesInclusive(num_data, data, num_condition, lower_bounds,
			upper_bounds, result);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntInRangesInclusive)(
		size_t num_data, int const data[], size_t num_condition,
		int const lower_bounds[], int const upper_bounds[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (lower_bounds == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (upper_bounds == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(lower_bounds)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(upper_bounds)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	// lower_bounds should be smaller or equals to corresponding upper_bounds.
	for (size_t i = 0; i < num_condition; ++i) {
		if (lower_bounds[i] > upper_bounds[i])
			return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}

	// Now actual operation
	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplInt();
	bfc->SetTrueInRangesInclusive(num_data, data, num_condition, lower_bounds,
			upper_bounds, result);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint8ToBool)(
		size_t num_data, uint8_t const data[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplUint8();
	bfc->ToBool(num_data, data, result);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint32ToBool)(
		size_t num_data, uint32_t const data[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplUint32();
	bfc->ToBool(num_data, data, result);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InvertBool)(
		size_t num_data,
		bool const data[], bool result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplUint8();
	bfc->InvertBool(num_data, data, result);

	return LIBSAKURA_SYMBOL(Status_kOK);
}
