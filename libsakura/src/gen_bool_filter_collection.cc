#include <cassert>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatInRangesInclusive)(
		size_t num_data, float const data[], size_t num_condition,
		float const lower_bounds[], float const upper_bounds[], bool result[]) {
	// Check parameter arguments.
	if (num_data == 0)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
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
	if (num_data == 0)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
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

	// Now actual operation
	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplInt();
	bfc->SetTrueInRangesInclusive(num_data, data, num_condition, lower_bounds,
			upper_bounds, result);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint8ToBool)(size_t num_in,
		uint8_t const in[], bool out[]) {
	// Check parameter arguments.
	if (num_in == 0)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (in == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplUint8();
	bfc->ToBool(num_in, in, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint32ToBool)(
		size_t num_in, uint32_t const in[], bool out[]) {
	// Check parameter arguments.
	if (num_in == 0)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (in == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplUint32();
	bfc->ToBool(num_in, in, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InvertBool)(size_t num_in,
		bool const in[], bool out[]) {
	// Check parameter arguments.
	if (num_in == 0)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (in == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (out == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(in)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(out)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bfc =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplUint8();
	bfc->InvertBool(num_in, in, out);

	return LIBSAKURA_SYMBOL(Status_kOK);
}
