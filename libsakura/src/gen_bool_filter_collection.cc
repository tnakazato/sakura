/*
 * @SAKURA_LICENSE_HEADER_START@
 * @SAKURA_LICENSE_HEADER_END@
 */
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
	try {
		bfc->SetTrueInRangesInclusive(num_data, data, num_condition,
				lower_bounds, upper_bounds, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
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
	try {
		bfc->SetTrueInRangesInclusive(num_data, data, num_condition,
				lower_bounds, upper_bounds, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatInRangesExclusive)(
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
	try {
		bfc->SetTrueInRangesExclusive(num_data, data, num_condition,
				lower_bounds, upper_bounds, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntInRangesExclusive)(
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
	try {
		bfc->SetTrueInRangesExclusive(num_data, data, num_condition,
				lower_bounds, upper_bounds, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatGreaterThan)(
		size_t num_data, float const data[], float threshold,
		bool result[]) {
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
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplFloat();
	try {
		bfc->SetTrueGreaterThan(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntGreaterThan)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) {
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
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplInt();
	try {
		bfc->SetTrueGreaterThan(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatGreaterThanOrEquals)(
		size_t num_data, float const data[], float threshold,
		bool result[]) {
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
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplFloat();
	try {
		bfc->SetTrueGreaterThanOrEquals(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntGreaterThanOrEquals)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) {
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
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplInt();
	try {
		bfc->SetTrueGreaterThanOrEquals(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatLessThan)(
		size_t num_data, float const data[], float threshold,
		bool result[]) {
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
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplFloat();
	try {
		bfc->SetTrueLessThan(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntLessThan)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) {
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
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplInt();
	try {
		bfc->SetTrueLessThan(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatLessThanOrEquals)(
		size_t num_data, float const data[], float threshold,
		bool result[]) {
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
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplFloat();
	try {
		bfc->SetTrueLessThanOrEquals(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntLessThanOrEquals)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) {
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
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplInt();
	try {
		bfc->SetTrueLessThanOrEquals(num_data, data, threshold, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetFalseFloatIfNanOrInf)(
		size_t num_data,
		float const data[], bool result[]) {
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
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBoolFilterCollectionImplFloat();
	try {
		bfc->SetFalseIfNanOrInf(num_data, data, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
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
	try {
		bfc->ToBool(num_data, data, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
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
	try {
		bfc->ToBool(num_data, data, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
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
	try {
		bfc->InvertBool(num_data, data, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
