/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2014
 * National Astronomical Observatory of Japan
 * 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
 * 
 * This file is part of Sakura.
 * 
 * Sakura is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * Sakura is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
 * @SAKURA_LICENSE_HEADER_END@
 */
#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

// Bit operation AND
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8And)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[],
		bool const edit_mask[], uint8_t result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (edit_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(edit_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImplUint8();
	try {
		bitop->OperateBitsAnd(bit_mask, num_data, data, edit_mask, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32And)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[],
		bool const edit_mask[], uint32_t result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (edit_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(edit_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImplUint32();
	try {
		bitop->OperateBitsAnd(bit_mask, num_data, data, edit_mask, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

// Bit operation, Converse Nonimplication
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8ConverseNonImplication)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[],
		bool const edit_mask[], uint8_t result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (edit_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(edit_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImplUint8();
	try {
		bitop->OperateBitsConverseNonImplication(bit_mask, num_data, data,
				edit_mask, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32ConverseNonImplication)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[],
		bool const edit_mask[], uint32_t result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (edit_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(edit_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImplUint32();
	try {
		bitop->OperateBitsConverseNonImplication(bit_mask, num_data, data,
				edit_mask, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

// Bit operation, Material Implication
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8Implication)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[],
		bool const edit_mask[], uint8_t result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (edit_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(edit_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImplUint8();
	try {
		bitop->OperateBitsImplication(bit_mask, num_data, data, edit_mask,
				result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32Implication)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[],
		bool const edit_mask[], uint32_t result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (edit_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(edit_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImplUint32();
	try {
		bitop->OperateBitsImplication(bit_mask, num_data, data, edit_mask,
				result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

// Bit operation NOT
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8Not)(
		size_t num_data, uint8_t const data[], bool const edit_mask[],
		uint8_t result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (edit_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(edit_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImplUint8();
	try {
		bitop->OperateBitsNot(num_data, data, edit_mask, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32Not)(
		size_t num_data, uint32_t const data[], bool const edit_mask[],
		uint32_t result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (edit_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(edit_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImplUint32();
	try {
		bitop->OperateBitsNot(num_data, data, edit_mask, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

// Bit operation OR
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8Or)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[],
		bool const edit_mask[], uint8_t result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (edit_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(edit_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImplUint8();
	try {
		bitop->OperateBitsOr(bit_mask, num_data, data, edit_mask, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32Or)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[],
		bool const edit_mask[], uint32_t result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (edit_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(edit_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImplUint32();
	try {
		bitop->OperateBitsOr(bit_mask, num_data, data, edit_mask, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

// Bit operation XOR
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8Xor)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[],
		bool const edit_mask[], uint8_t result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (edit_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(edit_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImplUint8();
	try {
		bitop->OperateBitsXor(bit_mask, num_data, data, edit_mask, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32Xor)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[],
		bool const edit_mask[], uint32_t result[]) {
	// Check parameter arguments.
	if (data == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (result == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (edit_mask == nullptr)
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(data)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(result)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	if (!( LIBSAKURA_SYMBOL(IsAligned)(edit_mask)))
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);

	// Now actual operation
	auto bitop =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetBitOperationImplUint32();
	try {
		bitop->OperateBitsXor(bit_mask, num_data, data, edit_mask, result);
	} catch (...) {
		// an exception is thrown during operation
		// abort if assertion is enabled. if not, return kUnknownError status.
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
