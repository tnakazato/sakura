/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2016
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
#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstddef>
#include <iostream>

#include "libsakura/localdef.h"
#include "libsakura/sakura.h"

// Vectorization by Compiler
namespace {

// Check parameter arguments.
template<typename DataType>
bool IsValidArguments(DataType const *data, bool const *edit_mask,
		DataType const *result) {
	if (data == nullptr || result == nullptr || edit_mask == nullptr) {
		return false;
	}
	if (!(LIBSAKURA_SYMBOL(IsAligned)(data))
			|| !(LIBSAKURA_SYMBOL(IsAligned)(result))
			|| !(LIBSAKURA_SYMBOL(IsAligned)(edit_mask))) {
		return false;
	}
	return true;
}

/* @brief A function to invoke various types of bit operation
 *
 *  This function loops over elements of an array, @a data, and call a function,
 *  @a operation , with each @a data element.
 *  In case an exception is thrown during operation, the function aborts
 *  if assertion is enabled. if not, it returns kUnknownError status.
 *
 *  @param[in] operation A function which defines bit operation.
 *  It should take a @a data and a @a edit_mask elements and return a boolean.
 */
template<typename Operation, typename DataType>
LIBSAKURA_SYMBOL(Status) DoBitOperation(Operation operation, size_t num_data,
		DataType const *data,
		bool const *edit_mask, DataType *result) {
	if (!IsValidArguments(data, edit_mask, result)) {
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument);
	}
	try {
		uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
		assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
		auto adata = AssumeAligned(data);
		auto aresult = AssumeAligned(result);
		auto amask = AssumeAligned(mask8);
		STATIC_ASSERT(sizeof(edit_mask[0]) == sizeof(amask[0]));
		STATIC_ASSERT(true == 1);
		STATIC_ASSERT(false == 0);
		for (size_t i = 0; i < num_data; ++i) {
			aresult[i] = operation(adata[i], amask[i]);
		}
	} catch (...) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

}
/* anonymous namespace */

// Bit operation AND
/* edit_mask = true: (mask8 - 1) = 00...0
 *                   -> (bit_mask | 00...0) = bit_mask,
 *           = false: (mask8 - 1) = 11...1
 *                   -> (bit_mask | 11...1) = 11...1
 *        where mask8 = (uint8_t) edit_mask */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseAndUint8)(
	uint8_t bit_mask, size_t num_data, uint8_t const data[],
	bool const edit_mask[], uint8_t result[]) noexcept {
return DoBitOperation(
		[bit_mask](decltype(data[0]) data_value, uint8_t mask) {
			return (data_value & (bit_mask | (static_cast<decltype(data_value)>(mask) - 1)));
		}, num_data, data, edit_mask, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseAndUint32)(
	uint32_t bit_mask, size_t num_data, uint32_t const data[],
	bool const edit_mask[], uint32_t result[]) noexcept {
return DoBitOperation(
		[bit_mask](decltype(data[0]) data_value, uint8_t mask) {
			return (data_value & (bit_mask | (static_cast<decltype(data_value)>(mask) - 1)));
		}, num_data, data, edit_mask, result);
}

// Bit operation, Converse Nonimplication
/* edit_mask = true: x = (mask8 - 1) = 00...0
 *                   -> (bit_mask | 00...0) = bit_mask,
 *                   ~x = 11...1 -> (data[i] ^ ~x) = ~data[i]
 *           = false: x = (mask8 - 1) = 11...1
 *                   -> (bit_mask | 11...1) = 11...1
 *                   ~x = 00...0 -> (data[i] ^ ~x) = data[i]
 *        where mask8 = (uint8_t) edit_mask  */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint8)(
	uint8_t bit_mask, size_t num_data, uint8_t const data[],
	bool const edit_mask[], uint8_t result[]) noexcept {
return DoBitOperation(
		[bit_mask](decltype(data[0]) data_value, uint8_t mask) {
			decltype(data_value) x = (static_cast<decltype(data_value)>(mask) - 1);
			return ((data_value ^ ~x) & (bit_mask | x));
		}, num_data, data, edit_mask, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint32)(
	uint32_t bit_mask, size_t num_data, uint32_t const data[],
	bool const edit_mask[], uint32_t result[]) noexcept {
return DoBitOperation(
		[bit_mask](decltype(data[0]) data_value, uint8_t mask) {
			decltype(data_value) x = (static_cast<decltype(data_value)>(mask) - 1);
			return ((data_value ^ ~x) & (bit_mask | x));
		}, num_data, data, edit_mask, result);
}

// Bit operation, Material Implication
/* edit_mask = true: x = ~(mask8 - 1) = 11...1
 *                   -> (bit_mask & 11...1) = bit_mask,
 *                   -> (data[i] ^ 11...1) = ~data[i]
 *           = false: x = ~(mask8 - 1) = 00...0
 *                   -> (bit_mask & 00...0) = 00...0
 *                   (data[i] ^ 00...0) = data[i]
 *        where mask8 = (uint8_t) edit_mask  */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint8)(
	uint8_t bit_mask, size_t num_data, uint8_t const data[],
	bool const edit_mask[], uint8_t result[]) noexcept {
return DoBitOperation(
		[bit_mask](decltype(data[0]) data_value, uint8_t mask) {
			decltype(data_value) x = ~(static_cast<decltype(data_value)>(mask) - 1);
			return ((data_value ^ x) | (bit_mask & x));
		}, num_data, data, edit_mask, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint32)(
	uint32_t bit_mask, size_t num_data, uint32_t const data[],
	bool const edit_mask[], uint32_t result[]) noexcept {
return DoBitOperation(
		[bit_mask](decltype(data[0]) data_value, uint8_t mask) {
			decltype(data_value) x = ~(static_cast<decltype(data_value)>(mask) - 1);
			return ((data_value ^ x) | (bit_mask & x));
		}, num_data, data, edit_mask, result);
}

// Bit operation NOT
/* edit_mask = true: ~(mask8 - 1) = 11...1
 *                   -> (data[i] ^ 11...1) = ~data[i],
 *           = false: ~(mask8 - 1) = 00...0
 *                   -> (data[i] ^ 00...0) = data[i]
 *        where mask8 = (uint8_t) edit_mask  */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseNotUint8)(
	size_t num_data, uint8_t const data[], bool const edit_mask[],
	uint8_t result[]) noexcept {
return DoBitOperation([](decltype(data[0]) data_value, uint8_t mask) {
			return (data_value ^ ~(static_cast<decltype(data_value)>(mask) - 1));
		}, num_data, data, edit_mask, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseNotUint32)(
	size_t num_data, uint32_t const data[], bool const edit_mask[],
	uint32_t result[]) noexcept {
return DoBitOperation([](decltype(data[0]) data_value, uint8_t mask) {
			return (data_value ^ ~(static_cast<decltype(data_value)>(mask) - 1));
		}, num_data, data, edit_mask, result);
}

// Bit operation OR
/* edit_mask = true: ~(mask8 - 1) = 11...1
 *                   -> (bit_mask & 11...1) = bit_mask,
 *           = false: ~(mask8 - 1) = 00...0
 *                   -> (bit_mask & 00...0) = 00...0
 *        where mask8 = (uint8_t) edit_mask  */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseOrUint8)(
	uint8_t bit_mask, size_t num_data, uint8_t const data[],
	bool const edit_mask[], uint8_t result[]) noexcept {
return DoBitOperation(
		[bit_mask](decltype(data[0]) data_value, uint8_t mask) {
			return (data_value | (bit_mask & ~(static_cast<decltype(data_value)>(mask) - 1)));
		}, num_data, data, edit_mask, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseOrUint32)(
	uint32_t bit_mask, size_t num_data, uint32_t const data[],
	bool const edit_mask[], uint32_t result[]) noexcept {
return DoBitOperation(
		[bit_mask](decltype(data[0]) data_value, uint8_t mask) {
			return (data_value | (bit_mask & ~(static_cast<decltype(data_value)>(mask) - 1)));
		}, num_data, data, edit_mask, result);
}

/* Bit operation XOR */
/* edit_mask = true: ~(mask8 - 1) = 11...1
 *                   -> (bit_mask & 11...1) = bit_mask,
 *           = false: ~(mask8 - 1) = 00...0
 *                   -> (bit_mask & 00...0) = 00...0
 *        where mask8 = (uint8_t) edit_mask
 */
extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseXorUint8)(
	uint8_t bit_mask, size_t num_data, uint8_t const data[],
	bool const edit_mask[], uint8_t result[]) noexcept {
return DoBitOperation(
		[bit_mask](decltype(data[0]) data_value, uint8_t mask) {
			return (data_value ^ (bit_mask & ~(static_cast<decltype(data_value)>(mask) - 1)));
		}, num_data, data, edit_mask, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseXorUint32)(
	uint32_t bit_mask, size_t num_data, uint32_t const data[],
	bool const edit_mask[], uint32_t result[]) noexcept {
return DoBitOperation(
		[bit_mask](decltype(data[0]) data_value, uint8_t mask) {
			return (data_value ^ (bit_mask & ~(static_cast<decltype(data_value)>(mask) - 1)));
		}, num_data, data, edit_mask, result);
}
