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
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

// Vectorization by Compiler
namespace {

template<typename DataType>
inline void OperateBitsAnd(DataType bit_mask, size_t num_data,
		DataType const *data,
		bool const *edit_mask, DataType *result) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	// cast bool array to uint8_t array
	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
	auto adata = AssumeAligned(data);
	auto aresult = AssumeAligned(result);
	auto amask = AssumeAligned(mask8);
	STATIC_ASSERT(sizeof(edit_mask[0]) == sizeof(amask[0]));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);

	/* edit_mask = true: (mask8 - 1) = 00...0
	 *                   -> (bit_mask | 00...0) = bit_mask,
	 *           = false: (mask8 - 1) = 11...1
	 *                   -> (bit_mask | 11...1) = 11...1 */
	for (size_t i = 0; i < num_data; ++i) {
		aresult[i] = adata[i]
				& (bit_mask | (static_cast<DataType>(amask[i]) - 1));
	}

}
//template<typename DataType>
//inline void OperateBitsAnd(DataType bit_mask, size_t num_data,
//		DataType const *data,
//		bool const *edit_mask, DataType *result) {
//
//	//std::cout << "Invoking OperateBitsAnd()" << std::endl;
//	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
//	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
//	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
//	// cast bool array to uint8_t array
//	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
//	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
//	STATIC_ASSERT(sizeof(edit_mask[0]) == sizeof(mask8[0]));
//	STATIC_ASSERT(true == 1);
//	STATIC_ASSERT(false == 0);
//
//	/* edit_mask = true: (mask8 - 1) = 00...0
//	 *                   -> (bit_mask | 00...0) = bit_mask,
//	 *           = false: (mask8 - 1) = 11...1
//	 *                   -> (bit_mask | 11...1) = 11...1 */
//	for (size_t i = 0; i < num_data; ++i) {
//		result[i] = data[i]
//				& (bit_mask | (static_cast<DataType>(mask8[i]) - 1));
//	}
//
//}

template<typename DataType>
inline void OperateBitsConverseNonImplication(DataType bit_mask,
		size_t num_data, DataType const *data, bool const *edit_mask,
		DataType *result) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	// cast bool array to uint8_t array
	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
	STATIC_ASSERT(sizeof(edit_mask[0]) == sizeof(mask8[0]));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);

	/* edit_mask = true: X = (mask8 - 1) = 00...0
	 *                   -> (bit_mask | 00...0) = bit_mask,
	 *                   ~X = 11...1 -> (data[i] ^ ~X) = ~data[i]
	 *           = false: X = (mask8 - 1) = 11...1
	 *                   -> (bit_mask | 11...1) = 11...1
	 *                   ~X = 00...0 -> (data[i] ^ ~X) = data[i] */
	for (size_t i = 0; i < num_data; ++i) {
		DataType X = (static_cast<DataType>(mask8[i]) - 1);
		result[i] = (data[i] ^ ~X) & (bit_mask | X);
	}

}

template<typename DataType>
inline void OperateBitsImplication(DataType bit_mask, size_t num_data,
		DataType const *data, bool const *edit_mask, DataType *result) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	// cast bool array to uint8_t array
	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
	STATIC_ASSERT(sizeof(edit_mask[0]) == sizeof(mask8[0]));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);

	/* edit_mask = true: X = ~(mask8 - 1) = 11...1
	 *                   -> (bit_mask & 11...1) = bit_mask,
	 *                   -> (data[i] ^ 11...1) = ~data[i]
	 *           = false: X = ~(mask8 - 1) = 00...0
	 *                   -> (bit_mask & 00...0) = 00...0
	 *                   (data[i] ^ 00...0) = data[i] */
	for (size_t i = 0; i < num_data; ++i) {
		DataType X = ~(static_cast<DataType>(mask8[i]) - 1);
		result[i] = (data[i] ^ X) | (bit_mask & X);
	}

}

template<typename DataType>
inline void OperateBitsNot(size_t num_data, DataType const *data,
bool const *edit_mask, DataType *result) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	// cast bool array to uint8_t array
	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
	STATIC_ASSERT(sizeof(edit_mask[0]) == sizeof(mask8[0]));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);

	/* edit_mask = true: ~(mask8 - 1) = 11...1
	 *                   -> (data[i] ^ 11...1) = ~data[i],
	 *           = false: ~(mask8 - 1) = 00...0
	 *                   -> (data[i] ^ 00...0) = data[i] */
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = (data[i] ^ ~(static_cast<DataType>(mask8[i]) - 1));
	}

}

template<typename DataType>
inline void OperateBitsOr(DataType bit_mask, size_t num_data,
		DataType const *data,
		bool const *edit_mask, DataType *result) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	// cast bool array to uint8_t array
	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
	STATIC_ASSERT(sizeof(edit_mask[0]) == sizeof(mask8[0]));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);

	/* edit_mask = true: ~(mask8 - 1) = 11...1
	 *                   -> (bit_mask & 11...1) = bit_mask,
	 *           = false: ~(mask8 - 1) = 00...0
	 *                   -> (bit_mask & 00...0) = 00...0 */
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = data[i]
				| (bit_mask & ~(static_cast<DataType>(mask8[i]) - 1));
	}

}

template<typename DataType>
inline void OperateBitsXor(DataType bit_mask, size_t num_data,
		DataType const *data,
		bool const *edit_mask, DataType *result) {

	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(result));
	assert(LIBSAKURA_SYMBOL(IsAligned)(edit_mask));
	// cast bool array to uint8_t array
	uint8_t const *mask8 = reinterpret_cast<uint8_t const *>(edit_mask);
	assert(LIBSAKURA_SYMBOL(IsAligned)(mask8));
	STATIC_ASSERT(sizeof(edit_mask[0]) == sizeof(mask8[0]));
	STATIC_ASSERT(true == 1);
	STATIC_ASSERT(false == 0);

	/* edit_mask = true: ~(mask8 - 1) = 11...1
	 *                   -> (bit_mask & 11...1) = bit_mask,
	 *           = false: ~(mask8 - 1) = 00...0
	 *                   -> (bit_mask & 00...0) = 00...0 */
	for (size_t i = 0; i < num_data; ++i) {
		result[i] = data[i]
				^ (bit_mask & ~(static_cast<DataType>(mask8[i]) - 1));
	}

}

} /* anonymous namespace */

namespace LIBSAKURA_PREFIX {
template<typename DataType>
void ADDSUFFIX(BitOperation, ARCH_SUFFIX)<DataType>::OperateBitsAnd(
		DataType bit_mask, size_t num_data, DataType const data[/*num_data*/],
		bool const edit_mask[/*num_data*/],
		DataType result[/*num_data*/]) const {
	::OperateBitsAnd(bit_mask, num_data, data, edit_mask, result);
}

template<typename DataType>
void ADDSUFFIX(BitOperation, ARCH_SUFFIX)<DataType>::OperateBitsConverseNonImplication(
		DataType bit_mask, size_t num_data, DataType const data[/*num_data*/],
		bool const edit_mask[/*num_data*/],
		DataType result[/*num_data*/]) const {
	::OperateBitsConverseNonImplication(bit_mask, num_data, data, edit_mask,
			result);
}

template<typename DataType>
void ADDSUFFIX(BitOperation, ARCH_SUFFIX)<DataType>::OperateBitsImplication(
		DataType bit_mask, size_t num_data, DataType const data[/*num_data*/],
		bool const edit_mask[/*num_data*/],
		DataType result[/*num_data*/]) const {
	::OperateBitsImplication(bit_mask, num_data, data, edit_mask, result);
}

template<typename DataType>
void ADDSUFFIX(BitOperation, ARCH_SUFFIX)<DataType>::OperateBitsNot(
		size_t num_data, DataType const data[/*num_data*/],
		bool const edit_mask[/*num_data*/],
		DataType result[/*num_data*/]) const {
	::OperateBitsNot(num_data, data, edit_mask, result);
}

template<typename DataType>
void ADDSUFFIX(BitOperation, ARCH_SUFFIX)<DataType>::OperateBitsOr(
		DataType bit_mask, size_t num_data, DataType const data[/*num_data*/],
		bool const edit_mask[/*num_data*/],
		DataType result[/*num_data*/]) const {
	::OperateBitsOr(bit_mask, num_data, data, edit_mask, result);
}

template<typename DataType>
void ADDSUFFIX(BitOperation, ARCH_SUFFIX)<DataType>::OperateBitsXor(
		DataType bit_mask, size_t num_data, DataType const data[/*num_data*/],
		bool const edit_mask[/*num_data*/],
		DataType result[/*num_data*/]) const {
	::OperateBitsXor(bit_mask, num_data, data, edit_mask, result);
}

template class ADDSUFFIX(BitOperation, ARCH_SUFFIX)<uint8_t> ;
template class ADDSUFFIX(BitOperation, ARCH_SUFFIX)<uint32_t> ;
}
