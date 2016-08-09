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
#include <cxxabi.h>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <tuple>
#include <typeinfo>

#include <libsakura/localdef.h>
#include <libsakura/sakura.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"
#include "testutil.h"

/* the number of elements in input/output array to test */
#define NUM_IN 8 // length of base data. DO NOT MODIFY!
#define NUM_IN_LONG (1 << 18) //2**18
#define UNALIGN_OFFSET 1 // should be != ALIGNMENT
using namespace std;

/*
 * Tests various bit operations of an uint8_t or uint32_t value and array
 * INPUTS (uint8_t example):
 * - bit_mask = 00000010
 * - in = [ 00000000, 00000001, 00000010, 00000011, 00000000, 00000001, 00000010, 00000011, ...repeated... ]
 * - edit_mask = [ F, F, F, F, T, T, T, T, ...repeated... ]
 *
 * RESULTS (uint8_t example):
 * - And:
 *   [ 00000000, 00000001, 00000010, 00000011, 00000000, 00000000, 00000010, 00000010, ...repeated... ]
 * - Nonimplication:
 *   [ 00000000, 00000001, 00000010, 00000011, 00000000, 00000001, 00000000, 00000001, ...repeated... ]
 * - ConverseNonImplication:
 *   [ 00000000, 00000001, 00000010, 00000011, 00000010, 00000010, 00000000, 00000000, ...repeated... ]
 * - Nor:
 *   [ 00000000, 00000001, 00000010, 00000011, 11111101, 11111100, 11111101, 11111100, ...repeated... ]
 * - Implication:
 *   [ 00000000, 00000001, 00000010, 00000011, 11111111, 11111110, 11111111, 11111110, ...repeated... ]
 * - Nand:
 *   [ 00000000, 00000001, 00000010, 00000011, 11111111, 11111111, 11111101, 11111101, ...repeated... ]
 * - Or:
 *   [ 00000000, 00000001, 00000010, 00000011, 00000010, 00000011, 00000010, 00000011, ...repeated... ]
 * - ConverseImplication:
 *   [ 00000000, 00000001, 00000010, 00000011, 11111101, 11111101, 11111111, 11111111, ...repeated... ]
 * - Xor:
 *   [ 00000000, 00000001, 00000010, 00000011, 00000010, 00000011, 00000000, 00000001, ...repeated... ]
 * - Xnor:
 *   [ 00000000, 00000001, 00000010, 00000011, 11111101, 11111100, 11111111, 11111110, ...repeated... ]
 * - Not:
 *   [ 00000000, 00000001, 00000010, 00000011, 11111111, 11111110, 11111101, 11111100, ...repeated... ]
 *
 */

// Wrapper functions of OperateBitwiseNot to align
// interface with the other bit operation functions.
static LIBSAKURA_SYMBOL(Status) NotWrapper32(uint32_t bit_mask, size_t num_data,
		uint32_t const *data, bool const *edit_mask, uint32_t *result) {
	return LIBSAKURA_SYMBOL(OperateBitwiseNotUint32)(num_data, data, edit_mask,
			result);
}

static LIBSAKURA_SYMBOL(Status) NotWrapper8(uint8_t bit_mask, size_t num_data,
		uint8_t const *data, bool const *edit_mask, uint8_t *result) {
	return LIBSAKURA_SYMBOL(OperateBitwiseNotUint8)(num_data, data, edit_mask,
			result);
}

// a struct to store the combination of a sakura function and answer array
template<typename DataType>
struct FuncAndAnswer {
	typedef LIBSAKURA_SYMBOL(Status) (*FuncType)(DataType, size_t,
			DataType const *, bool const *, DataType *);
	FuncType function;               // function
	DataType answer[NUM_IN];     // answer array
};

// a struct to store test parameters for both uint8 and uint32
struct TestComponent {
	string name;               // name of operation
	bool invert_mask;          // is operation needs inverting mask
	FuncAndAnswer<uint8_t> uint8_kit;
	FuncAndAnswer<uint32_t> uint32_kit;
};

/*
 * Define a list of test cases (TestComponent structure)
 * - name of bit operation
 * - whether or not invert bit pattern
 * - sakura function to test and the answer (FuncAndAnswer structure) for each data type.
 *
 * Operations and answers to be tested
 * INPUTS:
 * - bit_mask = 0...010
 * - in = [ 0...000, 0...001, 0...010, 0...011, 0...000, 0...001, 0...010, 0...011, ...repeated... ]
 * - edit_mask = [ F, F, F, F, T, T, T, T, ...repeated... ]
 *
 * RESULTS:
 * - And:
 *   [ 0...000, 0...001, 0...010, 0...011, 0...000, 0...000, 0...010, 0...010, ...repeated... ]
 * - Nonimplication:
 *   [ 0...000, 0...001, 0...010, 0...011, 0...000, 0...001, 0...000, 0...001, ...repeated... ]
 * - ConverseNonImplication:
 *   [ 0...000, 0...001, 0...010, 0...011, 0...010, 0...010, 0...000, 0...000, ...repeated... ]
 * - Nor:
 *   [ 0...000, 0...001, 0...010, 0...011, 1...101, 1...100, 1...101, 1...100, ...repeated... ]
 * - Implication:
 *   [ 0...000, 0...001, 0...010, 0...011, 1...111, 1...110, 1...111, 1...110, ...repeated... ]
 * - Nand:
 *   [ 0...000, 0...001, 0...010, 0...011, 1...111, 1...111, 1...101, 1...101, ...repeated... ]
 * - Or:
 *   [ 0...000, 0...001, 0...010, 0...011, 0...010, 0...011, 0...010, 0...011, ...repeated... ]
 * - ConverseImplication:
 *   [ 0...000, 0...001, 0...010, 0...011, 1...101, 1...101, 1...111, 1...111, ...repeated... ]
 * - Xor:
 *   [ 0...000, 0...001, 0...010, 0...011, 0...010, 0...011, 0...000, 0...001, ...repeated... ]
 * - Xnor:
 *   [ 0...000, 0...001, 0...010, 0...011, 1...101, 1...100, 1...111, 1...110, ...repeated... ]
 * - Not:
 *   [ 0...000, 0...001, 0...010, 0...011, 1...111, 1...110, 1...101, 1...100, ...repeated... ]
 *
 */
// 11...111
template<typename DataType>
constexpr DataType bit111() {
	return (~static_cast<DataType>(0));
}
;
// 11...110
template<typename DataType>
constexpr DataType bit110() {
	return (~static_cast<DataType>(1));
}
;
// 11...101
template<typename DataType>
constexpr DataType bit101() {
	return (~static_cast<DataType>(2));
}
;
// 11...100
template<typename DataType>
constexpr DataType bit100() {
	return (~static_cast<DataType>(3));
}
;

TestComponent StandardTestCase[] {
		{"AND", false,
				{LIBSAKURA_SYMBOL(OperateBitwiseAndUint8),
						{ 0, 1, 2, 3, 0, 0, 2, 2 }},
				{LIBSAKURA_SYMBOL(OperateBitwiseAndUint32),
						{ 0, 1, 2, 3, 0, 0, 2, 2 }}},
		{"Converse Nonimplication", false,
				{LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint8),
						{ 0, 1, 2, 3, 2, 2, 0, 0 }},
				{LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint32),
						{ 0, 1, 2, 3, 2, 2, 0, 0 }}},
		{"Material Implication", false,
				{LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint8),
						{ 0, 1, 2, 3,
								bit111<uint8_t>(), bit110<uint8_t>(),
								bit111<uint8_t>(), bit110<uint8_t>() }},
				{LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint32),
						{ 0, 1, 2, 3,
								bit111<uint32_t>(), bit110<uint32_t>(),
								bit111<uint32_t>(), bit110<uint32_t>()}}},
		{"OR", false,
				{LIBSAKURA_SYMBOL(OperateBitwiseOrUint8),
						{ 0, 1, 2, 3, 2, 3, 2, 3 }},
				{LIBSAKURA_SYMBOL(OperateBitwiseOrUint32),
						{ 0, 1, 2, 3, 2, 3, 2, 3 }}},
		{"XOR", false,
				{LIBSAKURA_SYMBOL(OperateBitwiseXorUint8),
						{ 0, 1, 2, 3, 2, 3, 0, 1 }},
				{LIBSAKURA_SYMBOL(OperateBitwiseXorUint32),
						{ 0, 1, 2, 3, 2, 3, 0, 1}}},
		{"NOT", false,
				{NotWrapper8,
						{ 0, 1, 2, 3,
								bit111<uint8_t>(),bit110<uint8_t>(),
								bit101<uint8_t>(), bit100<uint8_t>() }},
				{NotWrapper32,
						{ 0, 1, 2, 3,
								bit111<uint32_t>(), bit110<uint32_t>(),
								bit101<uint32_t>(), bit100<uint32_t>() }}}
};

TestComponent ExtendedTestCase[] {
		{"Material Nonimplication", true,
				{LIBSAKURA_SYMBOL(OperateBitwiseAndUint8),
						{ 0, 1, 2, 3, 0, 1, 0, 1 }},
				{LIBSAKURA_SYMBOL(OperateBitwiseAndUint32),
						{ 0, 1, 2, 3, 0, 1, 0, 1 }}},
		{"NOR", true,
				{LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint8),
						{ 0, 1, 2, 3,
								bit101<uint8_t>(), bit100<uint8_t>(),
								bit101<uint8_t>(), bit100<uint8_t>() }},
				{LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint32),
						{ 0, 1, 2, 3,
								bit101<uint32_t>(), bit100<uint32_t>(),
								bit101<uint32_t>(), bit100<uint32_t>() }}},
		{"NAND", true,
				{LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint8),
						{ 0, 1, 2, 3,
								bit111<uint8_t>(), bit111<uint8_t>(),
								bit101<uint8_t>(), bit101<uint8_t>() }},
				{LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint32),
						{ 0, 1, 2, 3,
								bit111<uint32_t>(), bit111<uint32_t>(),
								bit101<uint32_t>(), bit101<uint32_t>() }}},
		{"Converse Implication", true,
				{LIBSAKURA_SYMBOL(OperateBitwiseOrUint8),
						{ 0, 1, 2, 3,
								bit101<uint8_t>(), bit101<uint8_t>(),
								bit111<uint8_t>(), bit111<uint8_t>() }},
				{LIBSAKURA_SYMBOL(OperateBitwiseOrUint32),
						{ 0, 1, 2, 3,
								bit101<uint32_t>(), bit101<uint32_t>(),
								bit111<uint32_t>(), bit111<uint32_t>() }}},
		{"XNOR", true,
				{LIBSAKURA_SYMBOL(OperateBitwiseXorUint8),
						{ 0, 1, 2, 3,
								bit101<uint8_t>(), bit100<uint8_t>(),
								bit111<uint8_t>(), bit110<uint8_t>()}},
				{LIBSAKURA_SYMBOL(OperateBitwiseXorUint32),
						{ 0, 1, 2, 3,
								bit101<uint32_t>(), bit100<uint32_t>(),
								bit111<uint32_t>(), bit110<uint32_t>() }}},
};

/*
 * Helper functions to select proper (FuncAndAnswer) from
 *  TestComponents by data type and return test kit consists of
 * - name of bit operation
 * - whether or not invert bit pattern
 * - sakura function to test and the answer (FuncAndAnswer structure) for the given data type.
 */
template<typename DataType>
struct TestCase {
	typedef tuple<string, bool, FuncAndAnswer<DataType> > TestKit;
	static TestKit GetItem(size_t num_components,
			TestComponent const *test_components, size_t i);
};

template<>
TestCase<uint8_t>::TestKit TestCase<uint8_t>::GetItem(size_t num_components,
		TestComponent const *test_components, size_t i) {
	assert(i < num_components);
	auto testcase = test_components[i];
	return TestKit(testcase.name, testcase.invert_mask, testcase.uint8_kit);
}

template<>
TestCase<uint32_t>::TestKit TestCase<uint32_t>::GetItem(size_t num_components,
		TestComponent const *test_components, size_t i) {
	assert(i < num_components);
	auto testcase = test_components[i];
	return TestKit(testcase.name, testcase.invert_mask, testcase.uint32_kit);
}

/*
 * A super class to test various bit operation of an value and array
 * INPUTS:
 * - bit_mask = 0...010
 * - in = [ 0...000, 0...001, 0...010, 0...011, 0...000, 0...001, 0...010, 0...011 ]
 * - edit_mask = [F, F, F, F, T, T, T, T]
 */
template<typename DataType>
class BitOperation: public ::testing::Test {
protected:
	BitOperation() :
			verbose_(false) {
	}

	virtual void SetUp() {
		// Initialize sakura
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Initialize)(nullptr,
				nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}

	virtual void TearDown() {
		// Clean-up sakura
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	/*
	 * Actual test definitions
	 */
	// Test various bit operations with an Uint8 array of length 8, 11 and zero
	void RunVariousLengthTests() {
		size_t const array_length[] = { NUM_IN, 11, 0 };
		size_t const num_test(ELEMENTSOF(array_length));
		size_t const num_max(11);
		//num_max = max(array_length);
		SIMD_ALIGN
		DataType data[num_max];
		SIMD_ALIGN
		bool edit_mask[ELEMENTSOF(data)];
		SIMD_ALIGN
		DataType result[ELEMENTSOF(data)];

		// Loop over various array length
		for (size_t irun = 0; irun < num_test; ++irun) {
			size_t const num_data(array_length[irun]);
			// Loop over operation types  (ALL operations)
			cout << "[Tests with array length = " << num_data << "]" << endl;
			RunBitOperationTest<OutOfPlaceAction>(num_data, data, edit_mask,
					result, ELEMENTSOF(StandardTestCase), StandardTestCase,
					LIBSAKURA_SYMBOL(Status_kOK));
			RunBitOperationTest<OutOfPlaceAction>(num_data, data, edit_mask,
					result, ELEMENTSOF(ExtendedTestCase), ExtendedTestCase,
					LIBSAKURA_SYMBOL(Status_kOK));
		}
	}

	// Test various in-place (&out == &in) bit operations with an Uint8 array of length 10
	void RunMiscOpInPlaceTests() {
		size_t const num_data(10);
		SIMD_ALIGN
		DataType data[num_data];
		SIMD_ALIGN
		bool edit_mask[ELEMENTSOF(data)];

		// Loop over minimum set of operation types  (Standard tests)
		cout << "[In-place bit operations]" << endl;
		RunBitOperationTest<InPlaceAction>(num_data, data, edit_mask, data,
				ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kOK));
	}

	// Performance test of various bit operations with with a large array
	void RunMiscOpPerformanceTests(size_t const num_long, size_t const num_repeat) {
		assert(num_long > 0);
		assert(num_repeat > 0);
		size_t const num_large(num_long);
		SIMD_ALIGN
		DataType data[num_large];
		SIMD_ALIGN
		bool edit_mask[ELEMENTSOF(data)];
		SIMD_ALIGN
		DataType result[ELEMENTSOF(data)];
		//size_t const num_repeat = 20000;

		// Loop over num_repeat times for each operation type (Standard set)
		cout << "[Performance tests]" << endl;
		RunBitOperationTest<OutOfPlaceAction>(num_large, data, edit_mask,
				result, ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kOK), num_repeat);
	}

	// Test failure cases of various bit operations: NULL data, mask, or result array
	void RunMiscOpFailNullArrayTests() {
		size_t const num_data(NUM_IN);
		SIMD_ALIGN
		DataType data[num_data];
		SIMD_ALIGN
		bool edit_mask[ELEMENTSOF(data)];
		SIMD_ALIGN
		DataType result[ELEMENTSOF(data)];

		DataType *data_null = nullptr;
		bool *mask_null = nullptr;

		// Loop over operation types  (Standard tests)
		cout << "[Test NULL data array]" << endl;
		RunBitOperationTest<OutOfPlaceAction>(num_data, data_null, edit_mask,
				result, ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));

		cout << "[Test NULL mask array]" << endl;
		RunBitOperationTest<OutOfPlaceAction>(num_data, data, mask_null, result,
				ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));

		cout << "[Test NULL result array]" << endl;
		RunBitOperationTest<OutOfPlaceAction>(num_data, data, edit_mask,
				data_null, ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
	}

	// Test failure cases of various bit operations: Unaligned data array
	void RunMiscOpFailNotAlignedArrayTests() {
		size_t offset(UNALIGN_OFFSET);
		size_t const num_data(NUM_IN);
		size_t const num_elements(num_data + offset);
		SIMD_ALIGN
		DataType data[num_elements];
		SIMD_ALIGN
		bool edit_mask[ELEMENTSOF(data)];
		SIMD_ALIGN
		DataType result[ELEMENTSOF(data)];

		// Define unaligned arrays
		DataType *data_shift = &data[offset];
		bool *mask_shift = &edit_mask[offset];
		assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));
		assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

		// Loop over operation types  (Standard tests)
		cout << "[Test unaligned data array]" << endl;
		RunBitOperationTest<OutOfPlaceAction>(num_data, data_shift, edit_mask,
				result, ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));

		cout << "[Test unaligned mask array]" << endl;
		RunBitOperationTest<OutOfPlaceAction>(num_data, data, mask_shift,
				result, ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));

		cout << "[Test unaligned result array]" << endl;
		RunBitOperationTest<OutOfPlaceAction>(num_data, data, edit_mask,
				data_shift, ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
	}

private:
	// Set values to an arbitrary length of input data and edit mask.
	// Handle nullptr case properly.
	static void PrepareInputs(size_t num_data, DataType *data, bool *mask) {
		// Handling of nullptr array
		DataType dummy_data[num_data];
		bool dummy_mask[ELEMENTSOF(dummy_data)];
		DataType *data_ptr = (data != nullptr) ? data : dummy_data;
		bool *mask_ptr = (mask != nullptr) ? mask : dummy_mask;
		// Get array with length, num_data
		GetInputDataInLength(num_data, data_ptr, mask_ptr);
	}
	// Set values to an arbitrary length of input data and edit mask arrays
	// by repeating data_[] and edit_mask_[]
	static void GetInputDataInLength(size_t num_data, DataType *data,
	bool *mask) {
		/* repeat bit pattern of *00, *01, *10, *11,... */
		static constexpr DataType data_[] = { 0, 1, 2, 3, 0, 1, 2, 3 };
		/* boolean array of {F, F, F, F, T, T, T, T} */
		static constexpr bool edit_mask_[] = { false, false, false, false, true, true, true, true };
		STATIC_ASSERT(ELEMENTSOF(data_)==NUM_IN);
		STATIC_ASSERT(ELEMENTSOF(edit_mask_)==ELEMENTSOF(data_));

		size_t const num_in(ELEMENTSOF(data_));
		for (size_t i = 0; i < num_data; ++i) {
			data[i] = data_[i % num_in];
			mask[i] = edit_mask_[i % num_in];
		}
	}

	/* Converts an input value to a bit pattern.*/
	string ConvertToBitPattern(DataType value) {
		char buff[bit_size_ + 1];
		buff[bit_size_] = '\0';
		for (size_t i = 0; i < bit_size_; ++i) {
			if ((value >> i) % 2 == 0)
				buff[bit_size_ - 1 - i] = '0';
			else
				buff[bit_size_ - 1 - i] = '1';
		}
		return string(buff);
	}
//	/* Converts an bit pattern (char) to a value of DataType.*/
//	DataType ConvertBitPaternToNumber(char* bit_pattern) {
//		DataType result(0);
//		size_t i = 0;
//		while (bit_pattern[i] != '\0') {
//			//result = result * 2 + ((uint8_t) in_string[i]);
//			result <<= 1;
//			if (bit_pattern[i] == '1')
//				result += 1;
//			++i;
//		}
//		return result;
//	}

	// Convert an array to a string formatted for printing. Convert numeric values to bit patterns.
	void PrintArray(char const *name, size_t num_data, DataType *data_array) {
		constexpr size_t kMaxLength(20);
		cout << name << " = ";
		if (data_array == nullptr) { // array is nullptr
			cout << "NULL" << endl;
		} else if (num_data > kMaxLength) { // long array (just show the length)
			cout << num_data << " elements" << endl;
		} else { // normal array
			cout << "[ ";
			if (num_data > 0) {
				for (size_t i = 0; i < num_data - 1; ++i)
					cout << ConvertToBitPattern(data_array[i]) << ", ";
				cout << ConvertToBitPattern(data_array[num_data - 1]);
			}
			cout << " ]" << endl;
		}
	}
	// Convert an array to a string formatted for printing. Convert booleans to T/F strings.
	void PrintArray(char const *name, size_t num_data, bool const *data_array) {
		constexpr size_t kMaxLength(20);
		cout << name << " = ";
		if (data_array == nullptr) { // array is nullptr
			cout << "NULL" << endl;
		} else if (num_data > kMaxLength) { // long array (just show the length)
			cout << num_data << " elements" << endl;
		} else { // normal array
			cout << "[ ";
			if (num_data > 0) {
				for (size_t i = 0; i < num_data - 1; ++i)
					cout << (data_array[i] ? "T" : "F") << ", ";
				cout << (data_array[num_data - 1] ? "T" : "F");
			}
			cout << " ]" << endl;
		}
	}


	void string_replace(string &invalue, string const &from, string const &to) {
		assert(from.length()==to.length());
		size_t pos = invalue.find(from, 0);
		while (pos < string::npos) {
			invalue.replace(pos, to.length(), to);
			pos = invalue.find(from, pos);
		}
	}

	/*
	 * Compare data with reference array, and assert values of corresponding
	 * elements are the exact match.
	 * If num_data > num_reference, elements of reference_array
	 * are repeated from the beginning as many times as necessary.
	 */
	void ExactCompare(size_t num_data, DataType const *data_array,
			size_t num_reference, DataType const *reference_array) {
		for (size_t i = 0; i < num_data; ++i) {
			ASSERT_EQ(reference_array[i % num_reference], data_array[i]);
		}
	}

	/*
	 * A function to run a list of bit operations and compare result with expected answer.
	 */
	template<typename InitializeAction>
	void RunBitOperationTest(size_t num_data, DataType *in_data, bool *mask,
			DataType *out_data, size_t num_operation,
			TestComponent const *test_components,
			LIBSAKURA_SYMBOL(Status) return_value, size_t num_repeat = 1) {

		/* bit pattern of 0...010 */
		static constexpr DataType bit_mask_ = 2;

		PrepareInputs(num_data, in_data, mask);
		if (verbose_) {
			PrintArray("data", num_data, in_data);
			PrintArray("mask", num_data, mask);
		}
		int success = 0 ;
		string data_type_name = abi::__cxa_demangle(typeid(DataType).name(), nullptr, nullptr, &success);
		string_replace(data_type_name, " ", "_");

		if (num_repeat > 1)
			cout << "Iterating " << num_repeat
					<< " loops for each operation. The length of arrays is "
					<< num_data << endl;
		TestCase<DataType> my_testcase;
		for (size_t iop = 0; iop < num_operation; ++iop) {
			auto kit = my_testcase.GetItem(num_operation, test_components, iop);
			LIBSAKURA_SYMBOL(Status) status;
			cout << "Testing bit operation: " << std::get<0>(kit) << endl;
			double start = GetCurrentTime();
			for (size_t irun = 0; irun < num_repeat; ++irun) {
				// Need to refresh data for in-place operation
				InitializeAction::reinitialize(num_data, in_data, mask);
				// Actual execution of bit operation function
				status = (std::get<2>(kit).function)(
						(std::get<1>(kit) ? ~bit_mask_ : bit_mask_), num_data,
						in_data, mask, out_data);
			} // end of num_repeat loop
			double end = GetCurrentTime();
			if (num_repeat > 1){
				string test_name = std::get<0>(kit);
				string_replace(test_name, " ", "_");
				cout << "#x# benchmark Bit_" << test_name
				<< "_" << data_type_name
				<< " " << end - start << endl;
			}
			if (verbose_) {
				if (status == LIBSAKURA_SYMBOL(Status_kOK))
					PrintArray("result", num_data, out_data);
				else
					cout << "sakura_Status = " << status << endl;
			}
			// Verification
			EXPECT_EQ(return_value, status);
			if (status == LIBSAKURA_SYMBOL(Status_kOK))
				ExactCompare(num_data, out_data,
						ELEMENTSOF(std::get<2>(kit).answer),
						std::get<2>(kit).answer);
		} // end of bit operation loop
	}
	// Re-initialization of input data and mask for the case of in-place operation.
	struct InPlaceAction {
		static void reinitialize(size_t num_data, DataType *data, bool *mask) {
			PrepareInputs(num_data, data, mask);
		}
	};
	// Dummy function (does nothing) for re-initialization for the case of out-of-place operation.
	struct OutOfPlaceAction {
		static void reinitialize(size_t num_data, DataType *data, bool *mask) {
			// no need to initialize
		}
	};

	// Member variables
//	static constexpr DataType bit_mask_ = 2;
//	static DataType data_[]; //[NUM_IN];
//	static bool edit_mask_[]; //[NUM_IN];
	bool verbose_;
	static constexpr size_t bit_size_ = sizeof(DataType) * 8;

};

///*
// * Initialization of static input data templates
// */
///* bit pattern of 0...010 */
//template<typename DataType>
//DataType BitOperation<DataType>::bit_mask_ = 2;
///* repeat bit pattern of *00, *01, *10, *11,... */
//template<typename DataType>
//DataType BitOperation<DataType>::data_[] = { 0, 1, 2, 3, 0, 1, 2, 3 };
///* boolean array of {F, F, F, F, T, T, T, T} */
//template<typename DataType>
//bool BitOperation<DataType>::edit_mask_[] = { false, false, false, false, true,
//true, true, true };

/*
 * Tests various bit operations of an uint32_t value and array
 */
class BitOperation32: public BitOperation<uint32_t> {

};

/*
 * Tests various bit operations of an uint8_t value and array
 */
class BitOperation8: public BitOperation<uint8_t> {

};

/////////////////////////////////////////////////////////////////////////////////////////
/*
 * Test various bit operations with an Uint8 array of length 8, 11 and zero
 */
TEST_F(BitOperation8, VariousLength) {
	RunVariousLengthTests();
}

TEST_F(BitOperation32, VariousLength) {
	RunVariousLengthTests();
}
///////////////////////////////////////////////////////////////////////////////////////////
/*
 * Test various in-place (&out == &in) bit operations with an Uint8 array of length 10
 */
TEST_F(BitOperation8, MiscOpInPlace) {
	RunMiscOpInPlaceTests();
}

TEST_F(BitOperation32, MiscOpInPlace) {
	RunMiscOpInPlaceTests();
}
/////////////////////////////////////////////////////////////////////////////////////////
/*
 * Test failure cases of various bit operations: NULL data, mask, or result array
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
TEST_F(BitOperation8, MiscOpFailNullArray) {
	RunMiscOpFailNullArrayTests();
}

TEST_F(BitOperation32, MiscOpFailNullArray) {
	RunMiscOpFailNullArrayTests();
}
/////////////////////////////////////////////////////////////////////////////////////////
/*
 * Test failure cases of various bit operations: Unaligned data array
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
TEST_F(BitOperation8, MiscOpFailNotAlignedArray) {
	RunMiscOpFailNotAlignedArrayTests();
}

TEST_F(BitOperation32, MiscOpFailNotAlignedArray) {
	RunMiscOpFailNotAlignedArrayTests();
}
/////////////////////////////////////////////////////////////////////////////////////////
/*
 * Performance test of various bit operations with with a large Uint8 array
 */
TEST_F(BitOperation8, MiscOpPerformance) {
	RunMiscOpPerformanceTests(NUM_IN_LONG, 20000); // array_length, repeat
}

TEST_F(BitOperation32, MiscOpPerformance) {
	RunMiscOpPerformanceTests(NUM_IN_LONG, 20000); // array_length, repeat
}
/////////////////////////////////////////////////////////////////////////////////////////
