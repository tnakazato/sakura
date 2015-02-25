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
#include <iostream>
#include <string>
#include <sys/time.h>
#include <algorithm>
#include <tuple>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"

/* the number of elements in input/output array to test */
#define NUM_IN 8 // length of base data. DO NOT MODIFY!
#define NUM_IN_LONG (1 << 18) //2**18
#define UNALIGN_OFFSET 1 // should be != ALIGNMENT
using namespace std;

// Wrapper functions of OperateBitwideNot
//template<uint32_t>
static LIBSAKURA_SYMBOL(Status) NotWrapper32(uint32_t bit_mask,
		size_t num_data, uint32_t const *data, bool const *edit_mask,
		uint32_t *result) {
	return LIBSAKURA_SYMBOL(OperateBitwiseNotUint32)(num_data, data,
			edit_mask, result);
}

//template<uint8_t>
static LIBSAKURA_SYMBOL(Status) NotWrapper8(uint8_t bit_mask,
		size_t num_data, uint8_t const *data, bool const *edit_mask,
		uint8_t *result) {
	return LIBSAKURA_SYMBOL(OperateBitwiseNotUint8)(num_data, data,
			edit_mask, result);
}

// a struct to store the combination of a sakura function and answer array
template<typename DataType>
struct FuncAndAnswer {
	typedef LIBSAKURA_SYMBOL(Status) (*FuncType)(DataType, size_t, DataType const *, bool const *, DataType *);
	FuncType function;               // function
	DataType answer[NUM_IN];     // answer array
};

// a struct to store test parameters for both uint8 and uint32
struct TestComponent {
		string name;               // name of operation
		bool invert_mask;          // is operation needs inverting mask
		FuncAndAnswer<uint8_t> uint8kit;
		FuncAndAnswer<uint32_t> uint32kit;
};

//DataType const base_pattern1 = (~static_cast<DataType>(0)); // 11...111
//DataType const base_pattern00 = (~static_cast<DataType>(0) << 2); // 11...100
uint8_t const bit_pattern1_8 = (~static_cast<uint8_t>(0)); // 11...111
uint32_t const bit_pattern1_32 = (~static_cast<uint32_t>(0)); // 11...111
uint8_t const bit_pattern00_8 = (~static_cast<uint8_t>(0) << 2); // 11...100
uint32_t const bit_pattern00_32 = (~static_cast<uint32_t>(0) << 2); // 11...100

TestComponent StandardTestCase[] {
		{"AND", false,
				{LIBSAKURA_SYMBOL(OperateBitwiseAndUint8), { 0, 1, 2, 3, 0, 0, 2, 2 }},
				{LIBSAKURA_SYMBOL(OperateBitwiseAndUint32), { 0, 1, 2, 3, 0, 0, 2, 2 }}},
		{"Converse Nonimplication", false,
				{LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint8), { 0, 1, 2, 3, 2, 2, 0, 0 }},
				{LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint32), { 0, 1, 2, 3, 2, 2, 0, 0 }}},
		{"Material Implication", false,
				{LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint8),
						{ 0, 1, 2, 3, bit_pattern1_8, static_cast<uint8_t>(bit_pattern1_8 << 1),
								bit_pattern1_8, static_cast<uint8_t>(bit_pattern1_8 << 1) }},
				{LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint32),
						{ 0, 1, 2, 3, bit_pattern1_32, static_cast<uint32_t>(bit_pattern1_32 << 1),
								bit_pattern1_32, static_cast<uint32_t>(bit_pattern1_32 << 1)}}},
		{"OR", false,
				{LIBSAKURA_SYMBOL(OperateBitwiseOrUint8), { 0, 1, 2, 3, 2, 3, 2, 3 }},
				{LIBSAKURA_SYMBOL(OperateBitwiseOrUint32), { 0, 1, 2, 3, 2, 3, 2, 3 }}},
		{"XOR", false,
				{LIBSAKURA_SYMBOL(OperateBitwiseXorUint8), { 0, 1, 2, 3, 2, 3, 0, 1 }},
				{LIBSAKURA_SYMBOL(OperateBitwiseXorUint32), { 0, 1, 2, 3, 2, 3, 0, 1}}},
		{"NOT", false,
				{NotWrapper8,
						{ 0, 1, 2, 3, static_cast<uint8_t>(bit_pattern00_8 + 3), static_cast<uint8_t>(bit_pattern00_8 + 2),
								static_cast<uint8_t>(bit_pattern00_8 + 1), bit_pattern00_8 }},
				{NotWrapper32,
						{ 0, 1, 2, 3, static_cast<uint32_t>(bit_pattern00_32 + 3), static_cast<uint32_t>(bit_pattern00_32 + 2),
								static_cast<uint32_t>(bit_pattern00_32 + 1), bit_pattern00_32 }}}
};

TestComponent ExtendedTestCase[] {
		{"Material Nonimplication", true,
				{LIBSAKURA_SYMBOL(OperateBitwiseAndUint8), { 0, 1, 2, 3, 0, 1, 0, 1 }},
				{LIBSAKURA_SYMBOL(OperateBitwiseAndUint32), { 0, 1, 2, 3, 0, 1, 0, 1 }}},
		{"NOR", true,
				{LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint8),
						{ 0, 1, 2, 3, static_cast<uint8_t>(bit_pattern00_8 + 1), bit_pattern00_8,
								static_cast<uint8_t>(bit_pattern00_8 + 1), bit_pattern00_8 }},
				{LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint32),
						{ 0, 1, 2, 3, static_cast<uint32_t>(bit_pattern00_32 + 1), bit_pattern00_32,
								static_cast<uint32_t>(bit_pattern00_32 + 1), bit_pattern00_32 }}},
		{"NAND", true,
				{LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint8),
						{ 0, 1, 2, 3, bit_pattern1_8, bit_pattern1_8,
								static_cast<uint8_t>(bit_pattern1_8 - 2), static_cast<uint8_t>(bit_pattern1_8 - 2) }},
				{LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint32),
						{ 0, 1, 2, 3, bit_pattern1_32, bit_pattern1_32,
								static_cast<uint32_t>(bit_pattern1_32 - 2), static_cast<uint32_t>(bit_pattern1_32 - 2) }}},
		{"Converse Implication", true,
				{LIBSAKURA_SYMBOL(OperateBitwiseOrUint8),
						{ 0, 1, 2, 3, static_cast<uint8_t>(bit_pattern1_8 - 2), static_cast<uint8_t>(bit_pattern1_8 - 2),
								bit_pattern1_8, bit_pattern1_8 }},
				{LIBSAKURA_SYMBOL(OperateBitwiseOrUint32),
						{ 0, 1, 2, 3, static_cast<uint32_t>(bit_pattern1_32 - 2), static_cast<uint32_t>(bit_pattern1_32 - 2),
								bit_pattern1_32, bit_pattern1_32 }}},
		{"XNOR", true,
				{LIBSAKURA_SYMBOL(OperateBitwiseXorUint8),
						{ 0, 1, 2, 3,static_cast<uint8_t>(bit_pattern1_8 - 2), static_cast<uint8_t>(bit_pattern1_8 << 2),
								bit_pattern1_8,static_cast<uint8_t>(bit_pattern1_8 << 1) }},
				{LIBSAKURA_SYMBOL(OperateBitwiseXorUint32),
						{ 0, 1, 2, 3,static_cast<uint32_t>(bit_pattern1_32 - 2), static_cast<uint32_t>(bit_pattern1_32 << 2),
								bit_pattern1_32,static_cast<uint32_t>(bit_pattern1_32 << 1) }}},
};


// Functions to select proper FuncAndAnswer from TestComponents by data type and return test kit
template<typename DataType>
struct TestCase {
	typedef tuple< string, bool, FuncAndAnswer<DataType> > TestKit;
	static TestKit GetItem(size_t num_components, TestComponent const *test_components, size_t i);
};

template<>
TestCase<uint8_t>::TestKit TestCase<uint8_t>::GetItem(size_t num_components, TestComponent const *test_components, size_t i) {
	assert(i < num_components);
	auto testcase = test_components[i];
	return TestKit(testcase.name, testcase.invert_mask, testcase.uint8kit);
}

template<>
TestCase<uint32_t>::TestKit TestCase<uint32_t>::GetItem(size_t num_components, TestComponent const *test_components, size_t i) {
	assert(i < num_components);
	auto testcase = test_components[i];
	return TestKit(testcase.name, testcase.invert_mask, testcase.uint32kit);
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
			verbose(false) {
		Initialize();
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

	void Initialize() {
		// Initialize data template
		size_t const ntype(4);
		bit_mask_ = 2; /* bit pattern of 0...010 */
		for (size_t i = 0; i < NUM_IN; ++i) {
			data_[i] = i % ntype; /* repeat bit pattern of *00, *01, *10, *11,... */
			edit_mask_[i] = (i / ntype % 2 == 1); /*{F, F, F, F, T, T, T, T, (repeated)};*/
		}
	}

// Create arbitrary length of input data and edit mask by repeating in_[] and edit_mask_[]
	void GetInputDataInLength(size_t num_data, DataType *data, bool *mask) {
		size_t const num_in(NUM_IN);
		for (size_t i = 0; i < num_data; ++i) {
			data[i] = data_[i % num_in];
			mask[i] = edit_mask_[i % num_in];
		}
	}

	/* Converts an input value to a bit pattern.*/
	string BToS(DataType value) {
		char buff[bit_size + 1];
		buff[bit_size] = '\0';
		for (size_t i = 0; i < bit_size; ++i) {
			if ((value >> i) % 2 == 0)
				buff[bit_size - 1 - i] = '0';
			else
				buff[bit_size - 1 - i] = '1';
		}
		return string(buff);
	}

	/* Converts an bit pattern (char) to a value of DataType.*/
	DataType SToB(char* bit_pattern) {
		DataType result(0);
		size_t i = 0;
		while (bit_pattern[i] != '\0') {
			//result = result * 2 + ((uint8_t) in_string[i]);
			result <<= 1;
			if (bit_pattern[i] == '1')
				result += 1;
			++i;
		}
		return result;
	}

	void PrintInputs() {
		cout << "bit_mask = " << BToS(bit_mask_);
		cout << endl;
		PrintArray("data", NUM_IN, data_);
	}

	void PrintArray(char const *name, size_t num_data, DataType *data_array) {
		size_t max_length(20);
		cout << name << " = ";
		if (data_array == nullptr) { // array is nullptr
			cout << "NULL" << endl;
		} else if (num_data > max_length) { // long array (just show the length)
			cout << num_data << " elements" << endl;
		} else { // normal array
			cout << "[ ";
			if (num_data > 0) {
				for (size_t i = 0; i < num_data - 1; ++i)
					cout << BToS(data_array[i]) << ", ";
				cout << BToS(data_array[num_data - 1]);
			}
			cout << " ]" << endl;
		}
	}

	void PrintArray(char const *name, size_t num_data, bool const *data_array) {
		size_t max_length(20);
		cout << name << " = ";
		if (data_array == nullptr) { // array is nullptr
			cout << "NULL" << endl;
		} else if (num_data > max_length) { // long array (just show the length)
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
	void RunBitOperationTest(size_t num_data, DataType *in_data, bool *mask,
			DataType *out_data, size_t num_operation, TestComponent const *test_components,
			LIBSAKURA_SYMBOL(Status) return_value, size_t num_repeat = 1) {

		bool in_place(in_data == out_data);
		PrepareInputs(num_data, in_data, mask);
		if (verbose) {
			PrintArray("in", num_data, in_data);
			PrintArray("mask", num_data, mask);
		}

		if (num_repeat > 1)
			cout << "Iterating " << num_repeat
					<< " loops for each operation. The length of arrays is "
					<< num_data << endl;
		TestCase<DataType> my_testcase;
		for (size_t iop = 0; iop < num_operation; ++iop) {
			auto kit = my_testcase.GetItem(num_operation, test_components, iop);
			LIBSAKURA_SYMBOL(Status) status;
			double start, end;
			cout << "Testing bit operation: " << std::get<0>(kit) << endl;
			start = LIBSAKURA_SYMBOL(GetCurrentTime)();
			for (size_t irun = 0; irun < num_repeat; ++irun) {
				// Need to refresh data for in-place operation
				if (in_place)
					PrepareInputs(num_data, in_data, mask);
				// Actual execution of bit operation function
				status = (std::get<2>(kit).function)(
						(std::get<1>(kit) ? ~bit_mask_ : bit_mask_), num_data,
						in_data, mask, out_data);
			} // end of num_repeat loop
			end = LIBSAKURA_SYMBOL(GetCurrentTime)();
			if (num_repeat > 1)
				cout << "Elapse time of operation: " << end - start << " sec"
						<< endl;

			if (verbose) {
				if (status == LIBSAKURA_SYMBOL(Status_kOK))
					PrintArray("result", num_data, out_data);
				else
					cout << "sakura_Status = " << status << endl;
			}
			// Verification
			EXPECT_EQ(return_value, status);
			if (status == LIBSAKURA_SYMBOL(Status_kOK))
				ExactCompare(num_data, out_data, ELEMENTSOF(std::get<2>(kit).answer),
						std::get<2>(kit).answer);
		} // end of bit operation loop
	}

	void PrepareInputs(size_t num_data, DataType *data, bool *mask) {
		// Handling of nullptr array
		DataType dummy_data[num_data];
		bool dummy_mask[ELEMENTSOF(dummy_data)];
		DataType *data_ptr = (data != nullptr) ? data : dummy_data;
		bool *mask_ptr = (mask != nullptr) ? mask : dummy_mask;
		// Get array with length, num_data
		GetInputDataInLength(num_data, data_ptr, mask_ptr);
	}

	DataType data_[NUM_IN];

	bool edit_mask_[NUM_IN];

	bool verbose;
	DataType bit_mask_;
	static size_t const bit_size = sizeof(DataType) * 8;


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

		for (size_t irun = 0; irun < num_test; ++irun) {
			size_t const num_data(array_length[irun]);
			// Loop over operation types  (ALL operations)
			cout << "[Tests with array length = " << num_data << "]" << endl;
			RunBitOperationTest(num_data, data, edit_mask, result,
					ELEMENTSOF(StandardTestCase), StandardTestCase, LIBSAKURA_SYMBOL(Status_kOK));
			RunBitOperationTest(num_data, data, edit_mask, result,
					ELEMENTSOF(ExtendedTestCase), ExtendedTestCase, LIBSAKURA_SYMBOL(Status_kOK));
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
		RunBitOperationTest(num_data, data, edit_mask, data, ELEMENTSOF(StandardTestCase),
				StandardTestCase, LIBSAKURA_SYMBOL(Status_kOK));
	}

	// Long test of various bit operations with with a large array
	void RunMiscOpLongTests(size_t const num_long, size_t const num_repeat) {
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
		cout << "[Long tests]" << endl;
		RunBitOperationTest(num_large, data, edit_mask, result,
				ELEMENTSOF(StandardTestCase), StandardTestCase, LIBSAKURA_SYMBOL(Status_kOK),
				num_repeat);
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
		RunBitOperationTest(num_data, data_null, edit_mask, result,
				ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));

		cout << "[Test NULL mask array]" << endl;
		RunBitOperationTest(num_data, data, mask_null, result,
				ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));

		cout << "[Test NULL result array]" << endl;
		RunBitOperationTest(num_data, data, edit_mask, data_null,
				ELEMENTSOF(StandardTestCase), StandardTestCase,
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

		// Define unaligned array
		DataType *data_shift = &data[offset];
		assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));
		// Define unaligned array
		bool *mask_shift = &edit_mask[offset];
		assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

		// Loop over operation types  (Standard tests)
		cout << "[Test unaligned data array]" << endl;
		RunBitOperationTest(num_data, data_shift, edit_mask, result,
				ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));

		cout << "[Test unaligned mask array]" << endl;
		RunBitOperationTest(num_data, data, mask_shift, result,
				ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));

		cout << "[Test unaligned result array]" << endl;
		RunBitOperationTest(num_data, data, edit_mask, data_shift,
				ELEMENTSOF(StandardTestCase), StandardTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
	}

};

/*
 * Tests various bit operations of an uint32_t value and array
 * INPUTS:
 * - bit_mask = 0...010
 * - in = [ 0...000, 0...001, 0...010, 0...011, 0...000, 0...001, 0...010, 0...011, ...repeated... ]
 * - edit_mask = [F, F, F, F, T, T, T, T]
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
class BitOperation32: public BitOperation<uint32_t> {

};

/*
 * Tests various bit operations of an uint8_t value and array
 * INPUTS:
 * - bit_mask = 00000010
 * - in = [ 00000000, 00000001, 00000010, 00000011, 00000000, 00000001, 00000010, 00000011, ...repeated... ]
 * - edit_mask = [F, F, F, F, T, T, T, T]
 *
 * RESULTS (uint8_t):
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
class BitOperation8: public BitOperation<uint8_t> {

};

/////////////////////////////////////////////////////////////////////////////////////////
/*
 * Test various bit operations with an Uint8 array of length 8, 11 and zero
 */
//TEST_F(BitOperation, VariousLength) {
//	BitOperation<uint8_t>::RunVariousLengthTests();
//	BitOperation<uint32_t>::RunVariousLengthTests();
//}
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
 * Long test of various bit operations with with a large Uint8 array
 */
TEST_F(BitOperation8, MiscOpLong) {
	RunMiscOpLongTests(NUM_IN_LONG, 20000); // array_length, repeat
}

TEST_F(BitOperation32, MiscOpLong) {
	RunMiscOpLongTests(NUM_IN_LONG, 20000); // array_length, repeat
}

/////////////////////////////////////////////////////////////////////////////////////////
