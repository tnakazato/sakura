#include <iostream>
#include <string>
#include <sys/time.h>
#include <algorithm>

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
		for (size_t iop=0; iop < NUM_OPERATION;++iop){
			operation_functions[iop] = nullptr;
		}
	}

	typedef LIBSAKURA_SYMBOL(Status) (*function_ptr_t)(DataType, size_t,
			DataType const*, bool const*, DataType*);

//	virtual LIBSAKURA_SYMBOL(Status) NotWrapper(DataType bit_mask, size_t num_data,
//			DataType const *data, bool const *edit_mask, DataType *result) = 0;

// Types of operation
	enum operation_type {
		And = 0,
		NonImplication,
		ConverseNonImplication,
		Nor,
		Implication,
		Nand,
		Or,
		ConverseImplication,
		Xor,
		Xnor,
		Not, /// Not has different interface
		NUM_OPERATION
	};

	struct testKit {
		string name;               // name of operation
		function_ptr_t function;   // pointer to function
		DataType answer[NUM_IN];   // answer array
		bool invert_mask;          // is operation needs inverting mask
	};

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

	testKit GetTestKit(operation_type const operation) {
		testKit ret_kit;
		// Answers for each operation
		DataType const base_pattern1 = (~static_cast<DataType>(0)); // 11...111
		DataType const base_pattern00 = (~static_cast<DataType>(0) << 2); // 11...100

		switch (operation) {
		case And: {
			DataType answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
			STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
			ret_kit.name = "AND";
			ret_kit.invert_mask = false;
			for (size_t i = 0; i < NUM_IN; ++i)
				ret_kit.answer[i] = answer[i];
			break;
		}
		case NonImplication: {
			DataType answer[] = { 0, 1, 2, 3, 0, 1, 0, 1 };
			STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
			ret_kit.name = "Material Nonimplication";
			ret_kit.invert_mask = true;
			for (size_t i = 0; i < NUM_IN; ++i)
				ret_kit.answer[i] = answer[i];
			break;
		}
		case ConverseNonImplication: {
			DataType answer[] = { 0, 1, 2, 3, 2, 2, 0, 0 };
			STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
			ret_kit.name = "Converse Nonimplication";
			ret_kit.invert_mask = false;
			for (size_t i = 0; i < NUM_IN; ++i)
				ret_kit.answer[i] = answer[i];
			break;
		}
		case Nor: {
			DataType answer[] = { 0, 1, 2, 3,
					static_cast<DataType>(base_pattern00 + 1), base_pattern00,
					static_cast<DataType>(base_pattern00 + 1), base_pattern00 };
			STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
			ret_kit.name = "NOR";
			ret_kit.invert_mask = true;
			for (size_t i = 0; i < NUM_IN; ++i)
				ret_kit.answer[i] = answer[i];
			break;
		}
		case Implication: {
			DataType answer[] = { 0, 1, 2, 3, base_pattern1,
					static_cast<DataType>(base_pattern1 << 1), base_pattern1,
					static_cast<DataType>(base_pattern1 << 1) };
			STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
			ret_kit.name = "Material Implication";
			ret_kit.invert_mask = false;
			for (size_t i = 0; i < NUM_IN; ++i)
				ret_kit.answer[i] = answer[i];
			break;
		}
		case Nand: {
			DataType answer[] = { 0, 1, 2, 3, base_pattern1, base_pattern1,
					static_cast<DataType>(base_pattern1 - 2),
					static_cast<DataType>(base_pattern1 - 2) };
			STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
			ret_kit.name = "NAND";
			ret_kit.invert_mask = true;
			for (size_t i = 0; i < NUM_IN; ++i)
				ret_kit.answer[i] = answer[i];
			break;
		}
		case Or: {
			DataType answer[] = { 0, 1, 2, 3, 2, 3, 2, 3 };
			STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
			ret_kit.name = "OR";
			ret_kit.invert_mask = false;
			for (size_t i = 0; i < NUM_IN; ++i)
				ret_kit.answer[i] = answer[i];
			break;
		}
		case ConverseImplication: {
			DataType answer[] = { 0, 1, 2, 3,
					static_cast<DataType>(base_pattern1 - 2),
					static_cast<DataType>(base_pattern1 - 2), base_pattern1,
					base_pattern1 };
			STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
			ret_kit.name = "Converse Implication";
			ret_kit.invert_mask = true;
			for (size_t i = 0; i < NUM_IN; ++i)
				ret_kit.answer[i] = answer[i];
			break;
		}
		case Xor: {
			DataType answer[] = { 0, 1, 2, 3, 2, 3, 0, 1 };
			STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
			ret_kit.name = "XOR";
			ret_kit.invert_mask = false;
			for (size_t i = 0; i < NUM_IN; ++i)
				ret_kit.answer[i] = answer[i];
			break;
		}
		case Xnor: {
			DataType answer[] = { 0, 1, 2, 3,
					static_cast<DataType>(base_pattern1 - 2),
					static_cast<DataType>(base_pattern1 << 2), base_pattern1,
					static_cast<DataType>(base_pattern1 << 1) };
			STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
			ret_kit.name = "XNOR";
			ret_kit.invert_mask = true;
			for (size_t i = 0; i < NUM_IN; ++i)
				ret_kit.answer[i] = answer[i];
			break;
		}
		case Not: {
			DataType answer[] = { 0, 1, 2, 3,
					static_cast<DataType>(base_pattern00 + 3),
					static_cast<DataType>(base_pattern00 + 2),
					static_cast<DataType>(base_pattern00 + 1), base_pattern00 };
			STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
			ret_kit.name = "NOT";
			ret_kit.invert_mask = false;
			for (size_t i = 0; i < NUM_IN; ++i)
				ret_kit.answer[i] = answer[i];
			break;
		}
		default:
			ADD_FAILURE()<< "Invalid operation type.";
		}

		ret_kit.function = operation_functions[operation];
		return ret_kit;
	}

	void Initialize() {
		// Initialize data template
		size_t const ntype(4);
		bit_mask_ = 2; /* bit pattern of 0...010 */
		for (size_t i = 0; i < NUM_IN; ++i) {
			data_[i] = i % ntype; /* repeat bit pattern of *00, *01, *10, *11,... */
			edit_mask_[i] = (i / ntype % 2 == 1); /*{F, F, F, F, T, T, T, T, (repeated)};*/
		}
		// Initialize TestKit
		InitKit();
	}
	void InitKit() {
		for (size_t iop = 0; iop < num_operation; ++iop) {
			all_kit[iop] = GetTestKit((operation_type) iop);
		}
		for (size_t iop = 0; iop < num_standard_test; ++iop) {
			misc_kit[iop] = all_kit[test_list_standard[iop]];
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
//	char* BToS(DataType in_value) {
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
	 * A function ot run a list of bit operations and compare result with expected answer.
	 */
	void RunBitOperationTest(size_t num_data, DataType *in_data, bool *mask,
			DataType *out_data, size_t num_operation, testKit const *kits,
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
		for (size_t iop = 0; iop < num_operation; ++iop) {
			testKit kit = kits[iop];
			LIBSAKURA_SYMBOL(Status) status;
			double start, end;
			cout << "Testing bit operation: " << kit.name << endl;
			start = LIBSAKURA_SYMBOL(GetCurrentTime)();
			for (size_t irun = 0; irun < num_repeat; ++irun) {
				// Need to refresh data for in-place operation
				if (in_place)
					PrepareInputs(num_data, in_data, mask);
				// Actual execution of bit operation function
				status = (*kit.function)(
						(kit.invert_mask ? ~bit_mask_ : bit_mask_), num_data,
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
				ExactCompare(num_data, out_data, ELEMENTSOF(kit.answer),
						kit.answer);
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

	operation_type const num_operation = NUM_OPERATION;
	function_ptr_t operation_functions[NUM_OPERATION];
	testKit all_kit[NUM_OPERATION];
	// The standard test list
	static size_t const num_standard_test = 6;
	operation_type test_list_standard[num_standard_test] = { And,
			ConverseNonImplication, Implication, Or, Xor, Not };
	testKit misc_kit[num_standard_test];
};

/*
 * Tests various bit operations (except for NOT) of an uint32_t value and array
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

protected:
	BitOperation32() {
		operation_functions[And] = LIBSAKURA_SYMBOL(OperateBitsUint32And);
		operation_functions[NonImplication] =
		LIBSAKURA_SYMBOL(OperateBitsUint32And);
		operation_functions[ConverseNonImplication] =
		LIBSAKURA_SYMBOL(OperateBitsUint32ConverseNonImplication);
		operation_functions[Nor] =
		LIBSAKURA_SYMBOL(OperateBitsUint32ConverseNonImplication);
		operation_functions[Implication] =
		LIBSAKURA_SYMBOL(OperateBitsUint32Implication);
		operation_functions[Nand] =
		LIBSAKURA_SYMBOL(OperateBitsUint32Implication);
		operation_functions[Or] = LIBSAKURA_SYMBOL(OperateBitsUint32Or);
		operation_functions[ConverseImplication] =
		LIBSAKURA_SYMBOL(OperateBitsUint32Or);
		operation_functions[Xor] = LIBSAKURA_SYMBOL(OperateBitsUint32Xor);
		operation_functions[Xnor] = LIBSAKURA_SYMBOL(OperateBitsUint32Xor);
		operation_functions[Not] = Not32Wrapper;

		Initialize();
	}

	static LIBSAKURA_SYMBOL(Status) Not32Wrapper(uint32_t bit_mask,
			size_t num_data, uint32_t const *data, bool const *edit_mask,
			uint32_t *result) {
		return LIBSAKURA_SYMBOL(OperateBitsUint32Not)(num_data, data, edit_mask,
				result);
	}
};

/*
 * Tests various bit operations of an uint8_t value and array
 * INPUTS:
 * - bit_mask = 00000010
 * - in = [ 00000000, 00000001, 00000010, 00000011, 00000000, 00000001, 00000010, 00000011, ...repeated... ]
 * - edit_mask = [F, F, F, F, T, T, T, T]
 *
 * RESULTS:
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
 */
class BitOperation8: public BitOperation<uint8_t> {

protected:
	BitOperation8() {
		operation_functions[And] = LIBSAKURA_SYMBOL(OperateBitsUint8And);
		operation_functions[NonImplication] =
		LIBSAKURA_SYMBOL(OperateBitsUint8And);
		operation_functions[ConverseNonImplication] =
		LIBSAKURA_SYMBOL(OperateBitsUint8ConverseNonImplication);
		operation_functions[Nor] =
		LIBSAKURA_SYMBOL(OperateBitsUint8ConverseNonImplication);
		operation_functions[Implication] =
		LIBSAKURA_SYMBOL(OperateBitsUint8Implication);
		operation_functions[Nand] =
		LIBSAKURA_SYMBOL(OperateBitsUint8Implication);
		operation_functions[Or] = LIBSAKURA_SYMBOL(OperateBitsUint8Or);
		operation_functions[ConverseImplication] =
		LIBSAKURA_SYMBOL(OperateBitsUint8Or);
		operation_functions[Xor] = LIBSAKURA_SYMBOL(OperateBitsUint8Xor);
		operation_functions[Xnor] = LIBSAKURA_SYMBOL(OperateBitsUint8Xor);
		operation_functions[Not] = Not8Wrapper;

		Initialize();
	}
	static LIBSAKURA_SYMBOL(Status) Not8Wrapper(uint8_t bit_mask,
			size_t num_data, uint8_t const *data, bool const *edit_mask,
			uint8_t *result) {
		return LIBSAKURA_SYMBOL(OperateBitsUint8Not)(num_data, data, edit_mask,
				result);
	}
};

/////////////////////////////////////////////////////////////////////////////////////////
/*
 * Test various bit operations with an Uint8 array of length 8, 11 and zero
 */
//TEST_F(BitOperation, VariousLength) {
//	BitOperation<uint8_t>::RunVariousLength();
//	BitOperation<uint32_t>::RunVariousLength();
//}
TEST_F(BitOperation8, VariousLength) {
	size_t const array_length[] = { NUM_IN, 11, 0 };
	size_t const num_test(ELEMENTSOF(array_length));
	size_t const num_max(11);
	//num_max = max(array_length);
	SIMD_ALIGN
	uint8_t data[num_max];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	for (size_t irun = 0; irun < num_test; ++irun) {
		size_t const num_data(array_length[irun]);
		// Loop over operation types  (ALL operations)
		cout << "[Tests with array length = " << num_data << "]" << endl;
		RunBitOperationTest(num_data, data, edit_mask, result, num_operation,
				all_kit, LIBSAKURA_SYMBOL(Status_kOK));
	}
}

/*
 * Test various bit operations with an Uint32 array of length 8, 11 and zero
 */
TEST_F(BitOperation32, VariousLength) {
	size_t const array_length[] = { NUM_IN, 11, 0 };
	size_t const num_test(ELEMENTSOF(array_length));
	size_t const num_max(11);
	//num_max = max(array_length);
	SIMD_ALIGN
	uint32_t data[num_max];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	for (size_t irun = 0; irun < num_test; ++irun) {
		size_t const num_data(array_length[irun]);
		// Loop over operation types  (ALL operations)
		cout << "[Tests with array length = " << num_data << "]" << endl;
		RunBitOperationTest(num_data, data, edit_mask, result, num_operation,
				all_kit, LIBSAKURA_SYMBOL(Status_kOK));
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
/*
 * Test various in-place (&out == &in) bit operations with an Uint8 array of length 10
 */
TEST_F(BitOperation8, MiscOpInPlace) {
	size_t const num_data(10);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];

	// Loop over operation types  (Standard tests)
	RunBitOperationTest(num_data, data, edit_mask, data, num_standard_test,
			misc_kit, LIBSAKURA_SYMBOL(Status_kOK));
}
/*
 * Test various in-place (&out == &in) bit operations with an Uint32 array of length 10
 */
TEST_F(BitOperation32, MiscOpInPlace) {
	size_t const num_data(10);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];

	// Loop over operation types  (Standard tests)
	RunBitOperationTest(num_data, data, edit_mask, data, num_standard_test,
			misc_kit, LIBSAKURA_SYMBOL(Status_kOK));
}
/////////////////////////////////////////////////////////////////////////////////////////
/*
 * Test failure cases of various bit operations: NULL data, mask, or result array
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
TEST_F(BitOperation8, MiscOpFailNullArray) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t *data_null = nullptr;
	bool *mask_null = nullptr;

	// Loop over operation types  (Standard tests)
	cout << "[Test NULL data array]" << endl;
	RunBitOperationTest(num_data, data_null, edit_mask, result,
			num_standard_test, misc_kit,
			LIBSAKURA_SYMBOL(Status_kInvalidArgument));

	cout << "[Test NULL mask array]" << endl;
	RunBitOperationTest(num_data, data, mask_null, result, num_standard_test,
			misc_kit, LIBSAKURA_SYMBOL(Status_kInvalidArgument));

	cout << "[Test NULL result array]" << endl;
	RunBitOperationTest(num_data, data, edit_mask, data_null,
			num_standard_test, misc_kit,
			LIBSAKURA_SYMBOL(Status_kInvalidArgument));
}

/*
 * Test failure cases of various bit operations: NULL data, mask, or result array
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
TEST_F(BitOperation32, MiscOpFailNullArray) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t *data_null = nullptr;
	bool *mask_null = nullptr;

	// Loop over operation types  (Standard tests)
	cout << "[Test NULL data array]" << endl;
	RunBitOperationTest(num_data, data_null, edit_mask, result,
			num_standard_test, misc_kit,
			LIBSAKURA_SYMBOL(Status_kInvalidArgument));

	cout << "[Test NULL mask array]" << endl;
	RunBitOperationTest(num_data, data, mask_null, result, num_standard_test,
			misc_kit, LIBSAKURA_SYMBOL(Status_kInvalidArgument));

	cout << "[Test NULL result array]" << endl;
	RunBitOperationTest(num_data, data, edit_mask, data_null,
			num_standard_test, misc_kit,
			LIBSAKURA_SYMBOL(Status_kInvalidArgument));
}

/////////////////////////////////////////////////////////////////////////////////////////
/*
 * Test failure cases of various bit operations: Unaligned data array
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
TEST_F(BitOperation8, MiscOpFailNotAlignedArray) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	// Define unaligned array
	uint8_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));
	// Define unaligned array
	bool *mask_shift = &edit_mask[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

	// Loop over operation types  (Standard tests)
	cout << "[Test unaligned data array]" << endl;
	RunBitOperationTest(num_data, data_shift, edit_mask, result,
			num_standard_test, misc_kit,
			LIBSAKURA_SYMBOL(Status_kInvalidArgument));

	cout << "[Test unaligned mask array]" << endl;
	RunBitOperationTest(num_data, data, mask_shift, result, num_standard_test,
			misc_kit,
			LIBSAKURA_SYMBOL(Status_kInvalidArgument));

	cout << "[Test unaligned result array]" << endl;
	RunBitOperationTest(num_data, data, edit_mask, data_shift,
			num_standard_test, misc_kit,
			LIBSAKURA_SYMBOL(Status_kInvalidArgument));
}

TEST_F(BitOperation32, MiscOpFailNotAlignedArray) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	// Define unaligned array
	uint32_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));
	// Define unaligned array
	bool *mask_shift = &edit_mask[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

	// Loop over operation types  (Standard tests)
	cout << "[Test unaligned data array]" << endl;
	RunBitOperationTest(num_data, data_shift, edit_mask, result,
			num_standard_test, misc_kit,
			LIBSAKURA_SYMBOL(Status_kInvalidArgument));

	cout << "[Test unaligned mask array]" << endl;
	RunBitOperationTest(num_data, data, mask_shift, result, num_standard_test,
			misc_kit,
			LIBSAKURA_SYMBOL(Status_kInvalidArgument));

	cout << "[Test unaligned result array]" << endl;
	RunBitOperationTest(num_data, data, edit_mask, data_shift,
			num_standard_test, misc_kit,
			LIBSAKURA_SYMBOL(Status_kInvalidArgument));
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*
 * Long test of various bit operations with with a large Uint8 array
 */
TEST_F(BitOperation8, MiscOpLong) {
	size_t const num_large(NUM_IN_LONG);
	SIMD_ALIGN
	uint8_t data[num_large];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];
	size_t const num_repeat = 20000;

	// Loop over num_repeat times for each operation type (Standard set)
	RunBitOperationTest(num_large, data, edit_mask, result, num_standard_test,
			misc_kit,
			LIBSAKURA_SYMBOL(Status_kOK), num_repeat);
}

/*
 * Long test of various bit operations with with a large Uint32 array
 */

TEST_F(BitOperation32, MiscOpLong) {
	size_t const num_large(NUM_IN_LONG);
	SIMD_ALIGN
	uint32_t data[num_large];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];
	size_t const num_repeat = 20000;

	// Loop over num_repeat times for each operation type (Standard set)
	RunBitOperationTest(num_large, data, edit_mask, result, num_standard_test,
			misc_kit,
			LIBSAKURA_SYMBOL(Status_kOK), num_repeat);
}
