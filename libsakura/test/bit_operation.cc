#include <iostream>
#include <string>
#include <sys/time.h>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
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
	}

	virtual void SetUp() {
		size_t const ntype(4);
		bit_mask_ = 2; /* bit pattern of 0...010 */
		for (size_t i = 0; i < NUM_IN; ++i) {
			data_[i] = i % ntype; /* repeat bit pattern of *00, *01, *10, *11,... */
			edit_mask_[i] = (i / ntype % 2 == 1); /*{F, F, F, F, T, T, T, T, (repeated)};*/
		}

		// Initialize sakura
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Initialize)(nullptr,
				nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}

	virtual void TearDown() {
		// Clean-up sakura
		LIBSAKURA_SYMBOL(CleanUp)();
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
		cout << name << " = [";
		for (size_t i = 0; i < num_data - 1; ++i)
			cout << BToS(data_array[i]) << ", ";
		cout << BToS(data_array[num_data - 1]);
		cout << " ]" << endl;
	}

	void PrintArray(char const *name, size_t num_data, bool const *data_array) {
		cout << name << " = [ ";
		for (size_t i = 0; i < num_data - 1; ++i)
			cout << (data_array[i] ? "T" : "F") << ", ";
		cout << (data_array[num_data - 1] ? "T" : "F");
		cout << " ]" << endl;
	}

	DataType data_[NUM_IN];bool edit_mask_[NUM_IN];

	bool verbose;
	DataType bit_mask_;
	static size_t const bit_size = sizeof(DataType) * 8;

};

/*
 * Tests various bit operation of an uint32_t value and array
 * INPUTS:
 * - bit_mask = 0...010
 * - in = [ 0...000, 0...001, 0...010, 0...011, 0...000, 0...001, 0...010, 0...011 ]
 * - edit_mask = [F, F, F, F, T, T, T, T]
 */
class BitOperation32: public BitOperation<uint32_t> {
};

/*
 * Tests various bit operation of an uint8_t value and array
 * INPUTS:
 * - bit_mask = 00000010
 * - in = [ 00000000, 00000001, 00000010, 00000011, 00000000, 00000001, 00000010, 00000011 ]
 * - edit_mask = [F, F, F, F, T, T, T, T]
 */
class BitOperation8: public BitOperation<uint8_t> {
};

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bit operation AND by sakura_OperateBitsUint8And
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000000, 00000010, 00000010 ]
 */
TEST_F(BitOperation8, And) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation AND by sakura_OperateBitsUint8And in-place operation (&out == &in)
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000000, 00000010, 00000010 ]
 */
TEST_F(BitOperation8, AndInPlace) {
	size_t const num_data(12);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	uint8_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_data, data, edit_mask, data);

	if (verbose)
		PrintArray("in (after)", num_data, data);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], data[i]);
	}
}

/*
 * Test bit operation AND by sakura_OperateBitsUint8And with an array of length 11
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000000, 00000010, 00000010,
 *        00000000, 00000001, 00000010 ]
 */
TEST_F(BitOperation8, AndLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_data, data, edit_mask);

	if (verbose) {
		PrintArray("data", num_data, data);
//		PrintArray("edit_mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation AND by sakura_OperateBitsUint8And
 * with an array of length zero
 *
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kOK)
 */
TEST_F(BitOperation8, AndLengthZero) {
	SIMD_ALIGN
	uint8_t data[0];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	size_t const num_in(ELEMENTSOF(data));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_in, data, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Long test of bit operation AND by sakura_OperateBitsUint8And with a large array
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000000, 00000010, 00000010, .... repeated... ]
 */
TEST_F(BitOperation8, AndLong) {
	SIMD_ALIGN
	uint8_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	uint8_t result_long[ELEMENTSOF(data_long)];

	uint8_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(ELEMENTSOF(data_long));

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint8And(bit_mask_, num_large, data_long,
				edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test bit operation AND by sakura_OperateBitsUint32And
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...000, 0...000, 0...010, 0...010 ]
 */
TEST_F(BitOperation32, And) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("data", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation AND by sakura_OperateBitsUint32And in-place operation (&out == &in)
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000000, 00000010, 00000010 ]
 */
TEST_F(BitOperation32, AndInPlace) {
	size_t const num_data(10);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	uint32_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose)
		PrintArray("data (before)", num_data, data);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_data, data, edit_mask, data);

	if (verbose)
		PrintArray("data (after)", num_data, data);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], data[i]);
	}
}

/*
 * Test bit operation AND by sakura_OperateBitsUint32And with an array of length 11
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...000, 0...000, 0...010, 0...010
 *        0...000, 0...001, 0...010 ]
 */
TEST_F(BitOperation32, AndLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_data, data, edit_mask);

	if (verbose) {
		PrintArray("data", num_data, data);
//		PrintArray("edit_mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Long test of bit operation AND by sakura_OperateBitsUint32And with a large array
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...000, 0...000, 0...010, 0...010, .... repeated... ]
 */
TEST_F(BitOperation32, AndLong) {
	SIMD_ALIGN
	uint32_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	uint32_t result_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];

	uint32_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(NUM_IN_LONG);
	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint32And(bit_mask_, num_large, data_long,
				edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test  bit operation AND by sakura_OperateBitsUint32And
 * with an array of length zero
 *
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kOK)
 */
TEST_F(BitOperation32, AndLengthZero) {
	SIMD_ALIGN
	uint32_t data[0];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	size_t const num_data(ELEMENTSOF(data));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_data, data, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test failure cases of sakura_OperateBitsUint8And
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BitOperation8, AndFailNullData) {
	size_t const num_data(NUM_IN);
	uint8_t dummy[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(dummy)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(dummy)];

	uint8_t *data_null = nullptr;
	// assert(data_null == nullptr);

	GetInputDataInLength(num_data, dummy, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_data, data_null, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, AndFailNullMask) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];

	bool dummy[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t out[ELEMENTSOF(data)];

	bool *mask_null = nullptr;

	GetInputDataInLength(num_data, data, dummy);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_data, data, mask_null, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, AndFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];

	uint8_t *result_null = nullptr;

	GetInputDataInLength(num_data, data, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_data, data, edit_mask, result_null);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BitOperation8, AndFailNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint8_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_data, data_shift, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, AndFailNotAlignedMask) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	bool *mask_shift = &edit_mask[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_data, data, mask_shift, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, AndFailNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint8_t *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_data, data, edit_mask, result_shift);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_OperateBitsUint32And
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BitOperation32, AndFailNullData) {
	size_t const num_data(NUM_IN);
	uint32_t dummy[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(dummy)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(dummy)];

	uint32_t *data_null = nullptr;

	GetInputDataInLength(num_data, dummy, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_data, data_null, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, AndFailNullMask) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];

	bool dummy[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t out[ELEMENTSOF(data)];

	bool *mask_null = nullptr;

	GetInputDataInLength(num_data, data, dummy);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_data, data, mask_null, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, AndFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];

	uint32_t *result_null = nullptr;

	GetInputDataInLength(num_data, data, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_data, data, edit_mask, result_null);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BitOperation32, AndFailNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint32_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_data, data_shift, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, AndFailNotAlignedMask) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	bool *mask_shift = &edit_mask[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_data, data, mask_shift, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, AndFailNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint32_t *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_data, data, edit_mask, result_shift);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test Material Nonimplication bit operation using sakura_OperateBitsUint8And
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000001, 00000000, 00000001 ]
 */
TEST_F(BitOperation8, NonImplication) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t answer[] = { 0, 1, 2, 3, 0, 1, 0, 1 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(~bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Long test of Material Nonimplication bit operation of a large array
 * using sakura_OperateBitsUint8And
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000001, 00000000, 00000001, .... repeated... ]
 */
TEST_F(BitOperation8, NonImplicationLong) {
	SIMD_ALIGN
	uint8_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	uint8_t result_long[ELEMENTSOF(data_long)];

	uint8_t answer[] = { 0, 1, 2, 3, 0, 1, 0, 1 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(ELEMENTSOF(data_long));

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint8And(~bit_mask_, num_large, data_long,
				edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test Material Nonimplication bit operation using sakura_OperateBitsUint32And
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...000, 0...001, 0...000, 0...001 ]
 */
TEST_F(BitOperation32, NonImplication) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t answer[] = { 0, 1, 2, 3, 0, 1, 0, 1 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("data", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(~bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Long test of Material Nonimplication bit operation of a large array
 * using sakura_OperateBitsUint32And
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...000, 0...001, 0...000, 0...001, .... repeated... ]
 */
TEST_F(BitOperation32, NonImplicationLong) {
	SIMD_ALIGN
	uint32_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	uint32_t result_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];

	uint32_t answer[] = { 0, 1, 2, 3, 0, 1, 0, 1 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(NUM_IN_LONG);
	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint32And(~bit_mask_, num_large, data_long,
				edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bit operation, Converse Nonimplication,
 * by sakura_OperateBitsUint8ConverseNonImplication
 * RESULT:
 * out = [ 00000000, 00000001, 00000010, 00000011, 00000010, 00000010, 00000000, 00000000 ]
 */
TEST_F(BitOperation8, ConverseNonImplication) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t answer[] = { 0, 1, 2, 3, 2, 2, 0, 0 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint8ConverseNonImplication(bit_mask_, num_data,
					data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation, Converse Nonimplication,
 * by sakura_OperateBitsUint8ConverseNonImplication in-place operation (&out == &in)
 * RESULT:
 * out = [ 00000000, 00000001, 00000010, 00000011, 00000010, 00000010, 00000000, 00000000 ]
 */
TEST_F(BitOperation8, ConverseNonImplicationInPlace) {
	size_t const num_data(12);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	uint8_t answer[] = { 0, 1, 2, 3, 2, 2, 0, 0 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint8ConverseNonImplication(bit_mask_, num_data,
					data, edit_mask, data);

	if (verbose)
		PrintArray("in (after)", num_data, data);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], data[i]);
	}
}

/*
 * Test bit operation, Converse Nonimplication,
 * by sakura_OperateBitsUint8ConverseNonImplication with an array of length 11
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000010, 00000010, 00000000, 00000000,
 *        00000000, 00000001, 00000010 ]
 */
TEST_F(BitOperation8, ConverseNonImplicationLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t answer[] = { 0, 1, 2, 3, 2, 2, 0, 0 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_data, data, edit_mask);

	if (verbose) {
		PrintArray("data", num_data, data);
//		PrintArray("edit_mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint8ConverseNonImplication(bit_mask_, num_data,
					data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation, Converse Nonimplication,
 * by sakura_OperateBitsUint8ConverseNonImplication with an array of length zero
 *
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kOK)
 */
TEST_F(BitOperation8, ConverseNonImplicationLengthZero) {
	SIMD_ALIGN
	uint8_t data[0];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	size_t const num_in(ELEMENTSOF(data));

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint8ConverseNonImplication(bit_mask_, num_in,
					data, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Long test of bit operation, Converse Nonimplication,
 * by sakura_OperateBitsUint8ConverseNonImplication with a large array
 * RESULT:
 * out = [ 00000000, 00000001, 00000010, 00000011, 00000010, 00000010, 00000000, 00000000, .... repeated... ]
 */
TEST_F(BitOperation8, ConverseNonImplicationLong) {
	SIMD_ALIGN
	uint8_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	uint8_t result_long[ELEMENTSOF(data_long)];

	uint8_t answer[] = { 0, 1, 2, 3, 2, 2, 0, 0 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(ELEMENTSOF(data_long));

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint8ConverseNonImplication(bit_mask_,
				num_large, data_long, edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test bit operation, Converse Nonimplication,
 * by sakura_OperateBitsUint32ConverseNonImplication
 * RESULT:
 * out = [ 0...000, 0...001, 0...010, 0...011, 0...010, 0...010, 0...000, 0...000 ]
 */
TEST_F(BitOperation32, ConverseNonImplication) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t answer[] = { 0, 1, 2, 3, 2, 2, 0, 0 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("data", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint32ConverseNonImplication(bit_mask_, num_data,
					data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation, Converse Nonimplication,
 * by sakura_OperateBitsUint32ConverseNonImplication in-place operation (&out == &in)
 * RESULT:
 * out = [ 0...000, 0...001, 0...010, 0...011, 0...010, 0...010, 0...000, 0...000 ]
 */
TEST_F(BitOperation32, ConverseNonImplicationInPlace) {
	size_t const num_data(10);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	uint32_t answer[] = { 0, 1, 2, 3, 2, 2, 0, 0 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose)
		PrintArray("data (before)", num_data, data);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint32ConverseNonImplication(bit_mask_, num_data,
					data, edit_mask, data);

	if (verbose)
		PrintArray("data (after)", num_data, data);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], data[i]);
	}
}

/*
 * Test bit operation, Converse Nonimplication,
 * by sakura_OperateBitsUint32ConverseNonImplication with an array of length 11
 * RESULT:
 * out = [ 0...000, 0...001, 0...010, 0...011, 0...010, 0...010, 0...000, 0...000,
 *        0...000, 0...001, 0...010 ]
 */
TEST_F(BitOperation32, ConverseNonImplicationLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t answer[] = { 0, 1, 2, 3, 2, 2, 0, 0 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_data, data, edit_mask);

	if (verbose) {
		PrintArray("data", num_data, data);
//		PrintArray("edit_mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint32ConverseNonImplication(bit_mask_, num_data,
					data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Long test of bit operation, Converse Nonimplication,
 * by sakura_OperateBitsUint32ConverseNonImplication with a large array
 * RESULT:
 * out = [ 0...000, 0...001, 0...010, 0...011, 0...010, 0...010, 0...000, 0...000, .... repeated... ]
 */
TEST_F(BitOperation32, ConverseNonImplicationLong) {
	SIMD_ALIGN
	uint32_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	uint32_t result_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];

	uint32_t answer[] = { 0, 1, 2, 3, 2, 2, 0, 0 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(NUM_IN_LONG);
	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint32ConverseNonImplication(bit_mask_,
				num_large, data_long, edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test  bit operation, Converse Nonimplication,
 * by sakura_OperateBitsUint32ConverseNonImplication with an array of length zero
 *
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kOK)
 */
TEST_F(BitOperation32, ConverseNonImplicationLengthZero) {
	SIMD_ALIGN
	uint32_t data[0];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	size_t const num_data(ELEMENTSOF(data));

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint32ConverseNonImplication(bit_mask_, num_data,
					data, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test failure cases of sakura_OperateBitsUint8ConverseNonImplication
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BitOperation8, ConverseNonImplicationFailNullData) {
	size_t const num_data(NUM_IN);
	uint8_t dummy[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(dummy)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(dummy)];

	uint8_t *data_null = nullptr;
	// assert(data_null == nullptr);

	GetInputDataInLength(num_data, dummy, edit_mask);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint8ConverseNonImplication(bit_mask_, num_data,
					data_null, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, ConverseNonImplicationFailNullMask) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];

	bool dummy[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t out[ELEMENTSOF(data)];

	bool *mask_null = nullptr;

	GetInputDataInLength(num_data, data, dummy);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint8ConverseNonImplication(bit_mask_, num_data,
					data, mask_null, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, ConverseNonImplicationFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];

	uint8_t *result_null = nullptr;

	GetInputDataInLength(num_data, data, edit_mask);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint8ConverseNonImplication(bit_mask_, num_data,
					data, edit_mask, result_null);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BitOperation8, ConverseNonImplicationFailNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint8_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint8ConverseNonImplication(bit_mask_, num_data,
					data_shift, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, ConverseNonImplicationFailNotAlignedMask) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	bool *mask_shift = &edit_mask[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint8ConverseNonImplication(bit_mask_, num_data,
					data, mask_shift, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, ConverseNonImplicationFailNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint8_t *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint8ConverseNonImplication(bit_mask_, num_data,
					data, edit_mask, result_shift);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_OperateBitsUint32ConverseNonImplication
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BitOperation32, ConverseNonImplicationFailNullData) {
	size_t const num_data(NUM_IN);
	uint32_t dummy[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(dummy)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(dummy)];

	uint32_t *data_null = nullptr;

	GetInputDataInLength(num_data, dummy, edit_mask);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint32ConverseNonImplication(bit_mask_, num_data,
					data_null, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, ConverseNonImplicationFailNullMask) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];

	bool dummy[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t out[ELEMENTSOF(data)];

	bool *mask_null = nullptr;

	GetInputDataInLength(num_data, data, dummy);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint32ConverseNonImplication(bit_mask_, num_data,
					data, mask_null, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, ConverseNonImplicationFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];

	uint32_t *result_null = nullptr;

	GetInputDataInLength(num_data, data, edit_mask);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint32ConverseNonImplication(bit_mask_, num_data,
					data, edit_mask, result_null);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BitOperation32, ConverseNonImplicationFailNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint32_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint32ConverseNonImplication(bit_mask_, num_data,
					data_shift, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, ConverseNonImplicationFailNotAlignedMask) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	bool *mask_shift = &edit_mask[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint32ConverseNonImplication(bit_mask_, num_data,
					data, mask_shift, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, ConverseNonImplicationFailNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint32_t *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint32ConverseNonImplication(bit_mask_, num_data,
					data, edit_mask, result_shift);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bit operation NOR using sakura_OperateBitsUint8ConverseNonImplication
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 11111101, 11111100, 11111101, 11111100 ]
 */
TEST_F(BitOperation8, Nor) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t const base_pattern = (~static_cast<uint8_t>(0) << 2); // 11111100
	uint8_t answer[] = { 0, 1, 2, 3, static_cast<uint8_t>(base_pattern + 1),
			base_pattern, static_cast<uint8_t>(base_pattern + 1), base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint8ConverseNonImplication(~bit_mask_, num_data,
					data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Long test of bit operation NOR with a large array
 * using sakura_OperateBitsUint8ConverseNonImplication
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 11111101, 11111100, 11111101, 11111100, .... repeated... ]
 */
TEST_F(BitOperation8, NorLong) {
	SIMD_ALIGN
	uint8_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	uint8_t result_long[ELEMENTSOF(data_long)];

	uint8_t const base_pattern = (~static_cast<uint8_t>(0) << 2); // 11111100
	uint8_t answer[] = { 0, 1, 2, 3, static_cast<uint8_t>(base_pattern + 1),
			base_pattern, static_cast<uint8_t>(base_pattern + 1), base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(ELEMENTSOF(data_long));

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint8ConverseNonImplication(~bit_mask_,
				num_large, data_long, edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test bit operation NOR using sakura_OperateBitsUint32ConverseNonImplication
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 1...101, 1...100, 1...101, 1...100 ]
 */
TEST_F(BitOperation32, Nor) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t const base_pattern = (~static_cast<uint32_t>(0) << 2); // 11...100
	uint32_t answer[] =
			{ 0, 1, 2, 3, static_cast<uint32_t>(base_pattern + 1), base_pattern,
					static_cast<uint32_t>(base_pattern + 1), base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("data", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateBitsUint32ConverseNonImplication(~bit_mask_, num_data,
					data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Long test of bit operation NOR of a large array
 * using sakura_OperateBitsUint32ConverseNonImplication
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...000, 0...001, 0...000, 0...001, .... repeated... ]
 */
TEST_F(BitOperation32, NorLong) {
	SIMD_ALIGN
	uint32_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	uint32_t result_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];

	uint32_t const base_pattern = (~static_cast<uint32_t>(0) << 2); // 11...100
	uint32_t answer[] =
			{ 0, 1, 2, 3, static_cast<uint32_t>(base_pattern + 1), base_pattern,
					static_cast<uint32_t>(base_pattern + 1), base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(NUM_IN_LONG);
	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint32ConverseNonImplication(~bit_mask_,
				num_large, data_long, edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bit operation OR by sakura_OperateBitsUint8Or
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000010, 00000011, 00000010, 00000011 ]
 */
TEST_F(BitOperation8, Or) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t answer[] = { 0, 1, 2, 3, 2, 3, 2, 3 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Or(bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation OR by sakura_OperateBitsUint8Or in-place operation (&out == &in)
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000010, 00000011, 00000010, 00000011 ]
 */
TEST_F(BitOperation8, OrInPlace) {
	size_t const num_data(12);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	uint8_t answer[] = { 0, 1, 2, 3, 2, 3, 2, 3 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Or(bit_mask_,
			num_data, data, edit_mask, data);

	if (verbose)
		PrintArray("in (after)", num_data, data);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], data[i]);
	}
}

/*
 * Test bit operation OR by sakura_OperateBitsUint8Or with an array of length 11
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000010, 00000011, 00000010, 00000011,
 *        00000000, 00000001, 00000010 ]
 */
TEST_F(BitOperation8, OrLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t answer[] = { 0, 1, 2, 3, 2, 3, 2, 3 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_data, data, edit_mask);

	if (verbose) {
		PrintArray("data", num_data, data);
//		PrintArray("edit_mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Or(bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation OR by sakura_OperateBitsUint8Or
 * with an array of length zero
 *
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kOK)
 */
TEST_F(BitOperation8, OrLengthZero) {
	SIMD_ALIGN
	uint8_t data[0];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	size_t const num_in(ELEMENTSOF(data));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Or(bit_mask_,
			num_in, data, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Long test of bit operation OR by sakura_OperateBitsUint8Or with a large array
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000010, 00000011, 00000010, 00000011, .... repeated... ]
 */
TEST_F(BitOperation8, OrLong) {
	SIMD_ALIGN
	uint8_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	uint8_t result_long[ELEMENTSOF(data_long)];

	uint8_t answer[] = { 0, 1, 2, 3, 2, 3, 2, 3 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(ELEMENTSOF(data_long));

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint8Or(bit_mask_, num_large, data_long,
				edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test bit operation OR by sakura_OperateBitsUint32Or
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...010, 0...011, 0...010, 0...011 ]
 */
TEST_F(BitOperation32, Or) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t answer[] = { 0, 1, 2, 3, 2, 3, 2, 3 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("data", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Or(bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation OR by sakura_OperateBitsUint32Or in-place operation (&out == &in)
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...010, 0...011, 0...010, 0...011 ]
 */
TEST_F(BitOperation32, OrInPlace) {
	size_t const num_data(10);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	uint32_t answer[] = { 0, 1, 2, 3, 2, 3, 2, 3 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose)
		PrintArray("data (before)", num_data, data);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Or(bit_mask_,
			num_data, data, edit_mask, data);

	if (verbose)
		PrintArray("data (after)", num_data, data);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], data[i]);
	}
}

/*
 * Test bit operation OR by sakura_OperateBitsUint32Or with an array of length 11
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...010, 0...011, 0...010, 0...011,
 *        0...000, 0...001, 0...010 ]
 */
TEST_F(BitOperation32, OrLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t answer[] = { 0, 1, 2, 3, 2, 3, 2, 3 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_data, data, edit_mask);

	if (verbose) {
		PrintArray("data", num_data, data);
//		PrintArray("edit_mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Or(bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Long test of bit operation OR by sakura_OperateBitsUint32Or with a large array
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...010, 0...011, 0...010, 0...011, .... repeated... ]
 */
TEST_F(BitOperation32, OrLong) {
	SIMD_ALIGN
	uint32_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	uint32_t result_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];

	uint32_t answer[] = { 0, 1, 2, 3, 2, 3, 2, 3 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(NUM_IN_LONG);
	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint32Or(bit_mask_, num_large, data_long,
				edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test  bit operation OR by sakura_OperateBitsUint32Or
 * with an array of length zero
 *
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kOK)
 */
TEST_F(BitOperation32, OrLengthZero) {
	SIMD_ALIGN
	uint32_t data[0];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	size_t const num_data(ELEMENTSOF(data));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Or(bit_mask_,
			num_data, data, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test failure cases of sakura_OperateBitsUint8Or
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BitOperation8, OrFailNullData) {
	size_t const num_data(NUM_IN);
	uint8_t dummy[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(dummy)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(dummy)];

	uint8_t *data_null = nullptr;
	// assert(data_null == nullptr);

	GetInputDataInLength(num_data, dummy, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Or(bit_mask_,
			num_data, data_null, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, OrFailNullMask) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];

	bool dummy[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t out[ELEMENTSOF(data)];

	bool *mask_null = nullptr;

	GetInputDataInLength(num_data, data, dummy);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Or(bit_mask_,
			num_data, data, mask_null, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, OrFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];

	uint8_t *result_null = nullptr;

	GetInputDataInLength(num_data, data, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Or(bit_mask_,
			num_data, data, edit_mask, result_null);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BitOperation8, OrFailNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint8_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Or(bit_mask_,
			num_data, data_shift, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, OrFailNotAlignedMask) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	bool *mask_shift = &edit_mask[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Or(bit_mask_,
			num_data, data, mask_shift, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, OrFailNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint8_t *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Or(bit_mask_,
			num_data, data, edit_mask, result_shift);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_OperateBitsUint32Or
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BitOperation32, OrFailNullData) {
	size_t const num_data(NUM_IN);
	uint32_t dummy[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(dummy)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(dummy)];

	uint32_t *data_null = nullptr;

	GetInputDataInLength(num_data, dummy, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Or(bit_mask_,
			num_data, data_null, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, OrFailNullMask) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];

	bool dummy[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t out[ELEMENTSOF(data)];

	bool *mask_null = nullptr;

	GetInputDataInLength(num_data, data, dummy);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Or(bit_mask_,
			num_data, data, mask_null, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, OrFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];

	uint32_t *result_null = nullptr;

	GetInputDataInLength(num_data, data, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Or(bit_mask_,
			num_data, data, edit_mask, result_null);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BitOperation32, OrFailNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint32_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Or(bit_mask_,
			num_data, data_shift, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, OrFailNotAlignedMask) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	bool *mask_shift = &edit_mask[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Or(bit_mask_,
			num_data, data, mask_shift, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, OrFailNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint32_t *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Or(bit_mask_,
			num_data, data, edit_mask, result_shift);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test Converse Implication bit operation using sakura_OperateBitsUint8Or
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 11111101, 11111101, 11111111, 11111111 ]
 */
TEST_F(BitOperation8, ConverseImplication) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t const base_pattern = (~static_cast<uint8_t>(0)); // 11111111
	uint8_t answer[] = { 0, 1, 2, 3, static_cast<uint8_t>(base_pattern - 2),
			static_cast<uint8_t>(base_pattern - 2), base_pattern, base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Or(~bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Long test of Converse Implication bit operation of a large array
 * using sakura_OperateBitsUint8Or
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 11111101, 11111101, 11111111, 11111111, .... repeated... ]
 */
TEST_F(BitOperation8, ConverseImplicationLong) {
	SIMD_ALIGN
	uint8_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	uint8_t result_long[ELEMENTSOF(data_long)];

	uint8_t const base_pattern = (~static_cast<uint8_t>(0)); // 11111111
	uint8_t answer[] = { 0, 1, 2, 3, static_cast<uint8_t>(base_pattern - 2),
			static_cast<uint8_t>(base_pattern - 2), base_pattern, base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(ELEMENTSOF(data_long));

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint8Or(~bit_mask_, num_large, data_long,
				edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test Converse Implication bit operation using sakura_OperateBitsUint32Or
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 1...101, 1...101, 1...111, 1...111 ]
 */
TEST_F(BitOperation32, ConverseImplication) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t const base_pattern = (~static_cast<uint32_t>(0)); // 1...111
	uint32_t answer[] =
			{ 0, 1, 2, 3, static_cast<uint32_t>(base_pattern - 2),
					static_cast<uint32_t>(base_pattern - 2), base_pattern,
					base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("data", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Or(~bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Long test of Converse Implication bit operation of a large array
 * using sakura_OperateBitsUint32Or
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 1...101, 1...101, 1...111, 1...111, .... repeated... ]
 */
TEST_F(BitOperation32, ConverseImplicationLong) {
	SIMD_ALIGN
	uint32_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	uint32_t result_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];

	uint32_t const base_pattern = (~static_cast<uint32_t>(0)); // 1...111
	uint32_t answer[] =
			{ 0, 1, 2, 3, static_cast<uint32_t>(base_pattern - 2),
					static_cast<uint32_t>(base_pattern - 2), base_pattern,
					base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(NUM_IN_LONG);
	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint32Or(~bit_mask_, num_large, data_long,
				edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bit operation XOR by sakura_OperateBitsUint8Xor
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000010, 00000011, 00000000, 00000001 ]
 */
TEST_F(BitOperation8, Xor) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t answer[] = { 0, 1, 2, 3, 2, 3, 0, 1 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Xor(bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation XOR by sakura_OperateBitsUint8Xor in-place operation (&out == &in)
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000010, 00000011, 00000000, 00000001 ]
 */
TEST_F(BitOperation8, XorInPlace) {
	size_t const num_data(12);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	uint8_t answer[] = { 0, 1, 2, 3, 2, 3, 0, 1 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Xor(bit_mask_,
			num_data, data, edit_mask, data);

	if (verbose)
		PrintArray("in (after)", num_data, data);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], data[i]);
	}
}

/*
 * Test bit operation XOR by sakura_OperateBitsUint8Xor with an array of length 11
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000010, 00000011, 00000000, 00000001,
 *        00000000, 00000001, 00000010 ]
 */
TEST_F(BitOperation8, XorLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t answer[] = { 0, 1, 2, 3, 2, 3, 0, 1 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_data, data, edit_mask);

	if (verbose) {
		PrintArray("data", num_data, data);
//		PrintArray("edit_mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Xor(bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation XOR by sakura_OperateBitsUint8Xor
 * with an array of length zero
 *
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kOK)
 */
TEST_F(BitOperation8, XorLengthZero) {
	SIMD_ALIGN
	uint8_t data[0];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	size_t const num_in(ELEMENTSOF(data));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Xor(bit_mask_,
			num_in, data, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Long test of bit operation XOR by sakura_OperateBitsUint8Xor with a large array
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000010, 00000011, 00000000, 00000001, .... repeated... ]
 */
TEST_F(BitOperation8, XorLong) {
	SIMD_ALIGN
	uint8_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	uint8_t result_long[ELEMENTSOF(data_long)];

	uint8_t answer[] = { 0, 1, 2, 3, 2, 3, 0, 1 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(ELEMENTSOF(data_long));

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint8Xor(bit_mask_, num_large, data_long,
				edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test bit operation XOR by sakura_OperateBitsUint32Xor
 * RESULT:
 * out = [ 0...000, 0...001, 0...010, 0...011, 0...010, 0...011, 0...000, 0...001 ]
 */
TEST_F(BitOperation32, Xor) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t answer[] = { 0, 1, 2, 3, 2, 3, 0, 1 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("data", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Xor(bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation XOR by sakura_OperateBitsUint32Xor in-place operation (&out == &in)
 * RESULT:
 * out = [ 0...000, 0...001, 0...010, 0...011, 0...010, 0...011, 0...000, 0...001 ]
 */
TEST_F(BitOperation32, XorInPlace) {
	size_t const num_data(10);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	uint32_t answer[] = { 0, 1, 2, 3, 2, 3, 0, 1 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose)
		PrintArray("data (before)", num_data, data);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Xor(bit_mask_,
			num_data, data, edit_mask, data);

	if (verbose)
		PrintArray("data (after)", num_data, data);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], data[i]);
	}
}

/*
 * Test bit operation XOR by sakura_OperateBitsUint32Xor with an array of length 11
 * RESULT:
 * out = [ 0...000, 0...001, 0...010, 0...011, 0...010, 0...011, 0...000, 0...001,
 *         0...000, 0...001, 0...010 ]
 */
TEST_F(BitOperation32, XorLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t answer[] = { 0, 1, 2, 3, 2, 3, 0, 1 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_data, data, edit_mask);

	if (verbose) {
		PrintArray("data", num_data, data);
//		PrintArray("edit_mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Xor(bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Long test of bit operation XOR by sakura_OperateBitsUint32Xor with a large array
 * RESULT:
 * out = [ 0...000, 0...001, 0...010, 0...011, 0...010, 0...011, 0...000, 0...001, .... repeated... ]
 */
TEST_F(BitOperation32, XorLong) {
	SIMD_ALIGN
	uint32_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	uint32_t result_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];

	uint32_t answer[] = { 0, 1, 2, 3, 2, 3, 0, 1 };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(NUM_IN_LONG);
	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint32Xor(bit_mask_, num_large, data_long,
				edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test  bit operation XOR by sakura_OperateBitsUint32Xor
 * with an array of length zero
 *
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kOK)
 */
TEST_F(BitOperation32, XorLengthZero) {
	SIMD_ALIGN
	uint32_t data[0];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	size_t const num_data(ELEMENTSOF(data));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Xor(bit_mask_,
			num_data, data, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test failure cases of sakura_OperateBitsUint8Xor
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BitOperation8, XorFailNullData) {
	size_t const num_data(NUM_IN);
	uint8_t dummy[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(dummy)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(dummy)];

	uint8_t *data_null = nullptr;
	// assert(data_null == nullptr);

	GetInputDataInLength(num_data, dummy, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Xor(bit_mask_,
			num_data, data_null, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, XorFailNullMask) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];

	bool dummy[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t out[ELEMENTSOF(data)];

	bool *mask_null = nullptr;

	GetInputDataInLength(num_data, data, dummy);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Xor(bit_mask_,
			num_data, data, mask_null, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, XorFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];

	uint8_t *result_null = nullptr;

	GetInputDataInLength(num_data, data, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Xor(bit_mask_,
			num_data, data, edit_mask, result_null);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BitOperation8, XorFailNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint8_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Xor(bit_mask_,
			num_data, data_shift, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, XorFailNotAlignedMask) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	bool *mask_shift = &edit_mask[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Xor(bit_mask_,
			num_data, data, mask_shift, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, XorFailNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint8_t *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Xor(bit_mask_,
			num_data, data, edit_mask, result_shift);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_OperateBitsUint32Xor
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BitOperation32, XorFailNullData) {
	size_t const num_data(NUM_IN);
	uint32_t dummy[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(dummy)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(dummy)];

	uint32_t *data_null = nullptr;

	GetInputDataInLength(num_data, dummy, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Xor(bit_mask_,
			num_data, data_null, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, XorFailNullMask) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];

	bool dummy[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t out[ELEMENTSOF(data)];

	bool *mask_null = nullptr;

	GetInputDataInLength(num_data, data, dummy);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Xor(bit_mask_,
			num_data, data, mask_null, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, XorFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];

	uint32_t *result_null = nullptr;

	GetInputDataInLength(num_data, data, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Xor(bit_mask_,
			num_data, data, edit_mask, result_null);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BitOperation32, XorFailNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint32_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Xor(bit_mask_,
			num_data, data_shift, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, XorFailNotAlignedMask) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	bool *mask_shift = &edit_mask[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Xor(bit_mask_,
			num_data, data, mask_shift, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, XorFailNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint32_t *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Xor(bit_mask_,
			num_data, data, edit_mask, result_shift);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bit operation XNOR using sakura_OperateBitsUint8Xor
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 11111101, 11111100, 11111111, 11111110 ]
 */
TEST_F(BitOperation8, Xnor) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t const base_pattern = (~static_cast<uint8_t>(0)); // 11111111
	uint8_t answer[] = { 0, 1, 2, 3, static_cast<uint8_t>(base_pattern - 2),
			static_cast<uint8_t>(base_pattern << 2), base_pattern,
			static_cast<uint8_t>(base_pattern << 1) };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Xor(~bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation XNOR using sakura_OperateBitsUint32Xor
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 1...101, 1...100, 1...111, 1...110 ]
 */
TEST_F(BitOperation32, Xnor) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t const base_pattern = (~static_cast<uint32_t>(0)); // 1...111
	uint32_t answer[] = { 0, 1, 2, 3, static_cast<uint32_t>(base_pattern - 2),
			static_cast<uint32_t>(base_pattern << 2), base_pattern,
			static_cast<uint32_t>(base_pattern << 1) };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("data", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Xor(~bit_mask_,
			num_data, data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bit operation NOT by sakura_OperateBitsUint8Not
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 11111111, 11111110, 11111101, 11111100 ]
 */
TEST_F(BitOperation8, Not) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t const base_pattern = (~static_cast<uint8_t>(0) << 2); // 11111100
	uint8_t answer[] = { 0, 1, 2, 3, static_cast<uint8_t>(base_pattern + 3),
			static_cast<uint8_t>(base_pattern + 2),
			static_cast<uint8_t>(base_pattern + 1), base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Not(num_data, data,
			edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation NOT by sakura_OperateBitsUint8Not in-place operation (&out == &in)
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 11111111, 11111110, 11111101, 11111100 ]
 */
TEST_F(BitOperation8, NotInPlace) {
	size_t const num_data(12);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	uint8_t const base_pattern = (~static_cast<uint8_t>(0) << 2); // 11111100
	uint8_t answer[] = { 0, 1, 2, 3, static_cast<uint8_t>(base_pattern + 3),
			static_cast<uint8_t>(base_pattern + 2),
			static_cast<uint8_t>(base_pattern + 1), base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Not(num_data, data,
			edit_mask, data);

	if (verbose)
		PrintArray("in (after)", num_data, data);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], data[i]);
	}
}

/*
 * Test bit operation NOT by sakura_OperateBitsUint8Not with an array of length 11
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 11111111, 11111110, 11111101, 11111100,
 *        00000000, 00000001, 00000010 ]
 */
TEST_F(BitOperation8, NotLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	uint8_t const base_pattern = (~static_cast<uint8_t>(0) << 2); // 11111100
	uint8_t answer[] = { 0, 1, 2, 3, static_cast<uint8_t>(base_pattern + 3),
			static_cast<uint8_t>(base_pattern + 2),
			static_cast<uint8_t>(base_pattern + 1), base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_data, data, edit_mask);

	if (verbose) {
		PrintArray("data", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Not(num_data, data,
			edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Long test of bit operation NOT by sakura_OperateBitsUint8Not with a large array
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 11111111, 11111110, 11111101, 11111100, .... repeated... ]
 */
TEST_F(BitOperation8, NotLong) {
	SIMD_ALIGN
	uint8_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	uint8_t result_long[ELEMENTSOF(data_long)];

	uint8_t const base_pattern = (~static_cast<uint8_t>(0) << 2); // 11111100
	uint8_t answer[] = { 0, 1, 2, 3, static_cast<uint8_t>(base_pattern + 3),
			static_cast<uint8_t>(base_pattern + 2),
			static_cast<uint8_t>(base_pattern + 1), base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(ELEMENTSOF(data_long));

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint8Not(num_large, data_long,
				edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test bit operation NOT by sakura_OperateBitsUint8Not
 * with an array of length zero
 *
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kOK)
 */
TEST_F(BitOperation8, NotLengthZero) {
	SIMD_ALIGN
	uint8_t data[0];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	size_t const num_in(ELEMENTSOF(data));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Not(num_in, data,
			edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bit operation NOT by sakura_OperateBitsUint32Not
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 1...111, 1...110, 1...101, 1...100 ]
 */
TEST_F(BitOperation32, Not) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t const base_pattern = (~static_cast<uint32_t>(0) << 2); // 1...100
	uint32_t answer[] = { 0, 1, 2, 3, static_cast<uint32_t>(base_pattern + 3),
			static_cast<uint32_t>(base_pattern + 2),
			static_cast<uint32_t>(base_pattern + 1), base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Not(num_data,
			data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bit operation NOT by sakura_OperateBitsUint32Not in-place operation (&out == &in)
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 1...111, 1...110, 1...101, 1...100 ]
 */
TEST_F(BitOperation32, NotInPlace) {
	size_t const num_data(12);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	uint32_t const base_pattern = (~static_cast<uint32_t>(0) << 2); // 1...100
	uint32_t answer[] = { 0, 1, 2, 3, static_cast<uint32_t>(base_pattern + 3),
			static_cast<uint32_t>(base_pattern + 2),
			static_cast<uint32_t>(base_pattern + 1), base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	GetInputDataInLength(num_data, data, edit_mask);
	if (verbose) {
		PrintArray("in (before)", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Not(num_data,
			data, edit_mask, data);

	if (verbose)
		PrintArray("in (after)", num_data, data);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], data[i]);
	}
}

/*
 * Test bit operation NOT by sakura_OperateBitsUint32Not with an array of length 11
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 1...111, 1...110, 1...101, 1...100,
 *        0...000, 0...001, 0...010 ]
 */
TEST_F(BitOperation32, NotLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	uint32_t const base_pattern = (~static_cast<uint32_t>(0) << 2); // 1...100
	uint32_t answer[] = { 0, 1, 2, 3, static_cast<uint32_t>(base_pattern + 3),
			static_cast<uint32_t>(base_pattern + 2),
			static_cast<uint32_t>(base_pattern + 1), base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_data, data, edit_mask);

	if (verbose) {
		PrintArray("data", num_data, data);
		PrintArray("mask", num_data, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Not(num_data,
			data, edit_mask, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Long test of bit operation NOT by sakura_OperateBitsUint32Not with a large array
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 1...111, 1...110, 1...101, 1...100, .... repeated... ]
 */
TEST_F(BitOperation32, NotLong) {
	SIMD_ALIGN
	uint32_t data_long[NUM_IN_LONG];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(data_long)];
	SIMD_ALIGN
	uint32_t result_long[ELEMENTSOF(data_long)];

	uint32_t const base_pattern = (~static_cast<uint32_t>(0) << 2); // 1...100
	uint32_t answer[] = { 0, 1, 2, 3, static_cast<uint32_t>(base_pattern + 3),
			static_cast<uint32_t>(base_pattern + 2),
			static_cast<uint32_t>(base_pattern + 1), base_pattern };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	size_t const num_large(ELEMENTSOF(data_long));

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, data_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint32Not(num_large, data_long,
				edit_mask_long, result_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result_long[i]);
	}
}

/*
 * Test bit operation NOT by sakura_OperateBitsUint32Not
 * with an array of length zero
 *
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kOK)
 */
TEST_F(BitOperation32, NotLengthZero) {
	SIMD_ALIGN
	uint32_t data[0];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	size_t const num_in(ELEMENTSOF(data));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Not(num_in, data,
			edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test failure cases of sakura_OperateBitsUint8Not
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BitOperation8, NotFailNullData) {
	size_t const num_data(NUM_IN);
	uint8_t dummy[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(dummy)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(dummy)];

	uint8_t *data_null = nullptr;
	// assert(data_null == nullptr);

	GetInputDataInLength(num_data, dummy, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Not(num_data,
			data_null, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, NotFailNullMask) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];

	bool dummy[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t out[ELEMENTSOF(data)];

	bool *mask_null = nullptr;

	GetInputDataInLength(num_data, data, dummy);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Not(num_data, data,
			mask_null, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, NotFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];

	uint8_t *result_null = nullptr;

	GetInputDataInLength(num_data, data, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Not(num_data, data,
			edit_mask, result_null);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BitOperation8, NotFailNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint8_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Not(num_data,
			data_shift, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, NotFailNotAlignedMask) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	bool *mask_shift = &edit_mask[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Not(num_data, data,
			mask_shift, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation8, NotFailNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint8_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint8_t *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8Not(num_data, data,
			edit_mask, result_shift);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_OperateBitsUint32Not
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BitOperation32, NotFailNullData) {
	size_t const num_data(NUM_IN);
	uint32_t dummy[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(dummy)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(dummy)];

	uint32_t *data_null = nullptr;

	GetInputDataInLength(num_data, dummy, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Not(num_data,
			data_null, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, NotFailNullMask) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];

	bool dummy[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t out[ELEMENTSOF(data)];

	bool *mask_null = nullptr;

	GetInputDataInLength(num_data, data, dummy);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Not(num_data,
			data, mask_null, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, NotFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];

	uint32_t *result_null = nullptr;

	GetInputDataInLength(num_data, data, edit_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Not(num_data,
			data, edit_mask, result_null);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BitOperation32, NotFailNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint32_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Not(num_data,
			data_shift, edit_mask, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, NotFailNotAlignedMask) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	bool *mask_shift = &edit_mask[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(mask_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Not(num_data,
			data, mask_shift, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BitOperation32, NotFailNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	uint32_t result[ELEMENTSOF(data)];

	GetInputDataInLength(num_elements, data, edit_mask);

	// Define unaligned array
	uint32_t *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32Not(num_data,
			data, edit_mask, result_shift);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

