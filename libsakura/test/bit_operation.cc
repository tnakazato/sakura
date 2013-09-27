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

