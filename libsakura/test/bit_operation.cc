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
			in_[i] = i % ntype; /* repeat bit pattern of *00, *01, *10, *11,... */
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
			data[i] = in_[i % num_in];
			mask[i] = edit_mask_[i % num_in];
		}
	}

	/* Converts an input value to a bit pattern.*/
//	char* BToS(DataType in_value) {
	string BToS(DataType in_value) {
		char buff[bit_size + 1];
		buff[bit_size] = '\0';
		for (size_t i = 0; i < bit_size; ++i) {
			if ((in_value >> i) % 2 == 0)
				buff[bit_size - 1 - i] = '0';
			else
				buff[bit_size - 1 - i] = '1';
		}
		return string(buff);
	}

	/* Converts an bit pattern (char) to a value of DataType.*/
	DataType SToB(char* in_bit) {
		DataType result(0);
		size_t i = 0;
		while (in_bit[i] != '\0') {
			//result = result * 2 + ((uint8_t) in_string[i]);
			result <<= 1;
			if (in_bit[i] == '1')
				result += 1;
			++i;
		}
		return result;
	}

	void PrintInputs() {
		cout << "bit_mask = " << BToS(bit_mask_);
		cout << endl;
		PrintArray("in", NUM_IN, in_);
	}

	void PrintArray(char const *name, size_t num_in, DataType *in) {
		cout << name << " = [";
		for (size_t i = 0; i < num_in - 1; ++i)
			cout << BToS(in[i]) << ", ";
		cout << BToS(in[num_in - 1]);
		cout << " ]" << endl;
	}

	void PrintArray(char const *name, size_t num_in, bool const *in) {
		cout << name << " = [ ";
		for (size_t i = 0; i < num_in - 1; ++i)
			cout << (in[i] ? "T" : "F") << ", ";
		cout << (in[num_in - 1] ? "T" : "F");
		cout << " ]" << endl;
	}

	DataType in_[NUM_IN];
	bool edit_mask_[NUM_IN];

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
	size_t const num_in(NUM_IN);
	SIMD_ALIGN
	uint8_t in_data[num_in];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	uint8_t out[ELEMENTSOF(in_data)];

	uint8_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	ASSERT_EQ(NUM_IN, ELEMENTSOF(answer));

	GetInputDataInLength(num_in, in_data, edit_mask);
	if (verbose){
		PrintArray("in (before)", num_in, in_data);
		PrintArray("mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_in, in_data, edit_mask, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], out[i]);
	}
}

/*
 * Test bit operation AND by sakura_OperateBitsUint8And in-place operation (&out == &in)
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000000, 00000010, 00000010 ]
 */
TEST_F(BitOperation8, AndInPlace) {
	size_t const num_in(12);
	SIMD_ALIGN
	uint8_t in_data[num_in];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(in_data)];
	uint8_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	ASSERT_EQ(NUM_IN, ELEMENTSOF(answer));

	GetInputDataInLength(num_in, in_data, edit_mask);
	if (verbose){
		PrintArray("in (before)", num_in, in_data);
		PrintArray("mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_in, in_data, edit_mask, in_data);

	if (verbose)
		PrintArray("in (after)", num_in, in_data);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], in_data[i]);
	}
}

/*
 * Test bit operation AND by sakura_OperateBitsUint8And with an array of length 11
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000000, 00000010, 00000010,
 *        00000000, 00000001, 00000010 ]
 */
TEST_F(BitOperation8, AndLengthEleven) {
	size_t const num_in(11);
	SIMD_ALIGN
	uint8_t in[num_in];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(in)];
	SIMD_ALIGN
	uint8_t out[ELEMENTSOF(in)];

	uint8_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	ASSERT_EQ(NUM_IN, ELEMENTSOF(answer));

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_in, in, edit_mask);

	if (verbose){
		PrintArray("in", num_in, in);
//		PrintArray("edit_mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_in, in, edit_mask, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], out[i]);
	}
}

/*
 * Test bit operation AND by sakura_OperateBitsUint32And
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...000, 0...000, 0...010, 0...010 ]
 */
TEST_F(BitOperation32, And) {
	size_t const num_in(NUM_IN);
	SIMD_ALIGN
	uint32_t in_data[num_in];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	uint32_t out[ELEMENTSOF(in_data)];

	uint32_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	ASSERT_EQ(NUM_IN, ELEMENTSOF(answer));

	GetInputDataInLength(num_in, in_data, edit_mask);
	if (verbose){
		PrintArray("in (before)", num_in, in_data);
		PrintArray("mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_in, in_data, edit_mask, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], out[i]);
	}
}

/*
 * Test bit operation AND by sakura_OperateBitsUint32And in-place operation (&out == &in)
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000000, 00000010, 00000010 ]
 */
TEST_F(BitOperation32, AndInPlace) {
	size_t const num_in(10);
	SIMD_ALIGN
	uint32_t in[num_in];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(in)];
	uint32_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	ASSERT_EQ(NUM_IN, ELEMENTSOF(answer));

	GetInputDataInLength(num_in, in, edit_mask);
	if (verbose)
		PrintArray("in (before)", num_in, in);

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_in, in, edit_mask, in);

	if (verbose)
		PrintArray("in (after)", num_in, in);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], in[i]);
	}
}

/*
 * Test bit operation AND by sakura_OperateBitsUint32And with an array of length 11
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...000, 0...000, 0...010, 0...010
 *        0...000, 0...001, 0...010 ]
 */
TEST_F(BitOperation32, AndLengthEleven) {
	size_t const num_in(11);
	SIMD_ALIGN
	uint32_t in[num_in];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(in)];
	SIMD_ALIGN
	uint32_t out[ELEMENTSOF(in)];

	uint32_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	ASSERT_EQ(NUM_IN, ELEMENTSOF(answer));

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_in, in, edit_mask);

	if (verbose){
		PrintArray("in", num_in, in);
//		PrintArray("edit_mask", num_in, edit_mask);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_in, in, edit_mask, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], out[i]);
	}
}

/*
 * Long test of bit operation AND by sakura_OperateBitsUint8And with a large array
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000000, 00000010, 00000010, .... repeated... ]
 */
TEST_F(BitOperation8, AndLong) {
	SIMD_ALIGN
	uint8_t in_long[NUM_IN_LONG];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(in_long)];
	SIMD_ALIGN
	uint8_t out_long[ELEMENTSOF(in_long)];

	uint8_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	ASSERT_EQ(NUM_IN, ELEMENTSOF(answer));

	size_t const num_large(ELEMENTSOF(in_long));

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, in_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint8And(bit_mask_, num_large, in_long,
				edit_mask_long, out_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], out_long[i]);
	}
}

/*
 * Long test of bit operation AND by sakura_OperateBitsUint32And with a large array
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...000, 0...000, 0...010, 0...010, .... repeated... ]
 */
TEST_F(BitOperation32, AndLong) {
	SIMD_ALIGN
	uint32_t in_long[NUM_IN_LONG];
	SIMD_ALIGN
	uint32_t out_long[ELEMENTSOF(in_long)];
	SIMD_ALIGN
	bool edit_mask_long[ELEMENTSOF(in_long)];

	uint32_t answer[] = { 0, 1, 2, 3, 0, 0, 2, 2 };
	ASSERT_EQ(NUM_IN, ELEMENTSOF(answer));

	size_t const num_large(NUM_IN_LONG);
	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	// Create long input data by repeating in_data and edit_mask_
	GetInputDataInLength(num_large, in_long, edit_mask_long);

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_large << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_OperateBitsUint32And(bit_mask_, num_large, in_long,
				edit_mask_long, out_long);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_large; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], out_long[i]);
	}
}

/*
 * Test failure cases of sakura_OperateBitsUint8And
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Input array is zero length */
TEST_F(BitOperation8, AndLengthZero) {
	SIMD_ALIGN
	uint8_t in[0];
	SIMD_ALIGN
	uint8_t out[ELEMENTSOF(in)];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(in)];
	size_t const num_in(ELEMENTSOF(in));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_,
			num_in, in, edit_mask, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_OperateBitsUint32And
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Input array is zero length */
TEST_F(BitOperation32, AndLengthZero) {
	SIMD_ALIGN
	uint32_t in[0];
	SIMD_ALIGN
	uint32_t out[ELEMENTSOF(in)];
	SIMD_ALIGN
	bool edit_mask[ELEMENTSOF(in)];
	size_t const num_in(ELEMENTSOF(in));

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_,
			num_in, in, edit_mask, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

