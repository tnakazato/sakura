#include <iostream>
#include <string>
#include <sys/time.h>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "aligned_memory.h"
#include "gtest/gtest.h"

/* the number of elements in input/output array to test */
#define NUM_IN 8
#define NUM_RANGE 2 // DO NOT MODIFY THIS!
#define NUM_IN_LONG (1 << 18) //2**18
#define UNALIGN_OFFSET 1 // should be != ALIGNMENT
using namespace std;

/*
 * A super class to test various bit operation of an value and array
 */
template<typename DataType>
class BoolFilter: public ::testing::Test {
protected:

	BoolFilter() :
			verbose(false) {
	}

	virtual void SetUp() {
		// Initialize sakura
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Initialize)(nullptr,
				nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		PrepareInputs();
	}

	virtual void TearDown() {
		// Clean-up sakura
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	virtual void PrepareInputs() = 0;

	// Create arbitrary length of input data by repeating values of data_[]
	void GetDataInLength(size_t num_out, DataType *out_data) {
		size_t const num_data(ELEMENTSOF(data_));
		for (size_t i = 0; i < num_out; ++i) {
			out_data[i] = data_[i % num_data];
		}
	}

	void GetBounds(DataType *lower, DataType *upper) {
		EXPECT_EQ(NUM_RANGE, ELEMENTSOF(lower));
		EXPECT_EQ(NUM_RANGE, ELEMENTSOF(upper));
		for (size_t i = 0; i < NUM_RANGE; ++i) {
			lower[i] = lower_[i];
			upper[i] = upper_[i];
		}
	}

	void PrintArray(char const *name, size_t num_data,
			DataType const *data_array) {
		cout << name << " = [ ";
		for (size_t i = 0; i < num_data - 1; ++i)
			cout << data_array[i] << ", ";
		cout << data_array[num_data - 1];
		cout << " ]" << endl;
	}

	void PrintArray(char const *name, size_t num_data, bool const *data_array) {
		cout << name << " = [ ";
		for (size_t i = 0; i < num_data - 1; ++i)
			cout << (data_array[i] ? "T" : "F") << ", ";
		cout << (data_array[num_data - 1] ? "T" : "F");
		cout << " ]" << endl;
	}

	bool verbose;
	DataType data_[NUM_IN];
	DataType upper_[NUM_RANGE];
	DataType lower_[NUM_RANGE];
};

/*
 * Tests various bool filter generation using float array
 * INPUTS:
 * - data = [0., -0.5, -1.0, -0.5, 0., 0.5, 1.0, 0.5]
 * - lower_bound = [-0.75, 0.25]
 * - upper_bound = [-0.25, 0.75]
 */
class BoolFilterFloat: public BoolFilter<float> {
protected:
	virtual void PrepareInputs() {
		float const base_input[] = { 0., -0.5, -1.0, -0.5, 0., 0.5, 1.0, 0.5 };
		STATIC_ASSERT(ELEMENTSOF(base_input) == NUM_IN);
		STATIC_ASSERT(NUM_RANGE == 2);
		lower_[0] = -0.75;
		lower_[1] = 0.25;
		upper_[0] = -0.25;
		upper_[1] = 0.75;
		for (size_t i = 0; i < NUM_IN; ++i) {
			data_[i] = base_input[i];
		}
	}
//	SIMD_ALIGN
//	float upper_[NUM_RANGE];
//	SIMD_ALIGN
//	float lower_[NUM_RANGE];
//
};

/*
 * Tests various bool filter generation using int array
 * INPUTS:
 * - data = [0, -5, -10, -5, 0, 5, 10, 5]
 * - lower_bound = [-7, 3]
 * - upper_bound = [-3, 7]
 */
class BoolFilterInt: public BoolFilter<int> {
protected:
	virtual void PrepareInputs() {
		int const base_input[] = { 0, -5, -10, -5, 0, 5, 10, 5 };
		STATIC_ASSERT(ELEMENTSOF(base_input) == NUM_IN);
		STATIC_ASSERT(NUM_RANGE == 2);
		lower_[0] = -7;
		lower_[1] = 3;
		upper_[0] = -3;
		upper_[1] = 7;
		for (size_t i = 0; i < NUM_IN; ++i) {
			data_[i] = base_input[i];
		}
	}
//	SIMD_ALIGN
//	int upper_[NUM_RANGE];
//	SIMD_ALIGN
//	int lower_[NUM_RANGE];
};

/*
 * Tests various bool filter generation using int array
 * INPUTS:
 * - data = [0, -5, -10, -5, 0, 5, 10, 5]
 * - lower_bound = [-7, 3]
 * - upper_bound = [-3, 7]
 */
class BoolFilterOther: public BoolFilter<int> {
protected:
	virtual void PrepareInputs() {
	}
};

/*
 * Test bool filter generation sakura_SetTrueFloatInRangesInclusive
 * RESULT:
 * result = [F, T, F, T, F, T, F, T]
 */
TEST_F(BoolFilterFloat, RangesInclusive) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	bool answer[] = { false, true, false, true, false, true, false, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		PrintArray("lower_bound", NUM_RANGE, lower);
		PrintArray("upper_bound", NUM_RANGE, upper);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, in_data, num_range, lower, upper, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatInRangesInclusive
 * with an array of 11 elements (num_data=11).
 * RESULT:
 * result = [F, T, F, T, F, T, F, T
 *        F, T, F]
 */
TEST_F(BoolFilterFloat, RangesInclusiveLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	bool answer[] = { false, true, false, true, false, true, false, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		PrintArray("lower_bound", NUM_RANGE, lower);
		PrintArray("upper_bound", NUM_RANGE, upper);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, in_data, num_range, lower, upper, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatInRangesInclusive
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterFloat, RangesInclusiveLengthLenghZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, in_data, num_range, lower, upper, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_SetTrueFloatInRangesInclusive
 * without bounderies (num_condition=0)
 * RESULT:
 * result = [F, F, F, F, F, F, F, F]
 */
TEST_F(BoolFilterFloat, RangesInclusiveZeroCondition) {
	size_t const num_data(NUM_IN), num_range(0);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float lower[num_range];
	SIMD_ALIGN
	float upper[num_range];

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, in_data, num_range, lower, upper, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(false, result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntInRangesInclusive
 * RESULT:
 * result = [F, T, F, T, F, T, F, T]
 */
TEST_F(BoolFilterInt, RangesInclusive) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	int lower[NUM_RANGE];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	bool answer[] = { false, true, false, true, false, true, false, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		PrintArray("lower_bound", NUM_RANGE, lower);
		PrintArray("upper_bound", NUM_RANGE, upper);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, in_data, num_range, lower, upper, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntInRangesInclusive with an array of 11 elements
 * RESULT:
 * result = [F, T, F, T, F, T, F, T,
 *        F, T, F]
 */
TEST_F(BoolFilterInt, RangesInclusiveLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	int lower[NUM_RANGE];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	bool answer[] = { false, true, false, true, false, true, false, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		PrintArray("lower_bound", NUM_RANGE, lower);
		PrintArray("upper_bound", NUM_RANGE, upper);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, in_data, num_range, lower, upper, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntInRangesInclusive
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterInt, RangesInclusiveLenghZero) {
	size_t const num_data(0), num_range(NUM_RANGE);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	int lower[num_range];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];

	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, in_data, num_range, lower, upper, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_SetTrueIntInRangesInclusive
 * without bounderies (num_condition=0)
 * RESULT:
 * result = [F, F, F, F, F, F, F, F]
 */
TEST_F(BoolFilterInt, RangesInclusiveZeroCondition) {
	size_t const num_data(NUM_IN), num_range(0);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	int lower[num_range];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, in_data, num_range, lower, upper, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(false, result[i]);
	}
}

/*
 * Test bool filter generation sakura_InvertBool
 * INPUT:
 * in = [T, F, F, T]
 * RESULT:
 * result = [F, T, T, F]
 */
TEST_F(BoolFilterOther, InvertBool) {
	size_t const num_data(4);
	SIMD_ALIGN
	bool const data[] = { true, false, false, true };
	STATIC_ASSERT(ELEMENTSOF(data) == num_data);
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];
	bool answer[] = { false, true, true, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == num_data);

	if (verbose)
		PrintArray("data", num_data, data);

	LIBSAKURA_SYMBOL(Status) status = sakura_InvertBool(num_data, data, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i], result[i]);
	}
}

/*
 * Test bool filter generation sakura_InvertBool
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterOther, InvertBoolLenghZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	bool data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_InvertBool(num_data, data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_Uint8ToBool
 * INPUT:
 * in = [00000000, 00000001, 00000010, 00000100, 00001000, 00010000, 00100000, 01000000]
 * RESULT:
 * result = [F, T, T, T, T, T, T, T]
 */
TEST_F(BoolFilterOther, Uint8ToBool) {
	size_t const num_data(8);
	SIMD_ALIGN
	uint8_t const data[] = { 0, 1, 2, 4, 8, 16, 32, 64 };
	STATIC_ASSERT(ELEMENTSOF(data) == num_data);
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];
	bool answer[] = { false, true, true, true, true, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == num_data);

//	if (verbose) PrintArray("data", num_data, in);

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint8ToBool(num_data, data,
			result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i], result[i]);
	}
}

/*
 * Test bool filter generation sakura_Uint8ToBool
 * with an array of zero elements (num_data=0).
 *
 * INPUT: in = []
 * RESULT: result = []
 */
TEST_F(BoolFilterOther, Uint8ToBoolLenghZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint8ToBool(num_data, data,
			result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_Uint32ToBool with
 * INPUT:
 * in = [0...000, 0...001, 0...010, 0...01000, 0...0100000000]
 * RESULT:
 * result = [F, T, T, T, T]
 */
TEST_F(BoolFilterOther, Uint32ToBool) {
	size_t const num_data(5);
	SIMD_ALIGN
	uint32_t const data[] = { 0, 1, (1 << 1), (1 << 3), (1 << 8) };
	STATIC_ASSERT(ELEMENTSOF(data) == num_data);
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];
	bool answer[] = { false, true, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == num_data);

//	if (verbose) PrintArray("data", num_data, in);

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint32ToBool(num_data, data,
			result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i], result[i]);
	}
}

/*
 * Test bool filter generation sakura_Uint32ToBool
 * with an array of zero elements (num_data=0).
 *
 * INPUT: in = []
 * RESULT: result = []
 */
/* Input array is zero length */
TEST_F(BoolFilterOther, Uint32ToBoolLenghZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint32ToBool(num_data, data,
			result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test failure cases of sakura_SetTrueFloatInRangesInclusive
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* lower_bound > upper_bound */
TEST_F(BoolFilterFloat, RangesInclusiveFailExchangeBounds) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		// Exchange upper and lower bounds
		PrintArray("lower_bound", NUM_RANGE, upper);
		PrintArray("upper_bound", NUM_RANGE, lower);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, in_data, num_range, upper, lower, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Null pointer arrays */
TEST_F(BoolFilterFloat, RangesInclusiveFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	float *data_null = nullptr;

	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, data_null, num_range, lower, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, RangesInclusiveFailNullLower) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float upper[NUM_RANGE];

	float dummy[ELEMENTSOF(upper)];
	size_t const num_range(ELEMENTSOF(upper));

	float *lower_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// Copy bounds to aligned arrays
	GetBounds(dummy, upper);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, in_data, num_range, lower_null, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, RangesInclusiveFailNullUpper) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float lower[NUM_RANGE];

	float dummy[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	float *upper_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(in_data), in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, dummy);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, in_data, num_range, lower, upper_null, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, RangesInclusiveFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];

	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	bool *result_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(in_data), in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, in_data, num_range, lower, upper, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterFloat, RangesInclusiveNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	float data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	// Define unaligned array
	float *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, data_shift, num_range, lower, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, RangesInclusiveNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	float data[num_data];
	SIMD_ALIGN
	bool result[num_elements];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, data, num_range, lower, upper, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, RangesInclusiveNotAlignedLower) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];
	SIMD_ALIGN
	float lower[NUM_RANGE + offset];
	SIMD_ALIGN
	float upper[NUM_RANGE];
	size_t const num_range(ELEMENTSOF(upper));

	float dummy[ELEMENTSOF(upper)];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);
	// Copy bounds to aligned arrays
	GetBounds(dummy, upper);
	// Initialize lower
	for (size_t i = 0; i < ELEMENTSOF(dummy); ++i) {
		lower[i + offset] = dummy[i];
	}
	// Define unaligned array
	float *lower_shift = &lower[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(lower_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, data, num_range, lower_shift, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, RangesInclusiveNotAlignedUpper) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	size_t const num_range(ELEMENTSOF(lower));
	SIMD_ALIGN
	float upper[NUM_RANGE + offset];

	float dummy[ELEMENTSOF(lower)];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);
	// Copy bounds to aligned arrays
	GetBounds(lower, dummy);
	// Initialize upper
	for (size_t i = 0; i < ELEMENTSOF(dummy); ++i) {
		upper[i + offset] = dummy[i];
	}
	// Define unaligned array
	float *upper_shift = &upper[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(upper_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_data, data, num_range, lower, upper_shift, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_SetTrueIntInRangesInclusive
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* lower_bound > upper_bound */
TEST_F(BoolFilterInt, RangesInclusiveFailExchangeBounds) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	int lower[NUM_RANGE];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		// Exchange upper and lower bounds
		PrintArray("lower_bound", NUM_RANGE, upper);
		PrintArray("upper_bound", NUM_RANGE, lower);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, in_data, num_range, upper, lower, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Null pointer arrays */
TEST_F(BoolFilterInt, RangesInclusiveFailNullData) {
	size_t const num_data(NUM_IN);

	SIMD_ALIGN
	bool result[num_data];
	SIMD_ALIGN
	int lower[NUM_RANGE];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	int *data_null = nullptr;

	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, data_null, num_range, lower, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, RangesInclusiveFailNullLower) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	int upper[NUM_RANGE];

	int dummy[ELEMENTSOF(upper)];
	size_t const num_range(ELEMENTSOF(upper));

	int *lower_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// Copy bounds to aligned arrays
	GetBounds(dummy, upper);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, in_data, num_range, lower_null, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, RangesInclusiveFailNullUpper) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	int lower[NUM_RANGE];

	int dummy[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	int *upper_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, dummy);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, in_data, num_range, lower, upper_null, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, RangesInclusiveFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];

	SIMD_ALIGN
	int lower[NUM_RANGE];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	bool *result_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, in_data, num_range, lower, upper, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterInt, RangesInclusiveNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	int data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	SIMD_ALIGN
	int lower[NUM_RANGE];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	// Define unaligned array
	int *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, data_shift, num_range, lower, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, RangesInclusiveNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	int data[num_data];
	SIMD_ALIGN
	bool result[num_elements];
	SIMD_ALIGN
	int lower[NUM_RANGE];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, data, num_range, lower, upper, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, RangesInclusiveNotAlignedLower) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];
	SIMD_ALIGN
	int lower[NUM_RANGE + offset];
	SIMD_ALIGN
	int upper[NUM_RANGE];
	size_t const num_range(ELEMENTSOF(upper));

	int dummy[ELEMENTSOF(upper)];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);
	// Copy bounds to aligned arrays
	GetBounds(dummy, upper);
	// Initialize lower
	for (size_t i = 0; i < ELEMENTSOF(dummy); ++i) {
		lower[i + offset] = dummy[i];
	}
	// Define unaligned array
	int *lower_shift = &lower[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(lower_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, data, num_range, lower_shift, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, RangesInclusiveNotAlignedUpper) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];
	SIMD_ALIGN
	int lower[NUM_RANGE];
	size_t const num_range(ELEMENTSOF(lower));
	SIMD_ALIGN
	int upper[NUM_RANGE + offset];

	int dummy[ELEMENTSOF(lower)];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);
	// Copy bounds to aligned arrays
	GetBounds(lower, dummy);
	// Initialize upper
	for (size_t i = 0; i < ELEMENTSOF(dummy); ++i) {
		upper[i + offset] = dummy[i];
	}
	// Define unaligned array
	int *upper_shift = &upper[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(upper_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(
			num_data, data, num_range, lower, upper_shift, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_InvertBool
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterOther, InvertBoolFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	bool *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_InvertBool(num_data, data_null,
			result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterOther, InvertBoolFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool data[num_data];
	// Initialize data with whatever valid value
	for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
		data[i] = true;
	}

	bool *result_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_InvertBool(num_data, data,
			result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterOther, InvertBoolNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	bool data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	// Initialize data with whatever valid value
	for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
		data[i] = true;
	}

	// Define unaligned array
	bool *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_InvertBool(num_data, data_shift,
			result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterOther, InvertBoolNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	bool data[num_data];
	SIMD_ALIGN
	bool result[num_elements];
	// Initialize data with whatever valid value
	for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
		data[i] = true;
	}

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_InvertBool(num_data, data,
			result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_Uint8ToBool
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterOther, Uint8ToBoolFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	uint8_t *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint8ToBool(num_data, data_null,
			result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterOther, Uint8ToBoolFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	// Initialize data with whatever valid value
	uint8_t data[num_data];
	for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
		data[i] = 1;
	}

	bool *result_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint8ToBool(num_data, data,
			result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterOther, Uint8ToBoolNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(8);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	// Initialize data with whatever valid value
	for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
		data[i] = 1;
	}

	// Define unaligned array
	uint8_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint8ToBool(num_data, data_shift,
			result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterOther, Uint8ToBoolNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint8_t data[num_data];
	SIMD_ALIGN
	bool result[num_elements];
	// Initialize data with whatever valid value
	for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
		data[i] = 1;
	}

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint8ToBool(num_data, data,
			result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_Uint32ToBool
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterOther, Uint32ToBoolFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	uint32_t *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint32ToBool(num_data, data_null,
			result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterOther, Uint32ToBoolFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	uint32_t data[num_data];
	// Initialize data with whatever valid value
	for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
		data[i] = 1;
	}

	bool *result_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint32ToBool(num_data, data,
			result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterOther, Uint32ToBoolNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	// Initialize data with whatever valid value
	for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
		data[i] = 1;
	}

	// Define unaligned array
	uint32_t *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint32ToBool(num_data, data_shift,
			result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterOther, Uint32ToBoolNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	uint32_t data[num_data];
	SIMD_ALIGN
	bool result[num_elements];
	// Initialize data with whatever valid value
	for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
		data[i] = 1;
	}

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint32ToBool(num_data, data,
			result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

