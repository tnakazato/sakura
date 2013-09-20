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
	void GetDataInLength(size_t num_out, DataType *out) {
		size_t const num_in(NUM_IN);
		for (size_t i = 0; i < num_out; ++i) {
			out[i] = data_[i % num_in];
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

	void PrintArray(char const *name, size_t num_in, DataType const *in) {
		cout << name << " = [ ";
		for (size_t i = 0; i < num_in - 1; ++i)
			cout << in[i] << ", ";
		cout << in[num_in - 1];
		cout << " ]" << endl;
	}

	void PrintArray(char const *name, size_t num_in, bool const *in) {
		cout << name << " = [ ";
		for (size_t i = 0; i < num_in - 1; ++i)
			cout << (in[i] ? "T" : "F") << ", ";
		cout << (in[num_in - 1] ? "T" : "F");
		cout << " ]" << endl;
	}

	DataType data_[NUM_IN];
	DataType upper_[NUM_RANGE];
	DataType lower_[NUM_RANGE];bool verbose;
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
		ASSERT_EQ(NUM_IN, ELEMENTSOF(base_input));
		ASSERT_EQ(2, NUM_RANGE);
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
		ASSERT_EQ(NUM_IN, ELEMENTSOF(base_input));
		ASSERT_EQ(2, NUM_RANGE);
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
 * out = [F, T, F, T, F, T, F, T]
 */
TEST_F(BoolFilterFloat, RangesInclusive) {
	size_t const num_in(NUM_IN);
	SIMD_ALIGN
	float in_data[num_in];
	SIMD_ALIGN
	bool out[ELEMENTSOF(data_)];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	bool answer[] = { false, true, false, true, false, true, false, true };
	ASSERT_EQ(NUM_IN, ELEMENTSOF(answer));

	// Create long input data by repeating data_
	GetDataInLength(num_in, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	if (verbose) {
		PrintArray("in", num_in, in_data);
		PrintArray("lower_bound", NUM_RANGE, lower);
		PrintArray("upper_bound", NUM_RANGE, upper);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_in, in_data, num_range, lower, upper, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], out[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatInRangesInclusive
 * with an array of 11 elements (num_data=11).
 * RESULT:
 * out = [F, T, F, T, F, T, F, T
 *        F, T, F]
 */
TEST_F(BoolFilterFloat, RangesInclusiveLengthEleven) {
	size_t const num_in(11);
	SIMD_ALIGN
	float in_data[num_in];
	SIMD_ALIGN
	bool out[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	bool answer[] = { false, true, false, true, false, true, false, true };
	ASSERT_EQ(NUM_IN, ELEMENTSOF(answer));

	// Create long input data by repeating data_
	GetDataInLength(num_in, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	if (verbose) {
		PrintArray("in", num_in, in_data);
		PrintArray("lower_bound", NUM_RANGE, lower);
		PrintArray("upper_bound", NUM_RANGE, upper);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_in, in_data, num_range, lower, upper, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], out[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatInRangesInclusive
 * with an array of zero elements (num_data=0).
 * RESULT:
 * out = []
 */
TEST_F(BoolFilterFloat, RangesInclusiveLengthLenghZero) {
	size_t const num_in(0);
	SIMD_ALIGN
	float in_data[num_in];
	SIMD_ALIGN
	bool out[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_in, in_data, num_range, lower, upper, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_SetTrueIntInRangesInclusive
 * RESULT:
 * out = [F, T, F, T, F, T, F, T]
 */
TEST_F(BoolFilterInt, RangesInclusive) {
	size_t const num_in(NUM_IN);
	SIMD_ALIGN
	int in_data[num_in];
	SIMD_ALIGN
	bool out[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	int lower[NUM_RANGE];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	bool answer[] = { false, true, false, true, false, true, false, true };
	ASSERT_EQ(NUM_IN, ELEMENTSOF(answer));

	// Create long input data by repeating data_
	GetDataInLength(num_in, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	if (verbose) {
		PrintArray("in", num_in, in_data);
		PrintArray("lower_bound", NUM_RANGE, lower);
		PrintArray("upper_bound", NUM_RANGE, upper);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(num_in,
			in_data, num_range, lower, upper, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], out[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntInRangesInclusive with an array of 11 elements
 * RESULT:
 * out = [F, T, F, T, F, T, F, T,
 *        F, T, F]
 */
TEST_F(BoolFilterInt, RangesInclusiveLengthEleven) {
	size_t const num_in(11);
	SIMD_ALIGN
	int in_data[num_in];
	SIMD_ALIGN
	bool out[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	int lower[NUM_RANGE];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	bool answer[] = { false, true, false, true, false, true, false, true };
	ASSERT_EQ(NUM_IN, ELEMENTSOF(answer));

	// Create long input data by repeating data_
	GetDataInLength(num_in, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	if (verbose) {
		PrintArray("in", num_in, in_data);
		PrintArray("lower_bound", NUM_RANGE, lower);
		PrintArray("upper_bound", NUM_RANGE, upper);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(num_in,
			in_data, num_range, lower, upper, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], out[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntInRangesInclusive
 * with an array of zero elements (num_data=0).
 * RESULT:
 * out = []
 */
TEST_F(BoolFilterInt, RangesInclusiveLenghZero) {
	size_t const num_in(0), num_range(NUM_RANGE);
	SIMD_ALIGN
	int in_data[num_in];
	SIMD_ALIGN
	bool out[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	int lower[num_range];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];

	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(num_in,
			in_data, num_range, lower, upper, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_InvertBool
 * INPUT:
 * in = [T, F, F, T]
 * RESULT:
 * out = [F, T, T, F]
 */
TEST_F(BoolFilterOther, InvertBool) {
	size_t const num_in(4);
	SIMD_ALIGN
	bool const in[] = { true, false, false, true };
	ASSERT_EQ(num_in, ELEMENTSOF(in));
	SIMD_ALIGN
	bool out[num_in];
	bool answer[] = { false, true, true, false };
	ASSERT_EQ(num_in, ELEMENTSOF(answer));

	if (verbose)
		PrintArray("in", num_in, in);

	LIBSAKURA_SYMBOL(Status) status = sakura_InvertBool(num_in, in, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i], out[i]);
	}
}

/*
 * Test bool filter generation sakura_InvertBool
 * with an array of zero elements (num_data=0).
 * RESULT:
 * out = []
 */
TEST_F(BoolFilterOther, InvertBoolLenghZero) {
	size_t const num_in(0);
	SIMD_ALIGN
	bool in[num_in];
	SIMD_ALIGN
	bool out[ELEMENTSOF(in)];

	LIBSAKURA_SYMBOL(Status) status = sakura_InvertBool(num_in, in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_Uint8ToBool
 * INPUT:
 * in = [00000000, 00000001, 00000010, 00000100, 00001000, 00010000, 00100000, 01000000]
 * RESULT:
 * out = [F, T, T, T, T, T, T, T]
 */
TEST_F(BoolFilterOther, Uint8ToBool) {
	size_t const num_in(8);
	SIMD_ALIGN
	uint8_t const in[] = { 0, 1, 2, 4, 8, 16, 32, 64 };
	ASSERT_EQ(num_in, ELEMENTSOF(in));
	SIMD_ALIGN
	bool out[num_in];
	bool answer[] = { false, true, true, true, true, true, true, true };
	ASSERT_EQ(num_in, ELEMENTSOF(answer));

//	if (verbose) PrintArray("in", num_in, in);

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint8ToBool(num_in, in, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i], out[i]);
	}
}

/*
 * Test bool filter generation sakura_Uint8ToBool
 * with an array of zero elements (num_data=0).
 *
 * INPUT: in = []
 * RESULT: out = []
 */
TEST_F(BoolFilterOther, Uint8ToBoolLenghZero) {
	size_t const num_in(0);
	SIMD_ALIGN
	uint8_t in[num_in];
	SIMD_ALIGN
	bool out[ELEMENTSOF(in)];

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint8ToBool(num_in, in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_Uint32ToBool with
 * INPUT:
 * in = [0...000, 0...001, 0...010, 0...01000, 0...0100000000]
 * RESULT:
 * out = [F, T, T, T, T]
 */
TEST_F(BoolFilterOther, Uint32ToBool) {
	size_t const num_in(5);
	SIMD_ALIGN
	uint32_t const in[] = { 0, 1, (1 << 1), (1 << 3), (1 << 8) };
	ASSERT_EQ(num_in, ELEMENTSOF(in));
	SIMD_ALIGN
	bool out[ELEMENTSOF(in)];
	bool answer[] = { false, true, true, true, true };
	ASSERT_EQ(num_in, ELEMENTSOF(answer));

//	if (verbose) PrintArray("in", num_in, in);

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint32ToBool(num_in, in, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(answer[i], out[i]);
	}
}

/*
 * Test bool filter generation sakura_Uint32ToBool
 * with an array of zero elements (num_data=0).
 *
 * INPUT: in = []
 * RESULT: out = []
 */
/* Input array is zero length */
TEST_F(BoolFilterOther, Uint32ToBoolLenghZero) {
	size_t const num_in(0);
	SIMD_ALIGN
	uint32_t in[num_in];
	SIMD_ALIGN
	bool out[ELEMENTSOF(in)];

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint32ToBool(num_in, in, out);

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
	size_t const num_in(NUM_IN);
	SIMD_ALIGN
	float in_data[num_in];
	SIMD_ALIGN
	bool out[ELEMENTSOF(data_)];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	// Create long input data by repeating data_
	GetDataInLength(num_in, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	if (verbose) {
		PrintArray("in", num_in, in_data);
		// Exchange upper and lower bounds
		PrintArray("lower_bound", NUM_RANGE, upper);
		PrintArray("upper_bound", NUM_RANGE, lower);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_in, in_data, num_range, upper, lower, out);
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
	size_t const num_in(NUM_IN);
	SIMD_ALIGN
	int in_data[num_in];
	SIMD_ALIGN
	bool out[ELEMENTSOF(data_)];
	SIMD_ALIGN
	int lower[NUM_RANGE];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	// Create long input data by repeating data_
	GetDataInLength(num_in, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	if (verbose) {
		PrintArray("in", num_in, in_data);
		// Exchange upper and lower bounds
		PrintArray("lower_bound", NUM_RANGE, upper);
		PrintArray("upper_bound", NUM_RANGE, lower);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(num_in,
			in_data, num_range, upper, lower, out);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_InvertBool
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */

/*
 * Test failure cases of sakura_Uint8ToBool
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */

/*
 * Test failure cases of sakura_Uint32ToBool
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
