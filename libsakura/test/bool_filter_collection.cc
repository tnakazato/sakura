#include <iostream>
#include <string>
#include <sys/time.h>

#include <libsakura/sakura.h>
#include "aligned_memory.h"
#include "gtest/gtest.h"

/* the number of elements in input/output array to test */
#define NUM_IN 8
#define NUM_RANGE 2 // DO NOT MODIFY THIS!
#define NUM_IN_LONG 262144 //2**18
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
	}

	virtual void TearDown() {
		// Clean-up sakura
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	virtual void PrepareInputs() = 0;

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

	SIMD_ALIGN DataType data_[NUM_IN];SIMD_ALIGN DataType upper_[NUM_RANGE];SIMD_ALIGN DataType lower_[NUM_RANGE];bool verbose;
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
		float const base_input[8] = { 0., -0.5, -1.0, -0.5, 0., 0.5, 1.0, 0.5 };
		ASSERT_EQ(NUM_RANGE, 2);
		lower_[0] = -0.75;
		lower_[1] = 0.25;
		upper_[0] = -0.25;
		upper_[1] = 0.75;
		for (size_t i = 0; i < NUM_IN; ++i) {
			data_[i] = base_input[i % 8];
		}
		if (verbose) {
			PrintArray("data", NUM_IN, data_);
			PrintArray("lower_bound", NUM_RANGE, lower_);
			PrintArray("upper_bound", NUM_RANGE, upper_);
		}
	}
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
		int const base_input[8] = { 0, -5, -10, -5, 0, 5, 10, 5 };
		ASSERT_EQ(NUM_RANGE, 2);
		lower_[0] = -7;
		lower_[1] = 3;
		upper_[0] = -3;
		upper_[1] = 7;
		for (size_t i = 0; i < NUM_IN; ++i) {
			data_[i] = base_input[i % 8];
		}
		if (verbose) {
			PrintArray("data", NUM_IN, data_);
			PrintArray("lower_bound", NUM_RANGE, lower_);
			PrintArray("upper_bound", NUM_RANGE, upper_);
		}
	}
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
	SIMD_ALIGN
	bool out[NUM_IN];
	bool answer[8] = { false, true, false, true, false, true, false, true };
	size_t const num_in(NUM_IN), num_range(NUM_RANGE);

	PrepareInputs();

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(
			num_in, data_, num_range, lower_, upper_, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(out[i], answer[i % 8]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntInRangesInclusive
 * RESULT:
 * out = [F, T, F, T, F, T, F, T]
 */
TEST_F(BoolFilterInt, RangesInclusive) {
	SIMD_ALIGN
	bool out[NUM_IN];
	bool answer[8] = { false, true, false, true, false, true, false, true };
	size_t const num_in(NUM_IN), num_range(NUM_RANGE);

	PrepareInputs();

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(num_in,
			data_, num_range, lower_, upper_, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(out[i], answer[i % 8]);
	}
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
	bool const in[num_in] = { true, false, false, true };
	SIMD_ALIGN
	bool out[num_in];
	bool answer[num_in] = { false, true, true, false };

	if (verbose)
		PrintArray("in", num_in, in);

	LIBSAKURA_SYMBOL(Status) status = sakura_InvertBool(num_in, in, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(out[i], answer[i]);
	}
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
	uint8_t const in[num_in] = { 0, 1, 2, 4, 8, 16, 32, 64 };
	SIMD_ALIGN
	bool out[num_in];
	bool answer[num_in] = { false, true, true, true, true, true, true, true };

	//if (verbose) PrintArray("in", num_in, inint);

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint8ToBool(num_in, in, out);

	if (verbose)
		PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_in; ++i) {
		ASSERT_EQ(out[i], answer[i]);
	}
}

