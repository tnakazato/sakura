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
#include <math.h>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "loginit.h"
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
	DataType threshold_;
};

/*
 * Tests various bool filter generation using float array
 * INPUTS:
 * - data = [0., -0.5, -1.0, -0.5, 0., 0.5, 1.0, 0.5]
 * - lower_bound = [-0.75, 0.5]
 * - upper_bound = [-0.25, 0.75]
 */
class BoolFilterFloat: public BoolFilter<float> {
protected:
	virtual void PrepareInputs() {
		threshold_ = 0.;
		float const base_input[] = { 0., -0.5, -1.0, -0.5, 0., 0.5, 1.0, 0.5 };
		STATIC_ASSERT(ELEMENTSOF(base_input) == NUM_IN);
		STATIC_ASSERT(NUM_RANGE == 2);
		lower_[0] = -0.75;
		lower_[1] = 0.5;
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
 * - lower_bound = [-7, 5]
 * - upper_bound = [-3, 7]
 */
class BoolFilterInt: public BoolFilter<int> {
protected:
	virtual void PrepareInputs() {
		threshold_ = 0;
		int const base_input[] = { 0, -5, -10, -5, 0, 5, 10, 5 };
		STATIC_ASSERT(ELEMENTSOF(base_input) == NUM_IN);
		STATIC_ASSERT(NUM_RANGE == 2);
		lower_[0] = -7;
		lower_[1] = 5;
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
 * - lower_bound = [-7, 5]
 * - upper_bound = [-3, 7]
 */
class BoolFilterOther: public BoolFilter<int> {
protected:
	virtual void PrepareInputs() {
	}
};

/////////////////////////////////////////////////////////////////////////////////////////

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
 * with a long array
 * RESULT:
 * result = [F, T, F, T, F, T, F, T, ...]
 */
TEST_F(BoolFilterFloat, RangesInclusiveLong) {
	size_t const num_data(NUM_IN_LONG);
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

	double start, end;
	size_t const num_repeat = 2000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetTrueFloatInRangesInclusive(num_data, in_data,
				num_range, lower, upper, result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

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
TEST_F(BoolFilterFloat, RangesInclusiveLenghZero) {
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
		ASSERT_FALSE(result[i]);
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
 * with a long array
 * RESULT:
 * result = [F, T, F, T, F, T, F, T, ...]
 */
TEST_F(BoolFilterInt, RangesInclusiveLong) {
	size_t const num_data(NUM_IN_LONG);
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

	double start, end;
	size_t const num_repeat = 2000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetTrueIntInRangesInclusive(num_data, in_data,
				num_range, lower, upper, result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

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
		ASSERT_FALSE(result[i]);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bool filter generation sakura_SetTrueFloatInRangesExclusive
 * RESULT:
 * result = [F, T, F, T, F, F, F, F]
 */
TEST_F(BoolFilterFloat, RangesExclusive) {
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

	bool answer[] = { false, true, false, true, false, false, false, false };
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
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
 * Test bool filter generation sakura_SetTrueFloatInRangesExclusive
 * with an array of 11 elements (num_data=11).
 * RESULT:
 * result = [F, T, F, T, F, F, F, F,
 *           F, T, F]
 */
TEST_F(BoolFilterFloat, RangesExclusiveLengthEleven) {
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

	bool answer[] = { false, true, false, true, false, false, false, false };
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
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
 * Test bool filter generation sakura_SetTrueFloatInRangesExclusive
 * with a long array
 * RESULT:
 * result = [F, T, F, T, F, F, F, F, ...]
 */
TEST_F(BoolFilterFloat, RangesExclusiveLong) {
	size_t const num_data(NUM_IN_LONG);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float lower[NUM_RANGE];
	SIMD_ALIGN
	float upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	bool answer[] = { false, true, false, true, false, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	double start, end;
	size_t const num_repeat = 2000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetTrueFloatInRangesExclusive(num_data, in_data,
				num_range, lower, upper, result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatInRangesExclusive
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterFloat, RangesExclusiveLenghZero) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
			num_data, in_data, num_range, lower, upper, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_SetTrueFloatInRangesExclusive
 * without bounderies (num_condition=0)
 * RESULT:
 * result = [F, F, F, F, F, F, F, F]
 */
TEST_F(BoolFilterFloat, RangesExclusiveZeroCondition) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
			num_data, in_data, num_range, lower, upper, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_FALSE(result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntInRangesExclusive
 * RESULT:
 * result = [F, T, F, T, F, F, F, F]
 */
TEST_F(BoolFilterInt, RangesExclusive) {
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

	bool answer[] = { false, true, false, true, false, false, false, false };
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
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
 * Test bool filter generation sakura_SetTrueIntInRangesExclusive with an array of 11 elements
 * RESULT:
 * result = [F, T, F, T, F, F, F, F,
 *           F, T, F]
 */
TEST_F(BoolFilterInt, RangesExclusiveLengthEleven) {
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

	bool answer[] = { false, true, false, true, false, false, false, false };
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
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
 * Test bool filter generation sakura_SetTrueIntInRangesExclusive
 * with a long array
 * RESULT:
 * result = [F, T, F, T, F, F, F, F, ...]
 */
TEST_F(BoolFilterInt, RangesExclusiveLong) {
	size_t const num_data(NUM_IN_LONG);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	int lower[NUM_RANGE];
	SIMD_ALIGN
	int upper[ELEMENTSOF(lower)];
	size_t const num_range(ELEMENTSOF(lower));

	bool answer[] = { false, true, false, true, false, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// Copy bounds to aligned arrays
	GetBounds(lower, upper);

	double start, end;
	size_t const num_repeat = 2000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetTrueIntInRangesExclusive(num_data, in_data,
				num_range, lower, upper, result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntInRangesExclusive
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterInt, RangesExclusiveLenghZero) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
			num_data, in_data, num_range, lower, upper, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_SetTrueIntInRangesExclusive
 * without bounderies (num_condition=0)
 * RESULT:
 * result = [F, F, F, F, F, F, F, F]
 */
TEST_F(BoolFilterInt, RangesExclusiveZeroCondition) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
			num_data, in_data, num_range, lower, upper, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_FALSE(result[i]);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bool filter generation sakura_SetTrueFloatGreaterThan
 * RESULT:
 * result = [F, F, F, F, F, T, T, T]
 */
TEST_F(BoolFilterFloat, GreaterThan) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { false, false, false, false, false, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	//verbose = true;
	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThan(num_data,
			in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatGreaterThan
 * with an array of 11 elements (num_data=11).
 * RESULT:
 * result = [F, F, F, F, F, T, T, T,
 *           F, F, F]
 */
TEST_F(BoolFilterFloat, GreaterThanLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { false, false, false, false, false, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThan(num_data,
			in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatGreaterThan
 * with a long array
 * RESULT:
 * result = [F, F, F, F, F, T, T, T, ...]
 */
TEST_F(BoolFilterFloat, GreaterThanLong) {
	size_t const num_data(NUM_IN_LONG);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { false, false, false, false, false, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetTrueFloatGreaterThan(num_data, in_data, threshold_,
				result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatGreaterThan
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterFloat, GreaterThanLenghZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThan(num_data,
			in_data, threshold_, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_SetTrueIntGreaterThan
 * RESULT:
 * result = [F, F, F, F, F, T, T, T]
 */
TEST_F(BoolFilterInt, GreaterThan) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { false, false, false, false, false, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	//verbose = true;
	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThan(num_data,
			in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntGreaterThan
 * with an array of 11 elements (num_data=11).
 * RESULT:
 * result = [F, F, F, F, F, T, T, T,
 *           F, F, F]
 */
TEST_F(BoolFilterInt, GreaterThanLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { false, false, false, false, false, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThan(num_data,
			in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntGreaterThan
 * with a long array
 * RESULT:
 * result = [F, F, F, F, F, T, T, T, ...]
 */
TEST_F(BoolFilterInt, GreaterThanLong) {
	size_t const num_data(NUM_IN_LONG);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { false, false, false, false, false, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetTrueIntGreaterThan(num_data, in_data, threshold_,
				result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntGreaterThan
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterInt, GreaterThanLengthZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThan(num_data,
			in_data, threshold_, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bool filter generation sakura_SetTrueFloatGreaterThanOrEquals
 * RESULT:
 * result = [T, F, F, F, T, T, T, T]
 */
TEST_F(BoolFilterFloat, GreaterThanOrEquals) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, false, false, false, true, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	//verbose = true;
	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThanOrEquals(
			num_data, in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatGreaterThanOrEquals
 * with an array of 11 elements (num_data=11).
 * RESULT:
 * result = [T, F, F, F, T, T, T, T,
 *           T, F, F]
 */
TEST_F(BoolFilterFloat, GreaterThanOrEqualsLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, false, false, false, true, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThanOrEquals(
			num_data, in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatGreaterThanOrEquals
 * with a long array
 * RESULT:
 * result = [T, F, F, F, T, T, T, T, ...]
 */
TEST_F(BoolFilterFloat, GreaterThanOrEqualsLong) {
	size_t const num_data(NUM_IN_LONG);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, false, false, false, true, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetTrueFloatGreaterThanOrEquals(num_data, in_data,
				threshold_, result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatGreaterThanOrEquals
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterFloat, GreaterThanOrEqualsLenghZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThanOrEquals(
			num_data, in_data, threshold_, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_SetTrueIntGreaterThanOrEquals
 * RESULT:
 * result = [T, F, F, F, T, T, T, T]
 */
TEST_F(BoolFilterInt, GreaterThanOrEquals) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, false, false, false, true, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	//verbose = true;
	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThanOrEquals(
			num_data, in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntGreaterThanOrEquals
 * with an array of 11 elements (num_data=11).
 * RESULT:
 * result = [T, F, F, F, T, T, T, T,
 *           T, F, F]
 */
TEST_F(BoolFilterInt, GreaterThanOrEqualsLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, false, false, false, true, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThanOrEquals(
			num_data, in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntGreaterThanOrEquals
 * with a long array
 * RESULT:
 * result = [T, F, F, F, T, T, T, T, ...]
 */
TEST_F(BoolFilterInt, GreaterThanOrEqualsLong) {
	size_t const num_data(NUM_IN_LONG);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, false, false, false, true, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetTrueIntGreaterThanOrEquals(num_data, in_data,
				threshold_, result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntGreaterThanOrEquals
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterInt, GreaterThanOrEqualsLenghZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThanOrEquals(
			num_data, in_data, threshold_, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bool filter generation sakura_SetTrueFloatLessThan
 * RESULT:
 * result = [F, T, T, T, F, F, F, F]
 */
TEST_F(BoolFilterFloat, LessThan) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { false, true, true, true, false, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	//verbose = true;
	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThan(num_data,
			in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatLessThan
 * with an array of 11 elements (num_data=11).
 * RESULT:
 * result = [F, T, T, T, F, F, F, F,
 *           F, T, T]
 */
TEST_F(BoolFilterFloat, LessThanLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { false, true, true, true, false, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThan(num_data,
			in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatLessThan
 * with a long array
 * RESULT:
 * result = [F, T, T, T, F, F, F, F, ...]
 */
TEST_F(BoolFilterFloat, LessThanLong) {
	size_t const num_data(NUM_IN_LONG);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { false, true, true, true, false, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetTrueFloatLessThan(num_data, in_data, threshold_,
				result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatLessThan
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterFloat, LessThanLenghZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThan(num_data,
			in_data, threshold_, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_SetTrueIntLessThan
 * RESULT:
 * result = [F, T, T, T, F, F, F, F]
 */
TEST_F(BoolFilterInt, LessThan) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { false, true, true, true, false, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	//verbose = true;
	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThan(num_data,
			in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntLessThan
 * with an array of 11 elements (num_data=11).
 * RESULT:
 * result = [F, T, T, T, F, F, F, F,
 *           F, T, T]
 */
TEST_F(BoolFilterInt, LessThanLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { false, true, true, true, false, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThan(num_data,
			in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntLessThan
 * with a long array
 * RESULT:
 * result = [F, T, T, T, F, F, F, F, ...]
 */
TEST_F(BoolFilterInt, LessThanLong) {
	size_t const num_data(NUM_IN_LONG);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { false, true, true, true, false, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetTrueIntLessThan(num_data, in_data, threshold_,
				result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntLessThan
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterInt, LessThanLengthZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThan(num_data,
			in_data, threshold_, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bool filter generation sakura_SetTrueFloatLessThanOrEquals
 * RESULT:
 * result = [T, T, T, T, T, F, F, F]
 */
TEST_F(BoolFilterFloat, LessThanOrEquals) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, true, true, true, true, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	//verbose = true;
	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThanOrEquals(
			num_data, in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatLessThanOrEquals
 * with an array of 11 elements (num_data=11).
 * RESULT:
 * result = [T, T, T, T, T, F, F, F,
 *           T, T, T]
 */
TEST_F(BoolFilterFloat, LessThanOrEqualsLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, true, true, true, true, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThanOrEquals(
			num_data, in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatLessThanOrEquals
 * with a long array
 * RESULT:
 * result = [T, T, T, T, T, F, F, F, ...]
 */
TEST_F(BoolFilterFloat, LessThanOrEqualsLong) {
	size_t const num_data(NUM_IN_LONG);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, true, true, true, true, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetTrueFloatLessThanOrEquals(num_data, in_data,
				threshold_, result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueFloatLessThanOrEquals
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterFloat, LessThanOrEqualsLenghZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThanOrEquals(
			num_data, in_data, threshold_, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test bool filter generation sakura_SetTrueIntLessThanOrEquals
 * RESULT:
 * result = [T, T, T, T, T, F, F, F]
 */
TEST_F(BoolFilterInt, LessThanOrEquals) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, true, true, true, true, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	//verbose = true;
	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThanOrEquals(
			num_data, in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntLessThanOrEquals
 * with an array of 11 elements (num_data=11).
 * RESULT:
 * result = [T, T, T, T, T, F, F, F,
 *           T, T, T]
 */
TEST_F(BoolFilterInt, LessThanOrEqualsLengthEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, true, true, true, true, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThanOrEquals(
			num_data, in_data, threshold_, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntLessThanOrEquals
 * with a long array
 * RESULT:
 * result = [T, T, T, T, T, F, F, F, ...]
 */
TEST_F(BoolFilterInt, LessThanOrEqualsLong) {
	size_t const num_data(NUM_IN_LONG);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, true, true, true, true, false, false, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetTrueIntLessThanOrEquals(num_data, in_data,
				threshold_, result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetTrueIntLessThanOrEquals
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterInt, LessThanOrEqualsLenghZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	int in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThanOrEquals(
			num_data, in_data, threshold_, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bool filter generation sakura_SetFalseFloatIfNanOrInf
 * RESULT:
 * result = [T, T, F, F, F, T, T, T]
 */
TEST_F(BoolFilterFloat, NanOrInf) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];
	bool answer[] = { true, true, false, false, false, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// modify in_data and insert NaN and Infs
	float posinf(1.0 / 0.0), neginf(-1.0 / 0.0), nanval(0.0 / 0.0);
	float infnans[] = { posinf, neginf, nanval };
	size_t j(0);
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		if (!answer[i % ELEMENTSOF(answer)]) {
			in_data[i] = infnans[j % ELEMENTSOF(infnans)];
			++j;
		}
	}

	if (verbose) {
		PrintArray("data", num_data, in_data);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseFloatIfNanOrInf(num_data,
			in_data, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetFalseFloatIfNanOrInf
 * with an array of 11 elements (num_data=11).
 * RESULT:
 * result = [T, T, F, F, F, T, T, T,
 *           T, T, F]
 */
TEST_F(BoolFilterFloat, NanOrInfEleven) {
	size_t const num_data(11);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, true, false, false, false, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// modify in_data and insert NaN and Infs
	float posinf(1.0 / 0.0), neginf(-1.0 / 0.0), nanval(0.0 / 0.0);
	float infnans[] = { posinf, neginf, nanval };
	size_t j(0);
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		if (!answer[i % ELEMENTSOF(answer)]) {
			in_data[i] = infnans[j % ELEMENTSOF(infnans)];
			++j;
		}
	}

	if (verbose) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseFloatIfNanOrInf(num_data,
			in_data, result);

	if (verbose)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetFalseFloatIfNanOrInf
 * with a long array
 * RESULT:
 * result = [T, T, F, F, F, T, T, T, ...]
 */
TEST_F(BoolFilterFloat, NanOrInfLong) {
	size_t const num_data(NUM_IN_LONG);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	bool answer[] = { true, true, false, false, false, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);

	// Create long input data by repeating data_
	GetDataInLength(num_data, in_data);
	// modify in_data and insert NaN and Infs
	float posinf(1.0 / 0.0), neginf(-1.0 / 0.0), nanval(0.0 / 0.0);
	float infnans[] = { posinf, neginf, nanval };
	size_t j(0);
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		if (!answer[i % ELEMENTSOF(answer)]) {
			in_data[i] = infnans[j % ELEMENTSOF(infnans)];
			++j;
		}
	}

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_data << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_SetFalseFloatIfNanOrInf(num_data, in_data, result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetFalseFloatIfNanOrInf
 * with an array of zero elements (num_data=0).
 * RESULT:
 * result = []
 */
TEST_F(BoolFilterFloat, NanOrInfZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool result[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseFloatIfNanOrInf(num_data,
			in_data, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

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
 * with a long array
 * INPUT:
 * in = [T, F, F, T]
 * RESULT:
 * result = [F, T, T, F]
 */
TEST_F(BoolFilterOther, InvertBoolLong) {
	size_t const num_base(4);
	size_t const num_long(NUM_IN_LONG);
	SIMD_ALIGN
	bool const data_base[] = { true, false, false, true };
	STATIC_ASSERT(ELEMENTSOF(data_base) == num_base);
	bool answer[] = { false, true, true, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == num_base);
	SIMD_ALIGN
	bool data_long[num_long];
	SIMD_ALIGN
	bool result[ELEMENTSOF(data_long)];

	for (size_t i = 0; i < num_long; ++i) {
		data_long[i] = data_base[i % num_base];
	}

	double start, end;
	size_t const num_repeat = 100000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_long << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_InvertBool(num_long, data_long, result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_long; ++i) {
		ASSERT_EQ(answer[i % num_base], result[i]);
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

/////////////////////////////////////////////////////////////////////////////////////////

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
 * Test bool filter generation sakura_Uint8ToBool  (Long Test)
 * INPUT:
 * in = [00000000, 00000001, 00000010, 00000100, 00001000, 00010000, 00100000, 01000000, ...repeated...]
 * RESULT:
 * result = [F, T, T, T, T, T, T, T, ... repeated...]
 */
TEST_F(BoolFilterOther, Uint8ToBoolLong) {
	size_t const num_base(NUM_IN);
	size_t const num_long(NUM_IN_LONG);
	uint8_t const data_base[] = { 0, 1, 2, 4, 8, 16, 32, 64 };
	STATIC_ASSERT(ELEMENTSOF(data_base) == num_base);
	bool answer[] = { false, true, true, true, true, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == ELEMENTSOF(data_base));
	SIMD_ALIGN
	uint8_t data_long[num_long];
	SIMD_ALIGN
	bool result[ELEMENTSOF(data_long)];

	for (size_t i = 0; i < num_long; ++i) {
		data_long[i] = data_base[i % num_base];
	}

	double start, end;
	size_t const num_repeat = 100000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_long << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_Uint8ToBool(num_long, data_long, result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_long; ++i) {
		ASSERT_EQ(answer[i % num_base], result[i]);
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
 * Test bool filter generation sakura_Uint32ToBool
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
 * Test bool filter generation sakura_Uint8ToBool  (Long Test)
 * INPUT:
 * in = [00000000, 00000001, 00000010, 00000100, 00001000, 00010000, 00100000, 01000000, ...repeated...]
 * RESULT:
 * result = [F, T, T, T, T, T, T, T, ... repeated...]
 */
TEST_F(BoolFilterOther, Uint32ToBoolLong) {
	size_t const num_base(5);
	size_t const num_long(NUM_IN_LONG);
	uint32_t const data_base[] = { 0, 1, (1 << 1), (1 << 3), (1 << 8) };
	STATIC_ASSERT(ELEMENTSOF(data_base) == num_base);
	bool answer[] = { false, true, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == ELEMENTSOF(data_base));
	SIMD_ALIGN
	uint32_t data_long[num_long];
	SIMD_ALIGN
	bool result[ELEMENTSOF(data_long)];

	for (size_t i = 0; i < num_long; ++i) {
		data_long[i] = data_base[i % num_base];
	}

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;

	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
			<< num_long << endl;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_Uint32ToBool(num_long, data_long, result);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_long; ++i) {
		ASSERT_EQ(answer[i % num_base], result[i]);
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

/////////////////////////////////////////////////////////////////////////////////////////

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

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test failure cases of sakura_SetTrueFloatInRangesExclusive
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* lower_bound > upper_bound */
TEST_F(BoolFilterFloat, RangesExclusiveFailExchangeBounds) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
			num_data, in_data, num_range, upper, lower, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Null pointer arrays */
TEST_F(BoolFilterFloat, RangesExclusiveFailNullData) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
			num_data, data_null, num_range, lower, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, RangesExclusiveFailNullLower) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
			num_data, in_data, num_range, lower_null, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, RangesExclusiveFailNullUpper) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
			num_data, in_data, num_range, lower, upper_null, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, RangesExclusiveFailNullResult) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
			num_data, in_data, num_range, lower, upper, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterFloat, RangesExclusiveNotAlignedData) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
			num_data, data_shift, num_range, lower, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, RangesExclusiveNotAlignedResult) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
			num_data, data, num_range, lower, upper, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, RangesExclusiveNotAlignedLower) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
			num_data, data, num_range, lower_shift, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, RangesExclusiveNotAlignedUpper) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesExclusive(
			num_data, data, num_range, lower, upper_shift, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_SetTrueIntInRangesExclusive
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* lower_bound > upper_bound */
TEST_F(BoolFilterInt, RangesExclusiveFailExchangeBounds) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
			num_data, in_data, num_range, upper, lower, result);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Null pointer arrays */
TEST_F(BoolFilterInt, RangesExclusiveFailNullData) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
			num_data, data_null, num_range, lower, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, RangesExclusiveFailNullLower) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
			num_data, in_data, num_range, lower_null, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, RangesExclusiveFailNullUpper) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
			num_data, in_data, num_range, lower, upper_null, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, RangesExclusiveFailNullResult) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
			num_data, in_data, num_range, lower, upper, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterInt, RangesExclusiveNotAlignedData) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
			num_data, data_shift, num_range, lower, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, RangesExclusiveNotAlignedResult) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
			num_data, data, num_range, lower, upper, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, RangesExclusiveNotAlignedLower) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
			num_data, data, num_range, lower_shift, upper, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, RangesExclusiveNotAlignedUpper) {
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesExclusive(
			num_data, data, num_range, lower, upper_shift, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test failure cases of sakura_SetTrueFloatGreaterThan
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterFloat, GreaterThanFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	float *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThan(num_data,
			data_null, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, GreaterThanFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];

	bool *result_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(in_data), in_data);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThan(num_data,
			in_data, threshold_, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterFloat, GreaterThanNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	float data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	float *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThan(num_data,
			data_shift, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, GreaterThanNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	float data[num_data];
	SIMD_ALIGN
	bool result[num_elements];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThan(num_data,
			data, threshold_, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_SetTrueIntGreaterThan
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterInt, GreaterThanFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	int *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThan(num_data,
			data_null, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, GreaterThanFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];

	bool *result_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(in_data), in_data);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThan(num_data,
			in_data, threshold_, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterInt, GreaterThanNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	int data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	int *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThan(num_data,
			data_shift, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, GreaterThanNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	int data[num_data];
	SIMD_ALIGN
	bool result[num_elements];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThan(num_data,
			data, threshold_, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test failure cases of sakura_SetTrueFloatGreaterThanOrEquals
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterFloat, GreaterThanOrEqualsFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	float *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThanOrEquals(
			num_data, data_null, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, GreaterThanOrEqualsFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];

	bool *result_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(in_data), in_data);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThanOrEquals(
			num_data, in_data, threshold_, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterFloat, GreaterThanOrEqualsNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	float data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	float *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThanOrEquals(
			num_data, data_shift, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, GreaterThanOrEqualsNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	float data[num_data];
	SIMD_ALIGN
	bool result[num_elements];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatGreaterThanOrEquals(
			num_data, data, threshold_, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_SetTrueIntGreaterThanOrEquals
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterInt, GreaterThanOrEqualsFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	int *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThanOrEquals(
			num_data, data_null, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, GreaterThanOrEqualsFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];

	bool *result_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(in_data), in_data);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThanOrEquals(
			num_data, in_data, threshold_, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterInt, GreaterThanOrEqualsNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	int data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	int *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThanOrEquals(
			num_data, data_shift, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, GreaterThanOrEqualsNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	int data[num_data];
	SIMD_ALIGN
	bool result[num_elements];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntGreaterThanOrEquals(
			num_data, data, threshold_, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test failure cases of sakura_SetTrueFloatLessThan
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterFloat, LessThanFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	float *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThan(num_data,
			data_null, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, LessThanFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];

	bool *result_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(in_data), in_data);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThan(num_data,
			in_data, threshold_, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterFloat, LessThanNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	float data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	float *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThan(num_data,
			data_shift, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, LessThanNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	float data[num_data];
	SIMD_ALIGN
	bool result[num_elements];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThan(num_data,
			data, threshold_, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_SetTrueIntLessThan
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterInt, LessThanFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	int *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThan(num_data,
			data_null, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, LessThanFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];

	bool *result_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(in_data), in_data);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThan(num_data,
			in_data, threshold_, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterInt, LessThanNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	int data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	int *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThan(num_data,
			data_shift, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, LessThanNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	int data[num_data];
	SIMD_ALIGN
	bool result[num_elements];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThan(num_data, data,
			threshold_, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test failure cases of sakura_SetTrueFloatLessThanOrEquals
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterFloat, LessThanOrEqualsFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	float *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThanOrEquals(
			num_data, data_null, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, LessThanOrEqualsFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];

	bool *result_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(in_data), in_data);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThanOrEquals(
			num_data, in_data, threshold_, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterFloat, LessThanOrEqualsNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	float data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	float *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThanOrEquals(
			num_data, data_shift, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, LessThanOrEqualsNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	float data[num_data];
	SIMD_ALIGN
	bool result[num_elements];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatLessThanOrEquals(
			num_data, data, threshold_, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test failure cases of sakura_SetTrueIntLessThanOrEquals
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterInt, LessThanOrEqualsFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	int *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThanOrEquals(
			num_data, data_null, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, LessThanOrEqualsFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	int in_data[num_data];

	bool *result_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(in_data), in_data);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThanOrEquals(
			num_data, in_data, threshold_, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterInt, LessThanOrEqualsNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	int data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	int *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThanOrEquals(
			num_data, data_shift, threshold_, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterInt, LessThanOrEqualsNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	int data[num_data];
	SIMD_ALIGN
	bool result[num_elements];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntLessThanOrEquals(
			num_data, data, threshold_, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test failure cases of sakura_SetFalseFloatIfNanOrInf
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterFloat, NanOrInfFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	float *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseFloatIfNanOrInf(num_data,
			data_null, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, NanOrInfFailNullResult) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	float in_data[num_data];

	bool *result_null = nullptr;

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(in_data), in_data);

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseFloatIfNanOrInf(num_data,
			in_data, result_null);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/* Unaligned arrays */
TEST_F(BoolFilterFloat, NanOrInfNotAlignedData) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	float data[num_elements];
	SIMD_ALIGN
	bool result[num_data];
	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	float *data_shift = &data[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseFloatIfNanOrInf(num_data,
			data_shift, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterFloat, NanOrInfNotAlignedResult) {
	size_t offset(UNALIGN_OFFSET);
	size_t const num_data(NUM_IN);
	size_t const num_elements(num_data + offset);
	SIMD_ALIGN
	float data[num_data];
	SIMD_ALIGN
	bool result[num_elements];

	// Create long input data by repeating data_
	GetDataInLength(ELEMENTSOF(data), data);

	// Define unaligned array
	bool *result_shift = &result[offset];
	assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseFloatIfNanOrInf(num_data,
			data, result_shift);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/////////////////////////////////////////////////////////////////////////////////////////

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

/////////////////////////////////////////////////////////////////////////////////////////

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

//////////////// TEMPORARY ADDED FOR TEST ////////////////////////
//class BoolFilterUint8: public BoolFilter<uint8_t> {
//protected:
//	virtual void PrepareInputs() {
//		uint8_t const base_input[] = { 10, 5, 0, 5, 10, 15, 20, 15 };
//		STATIC_ASSERT(ELEMENTSOF(base_input) == NUM_IN);
//		STATIC_ASSERT(ELEMENTSOF(lower_) == 2);
//		STATIC_ASSERT(ELEMENTSOF(upper_) == 2);
//		lower_[0] = 3;
//		lower_[1] = 13;
//		upper_[0] = 7;
//		upper_[1] = 17;
//		for (size_t i = 0; i < NUM_IN; ++i) {
//			data_[i] = base_input[i];
//		}
//	}
//};
//
//TEST_F(BoolFilterUint8, RangesInclusive) {
//	size_t const num_data(NUM_IN);
//	SIMD_ALIGN
//	uint8_t in_data[num_data];
//	SIMD_ALIGN
//	bool result[ELEMENTSOF(in_data)];
//	SIMD_ALIGN
//	uint8_t lower[NUM_RANGE];
//	SIMD_ALIGN
//	uint8_t upper[ELEMENTSOF(lower)];
//	size_t const num_range(ELEMENTSOF(lower));
//
//	bool answer[] = { false, true, false, true, false, true, false, true };
//	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
//
//	// Create long input data by repeating data_
//	GetDataInLength(num_data, in_data);
//	// bounds
//	lower[0] = lower_[1];
//	lower[1] = lower_[0];
//	upper[0] = upper_[1];
//	upper[1] = upper_[0];
//
//	verbose=true;
//	if (verbose) {
//		PrintArray("data", num_data, in_data);
//		PrintArray("lower_bound", NUM_RANGE, lower);
//		PrintArray("upper_bound", NUM_RANGE, upper);
//		PrintArray("result (before)", num_data, result);
//	}
//
//	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueUint8InRangesInclusive(
//				num_data, in_data, num_range, lower, upper, result);
//	if (verbose)
//		PrintArray("result", num_data, result);
//
//	// Verification
//	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
//	for (size_t i = 0; i < num_data; ++i) {
//		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
//	}
//}
//
//TEST_F(BoolFilterUint8, RangesInclusiveLong) {
//	size_t const num_data(NUM_IN_LONG);
//	SIMD_ALIGN
//	uint8_t in_data[num_data];
//	SIMD_ALIGN
//	bool result[ELEMENTSOF(in_data)];
//	SIMD_ALIGN
//	uint8_t lower[NUM_RANGE];
//	SIMD_ALIGN
//	uint8_t upper[ELEMENTSOF(lower)];
//	size_t const num_range(ELEMENTSOF(lower));
//
//	bool answer[] = { false, true, false, true, false, true, false, true };
//	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
//
//	// Create long input data by repeating data_
//	GetDataInLength(num_data, in_data);
//	// bounds
//	lower[0] = lower_[0];
//	lower[1] = lower_[1];
//	upper[0] = upper_[0];
//	upper[1] = upper_[1];
//
//	double start, end;
//	size_t const num_repeat = 2000;
//	LIBSAKURA_SYMBOL(Status) status;
//
//	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
//			<< num_data << endl;
//	start = sakura_GetCurrentTime();
//	for (size_t i = 0; i < num_repeat; ++i) {
//		status = sakura_SetTrueUint8InRangesInclusive(
//				num_data, in_data, num_range, lower, upper, result);
//	}
//	end = sakura_GetCurrentTime();
//	cout << "Elapse time of actual operation: " << end - start << " sec"
//			<< endl;
//
//	// Verification
//	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
//	for (size_t i = 0; i < num_data; ++i) {
//		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
//	}
//}
//
//TEST_F(BoolFilterUint8, RangesInclusive8conditions) {
//	size_t const num_data(NUM_IN_LONG);
//	size_t const num_range(8);
//	SIMD_ALIGN
//	uint8_t in_data[num_data];
//	SIMD_ALIGN
//	bool result[ELEMENTSOF(in_data)];
//	SIMD_ALIGN
//	uint8_t lower[num_range];
//	SIMD_ALIGN
//	uint8_t upper[ELEMENTSOF(lower)];
//
//	// Create long input data by repeating data_
//	GetDataInLength(num_data, in_data);
//	// set bounds
//	for (size_t i=0; i < num_range-1; ++i){
//		lower[i] = 20 + (uint8_t) i;
//		upper[i] = 50 + (uint8_t) i;
//	}
//	lower[num_range-1] = lower_[0];
//	upper[num_range-1] = upper_[0];
//	//GetBounds(lower, upper);
//
//
//	double start, end;
//	size_t const num_repeat = 2000;
//	LIBSAKURA_SYMBOL(Status) status;
//
//	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
//			<< num_data << endl;
//	start = sakura_GetCurrentTime();
//	for (size_t i = 0; i < num_repeat; ++i) {
//		status = sakura_SetTrueUint8InRangesInclusive(
//				num_data, in_data, num_range, lower, upper, result);
//	}
//	end = sakura_GetCurrentTime();
//	cout << "Elapse time of actual operation: " << end - start << " sec"
//			<< endl;
//
//	// Verification
//	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
//}
//
//class BoolFilterUint32: public BoolFilter<uint32_t> {
//protected:
//	virtual void PrepareInputs() {
//		uint32_t const base_input[] = { 10, 5, 0, 5, 10, 15, 20, 15 };
//		STATIC_ASSERT(ELEMENTSOF(base_input) == NUM_IN);
//		STATIC_ASSERT(ELEMENTSOF(lower_) == 2);
//		STATIC_ASSERT(ELEMENTSOF(upper_) == 2);
//		lower_[0] = 3;
//		lower_[1] = 13;
//		upper_[0] = 7;
//		upper_[1] = 17;
//		for (size_t i = 0; i < NUM_IN; ++i) {
//			data_[i] = base_input[i];
//		}
//	}
//};
//
//TEST_F(BoolFilterUint32, RangesInclusive) {
//	size_t const num_data(NUM_IN);
//	SIMD_ALIGN
//	uint32_t in_data[num_data];
//	SIMD_ALIGN
//	bool result[ELEMENTSOF(in_data)];
//	SIMD_ALIGN
//	uint32_t lower[NUM_RANGE];
//	SIMD_ALIGN
//	uint32_t upper[ELEMENTSOF(lower)];
//	size_t const num_range(ELEMENTSOF(lower));
//
//	bool answer[] = { false, true, false, true, false, true, false, true };
//	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
//
//	// Create long input data by repeating data_
//	GetDataInLength(num_data, in_data);
//	// bounds
//	lower[0] = lower_[1];
//	lower[1] = lower_[0];
//	upper[0] = upper_[1];
//	upper[1] = upper_[0];
//
//	verbose=true;
//	if (verbose) {
//		PrintArray("data", num_data, in_data);
//		PrintArray("lower_bound", NUM_RANGE, lower);
//		PrintArray("upper_bound", NUM_RANGE, upper);
//		PrintArray("result (before)", num_data, result);
//	}
//
//	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueUint32InRangesInclusive(
//				num_data, in_data, num_range, lower, upper, result);
//	if (verbose)
//		PrintArray("result", num_data, result);
//
//	// Verification
//	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
//	for (size_t i = 0; i < num_data; ++i) {
//		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
//	}
//}
//
//TEST_F(BoolFilterUint32, RangesInclusiveLong) {
//	size_t const num_data(NUM_IN_LONG);
//	SIMD_ALIGN
//	uint32_t in_data[num_data];
//	SIMD_ALIGN
//	bool result[ELEMENTSOF(in_data)];
//	SIMD_ALIGN
//	uint32_t lower[NUM_RANGE];
//	SIMD_ALIGN
//	uint32_t upper[ELEMENTSOF(lower)];
//	size_t const num_range(ELEMENTSOF(lower));
//
//	bool answer[] = { false, true, false, true, false, true, false, true };
//	STATIC_ASSERT(ELEMENTSOF(answer) == NUM_IN);
//
//	// Create long input data by repeating data_
//	GetDataInLength(num_data, in_data);
//	// bounds
//	lower[0] = lower_[0];
//	lower[1] = lower_[1];
//	upper[0] = upper_[0];
//	upper[1] = upper_[1];
//
//	double start, end;
//	size_t const num_repeat = 2000;
//	LIBSAKURA_SYMBOL(Status) status;
//
//	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
//			<< num_data << endl;
//	start = sakura_GetCurrentTime();
//	for (size_t i = 0; i < num_repeat; ++i) {
//		status = sakura_SetTrueUint32InRangesInclusive(
//				num_data, in_data, num_range, lower, upper, result);
//	}
//	end = sakura_GetCurrentTime();
//	cout << "Elapse time of actual operation: " << end - start << " sec"
//			<< endl;
//
//	// Verification
//	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
//	for (size_t i = 0; i < num_data; ++i) {
//		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
//	}
//}
//
//TEST_F(BoolFilterUint32, RangesInclusive8conditions) {
//	size_t const num_data(NUM_IN_LONG);
//	size_t const num_range(8);
//	SIMD_ALIGN
//	uint32_t in_data[num_data];
//	SIMD_ALIGN
//	bool result[ELEMENTSOF(in_data)];
//	SIMD_ALIGN
//	uint32_t lower[num_range];
//	SIMD_ALIGN
//	uint32_t upper[ELEMENTSOF(lower)];
//
//	// Create long input data by repeating data_
//	GetDataInLength(num_data, in_data);
//	// set bounds
//	for (size_t i=0; i < num_range-1; ++i){
//		lower[i] = 20 + (uint32_t) i;
//		upper[i] = 50 + (uint32_t) i;
//	}
//	lower[num_range-1] = lower_[0];
//	upper[num_range-1] = upper_[0];
//	//GetBounds(lower, upper);
//
//	double start, end;
//	size_t const num_repeat = 2000;
//	LIBSAKURA_SYMBOL(Status) status;
//
//	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
//			<< num_data << endl;
//	start = sakura_GetCurrentTime();
//	for (size_t i = 0; i < num_repeat; ++i) {
//		status = sakura_SetTrueUint32InRangesInclusive(
//				num_data, in_data, num_range, lower, upper, result);
//	}
//	end = sakura_GetCurrentTime();
//	cout << "Elapse time of actual operation: " << end - start << " sec"
//			<< endl;
//
//	// Verification
//	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
//}
//
///*
// * Test bool filter generation sakura_SetTrueFloatInRangesInclusive
// * with a long array
// */
//TEST_F(BoolFilterFloat, RangesInclusive8conditions) {
//	size_t const num_data(NUM_IN_LONG);
//	size_t const num_range(8);
//	SIMD_ALIGN
//	float in_data[num_data];
//	SIMD_ALIGN
//	bool result[ELEMENTSOF(in_data)];
//	SIMD_ALIGN
//	float lower[num_range];
//	SIMD_ALIGN
//	float upper[ELEMENTSOF(lower)];
//
//	// Create long input data by repeating data_
//	GetDataInLength(num_data, in_data);
//	// set bounds
//	for (size_t i=0; i < num_range-1; ++i){
//		lower[i] = 20. + (float) i;
//		upper[i] = 50. + (float) i;
//	}
//	lower[num_range-1] = lower_[0];
//	upper[num_range-1] = upper_[0];
//	//GetBounds(lower, upper);
//
//
//	double start, end;
//	size_t const num_repeat = 2000;
//	LIBSAKURA_SYMBOL(Status) status;
//
//	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
//			<< num_data << endl;
//	start = sakura_GetCurrentTime();
//	for (size_t i = 0; i < num_repeat; ++i) {
//		status = sakura_SetTrueFloatInRangesInclusive(
//				num_data, in_data, num_range, lower, upper, result);
//	}
//	end = sakura_GetCurrentTime();
//	cout << "Elapse time of actual operation: " << end - start << " sec"
//			<< endl;
//
//	// Verification
//	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
//}
//
///*
// * Test bool filter generation sakura_SetTrueIntInRangesInclusive
// * with a long array
// */
//TEST_F(BoolFilterInt, RangesInclusive8conditions) {
//	size_t const num_data(NUM_IN_LONG);
//	size_t const num_range(8);
//	SIMD_ALIGN
//	int in_data[num_data];
//	SIMD_ALIGN
//	bool result[ELEMENTSOF(in_data)];
//	SIMD_ALIGN
//	int lower[num_range];
//	SIMD_ALIGN
//	int upper[ELEMENTSOF(lower)];
//
//	// Create long input data by repeating data_
//	GetDataInLength(num_data, in_data);
//	// set bounds
//	for (size_t i=0; i < num_range-1; ++i){
//		lower[i] = 20 + (int) i;
//		upper[i] = 50 + (int) i;
//	}
//	lower[num_range-1] = lower_[0];
//	upper[num_range-1] = upper_[0];
//	//GetBounds(lower, upper);
//
//
//	double start, end;
//	size_t const num_repeat = 2000;
//	LIBSAKURA_SYMBOL(Status) status;
//
//	cout << "Iterating " << num_repeat << " loops. The length of arrays is "
//			<< num_data << endl;
//	start = sakura_GetCurrentTime();
//	for (size_t i = 0; i < num_repeat; ++i) {
//		status = sakura_SetTrueIntInRangesInclusive(
//				num_data, in_data, num_range, lower, upper, result);
//	}
//	end = sakura_GetCurrentTime();
//	cout << "Elapse time of actual operation: " << end - start << " sec"
//			<< endl;
//
//	// Verification
//	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
//}

//////////////// END TEMPORARY ADDED FOR ////////////////////////
