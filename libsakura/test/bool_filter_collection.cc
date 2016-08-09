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
#include <cxxabi.h>
#include <iostream>
#include <math.h>
#include <string>
#include <sys/time.h>
#include <typeinfo>

#include <libsakura/localdef.h>
#include <libsakura/sakura.h>
#include "aligned_memory.h"
#include "loginit.h"
#include "gtest/gtest.h"
#include "testutil.h"

/* the number of elements in input/output array to test */
#define NUM_IN 9
#define NUM_RANGE 2 // DO NOT MODIFY THIS!
#define MAX_NUM_RANGE 17 // DO NOT MODIFY THIS!
#define NUM_IN_LONG (1 << 18) //2**18
#define UNALIGN_OFFSET 1 // should be != ALIGNMENT
#define POS_INF (1.0f / 0.0f) // positive infinity in float
#define NEG_INF (-1.0f / 0.0f) // netative infinity in float
#define NOT_A_NUM (0.0f / 0.0f) // NaN (not a number) in float
using namespace std;

/*
 * Test of generating a boolean filter using ranges
 * Inputs (float):
 * - data = [0., -0.5, -1.0, -0.5, 0., 0.5, 1.0, 0.5, 0.]
 * - ranges: lower = [-0.75, 0.5], upper = [-0.25, 0.75]
 * Inputs (int):
 * - data = [0, -5, -10, -5, 0, 5, 10, 5, 0]
 * - ranges: lower = [-7, 5], upper = [-3, 7]
 * RESULT:
 * - RangesInclusive  = [F, T, F, T, F, T, T, T, F]
 * - RangesExclusive = [F, T, F, T, F, F, F, F, F]
 * results should be [F, F, F, F, F, F, F, F, F] if lower or upper is empty
 */
// Generalize function type for boundary tests
template<typename DataType>
struct RangesFunction {
	typedef LIBSAKURA_SYMBOL(Status) (*RangesFuncType)(size_t, DataType const*,
			size_t, DataType const*, DataType const*, bool*);
	RangesFuncType function;
};
struct NumConditionAndAnswer {
	size_t num_condition;bool answer[NUM_IN];     // answer array
};
// a struct to store test parameters for both uint8 and uint32
struct RangesTestComponent {
	string name;               // name of operationGetNumTest()
	RangesFunction<float> funcfloat;
	RangesFunction<int> funcint;
	NumConditionAndAnswer num_condition_answer[MAX_NUM_RANGE + 1];
	static size_t GetNumTest() {
		return ELEMENTSOF(num_condition_answer);
	}
	;
};

RangesTestComponent RangesTestCase[] { { "Inclusive Ranges", { LIBSAKURA_SYMBOL(
		SetTrueIfInRangesInclusiveFloat) }, { LIBSAKURA_SYMBOL(
		SetTrueIfInRangesInclusiveInt) }, { { NUM_RANGE, { false, true, false,
true, false, true, false, true, false } }, { 0, { false, false, false, false,
false, false, false, false, false } }, { 1, { false, true, false,
true,
false, false, false, false, false } }, { 3, { false, true, false,
true, false, true, false, true, false } }, { 4, { false, true, false,
true, false, true, false, true, false } }, { 5, { false, true, false,
true, false, true, false, true, false } }, { 6, { false, true, false,
true, false, true, false, true, false } }, { 7, { false, true, false,
true, false, true, false, true, false } }, { 8, { false, true, false,
true, false, true, false, true, false } }, { 9, { false, true, false,
true, false, true, false, true, false } }, { 10, { false, true, false,
true, false, true, false, true, false } }, { 11, { false, true, false,
true, false, true, false, true, false } }, { 12, { false, true, false,
true, false, true, false, true, false } }, { 13, { false, true, false,
true, false, true, false, true, false } }, { 14, { false, true, false,
true, false, true, false, true, false } }, { 15, { false, true, false,
true, false, true, false, true, false } }, { 16, { false, true, false,
true, false, true, false, true, false } }, { 17, { false, true, false,
true, false, true, false, true, false } } } }, { "Exclusive Ranges", {
		LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveFloat) }, { LIBSAKURA_SYMBOL(
		SetTrueIfInRangesExclusiveInt) }, { { NUM_RANGE, { false, true, false,
true, false, false, false, false, false } }, { 0, { false, false, false, false,
false, false, false, false, false } }, { 1, { false, true, false, true,
false, false, false, false, false } }, { 3, { false, true, false,
true, false, false, false, false, false } }, { 4, { false, true, false,
true, false, false, false, false, false } }, { 5, { false, true, false,
true, false, false, false, false, false } }, { 6, { false, true, false,
true, false, false, false, false, false } }, { 7, { false, true, false,
true, false, false, false, false, false } }, { 8, { false, true, false,
true, false, false, false, false, false } }, { 9, { false, true, false,
true, false, false, false, false, false } }, { 10, { false, true, false,
true, false, false, false, false, false } }, { 11, { false, true, false,
true, false, false, false, false, false } }, { 12, { false, true, false,
true, false, false, false, false, false } }, { 13, { false, true, false,
true, false, false, false, false, false } }, { 14, { false, true, false,
true, false, false, false, false, false } }, { 15, { false, true, false,
true, false, false, false, false, false } }, { 16, { false, true, false,
true, false, false, false, false, false } }, { 17, { false, true, false,
true, false, false, false, false, false } } } }, };

/*
 * Helper functions to select proper ranges_func_ptr_t from
 *  BoundaryTestComponent by data type and return test kit consists of
 * - name of bit operation
 * - sakura function to test, number of ranges, and the answer for the given data type.
 */
template<typename DataType>
struct RangesTestHelper {
	typedef tuple<string, RangesFunction<DataType>, NumConditionAndAnswer> RangesTestKit;
	static RangesTestKit GetItem(size_t num_components,
			RangesTestComponent const *test_components, size_t i, size_t j);
};

template<>
RangesTestHelper<float>::RangesTestKit RangesTestHelper<float>::GetItem(
		size_t num_components, RangesTestComponent const *test_components,
		size_t i, size_t j) {
	assert(i < num_components);
	auto testcase = test_components[i];
	assert(j < testcase.GetNumTest());
	return RangesTestKit(testcase.name, testcase.funcfloat,
			testcase.num_condition_answer[j]);
}

template<>
RangesTestHelper<int>::RangesTestKit RangesTestHelper<int>::GetItem(
		size_t num_components, RangesTestComponent const *test_components,
		size_t i, size_t j) {
	assert(i < num_components);
	auto testcase = test_components[i];
	assert(j < testcase.GetNumTest());
	return RangesTestKit(testcase.name, testcase.funcint,
			testcase.num_condition_answer[j]);
}

/*
 * Test of generating a boolean filter using a threshold value
 * Inputs:
 * -  data = [0., -0.5, -1.0, -0.5, 0., 0.5, 1.0, 0.5, 0.] (float)
 *        or [0, -5, -10, -5, 0, 5, 10, 5, 0] (int)
 * - threshold_ = 0
 * RESULT:
 * - GreaterThan         = [F, F, F, F, F, T, T, T, F]
 * - GreaterThanOrEquals = [T, F, F, F, T, T, T, T, T]
 * - LessThan            = [F, T, T, T, F, F, F, F, F]
 * - LessThanOrEquals    = [T, T, T, T, T, F, F, F, T]
 */
// Generalize function type for boundary tests
template<typename DataType>
struct BoundaryFunction {
	typedef LIBSAKURA_SYMBOL(Status) (*BoundaryFuncType)(size_t,
			DataType const*, DataType, bool*);
	BoundaryFuncType function;
};
// a struct to store test parameters for both uint8 and uint32
struct BoundaryTestComponent {
	string name;               // name of operation
	BoundaryFunction<float> funcfloat;
	BoundaryFunction<int> funcint;bool answer[NUM_IN];     // answer array
};

BoundaryTestComponent BoundaryTestCase[] { { "GreaterThan", { LIBSAKURA_SYMBOL(
		SetTrueIfGreaterThanFloat) }, { LIBSAKURA_SYMBOL(
		SetTrueIfGreaterThanInt) }, { false, false, false, false, false, true,
true, true, false } }, { "GreaterThanOrEquals", { LIBSAKURA_SYMBOL(
		SetTrueIfGreaterThanOrEqualsFloat) }, { LIBSAKURA_SYMBOL(
		SetTrueIfGreaterThanOrEqualsInt) }, { true, false, false, false, true,
true, true, true, true } }, { "LessThan",
		{ LIBSAKURA_SYMBOL(SetTrueIfLessThanFloat) }, { LIBSAKURA_SYMBOL(
				SetTrueIfLessThanInt) }, {
		false, true, true, true, false, false, false, false, false } }, {
		"LessThanOrEquals",
		{ LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsFloat) }, {
				LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsInt) }, { true, true,
		true, true, true, false, false, false, true } } };

/*
 * Helper functions to select proper boundary_func_ptr_t from
 *  BoundaryTestComponent by data type and return test kit consists of
 * - name of bit operation
 * - sakura function to test and the answer for the given data type.
 */
template<typename DataType>
struct BoundaryTestHelper {
	typedef tuple<string, BoundaryFunction<DataType>, bool*> BoundaryTestKit;
	static BoundaryTestKit GetItem(size_t num_components,
			BoundaryTestComponent const *test_components, size_t i);
};

template<>
BoundaryTestHelper<float>::BoundaryTestKit BoundaryTestHelper<float>::GetItem(
		size_t num_components, BoundaryTestComponent const *test_components,
		size_t i) {
	assert(i < num_components);
	auto testcase = test_components[i];
//	cout << "answer = ";
//	for (size_t i = 0; i < ELEMENTSOF(testcase.answer); ++i) cout << testcase.answer[i];
//	cout << endl;
	return BoundaryTestKit(testcase.name, testcase.funcfloat, testcase.answer);
}

template<>
BoundaryTestHelper<int>::BoundaryTestKit BoundaryTestHelper<int>::GetItem(
		size_t num_components, BoundaryTestComponent const *test_components,
		size_t i) {
	assert(i < num_components);
	auto testcase = test_components[i];
	return BoundaryTestKit(testcase.name, testcase.funcint, testcase.answer);
}

/*
 * A Utility functions to support testing
 */
// Create arbitrary length of input data by repeating values of data_[]
template<typename DataType>
void GetDataInLength(size_t num_in, DataType *in_data, size_t num_out,
		DataType *out_data) {
	assert(in_data!=nullptr);
	assert(num_in != 0);
	// Handling of nullptr array
	if (out_data == nullptr)
		return;
	for (size_t i = 0; i < num_out; ++i) {
		out_data[i] = in_data[i % num_in];
	}
}

// Convert an array to a string formatted for printing.
template<typename DataType>
void PrintArray(char const *name, size_t num_data, DataType const *data_array) {
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
				cout << data_array[i] << ", ";
			cout << data_array[num_data - 1];
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
 * A super class to test various bool filter functions
 */
template<typename DataType>
class BoolFilter: public ::testing::Test {
protected:
	BoolFilter() :
			verbose_(false) {
		STATIC_ASSERT(ELEMENTSOF(data_)==NUM_IN);
		STATIC_ASSERT(ELEMENTSOF(lower_)==NUM_RANGE);
		STATIC_ASSERT(ELEMENTSOF(upper_)==NUM_RANGE);
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
	void RunRangesVariousLengthTest() {
		size_t const array_length[] = { 5, NUM_IN, 11, 16, 0 };
		size_t const num_test(ELEMENTSOF(array_length));
		size_t const num_max(16);
		//num_max = max(array_length);
		SIMD_ALIGN
		DataType data[num_max];
		SIMD_ALIGN
		bool result[ELEMENTSOF(data)];
		SIMD_ALIGN
		DataType lower[MAX_NUM_RANGE]; //[NUM_RANGE];
		SIMD_ALIGN
		DataType upper[ELEMENTSOF(lower)];
		size_t const num_range(ELEMENTSOF(lower));
		for (size_t irun = 0; irun < num_test; ++irun) {
			size_t const num_data(array_length[irun]);
			// Loop over sakura functions and number of conditions
			cout << "[Tests with array length = " << num_data << "]" << endl;
			RunRangesTest(num_data, data, num_range, lower, upper, result,
					ELEMENTSOF(RangesTestCase), RangesTestCase,
					LIBSAKURA_SYMBOL(Status_kOK), true,
					num_data == num_max ? MAX_NUM_RANGE + 1 : 1);
		}
	}

	void RunRangesPerformanceTest(size_t const num_long,
			size_t const num_repeat) {
		size_t const num_large(num_long);
		SIMD_ALIGN
		DataType data[num_large];
		SIMD_ALIGN
		bool result[ELEMENTSOF(data)];
		SIMD_ALIGN
		DataType lower[MAX_NUM_RANGE]; //[NUM_RANGE];
		SIMD_ALIGN
		DataType upper[ELEMENTSOF(lower)];
		size_t const num_range(ELEMENTSOF(lower));

		// Loop over sakura functions and number of conditions
		cout << "[Performance tests]" << endl;
		RunRangesTest(num_long, data, num_range, lower, upper, result,
				ELEMENTSOF(RangesTestCase), RangesTestCase,
				LIBSAKURA_SYMBOL(Status_kOK), true, 1,//MAX_NUM_RANGE + 1,
				num_repeat);
	}

	/*
	 * Failure cases of various bool filter functions using ranges.
	 * Testing null pointer array and unaligned arrays.
	 * RESULT:
	 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
	 */
	void RunRangesFailTest() {
		size_t offset(UNALIGN_OFFSET);
		size_t const num_data(NUM_IN);
		size_t const num_elements(num_data + offset);
		size_t const num_range_elements(NUM_RANGE + offset);
		SIMD_ALIGN
		DataType data[num_elements];
		SIMD_ALIGN
		bool result[ELEMENTSOF(data)];
		SIMD_ALIGN
		DataType lower[num_range_elements];
		SIMD_ALIGN
		DataType upper[ELEMENTSOF(lower)];
		size_t const num_range(ELEMENTSOF(lower));
		// For null array tests
		DataType *datatype_null = nullptr;
		bool *result_null = nullptr;
		// Define unaligned array
		DataType *data_shift = &data[offset];
		assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));
		bool *result_shift = &result[offset];
		assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

		// Null pointer array
		cout << "[Test NULL data array]" << endl;
		RunRangesTest(num_data, datatype_null, num_range, upper, lower, result,
				ELEMENTSOF(RangesTestCase), RangesTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		cout << "[Test NULL lowerr_bounds array]" << endl;
		RunRangesTest(num_data, data, num_range, datatype_null, upper, result,
				ELEMENTSOF(RangesTestCase), RangesTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		cout << "[Test NULL uppe_bounds array]" << endl;
		RunRangesTest(num_data, data, num_range, lower, datatype_null, result,
				ELEMENTSOF(RangesTestCase), RangesTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		cout << "[Test NULL result array]" << endl;
		RunRangesTest(num_data, data, num_range, lower, upper, result_null,
				ELEMENTSOF(RangesTestCase), RangesTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		// Unaligned array
		cout << "[Test unaligned data array]" << endl;
		RunRangesTest(num_data, data_shift, num_range, lower, upper, result,
				ELEMENTSOF(RangesTestCase), RangesTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		cout << "[Test unaligned lower_bounds array]" << endl;
		RunRangesTest(num_data, data, num_range, data_shift, upper, result,
				ELEMENTSOF(RangesTestCase), RangesTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		cout << "[Test unaligned upper_bounds array]" << endl;
		RunRangesTest(num_data, data, num_range, lower, data_shift, result,
				ELEMENTSOF(RangesTestCase), RangesTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		cout << "[Test unaligned result array]" << endl;
		RunRangesTest(num_data, data, num_range, lower, upper, result_shift,
				ELEMENTSOF(RangesTestCase), RangesTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
#ifndef NDEBUG
		// lower > upper
		cout << "[Test lower_bounds > upper_bounds]" << endl;
		GetDataOfLength(num_data, data);
		RunRangesTest(num_data, data, num_range, upper, lower, result,
				ELEMENTSOF(RangesTestCase), RangesTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument), false);
#endif
	}

	void RunBoundaryVariousLengthTest() {
		size_t const array_length[] = { 5, NUM_IN, 11, 16, 0 };
		size_t const num_test(ELEMENTSOF(array_length));
		size_t const num_max(16);
		//num_max = max(array_length);
		SIMD_ALIGN
		DataType data[num_max];
		SIMD_ALIGN
		bool result[ELEMENTSOF(data)];

		for (size_t irun = 0; irun < num_test; ++irun) {
			size_t const num_data(array_length[irun]);
			// Loop over sakura functions with a threshold
			cout << "[Tests with array length = " << num_data << "]" << endl;
			RunBoundaryTest(num_data, data, result,
					ELEMENTSOF(BoundaryTestCase), BoundaryTestCase,
					LIBSAKURA_SYMBOL(Status_kOK));
		}
	}

	// Performance test of various bool filter functions using a threshold value with a large array
	void RunBoundaryPerformanceTest(size_t const num_long,
			size_t const num_repeat) {
		assert(num_long > 0);
		assert(num_repeat > 0);
		size_t const num_large(num_long);
		SIMD_ALIGN
		DataType data[num_large];
		SIMD_ALIGN
		bool result[ELEMENTSOF(data)];

		// Loop over sakura functions with a threshold
		cout << "[Performance tests]" << endl;
		RunBoundaryTest(num_large, data, result, ELEMENTSOF(BoundaryTestCase),
				BoundaryTestCase, LIBSAKURA_SYMBOL(Status_kOK), num_repeat);
	}
	/*
	 * Failure cases of various bool filter functions using a threshold value
	 * Testing null pointer array and unaligned arrays.
	 * RESULT:
	 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
	 */
	void RunBoundaryFailTest() {
		size_t offset(UNALIGN_OFFSET);
		size_t const num_data(NUM_IN);
		size_t const num_elements(num_data + offset);
		SIMD_ALIGN
		DataType data[num_elements];
		SIMD_ALIGN
		bool result[ELEMENTSOF(data)];
		// For null array tests
		DataType *data_null = nullptr;
		bool *result_null = nullptr;
		// Define unaligned array
		DataType *data_shift = &data[offset];
		assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));
		bool *result_shift = &result[offset];
		assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

		// Null pointer array
		cout << "[Test NULL data array]" << endl;
		RunBoundaryTest(num_data, data_null, result,
				ELEMENTSOF(BoundaryTestCase), BoundaryTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		cout << "[Test NULL result array]" << endl;
		RunBoundaryTest(num_data, data, result_null,
				ELEMENTSOF(BoundaryTestCase), BoundaryTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		// Unaligned array
		cout << "[Test unaligned data array]" << endl;
		RunBoundaryTest(num_data, data_shift, result,
				ELEMENTSOF(BoundaryTestCase), BoundaryTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		cout << "[Test unaligned result array]" << endl;
		RunBoundaryTest(num_data, data, result_shift,
				ELEMENTSOF(BoundaryTestCase), BoundaryTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
	}

	void RunRangesTest(size_t num_data, DataType *in_data, size_t num_condition,
			DataType *lower_bounds, DataType *upper_bounds,
			bool *result, size_t num_operation,
			RangesTestComponent const *test_components,
			LIBSAKURA_SYMBOL(Status) return_value, bool initialize = true,
			size_t num_num_condition = 1, size_t num_repeat = 1) {
		// initialize input data only if asked
		if (initialize) {
			GetDataOfLength(num_data, in_data);
		}
		if (verbose_) {
			PrintArray("data", num_data, in_data);
		}
		// get data type name to identify benchmark tests
		int success = 0;
		string data_type_name = abi::__cxa_demangle(typeid(DataType).name(),
				nullptr, nullptr, &success);
		string_replace(data_type_name, " ", "_");

		if (num_repeat > 1)
			cout << "Iterating " << num_repeat
					<< " loops for each operation. The length of arrays is "
					<< num_data << endl;

		RangesTestHelper<DataType> my_testcase;
		for (size_t iop = 0; iop < num_operation; ++iop) {
			// make sure number of elements for variation of conditions in RangesTestComponent
			// is sufficient for iteration
			assert(num_num_condition <= test_components[iop].GetNumTest());
			for (size_t jc = 0; jc < num_num_condition; ++jc) {
				auto kit = my_testcase.GetItem(num_operation, test_components,
						iop, jc);
				NumConditionAndAnswer num_and_ans = std::get<2>(kit);
				cout << "Testing: " << std::get<0>(kit)
						<< " (number of ranges = " << num_and_ans.num_condition
						<< ")" << endl;
				if (initialize) {
					GetBounds(num_and_ans.num_condition, lower_bounds,
							upper_bounds);
				}
				if (verbose_) {
					PrintArray("lower_bound", num_and_ans.num_condition,
							lower_bounds);
					PrintArray("upper_bound", num_and_ans.num_condition,
							upper_bounds);
				}
				LIBSAKURA_SYMBOL(Status) status;
				double start = GetCurrentTime();
				for (size_t irun = 0; irun < num_repeat; ++irun) {
					// Actual execution of bit operation function
					status = (std::get<1>(kit).function)(num_data, in_data,
							num_and_ans.num_condition, lower_bounds,
							upper_bounds, result);
				} // end of num_repeat loop
				double end = GetCurrentTime();
				if (num_repeat > 1) {
					string test_name = std::get<0>(kit);
					string_replace(test_name, " ", "_");
					cout << "#x# benchmark BoolFilter_" << test_name << "_"
							<< num_and_ans.num_condition << "_"
							<< data_type_name << " " << end - start << endl;
				}

				if (verbose_) {
					if (status == LIBSAKURA_SYMBOL(Status_kOK))
						PrintArray("result", num_data, result);
					else
						cout << "sakura_Status = " << status << endl;
				}
				// Verification
				EXPECT_EQ(return_value, status);
				if (status == LIBSAKURA_SYMBOL(Status_kOK)) {
					auto *answer = num_and_ans.answer;
					size_t const num_answer = ELEMENTSOF(num_and_ans.answer);
					EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
					for (size_t i = 0; i < num_data; ++i) {
						ASSERT_EQ(answer[i % num_answer], result[i]);
					}
				}
			}
		} // end of operation loop
	}

	void RunBoundaryTest(size_t num_data, DataType *in_data, bool *result,
			size_t num_operation, BoundaryTestComponent const *test_components,
			LIBSAKURA_SYMBOL(Status) return_value, size_t num_repeat = 1) {

		GetDataOfLength(num_data, in_data);
		if (verbose_) {
			PrintArray("data", num_data, in_data);
			cout << "threshold = " << threshold_ << endl;
		}
		// get data type name to identify benchmark tests
		int success = 0;
		string data_type_name = abi::__cxa_demangle(typeid(DataType).name(),
				nullptr, nullptr, &success);
		string_replace(data_type_name, " ", "_");

		if (num_repeat > 1)
			cout << "Iterating " << num_repeat
					<< " loops for each operation. The length of arrays is "
					<< num_data << endl;

		BoundaryTestHelper<DataType> my_testcase;
		for (size_t iop = 0; iop < num_operation; ++iop) {
			auto kit = my_testcase.GetItem(num_operation, test_components, iop);
			LIBSAKURA_SYMBOL(Status) status;
			cout << "Testing: " << std::get<0>(kit) << endl;
			double start = GetCurrentTime();
			for (size_t irun = 0; irun < num_repeat; ++irun) {
				// Actual execution of bit operation function
				status = (std::get<1>(kit).function)(num_data, in_data,
						threshold_, result);
			} // end of num_repeat loop
			double end = GetCurrentTime();
			if (num_repeat > 1) {
				string test_name = std::get<0>(kit);
				string_replace(test_name, " ", "_");
				cout << "#x# benchmark BoolFilter_" << test_name << "_"
						<< data_type_name << " " << end - start << endl;
			}

			if (verbose_) {
				if (status == LIBSAKURA_SYMBOL(Status_kOK))
					PrintArray("result", num_data, result);
				else
					cout << "sakura_Status = " << status << endl;
			}
			// Verification
			EXPECT_EQ(return_value, status);
			if (status == LIBSAKURA_SYMBOL(Status_kOK)) {
				auto *answer = test_components[iop].answer; //std::get<2>(kit);
				size_t const num_answer = ELEMENTSOF(test_components[iop].answer);
//				PrintArray("answer", num_answer, answer);
//				PrintArray("kit answer", num_answer, std::get<2>(kit));
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
				for (size_t i = 0; i < num_data; ++i) {
					ASSERT_EQ(answer[i % num_answer], result[i]);
				}
			}
		} // end of operation loop
	}

	// Create arbitrary length of input data by repeating values of data_[]
	void GetDataOfLength(size_t num_out, DataType *out_data) {
		GetDataInLength(ELEMENTSOF(data_), data_, num_out, out_data);
	}
	// Copy values of lower and upper bounds to given arrays from lower_[] and upper_[]
	void GetBounds(size_t num_condition, DataType *lower, DataType *upper) {
		// Handling of nullptr array
		DataType dummy_bound[num_condition];
		DataType *lower_ptr = (lower != nullptr) ? lower : dummy_bound;
		DataType *upper_ptr = (upper != nullptr) ? upper : dummy_bound;
//		EXPECT_LE(num_condition, ELEMENTSOF(lower_ptr));
//		EXPECT_LE(num_condition, ELEMENTSOF(upper_ptr));
		size_t iswitch =
				num_condition > NUM_RANGE ? num_condition - NUM_RANGE : 0;
		size_t iend =
				num_condition > NUM_RANGE ?
						iswitch + NUM_RANGE : iswitch + num_condition;
		for (size_t i = 0; i < iswitch; ++i) {
			lower_ptr[i] = static_cast<DataType>(1000 + i * 3);
			upper_ptr[i] = static_cast<DataType>(1000 + i * 3 + 2);
		}
		for (size_t i = iswitch; i < iend; ++i) {
			lower_ptr[i] = lower_[i - iswitch];
			upper_ptr[i] = upper_[i - iswitch];
		}
	}

	// Member variables
	bool verbose_;
	static DataType data_[];
	static DataType upper_[];
	static DataType lower_[];
	static DataType threshold_;

};

/*
 * Tests various bool filter generation using float array
 * INPUTS:
 * - data = [0., -0.5, -1.0, -0.5, 0., 0.5, 1.0, 0.5]
 * - lower_bound = [-0.75, 0.5]
 * - upper_bound = [-0.25, 0.75]
 */
class BoolFilterFloat: public BoolFilter<float> {

};
template<>
float BoolFilter<float>::threshold_ = 0.0;
template<>
float BoolFilter<float>::data_[] = { 0.0, -0.5, -1.0, -0.5, 0., 0.5, 1.0, 0.5, 0.0 };
template<>
float BoolFilter<float>::lower_[] = { -0.75, 0.5 };
template<>
float BoolFilter<float>::upper_[] = { -0.25, 0.75 };

/*
 * Tests various bool filter generation using int array
 * INPUTS:
 * - data = [0, -5, -10, -5, 0, 5, 10, 5]
 * - lower_bound = [-7, 5]
 * - upper_bound = [-3, 7]
 */
class BoolFilterInt: public BoolFilter<int> {

};
template<>
int BoolFilter<int>::threshold_ = 0;
template<>
int BoolFilter<int>::data_[] = { 0, -5, -10, -5, 0, 5, 10, 5, 0 };
template<>
int BoolFilter<int>::lower_[] = { -7, 5 };
template<>
int BoolFilter<int>::upper_[] = { -3, 7 };

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Tests various simple bool filter generation based on input data array
 */
template<typename DataType>
struct SimpleTestComponent {
	typedef LIBSAKURA_SYMBOL(Status) (*FuncType)(size_t, DataType const*,
	bool*);
	string name;
	FuncType function;
	DataType data[NUM_IN];bool answer[NUM_IN];
};

template<typename DataType>
class BoolFilterSimple: public ::testing::Test {
protected:
	BoolFilterSimple() :
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

	void RunVariousLengthTest() {
		size_t const array_length[] = { 5, NUM_IN, 11, 16, 0 };
		size_t const num_test(ELEMENTSOF(array_length));
		size_t const num_max(16);
		//num_max = max(array_length);
		SIMD_ALIGN
		DataType data[num_max];
		SIMD_ALIGN
		bool result[ELEMENTSOF(data)];

		for (size_t irun = 0; irun < num_test; ++irun) {
			size_t const num_data(array_length[irun]);
			// Loop over sakura functions with a threshold
			cout << "[Tests with array length = " << num_data << "]" << endl;
			RunSimpleFilterTest<OutOfPlaceAction>(num_data, data, result,
					ELEMENTSOF(test_components), test_components,
					LIBSAKURA_SYMBOL(Status_kOK));
		}
	}

	// Performance test of various bool filter functions using a threshold value with a large array
	void RunPerformanceTest(size_t const num_long, size_t const num_repeat) {
		assert(num_long > 0);
		assert(num_repeat > 0);
		size_t const num_large(num_long);
		SIMD_ALIGN
		DataType data[num_large];
		SIMD_ALIGN
		bool result[ELEMENTSOF(data)];

		// Loop over sakura functions with a threshold
		cout << "[Performance tests]" << endl;
		RunSimpleFilterTest<OutOfPlaceAction>(num_large, data, result,
				ELEMENTSOF(test_components), test_components,
				LIBSAKURA_SYMBOL(Status_kOK), num_repeat);
	}

	// Test in-place (&out == &in) using an array of length 10
	// This is only supported for DataType == bool, defined as a specialization
	// of template later in this code. All the other data types should fail.
	void RunInPlaceTests() {
		assert(false);
	}

	/*
	 * Failure cases of various simple bool filter functions
	 * Testing null pointer array and unaligned arrays.
	 * RESULT:
	 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
	 */
	void RunFailTest() {
		size_t offset(UNALIGN_OFFSET);
		size_t const num_data(NUM_IN);
		size_t const num_elements(num_data + offset);
		SIMD_ALIGN
		DataType data[num_elements];
		SIMD_ALIGN
		bool result[ELEMENTSOF(data)];
		// For null array tests
		DataType *data_null = nullptr;
		bool *result_null = nullptr;
		// Define unaligned array
		DataType *data_shift = &data[offset];
		assert(! LIBSAKURA_SYMBOL(IsAligned)(data_shift));
		bool *result_shift = &result[offset];
		assert(! LIBSAKURA_SYMBOL(IsAligned)(result_shift));

		// Null pointer array
		cout << "[Test NULL data array]" << endl;
		RunSimpleFilterTest<OutOfPlaceAction>(num_data, data_null, result,
				ELEMENTSOF(test_components), test_components,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		cout << "[Test NULL result array]" << endl;
		RunSimpleFilterTest<OutOfPlaceAction>(num_data, data, result_null,
				ELEMENTSOF(test_components), test_components,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		// Unaligned array
		cout << "[Test unaligned data array]" << endl;
		RunSimpleFilterTest<OutOfPlaceAction>(num_data, data_shift, result,
				ELEMENTSOF(test_components), test_components,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		cout << "[Test unaligned result array]" << endl;
		RunSimpleFilterTest<OutOfPlaceAction>(num_data, data, result_shift,
				ELEMENTSOF(test_components), test_components,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
	}

	template<typename InitializeAction>
	void RunSimpleFilterTest(size_t num_data, DataType *in_data, bool *result,
			size_t num_operation,
			SimpleTestComponent<DataType> const *test_components,
			LIBSAKURA_SYMBOL(Status) return_value, size_t num_repeat = 1) {

		// get data type name to identify benchmark tests
		int success = 0;
		string data_type_name = abi::__cxa_demangle(typeid(DataType).name(),
				nullptr, nullptr, &success);
		string_replace(data_type_name, " ", "_");

		if (num_repeat > 1)
			cout << "Iterating " << num_repeat
					<< " loops for each operation. The length of arrays is "
					<< num_data << endl;

		SimpleTestComponent<DataType> test_kit;
		for (size_t iop = 0; iop < num_operation; ++iop) {
			test_kit = test_components[iop];
			LIBSAKURA_SYMBOL(Status) status;
			cout << "Testing: " << test_kit.name << endl;
			GetDataInLength(ELEMENTSOF(test_kit.data), test_kit.data, num_data,
					in_data);
			if (verbose_) {
				PrintArray("data", num_data, in_data);
			}
			double start = GetCurrentTime();
			for (size_t irun = 0; irun < num_repeat; ++irun) {
				// Need to refresh data for in-place operation
				InitializeAction::reinitialize(ELEMENTSOF(test_kit.data),
						test_kit.data, num_data, in_data);
				// Actual execution of bit operation function
				status = (test_kit.function)(num_data, in_data, result);
			} // end of num_repeat loop
			double end = GetCurrentTime();
			if (num_repeat > 1) {
				string test_name = test_kit.name;
				string_replace(test_name, " ", "_");
				cout << "#x# benchmark BoolFilter_" << test_name << "_"
						<< data_type_name << " " << end - start << endl;
			}

			if (verbose_) {
				if (status == LIBSAKURA_SYMBOL(Status_kOK))
					PrintArray("result", num_data, result);
				else
					cout << "sakura_Status = " << status << endl;
			}
			// Verification
			EXPECT_EQ(return_value, status);
			if (status == LIBSAKURA_SYMBOL(Status_kOK)) {
				auto *answer = test_kit.answer;
				size_t const num_answer = ELEMENTSOF(test_kit.answer);
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
				for (size_t i = 0; i < num_data; ++i) {
					ASSERT_EQ(answer[i % num_answer], result[i]);
				}
			}
		} // end of operation loop
	}
	// Re-initialization of input data and mask for the case of in-place operation.
	struct InPlaceAction {
		static void reinitialize(size_t num_base, DataType *base_data,
				size_t num_data, DataType *data) {
			GetDataInLength(num_base, base_data, num_data, data);
		}
	};
	// Dummy function (does nothing) for re-initialization for the case of out-of-place operation.
	struct OutOfPlaceAction {
		static void reinitialize(size_t num_base, DataType *base_data,
				size_t num_data, DataType *data) {
			// no need to initialize
		}
	};

	bool verbose_;
	static SimpleTestComponent<DataType> test_components[];
};

/*
 * Test cases of a float array
 * NanOrInf
 * - INPUT: [ 0, inf, -1, -0.5, -inf, nan, 1, 0.5, 0, ... repeated ... ]
 * - RESULT: [T, F, T, T, F, F, T, T, T, ... repeated ... ]
 */
class BoolFilterSimpleFloat: public BoolFilterSimple<float> {
};
template<>
SimpleTestComponent<float> BoolFilterSimple<float>::test_components[] = { {
		"NanOrInf", LIBSAKURA_SYMBOL(SetFalseIfNanOrInfFloat), { 0.0, POS_INF,
				-1.0, -0.5, NEG_INF, NOT_A_NUM, 1.0, 0.5, 0.0 }, { true, false, true,
		true, false, false, true, true, true } } };

/*
 * Test cases of an uint8_t array
 * Uint8ToBool
 * - INPUT: [00000000, 00000001, 00000010, 00000100, 00001000, 00010000, 00100000, 01000000, 10000000, ... repeated ... ]
 * - RESULT: [F, T, T, T, T, T, T, T, F, ... repeated ... ]
 */
class BoolFilterSimpleUint8: public BoolFilterSimple<uint8_t> {
};
template<>
SimpleTestComponent<uint8_t> BoolFilterSimple<uint8_t>::test_components[] = { {
		"Uint8ToBool", LIBSAKURA_SYMBOL(Uint8ToBool), { 0, 1, 2, 4, 8, 16, 32,
				64, 128 }, { false, true, true, true, true, true, true, true, true } } };

/*
 * Test cases of an uint32_t array
 * Uint32ToBool
 * - INPUT: [0...000, 0...001, 0...010, 0...01000, (1<<7), (1<<11), (1<<16), (1<<24), (1<<30), ... repeated ... ]
 * - RESULT: [F, T, T, T, T, T, T, T, ... repeated ... ]
 */
class BoolFilterSimpleUint32: public BoolFilterSimple<uint32_t> {
};
template<>
SimpleTestComponent<uint32_t> BoolFilterSimple<uint32_t>::test_components[] = {
		{ "Uint32ToBool", LIBSAKURA_SYMBOL(Uint32ToBool), { 0, 1, (1 << 1), (1
				<< 3), (1 << 7), (1 << 11), (1 << 16), (1 << 24), (1 << 30) }, { false,
		true, true, true, true, true, true,
		true, true } } };

/*
 * Test cases of a bool array
 * InvertBool
 * - INPUT: [F, T, F, F, T, T, F, F, T, ... repeated ... ]
 * - RESULT: [T, F, T, T, F, F, T, T, F, ... repeated ... ]
 */
class BoolFilterSimpleBool: public BoolFilterSimple<bool> {
};
template<>
SimpleTestComponent<bool> BoolFilterSimple<bool>::test_components[] = { {
		"InvertBool", LIBSAKURA_SYMBOL(InvertBool), { false, true, false, false,
		true, true, false, false, true }, { true, false, true, true, false,
		false, true, true, false } } };
// Test in-place (&out == &in) using an array of length 10
template<>
void BoolFilterSimple<bool>::RunInPlaceTests() {
	size_t const num_data(10);
	SIMD_ALIGN
	bool data[num_data];

	// Loop over minimum set of operation types  (Standard tests)
	cout << "[In-place bit operations]" << endl;
	RunSimpleFilterTest<InPlaceAction>(num_data, data, data,
			ELEMENTSOF(test_components), test_components,
			LIBSAKURA_SYMBOL(Status_kOK));
}

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bool filter generation using ranges for an array of length 8, 11 and zero
 */
TEST_F(BoolFilterFloat, RangesVariousLength) {
	RunRangesVariousLengthTest();
}

TEST_F(BoolFilterInt, RangesVariousLength) {
	RunRangesVariousLengthTest();
}

/*
 * Test bool filter generation using ranges for an array for a large array
 */
TEST_F(BoolFilterFloat, RangesPerformance) {
//	RunRangesPerformanceTest(NUM_IN_LONG, 2000);
	RunRangesPerformanceTest(20007, 20000);
}

TEST_F(BoolFilterInt, RangesPerformance) {
//	RunRangesPerformanceTest(NUM_IN_LONG, 2000);
	RunRangesPerformanceTest(20007, 20000);
}
/*
 * Test failure cases of bool filter generation using ranges
 */
TEST_F(BoolFilterFloat, RangesFail) {
	RunRangesFailTest();
}

TEST_F(BoolFilterInt, RangesFail) {
	RunRangesFailTest();
}

/////////////////////////////////////////////////////////////////////////////////////////
/*
 * Test bool filter generation using a threshold value for an array of length 8, 11 and zero
 */
TEST_F(BoolFilterFloat, BoundaryVariousLength) {
	RunBoundaryVariousLengthTest();
}

TEST_F(BoolFilterInt, BoundaryVariousLength) {
	RunBoundaryVariousLengthTest();
}

/*
 * Test bool filter generation using a boundary value for a large array
 */
TEST_F(BoolFilterFloat, BoundaryPerformance) {
	RunBoundaryPerformanceTest(NUM_IN_LONG, 20000); // array_length, repeat
}

TEST_F(BoolFilterInt, BoundaryPerformance) {
	RunBoundaryPerformanceTest(NUM_IN_LONG, 20000); // array_length, repeat
}

/*
 * Test failure case of bool filter generation using a boundary value
 */
TEST_F(BoolFilterFloat, BoundaryFail) {
	RunBoundaryFailTest();
}

TEST_F(BoolFilterInt, BoundaryFail) {
	RunBoundaryFailTest();
}

///////////////////////////////////////////////////////////////////////////////////////
/* Simple filtering by NanOrInf */

/* Input array of length 8, 11 and zero */
TEST_F(BoolFilterSimpleFloat, NanOrInfVariousLength) {
	RunVariousLengthTest();
}

/* Performance test */
TEST_F(BoolFilterSimpleFloat, NanOrInfPerformance) {
	RunPerformanceTest(NUM_IN_LONG, 20000);
}

/* Failure cases */
TEST_F(BoolFilterSimpleFloat, NanOrInfFail) {
	RunFailTest();
}

/////////////////////////////////////////////////////////////////////////////////////////
/* Simple filtering by inverting a bool array */

/* Input array of length 8, 11 and zero */
TEST_F(BoolFilterSimpleBool, InvertBoolVariousLength) {
	RunVariousLengthTest();
}

/* Performance test */
TEST_F(BoolFilterSimpleBool, InvertBoolPerformance) {
	RunPerformanceTest(NUM_IN_LONG, 100000);
}

/* In place operation */
TEST_F(BoolFilterSimpleBool, InvertBoolInPlace) {
	RunInPlaceTests();
}

/* Failure cases */
TEST_F(BoolFilterSimpleBool, InvertBoolFail) {
	RunFailTest();
}
/////////////////////////////////////////////////////////////////////////////////////////
/* Simple filtering by inverting an uint8_t array */

/* Input array of length 8, 11 and zero */
TEST_F(BoolFilterSimpleUint8, Uint8ToBoolVariousLength) {
	RunVariousLengthTest();
}

/* Performance test */
TEST_F(BoolFilterSimpleUint8, Uint8ToBoolPerformance) {
	RunPerformanceTest(NUM_IN_LONG, 100000);
}

/* Failure cases */
TEST_F(BoolFilterSimpleUint8, Uint8ToBoolBoolFail) {
	RunFailTest();
}
/////////////////////////////////////////////////////////////////////////////////////////
/* Simple filtering by inverting an uint32_t array */

/* Input array of length 8, 11 and zero */
TEST_F(BoolFilterSimpleUint32, Uint32ToBoolVariousLength) {
	RunVariousLengthTest();
}

/* Performance test */
TEST_F(BoolFilterSimpleUint32, Uint32ToBoolPerformance) {
	RunPerformanceTest(NUM_IN_LONG, 20000);
}

/* Failure cases */
TEST_F(BoolFilterSimpleUint32, Uint32ToBoolBoolFail) {
	RunFailTest();
}
