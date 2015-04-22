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
#include <libsakura/sakura.h>

#include <iostream>
#include <math.h>
#include <string>
#include <sys/time.h>

#include <libsakura/localdef.h>
#include "aligned_memory.h"
#include "loginit.h"
#include "gtest/gtest.h"

/* the number of elements in input/output array to test */
#define NUM_IN 8
#define NUM_RANGE 2 // DO NOT MODIFY THIS!
#define NUM_IN_LONG (1 << 18) //2**18
#define UNALIGN_OFFSET 1 // should be != ALIGNMENT
using namespace std;

/*
 * Test of generating a boolean filter using ranges
 * Inputs (float):
 * - data = [0., -0.5, -1.0, -0.5, 0., 0.5, 1.0, 0.5]
 * - ranges: lower = [-0.75, 0.5], upper = [-0.25, 0.75]
 * Inputs (int):
 * - data = [0, -5, -10, -5, 0, 5, 10, 5]
 * - ranges: lower = [-7, 5], upper = [-3, 7]
 * RESULT:
 * - RangesInclusive  = [F, T, F, T, F, T, T, T]
 * - RangesExclusive = [F, T, F, T, F, F, F, F]
 * results should be [F, F, F, F, F, F, F, F] if lower or upper is empty
 */
// Generalize function type for boundary tests
template<typename DataType>
struct RangesFunction {
	typedef LIBSAKURA_SYMBOL(Status) (*RangesFuncType)(size_t,
			DataType const*, size_t, DataType const*, DataType const*, bool*);
	RangesFuncType function;
};
struct NumConditionAndAnswer {
	size_t num_condition;
	bool answer[NUM_IN];     // answer array
};
// a struct to store test parameters for both uint8 and uint32
struct RangesTestComponent {
	string name;               // name of operation
	RangesFunction<float> funcfloat;
	RangesFunction<int> funcint;
	NumConditionAndAnswer num_condition_answer[2];
	static size_t GetNumTest() {
		return ELEMENTSOF(num_condition_answer);
	}
	;
};

RangesTestComponent RangesTestCase[] {
		{"Inclusive Ranges",
				{ LIBSAKURA_SYMBOL(SetTrueIfInRangesInclusiveFloat) },
				{ LIBSAKURA_SYMBOL(SetTrueIfInRangesInclusiveInt) },
				{ { NUM_RANGE, { false, true, false, true, false, true, false, true } },
					{ 0, { false, false, false, false, false, false, false, false } } } },
		{"Exclusive Ranges",
				{ LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveFloat) },
						{ LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveInt) },
						{ {NUM_RANGE, { false, true, false, true, false, false, false, false }},
							{0, { false, false, false, false, false, false, false, false } } }},
};

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
 * -  data = [0., -0.5, -1.0, -0.5, 0., 0.5, 1.0, 0.5] (float)
 *        or [0, -5, -10, -5, 0, 5, 10, 5] (int)
 * - threshold_ = 0
 * RESULT:
 * - GreaterThan         = [F, F, F, F, F, T, T, T]
 * - GreaterThanOrEquals = [T, F, F, F, T, T, T, T]
 * - LessThan            = [F, T, T, T, F, F, F, F]
 * - LessThanOrEquals    = [T, T, T, T, T, F, F, F]
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
	BoundaryFunction<int> funcint;
	bool answer[NUM_IN];     // answer array
};

BoundaryTestComponent BoundaryTestCase[] {
		{"GreaterThan",
				{ LIBSAKURA_SYMBOL(SetTrueIfGreaterThanFloat) },
				{ LIBSAKURA_SYMBOL(SetTrueIfGreaterThanInt) },
				{ false, false, false, false, false, true, true, true } },
		{"GreaterThanOrEquals",
				{ LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsFloat) },
						{ LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsInt) },
				{ true, false, false, false, true, true, true, true } },
		{"LessThan",
				{ LIBSAKURA_SYMBOL(SetTrueIfLessThanFloat) },
				{ LIBSAKURA_SYMBOL(SetTrueIfLessThanInt) },
				{ false, true, true, true, false, false, false, false } },
		{"LessThanOrEquals",
				{ LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsFloat) },
				{ LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsInt) },
				{ true, true, true, true, true, false, false, false } }
};

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
 * A super class to test various bit operation of an value and array
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
		size_t const array_length[] = { NUM_IN, 11, 0 };
		size_t const num_test(ELEMENTSOF(array_length));
		size_t const num_max(11);
		//num_max = max(array_length);
		SIMD_ALIGN
		DataType data[num_max];
		SIMD_ALIGN
		bool result[ELEMENTSOF(data)];
		SIMD_ALIGN
		DataType lower[NUM_RANGE];
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
					num_data == NUM_IN ? 2 : 1);
		}
	}

	void RunRangesPerformanceTest(size_t const num_long, size_t const num_repeat) {
		size_t const num_large(num_long);
		SIMD_ALIGN
		DataType data[num_large];
		SIMD_ALIGN
		bool result[ELEMENTSOF(data)];
		SIMD_ALIGN
		DataType lower[NUM_RANGE];
		SIMD_ALIGN
		DataType upper[ELEMENTSOF(lower)];
		size_t const num_range(ELEMENTSOF(lower));

		// Loop over sakura functions and number of conditions
		cout << "[Performance tests]" << endl;
		RunRangesTest(num_long, data, num_range, lower, upper, result,
				ELEMENTSOF(RangesTestCase), RangesTestCase,
				LIBSAKURA_SYMBOL(Status_kOK), true, 1, num_repeat);
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
		GetDataInLength(num_data, data);
		RunRangesTest(num_data, data, num_range, upper, lower, result,
				ELEMENTSOF(RangesTestCase), RangesTestCase,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument), false);
#endif
	}

	void RunBoundaryVariousLengthTest() {
		size_t const array_length[] = { NUM_IN, 11, 0 };
		size_t const num_test(ELEMENTSOF(array_length));
		size_t const num_max(11);
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
	void RunBoundaryPerformanceTest(size_t const num_long, size_t const num_repeat) {
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
			GetDataInLength(num_data, in_data);
			GetBounds(lower_bounds, upper_bounds);
		}
		if (verbose_) {
			PrintArray("data", num_data, in_data);
			PrintArray("lower_bound", num_condition, lower_bounds);
			PrintArray("upper_bound", num_condition, upper_bounds);
		}

		if (num_repeat > 1)
			cout << "Iterating " << num_repeat
					<< " loops for each operation. The length of arrays is "
					<< num_data << endl;

		RangesTestHelper<DataType> my_testcase;
		for (size_t iop = 0; iop < num_operation; ++iop) {
			size_t const num_test = test_components[iop].GetNumTest();
			// make sure number of elements for variation of conditions in RangesTestComponent
			// is sufficient for iteration
			assert(num_num_condition <= num_test);
			for (size_t jc = 0; jc < num_num_condition; ++jc) {
				auto kit = my_testcase.GetItem(num_operation, test_components,
						iop, jc);
				NumConditionAndAnswer num_and_ans = std::get<2>(kit);
				cout << "Testing: " << std::get<0>(kit)
						<< " (number of ranges = " << num_and_ans.num_condition
						<< ")" << endl;
				LIBSAKURA_SYMBOL(Status) status;
				double start = LIBSAKURA_SYMBOL(GetCurrentTime)();
				for (size_t irun = 0; irun < num_repeat; ++irun) {
					// Actual execution of bit operation function
					status = (std::get<1>(kit).function)(num_data, in_data,
							num_and_ans.num_condition, lower_bounds,
							upper_bounds, result);
				} // end of num_repeat loop
				double end = LIBSAKURA_SYMBOL(GetCurrentTime)();
				if (num_repeat > 1)
					cout << "Elapse time of operation: " << end - start
							<< " sec" << endl;

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

		GetDataInLength(num_data, in_data);
		if (verbose_) {
			PrintArray("data", num_data, in_data);
			cout << "threshold = " << threshold_ << endl;
		}

		if (num_repeat > 1)
			cout << "Iterating " << num_repeat
					<< " loops for each operation. The length of arrays is "
					<< num_data << endl;

		BoundaryTestHelper<DataType> my_testcase;
		for (size_t iop = 0; iop < num_operation; ++iop) {
			auto kit = my_testcase.GetItem(num_operation, test_components, iop);
			LIBSAKURA_SYMBOL(Status) status;
			cout << "Testing: " << std::get<0>(kit) << endl;
			double start = LIBSAKURA_SYMBOL(GetCurrentTime)();
			for (size_t irun = 0; irun < num_repeat; ++irun) {
				// Actual execution of bit operation function
				status = (std::get<1>(kit).function)(num_data, in_data,
						threshold_, result);
			} // end of num_repeat loop
			double end = LIBSAKURA_SYMBOL(GetCurrentTime)();
			if (num_repeat > 1)
				cout << "Elapse time of operation: " << end - start << " sec"
						<< endl;

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
				size_t const num_answer = ELEMENTSOF(std::get<2>(kit));
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
	void GetDataInLength(size_t num_out, DataType *out_data) {
		// Handling of nullptr array
		if (out_data==nullptr) return;
		size_t const num_data(ELEMENTSOF(data_));
		for (size_t i = 0; i < num_out; ++i) {
			out_data[i] = data_[i % num_data];
		}
	}
	// Copy values of lower and upper bounds to given arrays from lower_[] and upper_[]
	void GetBounds(DataType *lower, DataType *upper) {
		// Handling of nullptr array
		DataType dummy_bound[NUM_RANGE];
		DataType *lower_ptr = (lower != nullptr) ? lower : dummy_bound;
		DataType *upper_ptr = (upper != nullptr) ? upper : dummy_bound;
		EXPECT_EQ(NUM_RANGE, ELEMENTSOF(lower_ptr));
		EXPECT_EQ(NUM_RANGE, ELEMENTSOF(upper_ptr));
		for (size_t i = 0; i < NUM_RANGE; ++i) {
			lower_ptr[i] = lower_[i];
			upper_ptr[i] = upper_[i];
		}
	}

	// Convert an array to a string formatted for printing.
	void PrintArray(char const *name, size_t num_data,
			DataType const *data_array) {
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
float BoolFilter<float>::data_[] = { 0.0, -0.5, -1.0, -0.5, 0., 0.5, 1.0, 0.5 };
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
int BoolFilter<int>::data_[] = { 0, -5, -10, -5, 0, 5, 10, 5 };
template<>
int BoolFilter<int>::lower_[] = { -7, 5 };
template<>
int BoolFilter<int>::upper_[] = { -3, 7 };

/*
 * Tests various simple bool filter generation based on input data array
 */
class BoolFilterSimple: public BoolFilter<float> {

};

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
	RunRangesPerformanceTest(NUM_IN_LONG, 2000);
}

TEST_F(BoolFilterInt, RangesPerformance) {
	RunRangesPerformanceTest(NUM_IN_LONG, 2000);
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

/////////////////////////////////////////////////////////////////////////////////////////

/*
 * Test bool filter generation sakura_SetFalseIfNanOrInfFloat
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

	if (verbose_) {
		PrintArray("data", num_data, in_data);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseIfNanOrInfFloat(num_data,
			in_data, result);

	if (verbose_)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetFalseIfNanOrInfFloat
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

	if (verbose_) {
		PrintArray("data", num_data, in_data);
		cout << "threshold = " << threshold_ << endl;
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseIfNanOrInfFloat(num_data,
			in_data, result);

	if (verbose_)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i % ELEMENTSOF(answer)], result[i]);
	}
}

/*
 * Test bool filter generation sakura_SetFalseIfNanOrInfFloat
 * with a long array
 * RESULT:
 * result = [T, T, F, F, F, T, T, T, ...]
 */
TEST_F(BoolFilterFloat, NanOrInfPerformance) {
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
		status = sakura_SetFalseIfNanOrInfFloat(num_data, in_data, result);
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
 * Test bool filter generation sakura_SetFalseIfNanOrInfFloat
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseIfNanOrInfFloat(num_data,
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
TEST_F(BoolFilterSimple, InvertBool) {
	size_t const num_data(4);
	SIMD_ALIGN
	bool const data[] = { true, false, false, true };
	STATIC_ASSERT(ELEMENTSOF(data) == num_data);
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];
	bool answer[] = { false, true, true, false };
	STATIC_ASSERT(ELEMENTSOF(answer) == num_data);

	if (verbose_)
		PrintArray("data", num_data, data);

	LIBSAKURA_SYMBOL(Status) status = sakura_InvertBool(num_data, data, result);

	if (verbose_)
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
TEST_F(BoolFilterSimple, InvertBoolPerformance) {
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
TEST_F(BoolFilterSimple, InvertBoolLenghZero) {
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
TEST_F(BoolFilterSimple, Uint8ToBool) {
	size_t const num_data(8);
	SIMD_ALIGN
	uint8_t const data[] = { 0, 1, 2, 4, 8, 16, 32, 64 };
	STATIC_ASSERT(ELEMENTSOF(data) == num_data);
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];
	bool answer[] = { false, true, true, true, true, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == num_data);

//	if (verbose_) PrintArray("data", num_data, in);

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint8ToBool(num_data, data,
			result);

	if (verbose_)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i], result[i]);
	}
}

/*
 * Test bool filter generation sakura_Uint8ToBool  (Performance Test)
 * INPUT:
 * in = [00000000, 00000001, 00000010, 00000100, 00001000, 00010000, 00100000, 01000000, ...repeated...]
 * RESULT:
 * result = [F, T, T, T, T, T, T, T, ... repeated...]
 */
TEST_F(BoolFilterSimple, Uint8ToBoolPerformance) {
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
TEST_F(BoolFilterSimple, Uint8ToBoolLenghZero) {
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
TEST_F(BoolFilterSimple, Uint32ToBool) {
	size_t const num_data(5);
	SIMD_ALIGN
	uint32_t const data[] = { 0, 1, (1 << 1), (1 << 3), (1 << 8) };
	STATIC_ASSERT(ELEMENTSOF(data) == num_data);
	SIMD_ALIGN
	bool result[ELEMENTSOF(data)];
	bool answer[] = { false, true, true, true, true };
	STATIC_ASSERT(ELEMENTSOF(answer) == num_data);

//	if (verbose_) PrintArray("data", num_data, in);

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint32ToBool(num_data, data,
			result);

	if (verbose_)
		PrintArray("result", num_data, result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i], result[i]);
	}
}

/*
 * Test bool filter generation sakura_Uint8ToBool  (Performance Test)
 * INPUT:
 * in = [00000000, 00000001, 00000010, 00000100, 00001000, 00010000, 00100000, 01000000, ...repeated...]
 * RESULT:
 * result = [F, T, T, T, T, T, T, T, ... repeated...]
 */
TEST_F(BoolFilterSimple, Uint32ToBoolPerformance) {
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
TEST_F(BoolFilterSimple, Uint32ToBoolLenghZero) {
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
 * Test failure cases of sakura_SetFalseIfNanOrInfFloat
 * RESULT:
 *   LIBSAKURA_SYMBOL(Status_kInvalidArgument)
 */
/* Null pointer arrays */
TEST_F(BoolFilterFloat, NanOrInfFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	float *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseIfNanOrInfFloat(num_data,
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseIfNanOrInfFloat(num_data,
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseIfNanOrInfFloat(num_data,
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

	LIBSAKURA_SYMBOL(Status) status = sakura_SetFalseIfNanOrInfFloat(num_data,
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
TEST_F(BoolFilterSimple, InvertBoolFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	bool *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_InvertBool(num_data, data_null,
			result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterSimple, InvertBoolFailNullResult) {
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
TEST_F(BoolFilterSimple, InvertBoolNotAlignedData) {
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

TEST_F(BoolFilterSimple, InvertBoolNotAlignedResult) {
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
TEST_F(BoolFilterSimple, Uint8ToBoolFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	uint8_t *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint8ToBool(num_data, data_null,
			result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterSimple, Uint8ToBoolFailNullResult) {
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
TEST_F(BoolFilterSimple, Uint8ToBoolNotAlignedData) {
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

TEST_F(BoolFilterSimple, Uint8ToBoolNotAlignedResult) {
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
TEST_F(BoolFilterSimple, Uint32ToBoolFailNullData) {
	size_t const num_data(NUM_IN);
	SIMD_ALIGN
	bool result[num_data];

	uint32_t *data_null = nullptr;

	LIBSAKURA_SYMBOL(Status) status = sakura_Uint32ToBool(num_data, data_null,
			result);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(BoolFilterSimple, Uint32ToBoolFailNullResult) {
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
TEST_F(BoolFilterSimple, Uint32ToBoolNotAlignedData) {
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

TEST_F(BoolFilterSimple, Uint32ToBoolNotAlignedResult) {
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

// TODO: performance measurements that defines which of scalar or vector version is fast in range functions.
// It is because scalar version can be fast if the data is in the first range (stops loop)
// It should not necessary be unit test. Running manuall tests will be fine to determine in which case which is faster.
// 1. Use variety of num_bounds from 1~16 (templated range)
// 2. Set a data range that will be true to either, idx=0, num_range/2, or num_range-1.
//    The other ranges should be out of range of data distribution.
// 3. run scalar and vector version to compare the speed
//
