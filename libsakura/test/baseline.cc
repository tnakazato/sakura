/*
 * baseline.cc
 *
 *  Created on: 2013/11/11
 *      Author: wataru
 */

#include <iostream>
#include <string>
#include <sys/time.h>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"
#include "baseline.h"

/* the number of elements in input/output array to test */
#define NUM_DATA 5
#define NUM_DATA2 15
#define NUM_DATA3 4096
#define NUM_MODEL 3
#define NUM_REPEAT 15000
#define NUM_REPEAT2 20000

using namespace std;

/*
 * A super class to test baseline functions
 */
class Baseline: public ::testing::Test {
protected:

	Baseline() :
			verbose(false) {
	}

	virtual void SetUp() {
		// Initialize sakura
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(Initialize)(nullptr, nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}

	virtual void TearDown() {
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	//1D float array
	void PrintArray(char const *name, size_t print_length, float const *data,
			size_t start_idx = 0, bool print_name = true, bool newline = true) {
		if (print_name)
			cout << name << " = ";
		cout << "[";
		for (size_t i = start_idx; i < start_idx + print_length - 1; ++i)
			cout << data[i] << ", ";
		cout << data[start_idx + print_length - 1];
		cout << " ]";
		if (newline)
			cout << endl;
	}
	//1D double array
	void PrintArray(char const *name, size_t print_length, double const *data,
			size_t start_idx = 0, bool print_name = true, bool newline = true) {
		if (print_name)
			cout << name << " = ";
		cout << "[";
		for (size_t i = start_idx; i < start_idx + print_length - 1; ++i)
			cout << data[i] << ", ";
		cout << data[start_idx + print_length - 1];
		cout << " ]";
		if (newline)
			cout << endl;
	}
	//1D bool array
	void PrintArray(char const *name, size_t print_length, bool const *data,
			size_t start_idx = 0, bool print_name = true, bool newline = true) {
		if (print_name)
			cout << name << " = ";
		cout << "[";
		for (size_t i = start_idx; i < start_idx + print_length - 1; ++i)
			cout << (data[i] ? "T" : "F") << ", ";
		cout << (data[start_idx + print_length - 1] ? "T" : "F");
		cout << " ]";
		if (newline)
			cout << endl;
	}
	//given as 1D float array but actually stores (num_row * num_column) 2D data
	//for which column loop comes inside row loop.
	void PrintArray(char const *name, size_t num_row, size_t num_column,
			float const *data) {
		cout << name << " = [";
		for (size_t i = 0; i < num_row; ++i) {
			PrintArray(name, num_column, data, num_column * i, false, false);
			if (i < num_row - 1)
				cout << ", " << endl;
		}
		cout << " ]" << endl;
	}
	//given as 1D double array but actually stores (num_row * num_column) 2D data
	//for which column loop comes inside row loop.
	void PrintArray(char const *name, size_t num_row, size_t num_column,
			double const *data) {
		cout << name << " = [";
		for (size_t i = 0; i < num_row; ++i) {
			PrintArray(name, num_column, data, num_column * i, false, false);
			if (i < num_row - 1)
				cout << ", " << endl;
		}
		cout << " ]" << endl;
	}

	bool verbose;

};

/*
 * Test sakura_CreateBaselineContext
 * RESULT:
 * out = []
 */
TEST_F(Baseline, CreateBaselineContextForPolynomial) {
	uint16_t const order(20);
	size_t const num_chan(4096);

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;

	size_t const num_repeat = NUM_REPEAT;
	LIBSAKURA_SYMBOL(Status) create_status;
	for (size_t i = 0; i < num_repeat; ++i) {
		create_status = sakura_CreateBaselineContext(
		LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order, num_chan, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

/*
 * Test sakura_CreateBaselineContext
 * RESULT:
 * out = []
 */
TEST_F(Baseline, CreateBaselineContextForChebyshevPolynomial) {
	uint16_t const order(20);
	size_t const num_chan(4096);

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;

	size_t const num_repeat = NUM_REPEAT;
	LIBSAKURA_SYMBOL(Status) create_status;
	for (size_t i = 0; i < num_repeat; ++i) {
		create_status = sakura_CreateBaselineContext(
		LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_chan, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

/*
 * Test sakura_CreateBaselineContext: failure case for too large order
 * RESULT:
 * out = []
 */
TEST_F(Baseline, CreateBaselineContextWithOrderLargerThanNumData) {
	uint16_t const order(20);
	size_t const num_chan(10);

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;

	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_chan, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), create_status);

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), destroy_status);
}

/*
 * Test sakura_GetBaselineModelPolynomial
 * RESULT:
 * out = []
 */
TEST_F(Baseline, GetBaselineModelPolynomial) {
	uint16_t const order(20);
	size_t const num_chan(4096);
	size_t const num_model(order + 1);
	double out[num_chan * num_model];
	double answer[ELEMENTSOF(out)];

	// Setup answer----------------------------
	size_t idx = 0;
	for (size_t i = 0; i < num_chan; ++i) {
		for (size_t j = 0; j < num_model; ++j) {
			double value = 1.0;
			for (size_t k = 0; k < j; ++k) {
				value *= (double) i;
			}

			answer[idx] = value;
			idx++;
		}
	}
	//---------------------------------

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order, num_chan, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	for (size_t i = 0; i < ELEMENTSOF(out); ++i) {
		out[i] = context->basis_data[i];
		EXPECT_EQ(out[i], answer[i]);
	}
	LIBSAKURA_SYMBOL(BaselineType) type = context->baseline_type;
	EXPECT_EQ(type, LIBSAKURA_SYMBOL(BaselineType_kPolynomial));
	size_t num_bases = context->num_bases;
	EXPECT_EQ(num_bases, num_model);
	size_t num_basis_data = context->num_basis_data;
	EXPECT_EQ(num_basis_data, num_chan);
	/*
	 cout.precision(20);
	 cout << "-------------------------" << endl;
	 cout << "power0 : [" << out[0] << "], ..., [" << out[5994] << "]" << endl;
	 cout << "power1 : [" << out[1] << "], ..., [" << out[5995] << "]" << endl;
	 cout << "power2 : [" << out[2] << "], ..., [" << out[5996] << "]" << endl;
	 cout << "power3 : [" << out[3] << "], ..., [" << out[5997] << "]" << endl;
	 cout << "power4 : [" << out[4] << "], ..., [" << out[5998] << "]" << endl;
	 cout << "power5 : [" << out[5] << "], ..., [" << out[5999] << "]" << endl;
	 cout << "-------------------------" << endl;
	 */

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_GetBaselineModelChebyshev
 * RESULT:
 * out = []
 */
TEST_F(Baseline, GetBaselineModelChebyshev) {
	uint16_t const order(20);
	size_t const num_chan(4096);
	size_t const num_model(order + 1);
	double out[num_chan * num_model];
	double answer[ELEMENTSOF(out)];

	// Setup answer----------------------------
	size_t idx = 0;
	for (size_t i = 0; i < num_chan; ++i) {
		for (size_t j = 0; j < num_model; ++j) {
			double value = 1.0;
			double x = 2.0 * (double) i / (double) (num_chan - 1) - 1.0;
			if (j == 0) {
				value = 1.0;
			} else if (j == 1) {
				value = x;
			} else {
				value = 2.0 * x * answer[idx - 1] - answer[idx - 2];
			}

			answer[idx] = value;
			idx++;
		}
	}
	//---------------------------------

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_chan, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	for (size_t i = 0; i < ELEMENTSOF(out); ++i) {
		out[i] = context->basis_data[i];
		EXPECT_EQ(out[i], answer[i]);
	}

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_DestroyBaselineContext
 * RESULT:
 * out = []
 */
TEST_F(Baseline, DestroyBaselineContext) {
	uint16_t const order(20);
	size_t const num_chan(4096);

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order, num_chan, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_GetBestFitBaseline
 * RESULT:
 * out = []
 */
TEST_F(Baseline, GetBestFitBaseline) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN float in_data[num_data] = { 1.0, 3.0, 7.0, 130.0, 21.0 };
	SIMD_ALIGN bool in_mask[ELEMENTSOF(in_data)] = { true, true, true, false, true };
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)] = { 1.0, 3.0, 7.0, 13.0, 21.0 };

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(Status) create_status =
	LIBSAKURA_SYMBOL(CreateBaselineContext)(
	LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	LIBSAKURA_SYMBOL(Status) getbl_status =
	LIBSAKURA_SYMBOL(GetBestFitBaseline)(num_data, in_data, in_mask, context,
			out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), getbl_status);
	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(out[i], answer[i]);
	}
	if (verbose) {
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	LIBSAKURA_SYMBOL(Status) destroy_status =
	LIBSAKURA_SYMBOL(DestroyBaselineContext)(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_GetBestFitBaseline: failure case of too many masked data
 * RESULT:
 * out = []
 */
TEST_F(Baseline, GetBestFitBaselineWithTooManyMaskedData) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN float in_data[num_data] = { 1.0, 3.0, 7.0, 130.0, 21.0 };
	SIMD_ALIGN bool in_mask[ELEMENTSOF(in_data)] = { true, false, false, false, true };
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN float out[ELEMENTSOF(in_data)];

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(Status) create_status =
	LIBSAKURA_SYMBOL(CreateBaselineContext)(
	LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	LIBSAKURA_SYMBOL(Status) getbl_status =
	LIBSAKURA_SYMBOL(GetBestFitBaseline)(num_data, in_data, in_mask, context,
			out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kNG), getbl_status);

	LIBSAKURA_SYMBOL(Status) destroy_status =
	LIBSAKURA_SYMBOL(DestroyBaselineContext)(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaseline
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaseline) {
	size_t const num_data(NUM_DATA3);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float)(1.0 + x + x * x);
	}
	SIMD_ALIGN bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 3;
	bool get_residual = true;

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	double start, end;
	size_t const num_repeat = NUM_REPEAT2;
	LIBSAKURA_SYMBOL(Status) subbl_status;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		subbl_status =
		LIBSAKURA_SYMBOL(SubtractBaseline)(num_data, in_data, in_mask, context,
				clipping_threshold_sigma, num_fitting_max, get_residual,
				final_mask, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	}
	end = sakura_GetCurrentTime();

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(out[i], answer[i]);
	}

	if (verbose) {
		cout << "Elapse time of " << num_repeat << " repetition: "
				<< end - start << " sec." << endl;
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaseline: failure case of too many masked data
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselineWithTooManyMaskedData) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float)(1.0 + x + x * x);
		if (i == 3) in_data[i] += 1000.0f;
		if (i == 8) in_data[i] += 100.0f;
	}
	SIMD_ALIGN bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = ((i == 0) || (i == ELEMENTSOF(in_data)-1));
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 3;
	bool get_residual = true;
	SIMD_ALIGN bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN float out[ELEMENTSOF(in_data)];

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	LIBSAKURA_SYMBOL(Status) subbl_status =
	LIBSAKURA_SYMBOL(SubtractBaseline)(num_data, in_data, in_mask, context,
			clipping_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kNG), subbl_status);

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaselinePolynomialWithClipping
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselinePolynomial) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float)(1.0 + x + x * x);
		if (i == 3) in_data[i] += 1000.0f;
		if (i == 8) in_data[i] += 100.0f;
	}
	SIMD_ALIGN bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
		if (i == 3) answer[i] += 1000.0f;
		if (i == 8) answer[i] += 100.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	SIMD_ALIGN bool final_mask[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(num_data, in_data, in_mask,
			order, 3.0, 3, true, final_mask, out);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(out[i], answer[i]);
	}
	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}
}

/*
 * Test sakura_SubtractBaselinePolynomial: failure case of too many masked data
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselinePolynomialWithTooManyMaskedData) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN float in_data[num_data] = {
			1.0, 3.0, 7.0, 1013.0, 21.0, 31.0, 43.0, 57.0,
			173.0, 91.0, 111.0, 133.0, 157.0, 183.0, 211.0 };
	SIMD_ALIGN bool in_mask[ELEMENTSOF(in_data)] = {
	true, false, false, false, false,
	true, false, false, false, false,
	false, false, false, false, false };
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN float out[ELEMENTSOF(in_data)];

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	SIMD_ALIGN bool final_mask[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(num_data, in_data, in_mask,
			order, 3.0, 3, true, final_mask, out);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kNG), status);
}
