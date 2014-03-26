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
#define NUM_REPEAT 20000

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

	size_t const num_repeat(NUM_REPEAT);
	for (size_t i = 0; i < num_repeat; ++i) {
		LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
		LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order, num_chan, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
				context);
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

	size_t const num_repeat(NUM_REPEAT);
	for (size_t i = 0; i < num_repeat; ++i) {
		LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
		LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_chan, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
				context);
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
	SIMD_ALIGN
	float in_data[num_data] = { 1.0, 3.0, 7.0, 130.0, 21.0 };
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)] = { true, true, true, false, true };
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
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

	LIBSAKURA_SYMBOL(BaselineStatus) getbl_blstatus;
	LIBSAKURA_SYMBOL(Status) getbl_status =
	LIBSAKURA_SYMBOL(GetBestFitBaseline)(num_data, in_data, in_mask, context,
			out, &getbl_blstatus);
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
 * Test sakura_GetBestFitBaseline: failure case of too many masked
 * input data.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, GetBestFitBaselineWithTooManyMaskedData) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data] = { 1.0, 3.0, 7.0, 130.0, 21.0 };
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)] = { true, false, false, false, true };
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

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

	LIBSAKURA_SYMBOL(BaselineStatus) getbl_blstatus;
	LIBSAKURA_SYMBOL(Status) getbl_status =
	LIBSAKURA_SYMBOL(GetBestFitBaseline)(num_data, in_data, in_mask, context,
			out, &getbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kNG), getbl_status);

	LIBSAKURA_SYMBOL(Status) destroy_status =
	LIBSAKURA_SYMBOL(DestroyBaselineContext)(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaselineFromNormalDataWithoutClipping
 * the input data have smooth shape and no spiky feature, and
 * sakura_SubtractBaseline is executed without doing recursive
 * clipping.
 * the baseline-subtracted data should be zero throughout.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselineFromSmoothDataWithoutClipping) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 3;
	bool get_residual = true;

	LIBSAKURA_SYMBOL(BaselineStatus) subbl_blstatus;

	LIBSAKURA_SYMBOL(Status) subbl_status =
	LIBSAKURA_SYMBOL(SubtractBaseline)(num_data, in_data, in_mask, context,
			clipping_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaselineFromSmoothDataWithClipping
 * the input data have smooth shape and no spiky feature, and
 * execute sakura_SubtractBaseline where recursive baseline fitting,
 * 10 times at maximum, with clipping data outside 3 sigma level.
 * since the input data haven't any outliers the baseline-subtracted
 * data should be identical with those in case without clipping,
 * namely, zero throughout.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselineFromSmoothDataWithClipping) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 10;
	bool get_residual = true;

	LIBSAKURA_SYMBOL(BaselineStatus) subbl_blstatus;

	LIBSAKURA_SYMBOL(Status) subbl_status =
	LIBSAKURA_SYMBOL(SubtractBaseline)(num_data, in_data, in_mask, context,
			clipping_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaselineFromSpikyDataWithClipping
 * the input data have three outliers.
 * execute sakura_SubtractBaseline where recursive
 * baseline fitting, 5 times at maximum, with clipping data
 * outside 3 sigma level.
 * in the process of recursive clipping, the outliers will be
 * removed from baseline fitting procedure and the final result
 * should be zero throughout except for the outliers.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselineFromSpikyDataWithClipping) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
		if (i == 3)
			in_data[i] += 1000.0f;
		if (i == 6)
			in_data[i] -= 200.0f;
		if (i == 8)
			in_data[i] += 100.0f;
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
		if (i == 3)
			answer[i] += 1000.0f;
		if (i == 6)
			answer[i] -= 200.0f;
		if (i == 8)
			answer[i] += 100.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 5;
	bool get_residual = true;

	LIBSAKURA_SYMBOL(BaselineStatus) subbl_blstatus;

	LIBSAKURA_SYMBOL(Status) subbl_status =
	LIBSAKURA_SYMBOL(SubtractBaseline)(num_data, in_data, in_mask, context,
			clipping_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaselineWithZeroNumFittingMax
 * the input data have smooth shape and no spiky feature, and
 * execute sakura_SubtractBaseline with negative value of
 * num_fitting_max specifying the maximum number of fitting
 * procedure inside SubtractBaseline.
 * in case zero or negative value is given, num_fitting_max should
 * be changed to 1 inside SubtractBaseline so that baseline
 * subtraction IS executed.
 * hence the resulting data should be zero throughout.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselineWithZeroNumFittingMax) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 0;
	bool get_residual = true;

	LIBSAKURA_SYMBOL(BaselineStatus) subbl_blstatus;

	LIBSAKURA_SYMBOL(Status) subbl_status =
	LIBSAKURA_SYMBOL(SubtractBaseline)(num_data, in_data, in_mask, context,
			clipping_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaselineWithNegativeNumFittingMax
 * the input data have smooth shape and no spiky feature, and
 * execute sakura_SubtractBaseline with negative value
 * of num_fitting_max specifying the maximum number of fitting
 * procedure inside SubtractBaseline.
 * in case zero or negative value is given, num_fitting_max should
 * be changed to 1 inside SubtractBaseline so that baseline
 * subtraction IS executed.
 * hence the resulting data should be zero throughout.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselineWithNegativeNumFittingMax) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = -5;
	bool get_residual = true;

	LIBSAKURA_SYMBOL(BaselineStatus) subbl_blstatus;

	LIBSAKURA_SYMBOL(Status) subbl_status =
	LIBSAKURA_SYMBOL(SubtractBaseline)(num_data, in_data, in_mask, context,
			clipping_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaselineWithZeroClipThreshold
 * the input data have smooth shape and no spiky feature, and
 * execute sakura_SubtractBaseline with zero value of
 * clip_threshold_sigma.
 * in case zero is given for clipping threshold, SubtractBaseline()
 * should fail and return Status_kInvalidArgument.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselineWithZeroClipThreshold) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 0.0;
	uint16_t num_fitting_max = -5;
	bool get_residual = true;

	LIBSAKURA_SYMBOL(BaselineStatus) subbl_blstatus;

	LIBSAKURA_SYMBOL(Status) subbl_status =
	LIBSAKURA_SYMBOL(SubtractBaseline)(num_data, in_data, in_mask, context,
			clipping_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
	}

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaselineWithNegativeClipThreshold
 * the input data have smooth shape and no spiky feature, and
 * execute sakura_SubtractBaseline with negative value
 * of clip_threshold_sigma.
 * even if negative value is given for clipping threshold, there
 * should be no problem because the absolute value will be internally
 * used. hence the resulting data should be zero throughout.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselineWithNegativeClipThreshold) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = -3.0;
	uint16_t num_fitting_max = 5;
	bool get_residual = true;

	LIBSAKURA_SYMBOL(BaselineStatus) subbl_blstatus;

	LIBSAKURA_SYMBOL(Status) subbl_status =
	LIBSAKURA_SYMBOL(SubtractBaseline)(num_data, in_data, in_mask, context,
			clipping_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaseline: failure case of too many masked
 * input data.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselineWithTooManyMaskedData) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
		if (i == 3)
			in_data[i] += 1000.0f;
		if (i == 8)
			in_data[i] += 100.0f;
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = ((i == 0) || (i == ELEMENTSOF(in_data) - 1));
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 3;
	bool get_residual = true;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	LIBSAKURA_SYMBOL(BaselineStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL(Status) subbl_status =
	LIBSAKURA_SYMBOL(SubtractBaseline)(num_data, in_data, in_mask, context,
			clipping_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kNG), subbl_status);

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaseline: failure case of too many data
 * clipped in the process of recursive baseline fitting.
 * for this case, unmasked data should be less than the minimum needed
 * for fitting with the given baseline model just before the sixth
 * recursive baseline fitting.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselineTooManyDataClipped) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data] = { 1.0, 2.5, 8.0, 1013.0, 21.5, 31.5, 42.5, 57.5,
			173.5, 90.5, 111.5, 132.5, 157.0, 182.5, 211.5 };
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)] = { true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true };

	size_t order = 3;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status = sakura_CreateBaselineContext(
	LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 1.0;
	uint16_t num_fitting_max = 10;
	bool get_residual = true;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	LIBSAKURA_SYMBOL(BaselineStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(SubtractBaseline)(num_data, in_data, in_mask, context,
			clipping_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kNG), status);
	EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kNotEnoughData), subbl_blstatus);

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaselinePolynomialFromSmoothDataWithoutClipping
 * the input data have smooth shape and no spiky feature, and
 * sakura_SubtractBaselinePolynomial is executed without doing
 * recursive clipping.
 * the baseline-subtracted data should be zero throughout.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselinePolynomialFromSmoothDataWithoutClipping) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(num_data, in_data, in_mask,
			order, 3.0, 0, true, final_mask, out, &baseline_status);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}
	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}
}

/*
 * Test sakura_SubtractBaselinePolynomialFromSmoothDataWithClipping
 * the input data have smooth shape and no spiky feature, and
 * execute sakura_SubtractBaselinePolynomial where recursive
 * baseline fitting, 10 times at maximum, with clipping data
 * outside 3 sigma level.
 * since the input data haven't any outliers the baseline-subtracted
 * data should be identical with those in case without clipping,
 * namely, zero throughout.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselinePolynomialFromSmoothDataWithClipping) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(num_data, in_data, in_mask,
			order, 3.0, 10, true, final_mask, out, &baseline_status);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}
	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}
}

/*
 * Test sakura_SubtractBaselinePolynomialFromSpikyDataWithClipping
 * the input data have three outliers.
 * execute sakura_SubtractBaselinePolynomial where recursive
 * baseline fitting, 5 times at maximum, with clipping data
 * outside 3 sigma level.
 * in the process of recursive clipping, the outliers will be
 * removed from baseline fitting procedure and the final result
 * should be zero throughout except for the outliers.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselinePolynomialFromSpikyDataWithClipping) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
		if (i == 3)
			in_data[i] += 1000.0f;
		if (i == 6)
			in_data[i] -= 200.0f;
		if (i == 8)
			in_data[i] += 100.0f;
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
		if (i == 3)
			answer[i] += 1000.0f;
		if (i == 6)
			answer[i] -= 200.0f;
		if (i == 8)
			answer[i] += 100.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(num_data, in_data, in_mask,
			order, 3.0, 5, true, final_mask, out, &baseline_status);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}
	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}
}

/*
 * Test sakura_SubtractBaselinePolynomialWithZeroNumFittingMax
 * the input data have smooth shape and no spiky feature, and
 * execute sakura_SubtractBaselinePolynomial with negative value
 * of num_fitting_max specifying the maximum number of fitting
 * procedure inside SubtractBaselinePolynomial.
 * in case zero or negative value is given, num_fitting_max should
 * be changed to 1 inside SubtractBaseline so that baseline
 * subtraction IS executed.
 * hence the resulting data should be zero throughout.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselinePolynomialWithZeroNumFittingMax) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(num_data, in_data, in_mask,
			order, 3.0, 0, true, final_mask, out, &baseline_status);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}
	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}
}

/*
 * Test sakura_SubtractBaselinePolynomialWithNegativeNumFittingMax
 * the input data have smooth shape and no spiky feature, and
 * execute sakura_SubtractBaselinePolynomial with negative value
 * of num_fitting_max specifying the maximum number of fitting
 * procedure inside SubtractBaselinePolynomial.
 * in case zero or negative value is given, num_fitting_max should
 * be changed to 1 inside SubtractBaseline so that baseline
 * subtraction IS executed.
 * hence the resulting data should be zero throughout.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselinePolynomialWithNegativeNumFittingMax) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(num_data, in_data, in_mask,
			order, 3.0, -5, true, final_mask, out, &baseline_status);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}
	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}
}

/*
 * Test sakura_SubtractBaselinePolynomialWithZeroClipThreshold
 * the input data have smooth shape and no spiky feature, and
 * execute sakura_SubtractBaselinePolynomial with zero value
 * of clip_threshold_sigma.
 * in case zero is given for clipping threshold,
 * SubtractBaselinePolynomial() should fail and return
 * Status_kInvalidArgument.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselinePolynomialWithZeroClipThreshold) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(num_data, in_data, in_mask,
			order, 0.0, 5, true, final_mask, out, &baseline_status);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);

	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
	}
}

/*
 * Test sakura_SubtractBaselinePolynomialWithNegativeClipThreshold
 * the input data have smooth shape and no spiky feature, and
 * execute sakura_SubtractBaselinePolynomial with negative value
 * of clip_threshold_sigma.
 * even if negative value is given for clipping threshold, there
 * should be no problem because the absolute value will be internally
 * used. hence the resulting data should be zero throughout.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselinePolynomialWithNegativeClipThreshold) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		in_data[i] = (float) (1.0 + x + x * x);
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < num_data; ++i) {
		answer[i] = 0.0f;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(num_data, in_data, in_mask,
			order, -3.0, 5, true, final_mask, out, &baseline_status);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);

	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}
	if (verbose) {
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}
}

/*
 * Test sakura_SubtractBaselinePolynomial: failure case of too many
 * masked input data.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselinePolynomialWithTooManyMaskedData) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data] = { 1.0, 3.0, 7.0, 1013.0, 21.0, 31.0, 43.0, 57.0,
			173.0, 91.0, 111.0, 133.0, 157.0, 183.0, 211.0 };
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)] = {
	true, false, false, false, false,
	true, false, false, false, false,
	false, false, false, false, false };
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(BaselineStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(num_data, in_data, in_mask,
			order, 3.0, 3, true, final_mask, out, &subbl_blstatus);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kNG), status);
}

/*
 * Test sakura_SubtractBaselinePolynomial: failure case of too many data
 * clipped in the process of recursive baseline fitting.
 * for this case, unmasked data should be less than the minimum needed
 * for fitting with the given baseline model just before the sixth
 * recursive baseline fitting.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselinePolynomialTooManyDataClipped) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data] = { 1.0, 2.5, 8.0, 1013.0, 21.5, 31.5, 42.5, 57.5,
			173.5, 90.5, 111.5, 132.5, 157.0, 182.5, 211.5 };
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)] = { true, true, true, true, true, true,
	true, true, true, true, true, true, true, true, true };
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = 3;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(BaselineStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(num_data, in_data, in_mask,
			order, 1.0, 10, true, final_mask, out, &subbl_blstatus);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kNG), status);
	ASSERT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kNotEnoughData), subbl_blstatus);
}

/*
 * Test sakura_SubtractBaseline : subtract baseline from a realistic
 * spectrum data with 3840 channels. the data used here are taken from
 * the first spectrum of the test data used in the e2e test script for Sakura.
 * clipping conditions also are identical with those used in e2e test.
 * the baseline model here is a big set of Chebyshev polynomials with
 * order up to 400.
 * RESULT:
 * out = []
 */
TEST_F(Baseline, SubtractBaselineFromBigDataUsingBigChebyshevModel) {
	size_t const num_data(3840);
	size_t const num_model(400);
	float in_data_orig[3840] = { 8.12075996399, 8.14013004303, 8.16000080109,
			8.1760559082, 8.18900680542, 8.19975852966, 8.21892738342,
			8.2372674942, 8.2509469986, 8.26165103912, 8.27385234833,
			8.29015827179, 8.29803085327, 8.30214595795, 8.30259227753,
			8.29682731628, 8.28851222992, 8.27629184723, 8.26773643494,
			8.26203632355, 8.24806213379, 8.23066329956, 8.21049118042,
			8.17866611481, 8.14826965332, 8.13268280029, 8.12299346924,
			8.11022567749, 8.08844947815, 8.05942630768, 8.04142189026,
			8.03867053986, 8.03557300568, 8.03134918213, 8.03248786926,
			8.03361225128, 8.02828979492, 8.02212429047, 8.02766895294,
			8.04063129425, 8.04559326172, 8.04049396515, 8.03954029083,
			8.05141353607, 8.06872272491, 8.080037117, 8.08463096619,
			8.09135723114, 8.09775352478, 8.10004901886, 8.10222911835,
			8.10590934753, 8.11696147919, 8.12734127045, 8.12966823578,
			8.13006401062, 8.1245059967, 8.1180562973, 8.12057113647,
			8.12496185303, 8.13310050964, 8.14472866058, 8.1450548172,
			8.14438152313, 8.153424263, 8.16324996948, 8.16886043549,
			8.17666530609, 8.19106483459, 8.20676517487, 8.21819114685,
			8.22779941559, 8.24337100983, 8.25886058807, 8.27098369598,
			8.28762722015, 8.30189228058, 8.31071281433, 8.31874847412,
			8.3185749054, 8.31304264069, 8.31508636475, 8.31225681305,
			8.29509735107, 8.28453350067, 8.27665042877, 8.25447559357,
			8.2310628891, 8.21114349365, 8.18120098114, 8.14653778076,
			8.12124824524, 8.0987071991, 8.07039451599, 8.04445934296,
			8.01870727539, 7.98952388763, 7.95964288712, 7.93109512329,
			7.91236972809, 7.90481853485, 7.90168142319, 7.89691448212,
			7.90011692047, 7.90786027908, 7.90864229202, 7.91699743271,
			7.93609905243, 7.95349359512, 7.97292280197, 7.99950027466,
			8.0240240097, 8.05040836334, 8.08612346649, 8.12047195435,
			8.14892578125, 8.1740436554, 8.19662189484, 8.21628189087,
			8.23430156708, 8.25967884064, 8.27651977539, 8.28192329407,
			8.27618312836, 8.27581501007, 8.27398777008, 8.25656318665,
			8.23871040344, 8.21976566315, 8.19635677338, 8.17358398438,
			8.15049171448, 8.12299442291, 8.08843994141, 8.05372238159,
			8.01842212677, 7.99234056473, 7.98046445847, 7.96460819244,
			7.94320964813, 7.92727565765, 7.91835594177, 7.9087433815,
			7.8967962265, 7.89575338364, 7.90544891357, 7.91350221634,
			7.9277176857, 7.94727611542, 7.95976781845, 7.97664070129,
			7.99883317947, 8.01224994659, 8.02389621735, 8.0500535965,
			8.080991745, 8.09771633148, 8.10599327087, 8.11374092102,
			8.11828327179, 8.11477184296, 8.10269641876, 8.09242343903,
			8.08263587952, 8.06323242188, 8.04398441315, 8.02875900269,
			8.0060749054, 7.97612762451, 7.94669389725, 7.91384458542,
			7.87895441055, 7.8496761322, 7.82255077362, 7.79739236832,
			7.77262115479, 7.75328874588, 7.74347400665, 7.73239803314,
			7.72198581696, 7.7215681076, 7.72915935516, 7.73534631729,
			7.74204874039, 7.75922203064, 7.78162240982, 7.80480861664,
			7.83340978622, 7.8632516861, 7.89774799347, 7.93946218491,
			7.978495121, 8.01029205322, 8.03232097626, 8.05025959015,
			8.07281398773, 8.0999174118, 8.11682033539, 8.12128925323,
			8.12583637238, 8.12662029266, 8.11871528625, 8.10320949554,
			8.08525848389, 8.06492328644, 8.03468799591, 7.99826860428,
			7.96393585205, 7.93222475052, 7.90276813507, 7.86948394775,
			7.82638502121, 7.7907371521, 7.76683282852, 7.73626089096,
			7.7063331604, 7.69181728363, 7.68838119507, 7.68442630768,
			7.67272853851, 7.66588783264, 7.6700630188, 7.67634534836,
			7.68765449524, 7.70694446564, 7.73110294342, 7.75482130051,
			7.77898836136, 7.80154085159, 7.82368993759, 7.85228443146,
			7.87838172913, 7.89724254608, 7.91253423691, 7.92363548279,
			7.93133831024, 7.94097709656, 7.94649457932, 7.93646812439,
			7.91858291626, 7.90197610855, 7.87806987762, 7.85250663757,
			7.82348489761, 7.78578090668, 7.75032043457, 7.71190261841,
			7.66800308228, 7.6336812973, 7.61265945435, 7.58640003204,
			7.55141496658, 7.52179002762, 7.50370073318, 7.48872089386,
			7.47505664825, 7.47208452225, 7.47060346603, 7.47191476822,
			7.48727178574, 7.5068321228, 7.52751302719, 7.55151605606,
			7.57772445679, 7.6025800705, 7.62630462646, 7.65826129913,
			7.69212436676, 7.72396183014, 7.74967050552, 7.76435279846,
			7.7752571106, 7.78108263016, 7.78576469421, 7.78843402863,
			7.77717685699, 7.7523021698, 7.72194814682, 7.68650341034,
			7.64370393753, 7.59168767929, 7.53543615341, 7.48402070999,
			7.426445961, 7.35844564438, 7.29364871979, 7.23204660416,
			7.16787672043, 7.10465478897, 7.04803180695, 6.99645614624,
			6.95305919647, 6.92122793198, 6.89307785034, 6.86499595642,
			6.84217023849, 6.83409357071, 6.84502506256, 6.8603105545,
			6.87102127075, 6.89031648636, 6.91965961456, 6.95523691177,
			6.99712085724, 7.04324436188, 7.09282827377, 7.1440448761,
			7.19941234589, 7.25406122208, 7.30419874191, 7.35411453247,
			7.40337514877, 7.44373416901, 7.47509527206, 7.50509881973,
			7.51973152161, 7.51553583145, 7.50743579865, 7.49338150024,
			7.47181987762, 7.44305038452, 7.40384674072, 7.35607242584,
			7.29788351059, 7.23021125793, 7.1606631279, 7.09197759628,
			7.02676820755, 6.9671998024, 6.90535497665, 6.83601522446,
			6.76767063141, 6.70972967148, 6.66035652161, 6.61795949936,
			6.57948350906, 6.54442453384, 6.52096843719, 6.50720024109,
			6.50208091736, 6.50610542297, 6.51536083221, 6.53099870682,
			6.55338001251, 6.58043861389, 6.61492443085, 6.659283638,
			6.70566654205, 6.74889326096, 6.80102539062, 6.86321783066,
			6.91829776764, 6.9686961174, 7.02283477783, 7.07521677017,
			7.12642526627, 7.17720079422, 7.21603727341, 7.24156284332,
			7.26276016235, 7.29396724701, 7.31566810608, 7.32289218903,
			7.31900262833, 7.3077712059, 7.29758644104, 7.28303956985,
			7.26357746124, 7.24498128891, 7.2165927887, 7.18241262436,
			7.15189313889, 7.12148809433, 7.08942461014, 7.05474424362,
			7.02655172348, 7.00711584091, 6.98656511307, 6.96526622772,
			6.94593811035, 6.92827272415, 6.91369867325, 6.90319824219,
			6.89724254608, 6.89558076859, 6.89713811874, 6.90123224258,
			6.90536117554, 6.91017198563, 6.92112112045, 6.93430185318,
			6.94519710541, 6.95769500732, 6.97356843948, 6.99000883102,
			6.99649429321, 6.99665498734, 7.00787639618, 7.02662992477,
			7.03393602371, 7.0339179039, 7.04204130173, 7.05207014084,
			7.05915307999, 7.06457996368, 7.0679898262, 7.07726144791,
			7.08907318115, 7.09503412247, 7.09802532196, 7.10003614426,
			7.10211896896, 7.10813331604, 7.12300252914, 7.14033174515,
			7.15412759781, 7.16847038269, 7.18293142319, 7.19659996033,
			7.21498966217, 7.23845863342, 7.25808191299, 7.2738366127,
			7.28678846359, 7.29838228226, 7.31467294693, 7.33152008057,
			7.34107351303, 7.34662294388, 7.35745859146, 7.3679766655,
			7.37259244919, 7.37258529663, 7.37198877335, 7.37109899521,
			7.35597896576, 7.33397960663, 7.32441949844, 7.31617164612,
			7.29775285721, 7.27822971344, 7.25907945633, 7.23957252502,
			7.22060012817, 7.19554376602, 7.1694355011, 7.15178585052,
			7.13057613373, 7.10325908661, 7.08419466019, 7.07488679886,
			7.06893825531, 7.05994081497, 7.04905462265, 7.04426670074,
			7.04514694214, 7.04663848877, 7.05417299271, 7.06796360016,
			7.07979774475, 7.08876848221, 7.10159540176, 7.12155580521,
			7.13956403732, 7.15519142151, 7.1785068512, 7.20482683182,
			7.22220516205, 7.23671531677, 7.26199913025, 7.29081487656,
			7.30830907822, 7.32197618484, 7.33739042282, 7.3523182869,
			7.37109804153, 7.38440847397, 7.39513587952, 7.4084944725,
			7.41758298874, 7.43786096573, 7.44582748413, 7.44947099686,
			7.4576086998, 7.46680927277, 7.46546411514, 7.46395635605,
			7.46769332886, 7.46724510193, 7.47073030472, 7.47474193573,
			7.46950483322, 7.46989011765, 7.47301197052, 7.46851539612,
			7.46233177185, 7.45397424698, 7.44135284424, 7.4278254509,
			7.41576576233, 7.40061235428, 7.38312673569, 7.37171506882,
			7.35638618469, 7.33409023285, 7.31179189682, 7.28352308273,
			7.26162338257, 7.25456905365, 7.24107122421, 7.21767044067,
			7.19594192505, 7.18364477158, 7.17755126953, 7.16022157669,
			7.14126157761, 7.1329627037, 7.12560939789, 7.12444591522,
			7.13073253632, 7.13717794418, 7.14447975159, 7.15181446075,
			7.15719795227, 7.17126750946, 7.19828414917, 7.2257475853,
			7.2507610321, 7.27505588531, 7.30245685577, 7.33554077148,
			7.36972045898, 7.40125656128, 7.42739009857, 7.44836425781,
			7.46564149857, 7.48721694946, 7.51089286804, 7.52393388748,
			7.53575134277, 7.54749488831, 7.54917812347, 7.54609775543,
			7.53768825531, 7.52284002304, 7.50445890427, 7.48219966888,
			7.45161342621, 7.41650342941, 7.38719844818, 7.35738515854,
			7.31559085846, 7.2719912529, 7.23434734344, 7.19409227371,
			7.15120792389, 7.11495923996, 7.08636713028, 7.05603790283,
			7.02895593643, 7.01038360596, 6.9919090271, 6.97125387192,
			6.95658874512, 6.94792985916, 6.94391012192, 6.95299720764,
			6.96642827988, 6.97453641891, 6.98484754562, 7.00023460388,
			7.02356767654, 7.05664157867, 7.08934640884, 7.11670875549,
			7.14239692688, 7.16897106171, 7.20006370544, 7.23588371277,
			7.27000617981, 7.30117225647, 7.33173656464, 7.35954523087,
			7.38109731674, 7.40467500687, 7.42380571365, 7.43461036682,
			7.44584417343, 7.4462852478, 7.44770240784, 7.45388031006,
			7.45177268982, 7.45102071762, 7.45040941238, 7.44227838516,
			7.43293571472, 7.42309904099, 7.41221237183, 7.4028878212,
			7.39662122726, 7.39000654221, 7.38777256012, 7.38167524338,
			7.37587070465, 7.36949586868, 7.36725950241, 7.36857891083,
			7.36777496338, 7.37169694901, 7.37912225723, 7.38121318817,
			7.3858833313, 7.39857816696, 7.40633010864, 7.40869808197,
			7.41413593292, 7.41945028305, 7.4256644249, 7.43603610992,
			7.44831848145, 7.45369911194, 7.44784164429, 7.4362654686,
			7.42946624756, 7.43036842346, 7.42521429062, 7.41503000259,
			7.40727424622, 7.39666891098, 7.38546037674, 7.37762022018,
			7.36616897583, 7.35137701035, 7.34041595459, 7.32948064804,
			7.31983613968, 7.31468153, 7.30469846725, 7.29480361938,
			7.29581356049, 7.30042743683, 7.29973840714, 7.29743242264,
			7.29667377472, 7.29620361328, 7.300385952, 7.30842781067,
			7.31564617157, 7.32042551041, 7.32769012451, 7.33732461929,
			7.34823608398, 7.35955095291, 7.3635635376, 7.3736371994,
			7.38716697693, 7.38946914673, 7.39365053177, 7.4015545845,
			7.40407991409, 7.40293884277, 7.40679311752, 7.41059541702,
			7.40896320343, 7.40744495392, 7.39758348465, 7.38583135605,
			7.3802614212, 7.37234830856, 7.36265277863, 7.35017061234,
			7.33915185928, 7.33252668381, 7.32214927673, 7.31211948395,
			7.31156682968, 7.31147098541, 7.30255079269, 7.29064226151,
			7.27879953384, 7.27506637573, 7.27644062042, 7.26810407639,
			7.2610912323, 7.25880765915, 7.25329494476, 7.2477517128,
			7.24259185791, 7.23846006393, 7.23500108719, 7.23321962357,
			7.23261594772, 7.22985553741, 7.2273850441, 7.22069311142,
			7.21165513992, 7.21300411224, 7.21216917038, 7.20251178741,
			7.19823360443, 7.19941282272, 7.19721412659, 7.19054269791,
			7.18556451797, 7.17811155319, 7.16682577133, 7.16280412674,
			7.16801214218, 7.17212486267, 7.1760635376, 7.18280696869,
			7.18922042847, 7.19285821915, 7.19622612, 7.20518064499,
			7.21275901794, 7.22086763382, 7.23392534256, 7.24312591553,
			7.24937152863, 7.27102422714, 7.27957677841, 7.2900094986,
			7.30287265778, 7.31206274033, 7.32067489624, 7.32420444489,
			7.32544851303, 7.33018684387, 7.32965183258, 7.32411718369,
			7.32111167908, 7.31628990173, 7.30755662918, 7.29530620575,
			7.2789888382, 7.26533555984, 7.25621938705, 7.2435760498,
			7.2292766571, 7.22014570236, 7.20611143112, 7.19047689438,
			7.18460655212, 7.17795610428, 7.16532850266, 7.15320825577,
			7.14423084259, 7.1392788887, 7.1325211525, 7.12719011307,
			7.1325802803, 7.13932991028, 7.14376926422, 7.15123844147,
			7.16027975082, 7.16856098175, 7.17772483826, 7.19163990021,
			7.20455408096, 7.21709680557, 7.23159074783, 7.23876810074,
			7.24551773071, 7.25992250443, 7.27365350723, 7.28034257889,
			7.28389883041, 7.28806781769, 7.28486967087, 7.28060293198,
			7.27888011932, 7.27081251144, 7.26161384583, 7.25446510315,
			7.24624300003, 7.23374223709, 7.22269916534, 7.21406698227,
			7.19732141495, 7.17974233627, 7.17074108124, 7.16233634949,
			7.14715909958, 7.13699007034, 7.13875770569, 7.13605308533,
			7.1284403801, 7.12462615967, 7.1245303154, 7.13220453262,
			7.13844919205, 7.14393043518, 7.1594581604, 7.1772851944,
			7.19258737564, 7.20654821396, 7.22364807129, 7.24295186996,
			7.25549507141, 7.26450967789, 7.27773237228, 7.29530525208,
			7.3127617836, 7.32845211029, 7.33986091614, 7.34508943558,
			7.35282278061, 7.36022186279, 7.36034488678, 7.35868215561,
			7.36085224152, 7.36317110062, 7.35430526733, 7.33873796463,
			7.33117246628, 7.33199214935, 7.32635641098, 7.30455255508,
			7.28013849258, 7.27016639709, 7.2654004097, 7.26013088226,
			7.26109600067, 7.25931739807, 7.25124597549, 7.24667596817,
			7.24560785294, 7.2513461113, 7.26762104034, 7.28422451019,
			7.29740381241, 7.30861902237, 7.32161855698, 7.34473466873,
			7.37785339355, 7.4061870575, 7.42629861832, 7.44901561737,
			7.47250938416, 7.48797130585, 7.51387929916, 7.54278993607,
			7.56525039673, 7.58421373367, 7.60286951065, 7.61820650101,
			7.63612747192, 7.64692878723, 7.64927864075, 7.65635919571,
			7.66897773743, 7.68317651749, 7.6852850914, 7.67951202393,
			7.67772865295, 7.67118263245, 7.66951704025, 7.67873096466,
			7.67751216888, 7.67129993439, 7.67139625549, 7.67657804489,
			7.67688655853, 7.66998624802, 7.66837263107, 7.67410182953,
			7.68791055679, 7.69609451294, 7.70029973984, 7.71253108978,
			7.72021484375, 7.72721099854, 7.73458719254, 7.73880243301,
			7.74986362457, 7.764336586, 7.77630758286, 7.78498077393,
			7.79250049591, 7.79759788513, 7.80395746231, 7.8108625412,
			7.81357622147, 7.81834506989, 7.8182554245, 7.80784988403,
			7.7964220047, 7.7860622406, 7.7725892067, 7.75797367096,
			7.74567890167, 7.73742246628, 7.73258399963, 7.7244758606,
			7.71369934082, 7.70643424988, 7.69985294342, 7.6945271492,
			7.69061899185, 7.68529176712, 7.68505764008, 7.69140148163,
			7.69793367386, 7.70762968063, 7.7222700119, 7.7369017601,
			7.75099086761, 7.76789665222, 7.78630638123, 7.80368328094,
			7.81852912903, 7.82989835739, 7.84478855133, 7.85798597336,
			7.86314916611, 7.87223625183, 7.88794279099, 7.89679813385,
			7.89187002182, 7.88120508194, 7.87066602707, 7.85795688629,
			7.84157848358, 7.82203102112, 7.80018043518, 7.7726893425,
			7.73990917206, 7.70831251144, 7.68141031265, 7.65298748016,
			7.61882972717, 7.58770561218, 7.56066513062, 7.52926969528,
			7.49280834198, 7.46283388138, 7.4439034462, 7.42945957184,
			7.41876268387, 7.41082000732, 7.40243816376, 7.40181446075,
			7.41263723373, 7.42386102676, 7.43638277054, 7.45138311386,
			7.46871566772, 7.49596166611, 7.52716445923, 7.55870580673,
			7.58707714081, 7.61023569107, 7.64344930649, 7.68299055099,
			7.71091079712, 7.73693227768, 7.76369094849, 7.78177309036,
			7.79780483246, 7.80772829056, 7.81634664536, 7.821849823,
			7.8234205246, 7.81125736237, 7.79878616333, 7.79262161255,
			7.7728562355, 7.73766231537, 7.70658397675, 7.68517780304,
			7.65684223175, 7.6220164299, 7.59341096878, 7.5641579628,
			7.53332996368, 7.50884866714, 7.48821115494, 7.46977090836,
			7.45471096039, 7.44415044785, 7.43634605408, 7.4299030304,
			7.42841625214, 7.42968225479, 7.43168210983, 7.43833971024,
			7.45305299759, 7.47234201431, 7.49399757385, 7.51683521271,
			7.5370388031, 7.55874252319, 7.58504962921, 7.60699081421,
			7.62638187408, 7.65438890457, 7.6771903038, 7.69168567657,
			7.7034497261, 7.70588970184, 7.70719623566, 7.70922422409,
			7.70313739777, 7.69319725037, 7.68228578568, 7.6659450531,
			7.64343118668, 7.61889743805, 7.59564256668, 7.5702791214,
			7.54442930222, 7.52018070221, 7.48895978928, 7.4547586441,
			7.424223423, 7.39600372314, 7.37413644791, 7.35591030121,
			7.3371181488, 7.32147121429, 7.31439685822, 7.31173181534,
			7.30859422684, 7.31520652771, 7.32813453674, 7.34105396271,
			7.357026577, 7.36985397339, 7.38890218735, 7.4178442955,
			7.44713687897, 7.47846364975, 7.50843715668, 7.53969764709,
			7.57298517227, 7.60134601593, 7.62480258942, 7.64459514618,
			7.66711902618, 7.69338035583, 7.7134771347, 7.72321653366,
			7.73452234268, 7.74113988876, 7.73173522949, 7.7261390686,
			7.72366333008, 7.71203231812, 7.69356679916, 7.66852903366,
			7.64394617081, 7.62371587753, 7.59595537186, 7.56138277054,
			7.5275182724, 7.49813127518, 7.47464227676, 7.44927883148,
			7.42304325104, 7.40612459183, 7.39348554611, 7.36926412582,
			7.34794425964, 7.33885908127, 7.32804584503, 7.31911325455,
			7.31877279282, 7.32936382294, 7.34712028503, 7.35911035538,
			7.3720202446, 7.38709354401, 7.39404773712, 7.40526914597,
			7.42473363876, 7.44612550735, 7.46915388107, 7.49240016937,
			7.5139541626, 7.53750181198, 7.5583729744, 7.57932472229,
			7.58756303787, 7.59568166733, 7.61272573471, 7.62504768372,
			7.6348285675, 7.64352464676, 7.64865779877, 7.65363121033,
			7.65665864944, 7.65868663788, 7.66421175003, 7.67418050766,
			7.67996501923, 7.68048858643, 7.67968845367, 7.67782211304,
			7.68246746063, 7.69012975693, 7.6948928833, 7.7039642334,
			7.71075248718, 7.71052312851, 7.71944379807, 7.73800516129,
			7.74742841721, 7.74740934372, 7.75878000259, 7.77610445023,
			7.78703260422, 7.80537986755, 7.82770872116, 7.84540271759,
			7.8563451767, 7.85962772369, 7.86924123764, 7.8828792572,
			7.89242935181, 7.89809179306, 7.90184640884, 7.90954399109,
			7.91669368744, 7.91177225113, 7.90012025833, 7.89547872543,
			7.89494132996, 7.88976955414, 7.88169145584, 7.86871099472,
			7.85089111328, 7.83941364288, 7.83499288559, 7.82340431213,
			7.80694723129, 7.79700756073, 7.78599643707, 7.7789068222,
			7.78028249741, 7.7788028717, 7.76747369766, 7.75797796249,
			7.75942850113, 7.75746917725, 7.75511598587, 7.75481081009,
			7.75138282776, 7.74838495255, 7.7505941391, 7.76552295685,
			7.77750205994, 7.77929496765, 7.78556728363, 7.79563093185,
			7.80804395676, 7.81974554062, 7.82999801636, 7.83997964859,
			7.84459686279, 7.85568761826, 7.87321424484, 7.87802410126,
			7.87224388123, 7.8741440773, 7.88350439072, 7.89676856995,
			7.91217708588, 7.91396379471, 7.91644573212, 7.93265199661,
			7.94353342056, 7.9435839653, 7.94388103485, 7.94730520248,
			7.94927120209, 7.95259904861, 7.95686340332, 7.95935678482,
			7.96019697189, 7.96092891693, 7.96402263641, 7.96780061722,
			7.97185468674, 7.97772598267, 7.98070383072, 7.97977209091,
			7.97977972031, 7.9799323082, 7.97458028793, 7.96474075317,
			7.95917892456, 7.95623064041, 7.95267057419, 7.94477272034,
			7.93858957291, 7.94419288635, 7.9448633194, 7.93726348877,
			7.931640625, 7.92292976379, 7.91078281403, 7.89994764328,
			7.89950895309, 7.89612007141, 7.88330173492, 7.86849546432,
			7.8605298996, 7.86080789566, 7.86195707321, 7.85739135742,
			7.84625959396, 7.84491729736, 7.85140514374, 7.84696626663,
			7.84511852264, 7.85173988342, 7.8537197113, 7.85652828217,
			7.86043214798, 7.86306667328, 7.86864042282, 7.87668228149,
			7.88663244247, 7.89623832703, 7.90598201752, 7.9138250351,
			7.92060947418, 7.93022346497, 7.94087362289, 7.95249605179,
			7.96294355392, 7.96544694901, 7.96980571747, 7.98465538025,
			7.99530601501, 8.00395774841, 8.01284122467, 8.0133228302,
			8.01611328125, 8.02821540833, 8.03614616394, 8.03638362885,
			8.04499053955, 8.05343341827, 8.05628204346, 8.0597410202,
			8.05686092377, 8.05560779572, 8.06011009216, 8.06386470795,
			8.06881141663, 8.07204437256, 8.0724773407, 8.07656669617,
			8.08156681061, 8.08254146576, 8.08190727234, 8.08732032776,
			8.09604072571, 8.10279083252, 8.1077709198, 8.11094093323,
			8.11585998535, 8.12240028381, 8.12823104858, 8.13474464417,
			8.13929653168, 8.14160442352, 8.14708328247, 8.15375232697,
			8.15979003906, 8.16722583771, 8.17723464966, 8.18527793884,
			8.1900062561, 8.19507312775, 8.19888210297, 8.1975479126,
			8.19874191284, 8.20218086243, 8.20220756531, 8.20957660675,
			8.21345996857, 8.20805263519, 8.20603084564, 8.206199646,
			8.20200920105, 8.19400501251, 8.18927001953, 8.18848609924,
			8.18693733215, 8.17725753784, 8.16678237915, 8.16396808624,
			8.16009712219, 8.14932632446, 8.13758277893, 8.13101863861,
			8.12631320953, 8.11800003052, 8.10793113708, 8.10132694244,
			8.09561347961, 8.09071350098, 8.09251403809, 8.09502792358,
			8.09329986572, 8.09142971039, 8.08652019501, 8.07578372955,
			8.06736850739, 8.06982040405, 8.07739162445, 8.08055973053,
			8.08371543884, 8.09214115143, 8.10767173767, 8.11809253693,
			8.11334705353, 8.12094497681, 8.13580131531, 8.14248180389,
			8.15146064758, 8.16594028473, 8.17880916595, 8.1929693222,
			8.20161247253, 8.21062850952, 8.22759246826, 8.24045944214,
			8.245470047, 8.25302886963, 8.26477050781, 8.27261829376,
			8.27708053589, 8.27997875214, 8.28256511688, 8.28594017029,
			8.28858661652, 8.28808116913, 8.28802204132, 8.29207038879,
			8.28993225098, 8.28422355652, 8.27896785736, 8.27011489868,
			8.26236438751, 8.25346565247, 8.23919677734, 8.22921848297,
			8.22567844391, 8.21653747559, 8.20173740387, 8.18379688263,
			8.16135978699, 8.14494895935, 8.13179302216, 8.10928535461,
			8.08910942078, 8.07723140717, 8.05792808533, 8.02808952332,
			8.00094795227, 7.97582817078, 7.94917631149, 7.92945480347,
			7.9086151123, 7.87782096863, 7.85252475739, 7.83498477936,
			7.80735397339, 7.77524614334, 7.75733375549, 7.748067379,
			7.73451280594, 7.71885585785, 7.710252285, 7.70845031738,
			7.7015838623, 7.69658184052, 7.70042133331, 7.70156908035,
			7.70640420914, 7.71851921082, 7.72942399979, 7.74170160294,
			7.76014947891, 7.78483390808, 7.80711126328, 7.82922267914,
			7.86379289627, 7.90194272995, 7.93743944168, 7.97399377823,
			8.01144695282, 8.04484939575, 8.07136058807, 8.102850914,
			8.13945770264, 8.17034339905, 8.19402599335, 8.22113990784,
			8.25027561188, 8.26933288574, 8.28921031952, 8.30828666687,
			8.31445598602, 8.31519031525, 8.3120470047, 8.30115222931,
			8.29194450378, 8.28835010529, 8.2733001709, 8.2442741394,
			8.21689605713, 8.19259357452, 8.16825962067, 8.14166164398,
			8.10931015015, 8.0738735199, 8.03869247437, 8.00539493561,
			7.97433137894, 7.95084667206, 7.93400812149, 7.91096401215,
			7.88158178329, 7.85550165176, 7.83841705322, 7.82634782791,
			7.81282234192, 7.80706453323, 7.8120508194, 7.81428146362,
			7.81367063522, 7.81651878357, 7.81661510468, 7.81367778778,
			7.81571149826, 7.83448457718, 7.85028457642, 7.86972522736,
			7.88339042664, 7.89880657196, 7.91750240326, 7.92618751526,
			7.93582630157, 7.94816207886, 7.95617198944, 7.97035503387,
			7.99101066589, 7.9983625412, 8.00194740295, 8.01785469055,
			8.02729797363, 8.02806377411, 8.03111171722, 8.03461265564,
			8.03641796112, 8.03925704956, 8.04730319977, 8.05677223206,
			8.06223011017, 8.06377983093, 8.0672492981, 8.08069229126,
			8.09993743896, 8.11199760437, 8.11881351471, 8.1283454895,
			8.14508724213, 8.16801929474, 8.18829727173, 8.20136451721,
			8.21094036102, 8.22277736664, 8.23278331757, 8.2464799881,
			8.26721191406, 8.2781124115, 8.27765750885, 8.28029251099,
			8.28765487671, 8.28765583038, 8.27886581421, 8.26998329163,
			8.26150035858, 8.25570487976, 8.24668884277, 8.22724151611,
			8.20375537872, 8.17924308777, 8.16101646423, 8.15032577515,
			8.13052272797, 8.09796237946, 8.0692243576, 8.04964160919,
			8.0331993103, 8.01523303986, 7.98886871338, 7.96008682251,
			7.94431400299, 7.93288373947, 7.92188596725, 7.92204809189,
			7.92390441895, 7.91978025436, 7.91755819321, 7.92563295364,
			7.9391541481, 7.94779825211, 7.9506521225, 7.95870161057,
			7.97839736938, 7.99888181686, 8.01552772522, 8.02870178223,
			8.03908824921, 8.05135250092, 8.06275749207, 8.07809257507,
			8.0974187851, 8.10685920715, 8.11617660522, 8.13121509552,
			8.1361618042, 8.136510849, 8.13943576813, 8.14209842682,
			8.14547920227, 8.13940620422, 8.12221336365, 8.10440349579,
			8.0943775177, 8.08886432648, 8.07915592194, 8.06608200073,
			8.05385017395, 8.04351234436, 8.03125572205, 8.01826095581,
			8.00825214386, 8.00330924988, 8.00164604187, 7.99419784546,
			7.98777723312, 7.99309062958, 7.99785614014, 7.99850702286,
			8.00183200836, 8.00466632843, 8.00725364685, 8.01591777802,
			8.02947044373, 8.04274559021, 8.05018138885, 8.05762672424,
			8.07236099243, 8.08533287048, 8.08908462524, 8.1000585556,
			8.1081199646, 8.10737991333, 8.10841178894, 8.1076335907,
			8.10327243805, 8.0938615799, 8.07818126678, 8.06322002411,
			8.05473136902, 8.05430412292, 8.04357814789, 8.02090454102,
			8.01053905487, 8.00291061401, 7.98765897751, 7.97618293762,
			7.9688577652, 7.96003389359, 7.95253610611, 7.94756174088,
			7.93644046783, 7.92758512497, 7.93276882172, 7.94515752792,
			7.95648479462, 7.96745443344, 7.98261117935, 7.99899673462,
			8.02008342743, 8.05041313171, 8.07649707794, 8.09696006775,
			8.12817192078, 8.1654214859, 8.19467258453, 8.22054958344,
			8.25249290466, 8.27949714661, 8.29973125458, 8.32907390594,
			8.36018753052, 8.38244724274, 8.39462471008, 8.40434265137,
			8.42000865936, 8.43309020996, 8.44019699097, 8.44818210602,
			8.4523601532, 8.44953918457, 8.44711399078, 8.44277954102,
			8.43727779388, 8.43320178986, 8.42776679993, 8.41907405853,
			8.40760707855, 8.39582920074, 8.38512706757, 8.37880992889,
			8.37800121307, 8.37682342529, 8.373046875, 8.37230873108,
			8.37158107758, 8.36835384369, 8.36769771576, 8.37882232666,
			8.40246963501, 8.4202375412, 8.42630290985, 8.42285823822,
			8.42257785797, 8.43997764587, 8.46535110474, 8.48433017731,
			8.49385261536, 8.50587844849, 8.5204334259, 8.5278263092,
			8.53159713745, 8.53606510162, 8.53941440582, 8.54238891602,
			8.54304218292, 8.54106044769, 8.54251861572, 8.54559707642,
			8.54706192017, 8.53811740875, 8.52019786835, 8.51674461365,
			8.51939582825, 8.51021766663, 8.50270652771, 8.49911308289,
			8.48893642426, 8.47911071777, 8.47003078461, 8.46785163879,
			8.47099590302, 8.46349811554, 8.45697879791, 8.45515727997,
			8.453125, 8.45686340332, 8.46833896637, 8.48434352875,
			8.48797893524, 8.48212051392, 8.48927783966, 8.50193881989,
			8.51127147675, 8.52079963684, 8.52412319183, 8.53412246704,
			8.5527305603, 8.56171035767, 8.55998802185, 8.5654592514,
			8.5728187561, 8.57834339142, 8.58409690857, 8.58157634735,
			8.57708358765, 8.5900850296, 8.60315132141, 8.6060333252,
			8.60588932037, 8.59922981262, 8.60122871399, 8.60878372192,
			8.60922241211, 8.61146450043, 8.61116027832, 8.60934829712,
			8.6166601181, 8.62182331085, 8.62955474854, 8.63935375214,
			8.64517307281, 8.65904140472, 8.67726516724, 8.68775081635,
			8.6907043457, 8.69830036163, 8.71481227875, 8.73046112061,
			8.7350358963, 8.73908519745, 8.75208950043, 8.76104068756,
			8.76454639435, 8.76094341278, 8.75631523132, 8.76625347137,
			8.7700176239, 8.75919532776, 8.74896907806, 8.73895931244,
			8.72628974915, 8.71197605133, 8.69447898865, 8.67156410217,
			8.6522436142, 8.63641166687, 8.616938591, 8.5938501358,
			8.57089042664, 8.54966068268, 8.52265644073, 8.50081348419,
			8.48956680298, 8.479139328, 8.46874809265, 8.45402145386,
			8.43565273285, 8.42217826843, 8.41670894623, 8.41686916351,
			8.41648292542, 8.41208839417, 8.41033363342, 8.42042732239,
			8.43694496155, 8.45462417603, 8.47487545013, 8.49409866333,
			8.50907039642, 8.52535247803, 8.54294872284, 8.56055927277,
			8.58857440948, 8.62298774719, 8.65158462524, 8.67098903656,
			8.69236183167, 8.71860313416, 8.74315643311, 8.75912761688,
			8.76709461212, 8.78390979767, 8.80609416962, 8.81982135773,
			8.82165050507, 8.82336807251, 8.83252429962, 8.84101200104,
			8.84022045135, 8.82791900635, 8.81840229034, 8.81108856201,
			8.79709625244, 8.77909755707, 8.75783920288, 8.73830318451,
			8.72225284576, 8.70154762268, 8.67311573029, 8.64464855194,
			8.61246299744, 8.58212471008, 8.55392742157, 8.51890563965,
			8.48592185974, 8.45454883575, 8.42849349976, 8.40036773682,
			8.35582637787, 8.31654930115, 8.29163742065, 8.25892543793,
			8.2269153595, 8.2043762207, 8.17944717407, 8.15489292145,
			8.13393592834, 8.10920715332, 8.08725261688, 8.06917858124,
			8.05560112, 8.04668521881, 8.04496574402, 8.03890705109,
			8.02960205078, 8.03019809723, 8.04202461243, 8.06022834778,
			8.07244205475, 8.0830411911, 8.10368442535, 8.12589073181,
			8.15268802643, 8.18854522705, 8.21258735657, 8.23454761505,
			8.27639865875, 8.32245922089, 8.36564445496, 8.4053812027,
			8.44486141205, 8.48923110962, 8.5286655426, 8.56522750854,
			8.59933185577, 8.63386631012, 8.67156505585, 8.69537830353,
			8.71320533752, 8.73652172089, 8.76090049744, 8.77928447723,
			8.78597068787, 8.78885173798, 8.79120826721, 8.78296852112,
			8.76660442352, 8.75299167633, 8.74104595184, 8.72427272797,
			8.70394134521, 8.68170547485, 8.64847755432, 8.61767292023,
			8.59720802307, 8.56428432465, 8.52126979828, 8.48477458954,
			8.4527425766, 8.42115116119, 8.39040279388, 8.3608455658,
			8.3306388855, 8.30311870575, 8.27539157867, 8.24709129333,
			8.23041725159, 8.21603870392, 8.19637680054, 8.17853927612,
			8.16152191162, 8.14853954315, 8.14424228668, 8.14244747162,
			8.14183998108, 8.14499855042, 8.14484119415, 8.14966773987,
			8.16097450256, 8.16472244263, 8.17295837402, 8.18916893005,
			8.20871257782, 8.23186016083, 8.24964046478, 8.26715755463,
			8.28885173798, 8.31385040283, 8.33936023712, 8.36141967773,
			8.38468933105, 8.41045570374, 8.4335193634, 8.45475006104,
			8.481590271, 8.51040172577, 8.53642272949, 8.55854797363,
			8.57421875, 8.58964824677, 8.60865592957, 8.62957859039,
			8.64600372314, 8.6542930603, 8.66097831726, 8.67244911194,
			8.6796131134, 8.67263317108, 8.66524791718, 8.6613368988,
			8.65614414215, 8.65094184875, 8.63935661316, 8.62225151062,
			8.6027431488, 8.58055782318, 8.56503772736, 8.54844474792,
			8.51564311981, 8.47978210449, 8.45097160339, 8.42550468445,
			8.40064430237, 8.37327098846, 8.33735370636, 8.30659675598,
			8.29136276245, 8.27040290833, 8.24339771271, 8.21022510529,
			8.18413829803, 8.16638278961, 8.14966964722, 8.13397693634,
			8.11731243134, 8.10600471497, 8.09916591644, 8.08992481232,
			8.08347415924, 8.08285903931, 8.085521698, 8.09181880951,
			8.09746360779, 8.10210704803, 8.10733222961, 8.11334037781,
			8.12590694427, 8.14626121521, 8.15655231476, 8.15448093414,
			8.16539382935, 8.18401908875, 8.19837379456, 8.21731090546,
			8.23066711426, 8.24215126038, 8.26081943512, 8.27018737793,
			8.2725687027, 8.27757167816, 8.28408145905, 8.28814983368,
			8.29102897644, 8.29270267487, 8.29063034058, 8.2891330719,
			8.28383636475, 8.27914047241, 8.2814874649, 8.27918052673,
			8.26929950714, 8.25757598877, 8.24315834045, 8.23264598846,
			8.22807312012, 8.22198963165, 8.20866298676, 8.18967247009,
			8.17512702942, 8.16254043579, 8.15053939819, 8.14206981659,
			8.12976455688, 8.11797237396, 8.1081943512, 8.10056495667,
			8.09845066071, 8.09035682678, 8.0835609436, 8.08365058899,
			8.08175754547, 8.07740306854, 8.07583808899, 8.07528114319,
			8.07213115692, 8.07848453522, 8.09290790558, 8.10063266754,
			8.10398006439, 8.11578273773, 8.13127422333, 8.1387348175,
			8.14239025116, 8.14617443085, 8.14684104919, 8.14973831177,
			8.15855121613, 8.16462135315, 8.16651630402, 8.16943264008,
			8.1736240387, 8.17929935455, 8.17649841309, 8.16615486145,
			8.16076850891, 8.15847682953, 8.15906715393, 8.15587711334,
			8.14550018311, 8.13190937042, 8.11627483368, 8.10281276703,
			8.08499717712, 8.07106208801, 8.06450176239, 8.04638195038,
			8.0274181366, 8.0214881897, 8.0121717453, 7.99555110931,
			7.98378658295, 7.97527599335, 7.96333408356, 7.95803928375,
			7.9578742981, 7.95024776459, 7.94644117355, 7.94628477097,
			7.94229125977, 7.9434633255, 7.95550632477, 7.96650123596,
			7.97001266479, 7.97813129425, 7.98298072815, 7.98123168945,
			7.987844944, 8.00622749329, 8.01207065582, 8.01240825653,
			8.01420974731, 8.01412582397, 8.01056194305, 8.01174354553,
			8.01465702057, 8.01417541504, 8.01328659058, 8.00541591644,
			7.99083328247, 7.97608661652, 7.96217346191, 7.95311117172,
			7.93917703629, 7.914041996, 7.88911581039, 7.8722319603,
			7.86020517349, 7.84219312668, 7.82363510132, 7.80884170532,
			7.7899222374, 7.77101516724, 7.75485086441, 7.74352455139,
			7.73681592941, 7.73352527618, 7.73734617233, 7.74125623703,
			7.74720954895, 7.75708770752, 7.75978660583, 7.75976085663,
			7.77070236206, 7.79165697098, 7.81267738342, 7.83241462708,
			7.84935903549, 7.86441278458, 7.88396930695, 7.90109491348,
			7.91491365433, 7.93791389465, 7.95749473572, 7.96436786652,
			7.97738313675, 7.99677658081, 8.00944709778, 8.01682090759,
			8.02036762238, 8.0215177536, 8.02194023132, 8.01295089722,
			8.00166320801, 8.00256061554, 8.00358963013, 7.99451684952,
			7.98081636429, 7.96657943726, 7.95762252808, 7.94761323929,
			7.92990493774, 7.9192404747, 7.91975831985, 7.92061185837,
			7.9178276062, 7.91965913773, 7.92627763748, 7.92937469482,
			7.93645143509, 7.94256210327, 7.94594860077, 7.95875692368,
			7.97505950928, 7.99081802368, 8.01022338867, 8.02683258057,
			8.04180526733, 8.05618286133, 8.0687122345, 8.08584117889,
			8.10120677948, 8.11323738098, 8.11932373047, 8.12347507477,
			8.12703132629, 8.12467002869, 8.12807559967, 8.13009929657,
			8.12509918213, 8.11622905731, 8.09817504883, 8.07773303986,
			8.05582809448, 8.03408336639, 8.01940727234, 8.00762844086,
			7.99504470825, 7.9779176712, 7.95887708664, 7.94193124771,
			7.92559576035, 7.91250324249, 7.90504741669, 7.90154314041,
			7.89955282211, 7.90156888962, 7.90695333481, 7.91658401489,
			7.93647003174, 7.95165300369, 7.95865774155, 7.98552179337,
			8.02060222626, 8.04000473022, 8.05763816833, 8.08115577698,
			8.10701274872, 8.1279001236, 8.1558265686, 8.18471050262,
			8.20629405975, 8.22286510468, 8.23676586151, 8.24706459045,
			8.25141143799, 8.25315475464, 8.24898719788, 8.23819637299,
			8.22916030884, 8.21339797974, 8.20057106018, 8.18796348572,
			8.15643978119, 8.12331581116, 8.09493541718, 8.0592250824,
			8.02241230011, 7.98753547668, 7.95424604416, 7.92430877686,
			7.89721012115, 7.87378549576, 7.84861326218, 7.82665014267,
			7.81456041336, 7.80281257629, 7.78907251358, 7.77850246429,
			7.77373075485, 7.7762966156, 7.78176593781, 7.78808832169,
			7.80356216431, 7.82630777359, 7.84558677673, 7.86742734909,
			7.89493703842, 7.91988229752, 7.94623756409, 7.97864675522,
			8.00567531586, 8.02857017517, 8.05636978149, 8.07888317108,
			8.09322929382, 8.10212135315, 8.10752487183, 8.11847877502,
			8.12989425659, 8.13033103943, 8.12199783325, 8.11404514313,
			8.10856342316, 8.09463691711, 8.07516002655, 8.0606842041,
			8.04429626465, 8.02862167358, 8.01630592346, 7.99832201004,
			7.97359848022, 7.95045852661, 7.9403424263, 7.9374370575,
			7.93057823181, 7.92933797836, 7.93099784851, 7.92611932755,
			7.92421627045, 7.93072414398, 7.93897914886, 7.94601202011,
			7.95906591415, 7.97292804718, 7.98867893219, 8.01072597504,
			8.02522945404, 8.03672599792, 8.05204963684, 8.06488323212,
			8.08038330078, 8.09582901001, 8.10705089569, 8.11206531525,
			8.10976600647, 8.11283588409, 8.1147518158, 8.10827827454,
			8.0986623764, 8.08327102661, 8.06754684448, 8.05297470093,
			8.03747844696, 8.02433681488, 8.00309944153, 7.96787881851,
			7.93402242661, 7.90634059906, 7.88426017761, 7.86573696136,
			7.84505462646, 7.83006381989, 7.82144021988, 7.81271839142,
			7.8078417778, 7.80437469482, 7.80336380005, 7.8109369278,
			7.81934165955, 7.82671689987, 7.83588600159, 7.85297203064,
			7.8776550293, 7.89720344543, 7.91787433624, 7.94791173935,
			7.97441053391, 8.00080966949, 8.02242946625, 8.04664993286,
			8.06117248535, 8.07886791229, 8.09781169891, 8.11084651947,
			8.11897563934, 8.11888408661, 8.1235332489, 8.12145996094,
			8.10855674744, 8.0976524353, 8.08591938019, 8.07401180267,
			8.05179309845, 8.0274810791, 8.00279045105, 7.96401691437,
			7.93130397797, 7.91089677811, 7.88188505173, 7.84560966492,
			7.82038927078, 7.8028421402, 7.78096866608, 7.75903940201,
			7.74502658844, 7.73490953445, 7.71971988678, 7.71276140213,
			7.71540927887, 7.71359729767, 7.71482896805, 7.72317647934,
			7.73244047165, 7.74627494812, 7.7677989006, 7.79114770889,
			7.80957365036, 7.82397794724, 7.84686994553, 7.87697410583,
			7.90247058868, 7.92553853989, 7.94786787033, 7.97062110901,
			7.99542808533, 8.01482582092, 8.03016853333, 8.0438117981,
			8.04632282257, 8.04724407196, 8.05627441406, 8.0623626709,
			8.06123447418, 8.06007862091, 8.05867576599, 8.05533123016,
			8.04979991913, 8.0372133255, 8.02383899689, 8.01679801941,
			8.00748157501, 7.99249887466, 7.98589277267, 7.98120737076,
			7.97230672836, 7.96703720093, 7.95780277252, 7.94930648804,
			7.94986486435, 7.95316743851, 7.95547962189, 7.95929765701,
			7.95965862274, 7.9575510025, 7.9636888504, 7.9764919281,
			7.98884010315, 8.0000629425, 8.01032066345, 8.0184879303,
			8.02512264252, 8.02591514587, 8.02631092072, 8.03721141815,
			8.04961490631, 8.05681037903, 8.05956935883, 8.05864906311,
			8.05764293671, 8.05824565887, 8.06252861023, 8.05849933624,
			8.04572868347, 8.04095840454, 8.0378704071, 8.02921485901,
			8.0272064209, 8.03067207336, 8.02477169037, 8.01268005371,
			8.00638008118, 8.00598526001, 8.00617980957, 8.00500583649,
			8.00778675079, 8.01406288147, 8.01398277283, 8.01285648346,
			8.01827049255, 8.02656078339, 8.03614616394, 8.05166339874,
			8.06782627106, 8.07035255432, 8.07391929626, 8.08851718903,
			8.09603691101, 8.09920978546, 8.10350513458, 8.10431480408,
			8.10460281372, 8.10097026825, 8.09704589844, 8.10573101044,
			8.1165304184, 8.11258220673, 8.10056591034, 8.08940792084,
			8.07553577423, 8.06416988373, 8.06133747101, 8.05098438263,
			8.02857971191, 8.00943088531, 7.99775695801, 7.99024152756,
			7.97324895859, 7.95495557785, 7.94698286057, 7.93586158752,
			7.92636919022, 7.92343521118, 7.91314220428, 7.89678955078,
			7.88823795319, 7.88960170746, 7.89086341858, 7.88834190369,
			7.8857460022, 7.87930440903, 7.87472724915, 7.87606668472,
			7.87603664398, 7.87745046616, 7.88698339462, 7.89891576767,
			7.91248655319, 7.92502975464, 7.92958831787, 7.93137645721,
			7.9351940155, 7.94848775864, 7.96977043152, 7.98190784454,
			7.9828119278, 7.98692655563, 7.99910831451, 8.01429748535,
			8.0291519165, 8.04279232025, 8.05248641968, 8.06007862091,
			8.06725025177, 8.07155227661, 8.08245849609, 8.09367179871,
			8.09355163574, 8.0982093811, 8.10668182373, 8.10890197754,
			8.11395263672, 8.12080192566, 8.11918258667, 8.11107635498,
			8.10538196564, 8.09896087646, 8.09236812592, 8.09209728241,
			8.08216094971, 8.05889129639, 8.04368209839, 8.03299236298,
			8.01527023315, 8.00213909149, 7.98806619644, 7.96695947647,
			7.94608688354, 7.92065906525, 7.89285230637, 7.86815214157,
			7.84844350815, 7.83527755737, 7.82396030426, 7.80706262589,
			7.7871928215, 7.77477502823, 7.76968812943, 7.76201581955,
			7.75361585617, 7.75000143051, 7.74638652802, 7.74619436264,
			7.74815034866, 7.75573396683, 7.77593183517, 7.80075073242,
			7.82373094559, 7.84372234344, 7.86365652084, 7.88445186615,
			7.90669298172, 7.9334359169, 7.96518373489, 7.9974603653,
			8.028383255, 8.06354999542, 8.09663772583, 8.12056541443,
			8.1469783783, 8.17952060699, 8.2117099762, 8.24175071716,
			8.26283740997, 8.2753572464, 8.29265213013, 8.31352329254,
			8.32247543335, 8.31682300568, 8.31445503235, 8.30654621124,
			8.30077266693, 8.29605674744, 8.28337287903, 8.26080608368,
			8.23676776886, 8.21669197083, 8.18619632721, 8.1534910202,
			8.12752819061, 8.08954524994, 8.0504732132, 8.0152759552,
			7.97236871719, 7.9376411438, 7.90288066864, 7.85722780228,
			7.81747245789, 7.78164100647, 7.74490070343, 7.70848321915,
			7.67669010162, 7.64842939377, 7.6148724556, 7.57888936996,
			7.54172229767, 7.508664608, 7.48379707336, 7.46046972275,
			7.43728208542, 7.41721916199, 7.40450286865, 7.39292478561,
			7.37615871429, 7.36441230774, 7.36388206482, 7.36583375931,
			7.35910129547, 7.35236406326, 7.36095809937, 7.37675285339,
			7.38984632492, 7.40640926361, 7.42680454254, 7.44803667068,
			7.47155761719, 7.49774551392, 7.52938222885, 7.56322145462,
			7.59220981598, 7.62657022476, 7.66828870773, 7.70178794861,
			7.73767709732, 7.77957010269, 7.81648683548, 7.8530292511,
			7.88836240768, 7.92106199265, 7.9510216713, 7.98057699203,
			8.00653076172, 8.02073287964, 8.02902793884, 8.03838920593,
			8.04463100433, 8.03886508942, 8.02777671814, 8.01416015625,
			7.98697853088, 7.95278501511, 7.92303752899, 7.89523124695,
			7.86416149139, 7.8208565712, 7.77069759369, 7.73088788986,
			7.69000577927, 7.63652276993, 7.58099365234, 7.52870941162,
			7.47930145264, 7.43618965149, 7.39538288116, 7.35069799423,
			7.30868721008, 7.27143907547, 7.23760032654, 7.20584535599,
			7.17548036575, 7.1524720192, 7.12553691864, 7.09763813019,
			7.08527755737, 7.07927846909, 7.07491350174, 7.0706911087,
			7.06668806076, 7.06729650497, 7.06828832626, 7.08026456833,
			7.10380554199, 7.11961841583, 7.1340265274, 7.15830087662,
			7.18004417419, 7.19556760788, 7.2178311348, 7.24293994904,
			7.26498794556, 7.29029560089, 7.31816530228, 7.34862327576,
			7.37806367874, 7.40276861191, 7.42131614685, 7.43926429749,
			7.46066617966, 7.48698186874, 7.50341320038, 7.51798820496,
			7.53644514084, 7.55985736847, 7.57718086243, 7.58316802979,
			7.58795976639, 7.59943580627, 7.60762500763, 7.60337114334,
			7.59895324707, 7.60019207001, 7.59701967239, 7.58854055405,
			7.580057621, 7.57005214691, 7.55473804474, 7.53875398636,
			7.51663446426, 7.48870038986, 7.46195459366, 7.43736410141,
			7.41973495483, 7.40035152435, 7.37011194229, 7.33852243423,
			7.31049394608, 7.27769804001, 7.24490118027, 7.21624946594,
			7.184633255, 7.15414094925, 7.13201332092, 7.11711883545,
			7.10717821121, 7.09240865707, 7.07345819473, 7.06006574631,
			7.05254793167, 7.05556869507, 7.06222581863, 7.0647187233,
			7.07306861877, 7.08275651932, 7.09105396271, 7.11094808578,
			7.1412768364, 7.17603778839, 7.20678853989, 7.23275613785,
			7.26426839828, 7.30336904526, 7.3449549675, 7.38864660263,
			7.43411064148, 7.47198152542, 7.50999879837, 7.55424308777,
			7.59380817413, 7.63609552383, 7.68387651443, 7.72041416168,
			7.74715948105, 7.78106641769, 7.81638145447, 7.84249830246,
			7.86964654922, 7.89424800873, 7.90618562698, 7.91610527039,
			7.93164634705, 7.94595003128, 7.95048618317, 7.9488658905,
			7.94873762131, 7.9445734024, 7.93837070465, 7.93646526337,
			7.93930578232, 7.9431643486, 7.94124031067, 7.93673038483,
			7.93354511261, 7.93474149704, 7.93744754791, 7.93512487411,
			7.92979812622, 7.9325504303, 7.93791818619, 7.93598461151,
			7.94631719589, 7.96349239349, 7.96897697449, 7.97226285934,
			7.98111343384, 7.99856472015, 8.01602077484, 8.03056907654,
			8.04693508148, 8.06158733368, 8.07710647583, 8.0958480835,
			8.11620044708, 8.12934970856, 8.14094257355, 8.1618013382,
			8.17868232727, 8.19012260437, 8.21336078644, 8.24546909332,
			8.27228450775, 8.29243850708, 8.31325149536, 8.33903312683,
			8.36073684692, 8.38379001617, 8.4112739563, 8.43443107605,
			8.46226119995, 8.5097026825, 8.5420923233, 8.57272911072,
			8.59282016754, 8.61180305481, 8.63645362854, 8.66426944733,
			8.69983386993, 8.72837162018, 8.74774837494, 8.77154922485,
			8.79400253296, 8.81126594543, 8.83026790619, 8.84740924835,
			8.85813617706, 8.86283874512, 8.86305904388, 8.87088298798,
			8.88127231598, 8.8765745163, 8.86592388153, 8.85805225372,
			8.84579372406, 8.83888435364, 8.83751678467, 8.82890224457,
			8.81696605682, 8.80245304108, 8.7867641449, 8.77338218689,
			8.77211380005, 8.78382873535, 8.78965759277, 8.78530693054,
			8.78461456299, 8.79415416718, 8.80043888092, 8.80343532562,
			8.82076644897, 8.84230422974, 8.85735416412, 8.87958717346,
			8.90711116791, 8.93562984467, 8.96898460388, 9.00435733795,
			9.03919887543, 9.06735134125, 9.0891494751, 9.11720466614,
			9.15299606323, 9.17986774445, 9.19734477997, 9.21802902222,
			9.24070072174, 9.25604343414, 9.26302528381, 9.26618766785,
			9.26438617706, 9.25904560089, 9.25528335571, 9.2527179718,
			9.24340724945, 9.22673606873, 9.20702552795, 9.183842659,
			9.15718078613, 9.12649726868, 9.10314559937, 9.08885383606,
			9.06251716614, 9.03045845032, 9.00986099243, 8.98926734924,
			8.97032356262, 8.95202064514, 8.92939186096, 8.9147233963,
			8.90509319305, 8.89811992645, 8.89166259766, 8.88670730591,
			8.88739395142, 8.89155101776, 8.89944839478, 8.91003131866,
			8.924492836, 8.93530750275, 8.93812084198, 8.94413757324,
			8.95843219757, 8.9796333313, 8.99796295166, 9.00492572784,
			9.00811958313, 9.01684951782, 9.02856254578, 9.03725147247,
			9.04589080811, 9.05056285858, 9.05894470215, 9.07532691956,
			9.07605934143, 9.06739139557, 9.06939792633, 9.07129573822,
			9.06507396698, 9.06613731384, 9.07288742065, 9.06759166718,
			9.06147384644, 9.05836105347, 9.05193424225, 9.05603885651,
			9.05771160126, 9.04844665527, 9.05020523071, 9.05387401581,
			9.0497379303, 9.05460357666, 9.05587291718, 9.04615402222,
			9.03758811951, 9.03048038483, 9.03011035919, 9.03077793121,
			9.0185956955, 9.00527572632, 8.99473762512, 8.98477935791,
			8.9749956131, 8.95422649384, 8.92944908142, 8.90646839142,
			8.88236045837, 8.85961341858, 8.83289432526, 8.80227184296,
			8.77404880524, 8.7384519577, 8.69419765472, 8.65835762024,
			8.62932682037, 8.59500408173, 8.55800247192, 8.53522109985,
			8.5210647583, 8.50103569031, 8.47715568542, 8.45037174225,
			8.42890548706, 8.41841220856, 8.415599823, 8.416888237,
			8.42402935028, 8.43095684052, 8.43307590485, 8.43934440613,
			8.45232772827, 8.47066688538, 8.49050712585, 8.51458835602,
			8.54197216034, 8.56639957428, 8.59591674805, 8.62212848663,
			8.63752937317, 8.65888404846, 8.69217205048, 8.72128582001,
			8.74486732483, 8.76850509644, 8.7867641449, 8.80733299255,
			8.82584762573, 8.83770561218, 8.84750556946, 8.85550022125,
			8.86154556274, 8.85960769653, 8.8572883606, 8.85938930511,
			8.85383319855, 8.84424304962, 8.83831977844, 8.8278799057,
			8.8094625473, 8.79155254364, 8.77886295319, 8.76490497589,
			8.74396705627, 8.72332572937, 8.71000862122, 8.6964635849,
			8.67926025391, 8.66042423248, 8.63750171661, 8.61997890472,
			8.61191368103, 8.60063552856, 8.58731842041, 8.57721233368,
			8.56639194489, 8.55415821075, 8.54437828064, 8.54225635529,
			8.53579711914, 8.51976013184, 8.50736999512, 8.50066184998,
			8.50468444824, 8.51100826263, 8.50107288361, 8.48988342285,
			8.48544597626, 8.48881626129, 8.50056743622, 8.50863933563,
			8.51459407806, 8.51846313477, 8.53083515167, 8.55496692657,
			8.57367134094, 8.59433841705, 8.62432479858, 8.65167999268,
			8.67533493042, 8.70475959778, 8.73267650604, 8.76753616333,
			8.81079387665, 8.84567832947, 8.88494968414, 8.93697071075,
			8.98551559448, 9.02143478394, 9.0581035614, 9.1014623642,
			9.14610290527, 9.20051288605, 9.23249912262, 9.26359462738,
			9.29514408112, 9.32794857025, 9.35932826996, 9.38289165497,
			9.40723609924, 9.43443489075, 9.45656871796, 9.46492862701,
			9.46679782867, 9.47117996216, 9.47204780579, 9.4729385376,
			9.47546768188, 9.47761821747, 9.47433567047, 9.46589183807,
			9.47029018402, 9.48671340942, 9.49512195587, 9.49359703064,
			9.48854923248, 9.49119758606, 9.50436878204, 9.51592826843,
			9.53189277649, 9.54752731323, 9.55767822266, 9.57679176331,
			9.60015583038, 9.61745452881, 9.6329498291, 9.65721607208,
			9.68789958954, 9.70963096619, 9.73753166199, 9.77277755737,
			9.79416656494, 9.81324481964, 9.83187580109, 9.84925937653,
			9.87716388702, 9.90539836884, 9.92084598541, 9.9285736084,
			9.94198703766, 9.96739578247, 9.99095249176, 10.0006103516,
			10.010848999, 10.0229654312, 10.0322675705, 10.0451507568,
			10.0561819077, 10.0671215057, 10.0849123001, 10.1103906631,
			10.1305828094, 10.147939682, 10.1733150482, 10.1921777725,
			10.2092189789, 10.2362661362, 10.2682266235, 10.3022546768,
			10.330783844, 10.362036705, 10.3962440491, 10.4233169556,
			10.4524488449, 10.4834299088, 10.5140342712, 10.5362854004,
			10.5494661331, 10.5680608749, 10.5838441849, 10.5863695145,
			10.583562851, 10.5808067322, 10.5792798996, 10.5726327896,
			10.5591316223, 10.5445861816, 10.5213804245, 10.4935789108,
			10.4752807617, 10.4597196579, 10.4291238785, 10.3905715942,
			10.3654279709, 10.3433914185, 10.3173236847, 10.3054094315,
			10.3009614944, 10.2882242203, 10.276845932, 10.275894165,
			10.2868432999, 10.300283432, 10.3077297211, 10.3202962875,
			10.3381519318, 10.3657598495, 10.4030599594, 10.4395008087,
			10.475312233, 10.5159769058, 10.5604829788, 10.5976543427,
			10.6353721619, 10.6754112244, 10.7116441727, 10.748673439,
			10.78677845, 10.8144788742, 10.8260259628, 10.8337097168,
			10.8472442627, 10.8739833832, 10.8664312363, 10.849064827,
			10.8338956833, 10.8120260239, 10.7790575027, 10.7334012985,
			10.6938095093, 10.6627998352, 10.6242599487, 10.5810279846,
			10.5452213287, 10.5124750137, 10.4772901535, 10.4454317093,
			10.4233818054, 10.4057569504, 10.3888597488, 10.3769845963,
			10.364824295, 10.3628082275, 10.3821563721, 10.4039621353,
			10.4153766632, 10.4265089035, 10.451965332, 10.4849061966,
			10.5135040283, 10.5506944656, 10.5982646942, 10.6381797791,
			10.6749668121, 10.7137727737, 10.7421627045, 10.7712783813,
			10.8104801178, 10.8333053589, 10.8320426941, 10.830239296,
			10.8372955322, 10.8384752274, 10.8181533813, 10.7831249237,
			10.7455205917, 10.709274292, 10.6824960709, 10.6512813568,
			10.6035404205, 10.5519590378, 10.502202034, 10.4540309906,
			10.4039649963, 10.3548374176, 10.3171329498, 10.2861871719,
			10.2560052872, 10.2320766449, 10.215470314, 10.2045545578,
			10.1988277435, 10.1984949112, 10.2062797546, 10.2233791351,
			10.2500782013, 10.2846193314, 10.321852684, 10.3623361588,
			10.3965167999, 10.4210720062, 10.454416275, 10.4974813461,
			10.5331020355, 10.5538368225, 10.5773334503, 10.6005735397,
			10.6099796295, 10.6161413193, 10.615149498, 10.599770546,
			10.574930191, 10.5399389267, 10.4970788956, 10.4514331818,
			10.3991727829, 10.3403606415, 10.2685432434, 10.1871795654,
			10.117726326, 10.0569868088, 9.98915100098, 9.91470432281,
			9.84278869629, 9.78225898743, 9.73212718964, 9.68856239319,
			9.64961910248, 9.61619567871, 9.59454917908, 9.58640766144,
			9.5822019577, 9.57709407806, 9.58692264557, 9.61194133759,
			9.64644813538, 9.68899250031, 9.73058319092, 9.78229808807,
			9.84271812439, 9.89753437042, 9.95823669434, 10.022395134,
			10.0830984116, 10.1451435089, 10.2115259171, 10.2838926315,
			10.33614254, 10.3527460098, 10.3591089249, 10.3648223877,
			10.3715028763, 10.370308876, 10.352845192, 10.3227853775,
			10.2954845428, 10.2554769516, 10.1908693314, 10.1223812103,
			10.0570955276, 9.98553276062, 9.91130065918, 9.84477996826,
			9.79121112823, 9.73248386383, 9.66252422333, 9.61019229889,
			9.57046985626, 9.53252220154, 9.51004219055, 9.49350738525,
			9.47366428375, 9.46138763428, 9.47087097168, 9.49817657471,
			9.52549552917, 9.55762195587, 9.59518623352, 9.64422607422,
			9.70612049103, 9.76457023621, 9.82859802246, 9.89949893951,
			9.9686088562, 10.0289268494, 10.0798501968, 10.142334938,
			10.2126731873, 10.2615613937, 10.2965307236, 10.3284950256,
			10.3587169647, 10.3831310272, 10.382642746, 10.3657197952,
			10.3519096375, 10.3297138214, 10.2887096405, 10.2397994995,
			10.1875810623, 10.137471199, 10.0909824371, 10.0318975449,
			9.96884059906, 9.92136192322, 9.8743429184, 9.82722091675,
			9.78621864319, 9.74321269989, 9.71473693848, 9.70408344269,
			9.69111251831, 9.67748737335, 9.66719245911, 9.67398643494,
			9.71137428284, 9.75535678864, 9.7909412384, 9.82914161682,
			9.8795633316, 9.94044685364, 9.991979599, 10.0424823761,
			10.1040992737, 10.1642198563, 10.2247714996, 10.2862958908,
			10.3518323898, 10.4041891098, 10.4401168823, 10.4785575867,
			10.5050144196, 10.523516655, 10.5447511673, 10.557636261,
			10.55560112, 10.5472593307, 10.5467500687, 10.5373868942,
			10.5108633041, 10.490231514, 10.4742946625, 10.4531526566,
			10.431883812, 10.411945343, 10.3922452927, 10.368771553,
			10.3509893417, 10.3428916931, 10.3369607925, 10.3320827484,
			10.3260383606, 10.3326177597, 10.3539800644, 10.3721199036,
			10.391409874, 10.4168834686, 10.4375686646, 10.4620552063,
			10.4887657166, 10.5063476562, 10.5331287384, 10.5670051575,
			10.581202507, 10.5883340836, 10.6084346771, 10.6090431213,
			10.5895843506, 10.58273983, 10.5723342896, 10.5506238937,
			10.5288124084, 10.4985313416, 10.4543037415, 10.3913030624,
			10.3346424103, 10.2803840637, 10.2212209702, 10.1667985916,
			10.1126480103, 10.0584630966, 10.0132789612, 9.97293376923,
			9.93074893951, 9.89605712891, 9.88102054596, 9.8728055954,
			9.86134815216, 9.85831260681, 9.87147903442, 9.88943386078,
			9.9100522995, 9.95050239563, 10.0074138641, 10.0756025314,
			10.1458969116, 10.2094421387, 10.2813053131, 10.3571166992,
			10.4389629364, 10.5286636353, 10.6113615036, 10.6957960129,
			10.7838020325, 10.861946106, 10.9282102585, 10.9946041107,
			11.0486354828, 11.0855150223, 11.1089515686, 11.117975235,
			11.1259670258, 11.120057106, 11.091794014, 11.0542097092,
			11.0112905502, 10.9543571472, 10.8819093704, 10.8029384613,
			10.7260570526, 10.6503515244, 10.5673027039, 10.4840517044,
			10.400844574, 10.319984436, 10.2444143295, 10.1718215942,
			10.108669281, 10.0451917648, 9.97826194763, 9.92585659027,
			9.89709091187, 9.87382030487, 9.84582996368, 9.8252954483,
			9.80916309357, 9.80147266388, 9.81764030457, 9.83853816986,
			9.84668731689, 9.85850811005, 9.88575744629, 9.92117881775,
			9.95739650726, 9.98702812195, 10.0142555237, 10.0537528992,
			10.094004631, 10.1250658035, 10.1518793106, 10.174246788,
			10.2016010284, 10.2336158752, 10.2624111176, 10.2885246277,
			10.2977743149, 10.2945690155, 10.3022060394, 10.3057165146,
			10.3057670593, 10.3188123703, 10.316608429, 10.2990379333,
			10.2956075668, 10.2955102921, 10.2930459976, 10.2968626022,
			10.2962284088, 10.2840585709, 10.2687826157, 10.2582025528,
			10.2566928864, 10.257601738, 10.2526855469, 10.2460660934,
			10.2382965088, 10.2338132858, 10.2347669601, 10.228884697,
			10.2128019333, 10.1946430206, 10.1766700745, 10.157959938,
			10.1403055191, 10.1230764389, 10.0977487564, 10.069565773,
			10.0463666916, 10.0146121979, 9.97633934021, 9.94124317169,
			9.90774726868, 9.85868835449, 9.82887268066, 9.79410934448,
			9.7555141449, 9.71692371368, 9.67935180664, 9.64661026001,
			9.61658859253, 9.58608341217, 9.56206321716, 9.5441827774,
			9.52537536621, 9.51406860352, 9.50880908966, 9.50643062592,
			9.51738452911, 9.53403091431, 9.53822040558, 9.54361057281,
			9.57014560699, 9.59776973724, 9.61891555786, 9.64807605743,
			9.67604160309, 9.7021074295, 9.73254013062, 9.76813793182,
			9.80630207062, 9.83059310913, 9.84663677216, 9.87332057953,
			9.9007101059, 9.9166469574, 9.92750072479, 9.93980693817,
			9.95473575592, 9.96490383148, 9.96609401703, 9.96467208862,
			9.96201610565, 9.95935058594, 9.9537525177, 9.93957805634,
			9.92557621002, 9.91346549988, 9.8971862793, 9.88962078094,
			9.88484573364, 9.87056732178, 9.86117076874, 9.85793876648,
			9.85667610168, 9.85972690582, 9.86462211609, 9.86545085907,
			9.86494541168, 9.87139892578, 9.88342285156, 9.89833450317,
			9.90982532501, 9.92206096649, 9.94151210785, 9.95924091339,
			9.97497367859, 9.99114322662, 10.0041074753, 10.0178756714,
			10.0334329605, 10.0512580872, 10.0687351227, 10.077173233,
			10.0823774338, 10.0886878967, 10.1019058228, 10.1177434921,
			10.1250658035, 10.1308441162, 10.1387224197, 10.1484918594,
			10.1474866867, 10.1458721161, 10.1646547318, 10.1784181595,
			10.1765203476, 10.1835670471, 10.2094039917, 10.2375679016,
			10.2562980652, 10.2728738785, 10.2896375656, 10.3131637573,
			10.3432369232, 10.3676900864, 10.3929796219, 10.4255743027,
			10.4584655762, 10.4885311127, 10.5180006027, 10.5519552231,
			10.5824317932, 10.6012010574, 10.6253261566, 10.6523723602,
			10.66801548, 10.6753425598, 10.6798410416, 10.6815443039,
			10.6775169373, 10.6751976013, 10.6729888916, 10.6660642624,
			10.6566019058, 10.6464233398, 10.6466522217, 10.6476612091,
			10.6328201294, 10.6142473221, 10.6034898758, 10.6014413834,
			10.5993509293, 10.5849266052, 10.5935573578, 10.5977745056,
			10.6105461121, 10.6382865906, 10.6602754593, 10.6841630936,
			10.7143621445, 10.7549829483, 10.8072776794, 10.8508005142,
			10.8923082352, 10.945019722, 11.0106172562, 11.0731716156,
			11.1263189316, 11.1947088242, 11.2692184448, 11.3372020721,
			11.4038543701, 11.4680528641, 11.5332603455, 11.5933980942,
			11.6392364502, 11.6801681519, 11.7285795212, 11.7725982666,
			11.8039207458, 11.8332395554, 11.8595609665, 11.8783340454,
			11.8930301666, 11.8987379074, 11.8956098557, 11.8972883224,
			11.8994274139, 11.8966779709, 11.8988103867, 11.9010705948,
			11.8998022079, 11.9045782089, 11.9101285934, 11.9043636322,
			11.8992214203, 11.9055652618, 11.9179840088, 11.9351673126,
			11.9494838715, 11.9502801895, 11.9513835907, 11.9768733978,
			12.0135498047, 12.0384426117, 12.0564098358, 12.0809965134,
			12.1156215668, 12.1448459625, 12.1569576263, 12.1727600098,
			12.1994075775, 12.2198753357, 12.2398729324, 12.2541503906,
			12.2573423386, 12.2607393265, 12.2677783966, 12.2702188492,
			12.2733402252, 12.2837247849, 12.2924385071, 12.2955589294,
			12.2844619751, 12.2755699158, 12.2900733948, 12.303691864,
			12.3011665344, 12.3035774231, 12.317527771, 12.3393001556,
			12.3658847809, 12.4025354385, 12.4473972321, 12.484869957,
			12.5230941772, 12.5715284348, 12.6331996918, 12.7017803192,
			12.7598276138, 12.8194828033, 12.8962516785, 12.9710912704,
			13.028427124, 13.0805015564, 13.1471138, 13.225684166,
			13.2886447906, 13.3295269012, 13.3732414246, 13.4313383102,
			13.4700307846, 13.4780769348, 13.489153862, 13.4942417145,
			13.4872560501, 13.4718809128, 13.4347200394, 13.3928747177,
			13.3580551147, 13.3198719025, 13.26882267, 13.213300705,
			13.1648521423, 13.1103839874, 13.057518959, 13.0114660263,
			12.9473438263, 12.8772211075, 12.8259305954, 12.7896547318,
			12.7563295364, 12.7055749893, 12.6655073166, 12.6404180527,
			12.6228542328, 12.607011795, 12.6112031937, 12.6213436127,
			12.6221179962, 12.6270284653, 12.6375293732, 12.6537828445,
			12.6763906479, 12.7041921616, 12.7407226562, 12.7750911713,
			12.7975263596, 12.8185634613, 12.8411111832, 12.8636312485,
			12.89199543, 12.9131126404, 12.9186792374, 12.9296045303,
			12.9431848526, 12.9431238174, 12.9504308701, 12.9612197876,
			12.9551773071, 12.9445543289, 12.9239311218, 12.8956737518,
			12.8783979416, 12.8686733246, 12.8583021164, 12.8445129395,
			12.825381279, 12.8059406281, 12.7908630371, 12.7761821747,
			12.7570123672, 12.7351722717, 12.720246315, 12.7139196396,
			12.7058420181, 12.6850643158, 12.6590127945, 12.6457252502,
			12.6286916733, 12.6017274857, 12.5911655426, 12.5882234573,
			12.5655632019, 12.5268726349, 12.4842882156, 12.4493227005,
			12.4181127548, 12.3728084564, 12.3283309937, 12.2839870453,
			12.2221660614, 12.1540699005, 12.089345932, 12.0287389755,
			11.9666051865, 11.9028501511, 11.8308477402, 11.7472324371,
			11.6646556854, 11.5925798416, 11.5325641632, 11.4732933044,
			11.4066162109, 11.3339424133, 11.2624387741, 11.1985902786,
			11.1534366608, 11.1150188446, 11.0674667358, 11.0215091705,
			10.9771404266, 10.9372091293, 10.9092473984, 10.887793541,
			10.8566999435, 10.8229494095, 10.8008718491, 10.7783298492,
			10.7543392181, 10.7306413651, 10.698307991, 10.6690835953,
			10.6461753845, 10.6082105637, 10.5576534271, 10.512093544,
			10.4720449448, 10.4268341064, 10.3760900497, 10.3207149506,
			10.2612476349, 10.1956472397, 10.1244955063, 10.0522413254,
			9.96832466125, 9.88236618042, 9.80559253693, 9.72632598877,
			9.64264678955, 9.5582780838, 9.47072219849, 9.38729858398,
			9.31250572205, 9.23576545715, 9.15560626984, 9.07088851929,
			8.98687934875, 8.91624641418, 8.85789680481, 8.80196762085,
			8.74350452423 };

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		in_data[i] = in_data_orig[i];
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	for (size_t i = 0; i < ELEMENTSOF(in_data); ++i) {
		in_mask[i] = true;
	}
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
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
	uint16_t num_fitting_max = 2;
	bool get_residual = true;

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	LIBSAKURA_SYMBOL(BaselineStatus) subbl_blstatus;

	/*
	 size_t const num_repeat = 1;
	 double start = sakura_GetCurrentTime();
	 for (size_t i = 0; i < num_repeat; ++i) {
	 */
	LIBSAKURA_SYMBOL(Status) subbl_status =
	LIBSAKURA_SYMBOL(SubtractBaseline)(num_data, in_data, in_mask, context,
			clipping_threshold_sigma, num_fitting_max, get_residual, final_mask,
			out, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	/*
	 }
	 double end = sakura_GetCurrentTime();
	 */

	float limit_residual = 0.014f;
	//float max_residual = out[0];
	//float min_residual = out[0];
	for (size_t i = 0; i < num_model; ++i) {
		EXPECT_TRUE((out[i] > -limit_residual) && (out[i] < limit_residual));
		/*
		 if (out[i] > max_residual) {
		 max_residual = out[i];
		 }
		 if (out[i] < min_residual) {
		 min_residual = out[i];
		 }
		 */
	}
	//cout << "******************" << endl;
	//cout << "{residual: max = " << max_residual << ", min = " << min_residual << "}" << endl;

	if (verbose) {
		/*
		 cout << "Elapse time of " << num_repeat << " repetition: "
		 << end - start << " sec." << endl;
		 */
		PrintArray("fmask ", num_data, final_mask);
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	LIBSAKURA_SYMBOL(Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}
