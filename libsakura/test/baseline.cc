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
#include "aligned_memory.h"
#include "gtest/gtest.h"
#include "baseline.h"

/* the number of elements in input/output array to test */
#define NUM_DATA 5
#define NUM_MODEL 3

using namespace std;

/*
 * A super class to test baseline functions
 */
class Baseline : public ::testing::Test {
protected:

	Baseline() : verbose(false)
	{}

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
		if (print_name) cout << name << " = ";
		cout << "[";
		for (size_t i = start_idx; i < start_idx+print_length-1; ++i)
			cout << data[i] << ", " ;
		cout << data[start_idx+print_length-1];
		cout << " ]";
		if (newline) cout << endl;
	}
	//1D double array
	void PrintArray(char const *name, size_t print_length, double const *data,
			size_t start_idx = 0, bool print_name = true, bool newline = true) {
		if (print_name) cout << name << " = ";
		cout << "[";
		for (size_t i = start_idx; i < start_idx+print_length-1; ++i)
			cout << data[i] << ", " ;
		cout << data[start_idx+print_length-1];
		cout << " ]";
		if (newline) cout << endl;
	}
	//1D bool array
	void PrintArray(char const *name, size_t print_length, bool const *data,
			size_t start_idx = 0, bool print_name = true, bool newline = true) {
		if (print_name) cout << name << " = ";
		cout << "[";
		for (size_t i = start_idx; i < start_idx+print_length-1; ++i)
			cout << (data[i] ? "T" : "F") << ", " ;
		cout << (data[start_idx+print_length-1] ? "T" : "F");
		cout << " ]";
		if (newline) cout << endl;
	}
	//given as 1D float array but actually stores (num_row * num_column) 2D data
	//for which column loop comes inside row loop.
	void PrintArray(char const *name, size_t num_row, size_t num_column, float const *data) {
		cout << name << " = [";
		for (size_t i = 0; i < num_row; ++i) {
			PrintArray(name, num_column, data, num_column*i, false, false);
			if (i < num_row-1) cout << ", " << endl;
		}
		cout << " ]" << endl;
	}
	//given as 1D double array but actually stores (num_row * num_column) 2D data
	//for which column loop comes inside row loop.
	void PrintArray(char const *name, size_t num_row, size_t num_column, double const *data) {
		cout << name << " = [";
		for (size_t i = 0; i < num_row; ++i) {
			PrintArray(name, num_column, data, num_column*i, false, false);
			if (i < num_row-1) cout << ", " << endl;
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
TEST_F(Baseline, CreateBaselineContext) {
	uint16_t const order(20);
	size_t const num_chan(4096);

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL(Status) status =
			sakura_CreateBaselineContext(
					LIBSAKURA_SYMBOL(BaselineType_kChebyshev),
					order, num_chan, &context);

	sakura_DestroyBaselineContext(context);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
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
	LIBSAKURA_SYMBOL(Status) status =
			sakura_CreateBaselineContext(
					LIBSAKURA_SYMBOL(BaselineType_kChebyshev),
					order, num_chan, &context);
	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);

	sakura_DestroyBaselineContext(context);

}

/*
 * Test sakura_GetBaselineModelPolynomial
 * RESULT:
 * out = []
 */
TEST_F(Baseline, GetBaselineModelPolynomial) {
	uint16_t const order(20);
	size_t const num_chan(4096);
	size_t const num_model(order+1);
	double out[num_chan*num_model];
	double answer[ELEMENTSOF(out)];

	// Setup answer----------------------------
	size_t idx = 0;
	for (size_t i = 0; i < num_chan; ++i) {
		for (size_t j = 0; j < num_model; ++j) {
			double value = 1.0;
			for (size_t k = 0; k < j; ++k) {
				value *= (double)i;
			}

			answer[idx] = value;
			idx++;
		}
	}
	//---------------------------------

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kPolynomial),
			order, num_chan, &context);

	for (size_t i = 0; i < ELEMENTSOF(out); ++i) {
		out[i] = context->basis_data[i];
	}
	LIBSAKURA_SYMBOL(BaselineType) type = context->baseline_type;
	size_t num_bases = context->num_bases;
	size_t num_basis_data = context->num_basis_data;

	sakura_DestroyBaselineContext(context);

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

	// Verification
	for (size_t i = 0; i < ELEMENTSOF(out); ++i) {
		EXPECT_EQ(out[i], answer[i]);
	}
	EXPECT_EQ(type, LIBSAKURA_SYMBOL(BaselineType_kPolynomial));
	EXPECT_EQ(num_bases, num_model);
	EXPECT_EQ(num_basis_data, num_chan);
	}

/*
 * Test sakura_GetBaselineModelChebyshev
 * RESULT:
 * out = []
 */
TEST_F(Baseline, GetBaselineModelChebyshev) {
	uint16_t const order(20);
	size_t const num_chan(4096);
	size_t const num_model(order+1);
	double out[num_chan*num_model];
	double answer[ELEMENTSOF(out)];

	// Setup answer----------------------------
	size_t idx = 0;
	for (size_t i = 0; i < num_chan; ++i) {
		for (size_t j = 0; j < num_model; ++j) {
			double value = 1.0;
			double x = 2.0*(double)i/(double)(num_chan-1)-1.0;
			if (j == 0) {
				value = 1.0;
			} else if (j == 1) {
				value = x;
			} else {
				value = 2.0 * x * answer[idx-1] - answer[idx-2];
			}

			answer[idx] = value;
			idx++;
		}
	}
	//---------------------------------

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kChebyshev),
			order, num_chan, &context);
	for (size_t i = 0; i < ELEMENTSOF(out); ++i) {
		out[i] = context->basis_data[i];
	}
	sakura_DestroyBaselineContext(context);

	// Verification
	for (size_t i = 0; i < ELEMENTSOF(out); ++i) {
		EXPECT_EQ(out[i], answer[i]);
	}
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
	sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kPolynomial),
			order, num_chan, &context);
	LIBSAKURA_SYMBOL(Status) status = sakura_DestroyBaselineContext(context);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

/*
 * Test sakura_GetBestFitBaselineModel
 * RESULT:
 * out = []
 */
TEST_F(Baseline, GetBestFitBaseline) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN float in_data[num_data] = {1.0, 3.0, 7.0, 130.0, 21.0};
	SIMD_ALIGN bool in_mask[ELEMENTSOF(in_data)] = {true, true, true, false, true};
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)] = {1.0, 3.0, 7.0, 13.0, 21.0};

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(CreateBaselineContext)(
					LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order, num_data, &context);
	LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL(GetBestFitBaseline)(num_data, in_data, in_mask, context, out);
	LIBSAKURA_SYMBOL(DestroyBaselineContext)(context);

	if (verbose) {
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0 ; i < num_model; ++i){
		ASSERT_EQ(out[i], answer[i]);
	}
}

TEST_F(Baseline, SubtractBaselinePolynomial) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN float in_data[num_data] = {1.0, 3.0, 7.0, 130.0, 21.0};
	SIMD_ALIGN bool in_mask[ELEMENTSOF(in_data)] = {true, true, true, false, true};
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)] = {0.0, 0.0, 0.0, 117.0, 0.0};

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	SIMD_ALIGN bool final_mask[ELEMENTSOF(in_data)];
	LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(
					num_data, in_data, in_mask, order,
					3.0, 1, true, final_mask, out);

	if (verbose) {
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0 ; i < num_model; ++i){
		ASSERT_EQ(out[i], answer[i]);
	}
}
