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
			if (i < num_row-1) cout << ", ";
		}
		cout << " ]" << endl;
	}
	//given as 1D double array but actually stores (num_row * num_column) 2D data
	//for which column loop comes inside row loop.
	void PrintArray(char const *name, size_t num_row, size_t num_column, double const *data) {
		cout << name << " = [";
		for (size_t i = 0; i < num_row; ++i) {
			PrintArray(name, num_column, data, num_column*i, false, false);
			if (i < num_row-1) cout << ", ";
		}
		cout << " ]" << endl;
	}

	bool verbose;

};

/*
 * Test sakura_SolveSimultaneousEquationsByLU
 * RESULT:
 * out = []
 */
TEST_F(Baseline, GetBaselineModelPolynomial) {
	size_t const num_chan(1000);
	uint16_t const order(5);
	SIMD_ALIGN double out[num_chan*(order+1)];

	double start, end;
	size_t const num_repeat = 20000;
	LIBSAKURA_SYMBOL(Status) status;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_GetBaselineModelPolynomial(num_chan, order, out);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec" << endl;
	cout.precision(20);
	cout << "-------------------------" << endl;
	cout << "power0 : [" << out[0] << "], ..., [" << out[999] << "]" << endl;
	cout << "power1 : [" << out[1000] << "], ..., [" << out[1999] << "]" << endl;
	cout << "power2 : [" << out[2000] << "], ..., [" << out[2999] << "]" << endl;
	cout << "power3 : [" << out[3000] << "], ..., [" << out[3999] << "]" << endl;
	cout << "power4 : [" << out[4000] << "], ..., [" << out[4999] << "]" << endl;
	cout << "power5 : [" << out[5000] << "], ..., [" << out[5999] << "]" << endl;
	cout << "-------------------------" << endl;
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
	SIMD_ALIGN double model[num_model*ELEMENTSOF(in_data)] = {1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 4.0, 9.0, 16.0};
	SIMD_ALIGN float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)] = {1.0, 3.0, 7.0, 13.0, 21.0};

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
	}

	LIBSAKURA_SYMBOL(Status) status =
			sakura_GetBestFitBaseline(num_data, in_data, in_mask, num_model, model, out);

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
