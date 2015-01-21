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
/*
 * Test calses to be implemented.
 * test order parameters in
 * - sakura_GetNumberOfCoefficients
 * - sakura_SubtractBaselineFloat
 * - sakura_GetBestFitBaselineFloat
 *
 */

#include <cmath>
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
#define NUM_MODEL2 4
#define NUM_REPEAT2 2000000

using namespace std;

/*
 * A super class to test baseline functions
 */
class BaselineSK: public ::testing::Test {
protected:

	BaselineSK() :
			verbose(false) {
	}

	virtual void SetUp() {
		// Initialize sakura
		LIBSAKURA_SYMBOL (Status)
		status = LIBSAKURA_SYMBOL(Initialize)(nullptr, nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}

	virtual void TearDown() {
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	// Set (1+x+x*x) float values into an array
	void SetFloatPolynomial(size_t const num_data, float *data) {
		for (size_t i = 0; i < num_data; ++i) {
			double x = (double) i;
			data[i] = (float) (1.0 + x + x * x);
		}
	}

	//Set (A[0]+A[1]*x+A[2]*x*x+A[3]*x*x*x) float values into an array
	void SetFloatPolynomial(size_t num_data, float *data,
			double *coeff_answer) {
		for (size_t i = 0; i < num_data; ++i) {
			double x = (double) i;
			data[i] = (float) (coeff_answer[0] + coeff_answer[1] * x
					+ coeff_answer[2] * x * x + coeff_answer[3] * x * x * x);
		}
	}

	// Set constant float values into an array
	void SetFloatConstant(float value, size_t const num_data, float *data) {
		for (size_t i = 0; i < num_data; ++i) {
			data[i] = value;
		}
	}

	// Set constant boolean values into an array
	void SetBoolConstant(bool value, size_t const num_data, bool *data) {
		for (size_t i = 0; i < num_data; ++i) {
			data[i] = value;
		}
	}

	// Check if the expected and actual values are enough close to each other
	void CheckAlmostEqual(double expected, double actual, double tolerance) {
		double deviation = fabs(actual - expected);
		double val = max(fabs(actual), fabs(expected)) * tolerance + tolerance;
		ASSERT_LE(deviation, val);
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
 * Test sakura_SubtractBaselineCubicSplineUsingCoefficientsFloat
 * successful case : out == 0
 * subtract best fit model from input data using input coeff
 */
TEST_F(BaselineSK, SubtractBaselineCubicSplineUsingCoefficientsFloat) {
	size_t const num_data(8);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data] = { 2, 2.838888889, 4.144444444, 5, 3.425925926,
			5.907407407, 11.35, 20.65926 };

	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float answer[ELEMENTSOF(in_data)];
	SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
	}

	size_t order = num_model;
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status)
	create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
			&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	SIMD_ALIGN
	size_t num_pieces = 2;
	double coeff[8] = { 2, 0.3, 0.691666667, -0.152777778, 5, 0.325,
			-0.783333333, 0.150925926 };

	double boundary[3] = { 0.0f, 4.0f };

	LIBSAKURA_SYMBOL (Status)
	subbl_status = LIBSAKURA_SYMBOL(
			SubtractBaselineCubicSplineUsingCoefficientsFloat)(num_data,
			in_data, context, coeff, num_pieces, boundary, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);

	for (size_t i = 0; i < num_data; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	LIBSAKURA_SYMBOL (Status)
	destroy_status = sakura_DestroyBaselineContext(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}
