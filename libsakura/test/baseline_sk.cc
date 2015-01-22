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
#include <cstddef>
#include <climits>
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
		LIBSAKURA_SYMBOL (Status) status = LIBSAKURA_SYMBOL(Initialize)(nullptr,
				nullptr);
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
	size_t const num_pieces(2);
	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		float x = (float) i;
		in_data[i] = 2 + x - 1 * x * x + 0.2 * x * x * x;
	}
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float answer[ELEMENTSOF(in_data)];
	SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);
	if (verbose) {
		PrintArray("in_data", num_data, in_data);
	}
	size_t order = 2;
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
			&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	double coeff[num_data] = { 2.0, 1.0, -1.0, 0.2, 2.0, 1.0, -1.0, 0.2 };
	SIMD_ALIGN
	double boundary[num_pieces] = { 0.0f, 4.0f };
	LIBSAKURA_SYMBOL (Status) subbl_status =
			LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
					context, num_data, in_data, num_pieces, coeff, boundary,
					out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);

	for (size_t i = 0; i < num_data; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaselineCubicSplineUsingCoefficientsFloatWithNotAligned
 * failure case: data/coeff/out is not aligned
 * subtract best fit model from input data using input coeff
 */

TEST_F(BaselineSK, SubtractBaselineCubicSplineUsingCoefficientsFloatWithNotAligned) {
	{ // in_data is not aligned
		size_t const num_data(8);
		SIMD_ALIGN
		float in_data[num_data + 1];
		float *in_data_unaligned = in_data + 1;
		assert(!LIBSAKURA_SYMBOL(IsAligned)(in_data_unaligned));
		for (size_t i = 0; i < num_data; ++i) {
			in_data_unaligned[i] = 1.0 + i;
		}
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

		if (verbose) {
			PrintArray("in_data", num_data, in_data);
		}

		size_t order = 2;
		LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
				&context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		size_t num_pieces = 2;
		SIMD_ALIGN
		double coeff[8] = { 1.0f, 1.0f, -1.0f, 0.5f, 1.0f, 1.0f, -1.0f, 0.5f };
		SIMD_ALIGN
		double boundary[2] = { 0.0f, 4.0f };
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, num_data, in_data_unaligned, num_pieces, coeff,
						boundary, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	{ // coeff is not aligned
		size_t const num_data(8);
		SIMD_ALIGN
		float in_data[num_data];
		for (size_t i = 0; i < num_data; ++i) {
			float x = (float) i;
			in_data[i] = 1 + x - 1 * x * x + 0.5 * x * x * x;
		}
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

		if (verbose) {
			PrintArray("in_data", num_data, in_data);
		}

		size_t order = 2;
		LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
				&context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		size_t num_pieces = 2;
		SIMD_ALIGN
		double coeff[4 * num_pieces + 1];
		double *coeff_unaligned = coeff + 1;
		assert(!LIBSAKURA_SYMBOL(IsAligned)(coeff_unaligned));
		SIMD_ALIGN
		double boundary[2] = { 0.0f, 4.0f };
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, num_data, in_data, num_pieces, coeff_unaligned,
						boundary, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	{ // boundary is not aligned
		size_t const num_data(8);
		SIMD_ALIGN
		float in_data[num_data];
		for (size_t i = 0; i < num_data; ++i) {
			float x = (float) i;
			in_data[i] = 1 + x - 1 * x * x + 0.5 * x * x * x;
		}
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

		if (verbose) {
			PrintArray("in_data", num_data, in_data);
		}

		size_t order = 2;
		LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
				&context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		size_t num_pieces = 2;
		SIMD_ALIGN
		double coeff[8] = { 1.0f, 1.0f, -1.0f, 0.5f, 1.0f, 1.0f, -1.0f, 0.5f };
		SIMD_ALIGN
		double boundary[num_pieces + 1];
		double *boundary_unaligned = boundary + 1;
		assert(!LIBSAKURA_SYMBOL(IsAligned)(boundary_unaligned));
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, num_data, in_data, num_pieces, coeff,
						boundary_unaligned, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	{ // out is not aligned
		size_t const num_data(8);
		SIMD_ALIGN
		float in_data[num_data];
		for (size_t i = 0; i < num_data; ++i) {
			float x = (float) i;
			in_data[i] = 1 + x - 1 * x * x + 0.5 * x * x * x;
		}
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data) + 1];
		float *out_unaligned = out + 1;
		assert(!LIBSAKURA_SYMBOL(IsAligned)(out_unaligned));
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

		if (verbose) {
			PrintArray("in_data", num_data, in_data);
		}

		size_t order = 2;
		LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
				&context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		size_t num_pieces = 2;
		SIMD_ALIGN
		double coeff[8] = { 1.0f, 1.0f, -1.0f, 0.5f, 1.0f, 1.0f, -1.0f, 0.5f };
		SIMD_ALIGN
		double boundary[2] = { 0.0f, 4.0f };
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, num_data, in_data, num_pieces, coeff, boundary,
						out_unaligned);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

/*
 * Test sakura_SubtractBaselineCubicSplineUsingCoefficientsFloatWithNullPointer
 * failure case : data/coeff/out is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(BaselineSK, SubtractBaselineCubicSplineUsingCoefficientsFloatWithNullPointer) {
	{ // data is nullpointer
		size_t const num_data(8);

		SIMD_ALIGN
		float *in_data = nullptr;
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

		if (verbose) {
			PrintArray("in_data", num_data, in_data);
		}

		size_t order = 2;
		LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
				&context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		SIMD_ALIGN
		size_t num_pieces = 2;
		double coeff[8] = { 1.0f, 1.0f, -1.0f, 0.5f, 1.0f, 1.0f, -1.0f, 0.5f };
		double boundary[2] = { 0.0f, 4.0f };
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, num_data, in_data, num_pieces, coeff, boundary,
						out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	{ // coeff is nullpointer
		size_t const num_data(8);

		SIMD_ALIGN
		float in_data[num_data];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), in_data);
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

		if (verbose) {
			PrintArray("in_data", num_data, in_data);
		}

		size_t order = 2;
		LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
				&context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		SIMD_ALIGN
		size_t num_pieces = 2;
		double *coeff = nullptr;
		double boundary[2] = { 0.0f, 4.0f };
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, num_data, in_data, num_pieces, coeff, boundary,
						out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	{ // boundary is nullpointer
		size_t const num_data(8);

		SIMD_ALIGN
		float in_data[num_data];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), in_data);
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

		size_t order = 2;
		LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
				&context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		SIMD_ALIGN
		size_t num_pieces = 2;
		double coeff[4 * num_pieces];
		//SetFloatConstant(0.0f, (4*num_pieces), coeff);
		double *boundary = nullptr;
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, num_data, in_data, num_pieces, coeff, boundary,
						out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	{ // out is nullpointer
		size_t const num_data(8);

		SIMD_ALIGN
		float in_data[num_data];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), in_data);
		SIMD_ALIGN
		float *out = nullptr;
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

		size_t order = 2;
		LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
				&context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		SIMD_ALIGN
		size_t num_pieces = 2;
		double coeff[4 * num_pieces];
		//SetFloatConstant(0.0f, 4*num_pieces, coeff);
		double boundary[2] = { 0.0f, 4.0f };
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, num_data, in_data, num_pieces, coeff, boundary,
						out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	{ // context is nullpointer
		size_t const num_data(8);

		SIMD_ALIGN
		float in_data[num_data];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), in_data);
		SIMD_ALIGN
		float *out = nullptr;
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

		LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;

		SIMD_ALIGN
		size_t num_pieces = 2;
		double coeff[4 * num_pieces];
		//SetFloatConstant(0.0f, 4*num_pieces, coeff);
		double boundary[2] = { 0.0f, 4.0f };
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, num_data, in_data, num_pieces, coeff, boundary,
						out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), destroy_status);
	}
}

/*
 * Test sakura_SubtractBaselineCubicSplineUsingCoefficientsFloatInvalidArguments
 * failure case :  num_data != context->num_basis_data
 * failure case :  num_data < context->num_bases
 * failure case :
 * returned value : Status_kInvalidArgument
 */
TEST_F(BaselineSK, SubtractBaselineCubicSplineUsingCoefficientsFloatInvalidArguments) {
	{ // num_data != context->num_basis_data
		size_t const num_data(8);
		size_t const num_pieces(size_t(INT_MAX) + 1);

		SIMD_ALIGN
		float in_data[num_data];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), in_data);
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

		size_t order = 2;
		LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
				&context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		SIMD_ALIGN
		double coeff[4];
		double boundary[2] = { 0.0f, 4.0f };
		size_t const bad_num_data(num_data + 1);
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, bad_num_data, in_data, num_pieces, coeff,
						boundary, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	{ // num_data < context->num_bases
		size_t const num_data(8);
		size_t const num_pieces(size_t(INT_MAX) + 1);

		SIMD_ALIGN
		float in_data[num_data];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), in_data);
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

		size_t order = 2;
		LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
				&context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		SIMD_ALIGN
		double coeff[4];
		double boundary[2] = { 0.0f, 4.0f };
		size_t bad_num_data = context->num_bases - 1;
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, bad_num_data, in_data, num_pieces, coeff,
						boundary, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	{ // num_pieces > INT_MAX
		size_t const num_data(8);
		size_t const num_pieces(size_t(INT_MAX) + 1);

		SIMD_ALIGN
		float in_data[num_data];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), in_data);
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

		size_t order = 2;
		LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
				&context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		SIMD_ALIGN
		//size_t num_pieces = 2;
		//double coeff[4*num_pieces];
		double coeff[4];
		//SetFloatConstant(0.0f, 4*num_pieces, coeff);
		double boundary[2] = { 0.0f, 4.0f };
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, num_data, in_data, num_pieces, coeff, boundary,
						out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

/*
 * Test sakura_SubtractBaselineCubicSplineUsingCoefficientsFloatPerformanceTest
 * successful case : out == 0
 * subtract best fit model from input data using input coeff
 */
TEST_F(BaselineSK, SubtractBaselineCubicSplineUsingCoefficientsFloatPerformanceTest) {
	size_t const num_data(300000);
	size_t const num_pieces(75000);
	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		float x = (float) i;
		in_data[i] = 0.00001 - 0.00001 * x + 0.00001 * x * x
				+ 0.00001 * x * x * x;
	}
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float answer[ELEMENTSOF(in_data)];
	SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);
	if (verbose) {
		PrintArray("in_data", num_data, in_data);
	}
	size_t order = 2;
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_data,
			&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	double coeff[num_data];
	SIMD_ALIGN
	double tmp_coeff[4] = { 0.00001, -0.00001, 0.00001, 0.00001 };
	for (size_t i = 0; i < num_data; ++i) {
		coeff[i] = tmp_coeff[i % 4];
	}
	if (verbose) {
		PrintArray("coeff", num_data, coeff);
	}
	SIMD_ALIGN
	double boundary[num_pieces];
	for (size_t i = 0; i < num_pieces; ++i) {
		double x = (double) i;
		boundary[i] = x * 4.0f;
	}
	if (verbose) {
		PrintArray("boundary", num_pieces, boundary);
	}
	size_t loop_max = 1000;
	double start_time = sakura_GetCurrentTime();
	for (size_t i = 0; i < loop_max; ++i) {
		LIBSAKURA_SYMBOL (Status) subbl_status =
				LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
						context, num_data, in_data, num_pieces, coeff, boundary,
						out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	}
	double end_time = sakura_GetCurrentTime();
	std::cout << "Elapsed Time: " << end_time - start_time << "sec\n";

	for (size_t i = 0; i < num_data; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}
	//verbose = true;
	if (verbose) {
		PrintArray("out   ", num_data, out);
		PrintArray("answer", num_data, answer);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

