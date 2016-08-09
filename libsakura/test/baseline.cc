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
/*
 * baseline.cc
 *
 *  Created on: 2013/11/11
 *      Author: wataru
 */

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <sys/time.h>
#include <random>
#include <stdio.h>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"
#include "baseline.h"
#include "testutil.h"

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
class Baseline: public ::testing::Test {
protected:

	Baseline() :
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

	// TODO Set constant XXX type values into an array
	template<typename T>
	void Set_XXX_Constant(T value, size_t const num_data, T *data) {
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

	template<typename T>
	void Print1DArray(char const *name, size_t print_length, T const *data,
			size_t start_idx = 0, bool print_name = true, bool newline = true) {
		if (print_name)
			cout << name << " = [";
		size_t index = start_idx + print_length - 1;
		if (typeid(float const*) == typeid(data)) {
			for (size_t i = start_idx; i < index; ++i)
				cout << data[i] << ", ";
		} else if (typeid(double const*) == typeid(data)) {
			for (size_t i = start_idx; i < index; ++i)
				cout << data[i] << ", ";
		} else if (typeid(bool const*) == typeid(data)) {
			for (size_t i = start_idx; i < index; ++i)
				cout << (data[i] ? "T" : "F") << ", ";
		}
		cout << data[index] << " ]";
		if (newline)
			cout << endl;
	}

	//given as 1D float array but actually stores (num_row * num_column) 2D data
	//for which column loop comes inside row loop.
	template<typename T>
	void PrintArray(char const *name, size_t num_row, size_t num_column,
			T const *data) {
		cout << name << " = [";
		for (size_t i = 0; i < num_row; ++i) {
			Print1DArray(name, num_column, data, num_column * i, false, false);
			if (i < num_row - 1)
				cout << ", " << endl;
		}
		cout << " ]" << endl;
	}

	bool verbose;
};

/*
 void Create(LIBSAKURA_SYMBOL(LSQFitContextFloat) * context, size_t status,
 uint16_t const order, uint16_t const npiece, uint16_t const nwave,
 size_t const num_chan) {
 LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateLSQFitContextFloat(
 LIBSAKURA_SYMBOL(LSQFitType_kNumElements), order, npiece, nwave,
 num_chan, &context);
 EXPECT_EQ(status, create_status);
 std::cout << "Create status : " << create_status << std::endl;
 }
 */
/*
 *
 //LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateLSQFitContextFloat(
 //		LIBSAKURA_SYMBOL(LSQFitType_kNumElements), order, num_chan,
 //		&context);
 //EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), create_status);



 LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateLSQFitContextFloat(
 LIBSAKURA_SYMBOL(LSQFitType_kNumElements), order, num_chan,
 &context);
 EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), create_status);


 LIBSAKURA_SYMBOL (Status) create_status =
 LIBSAKURA_SYMBOL(CreateLSQFitContextFloat)(
 LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
 &context);
 EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
 */

//Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
void Destroy(LIBSAKURA_SYMBOL(LSQFitContextFloat) *context, size_t status) {
	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(status, destroy_status);
	std::cout << "Destroy Status : " << destroy_status << std::endl;
}

/*
 * Test sakura_CreateLSQFitContextFloatWithPolynomialPerformanceTest
 * successful case (with normal polynomial model)
 */
TEST_F(Baseline, CreateLSQFitContextFloatWithPolynomialPerformanceTest) {
	uint16_t const order(8000);
	size_t const num_chan(65535);

	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;

	size_t const num_repeat(1);
	double start, end;
	double elapsed_time = 0.0;
	for (size_t i = 0; i < num_repeat; ++i) {
		start = GetCurrentTime();
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(
						LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order,
						num_chan, &context);
		end = GetCurrentTime();
		elapsed_time += (end - start);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
	}
	std::cout << std::setprecision(5)
			<< "#x# benchmark Baseline_CreateLSQFitContextPolynomialFloatWithPolynomialPerformanceTest"
			<< " " << elapsed_time << std::endl;
}

/*
 * Test sakura_CreateLSQFitContextFloatWithChebyshevPolynomial
 * successful case (with Chebyshev polynomial model)
 */
TEST_F(Baseline, CreateLSQFitContextFloatWithChebyshevPolynomial) {
	uint16_t const order(7000);
	size_t const num_chan(65535);

	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	double start = GetCurrentTime();
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_chan,
					&context);
	double end = GetCurrentTime();
	cout << "Elapsed Time: " << (end - start) << " sec." << endl;
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_CreateLSQFitContextFloatWithZeroNumPieces
 * failure case : npiece = 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, CreateLSQFitContextFloatWithZeroNumPieces) {
	uint16_t const npiece(0);
	size_t const num_chan(4096);

	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextCubicSplineFloat(npiece, num_chan,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), create_status);
	EXPECT_EQ(nullptr, context);
}

/*
 * Test sakura_CreateLSQFitContextFloatWithInvalidLSQFitType
 * failure case : invalid baseline type
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, CreateLSQFitContextFloatWithInvalidLSQFitType) {
	uint16_t const order(20);
	size_t const num_chan(4096);

	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kNumElements), order,
					num_chan, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), create_status);
	EXPECT_EQ(nullptr, context);
}

/*
 * Test sakura_CreateLSQFitContextFloatWithContextNullPointer
 * failure case : context is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, CreateLSQFitContextFloatWithContextNullPointer) {
	uint16_t const order(20);
	size_t const num_chan(4096);

	LIBSAKURA_SYMBOL(LSQFitContextFloat) * *context_ptr_ptr = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_chan,
					context_ptr_ptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), create_status);
}

/*
 * Test sakura_CreateLSQFitContextFloatWithOrderLargerThanNumData
 * failure case : order is too large (> num_data)
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, CreateLSQFitContextFloatWithOrderLargerThanNumData) {
	uint16_t const order(20);
	size_t const num_chan(10);

	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_chan,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), create_status);
	EXPECT_EQ(nullptr, context);

	Destroy(context, LIBSAKURA_SYMBOL(Status_kInvalidArgument));
}

/*
 * Test sakura_CreateLSQFitContextFloatWithTooSmallNumData
 * failure case : num_data is smaller than minimum acceptable value
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, CreateLSQFitContextFloatWithTooSmallNumData) {
	uint16_t const order(20);
	uint16_t const npiece(1);

	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status;

	size_t const num_func_types = 3;
	LSQFitTypeInternal func_types[num_func_types] = {
			LSQFitTypeInternal_kPolynomial, LSQFitTypeInternal_kChebyshev,
			LSQFitTypeInternal_kCubicSpline };
	string func_names[num_func_types] = { "Polynomial", "Chebyshev",
			"CubicSpline" };
	for (size_t i = 0; i < num_func_types; ++i) {
		std::cout << "testing for " << func_names[i] << "..." << std::flush;
		size_t num_data_min = 0;
		switch (func_types[i]) {
		case (LSQFitTypeInternal_kPolynomial):
		case (LSQFitTypeInternal_kChebyshev):
			num_data_min = order + 1;
			break;
		case (LSQFitTypeInternal_kCubicSpline):
			num_data_min = 4;
			break;
		default:
			assert(false);
		};
		size_t num_data = num_data_min - 1;
		if (func_types[i] == LSQFitTypeInternal_kCubicSpline) {
			create_status = sakura_CreateLSQFitContextCubicSplineFloat(npiece,
					num_data, &context);
		} else {
			LIBSAKURA_SYMBOL(LSQFitType) type_ext = LIBSAKURA_SYMBOL(
					LSQFitType_kPolynomial);
			if (func_types[i] == LSQFitTypeInternal_kChebyshev) {
				type_ext = LIBSAKURA_SYMBOL(LSQFitType_kChebyshev);
			}
			create_status = sakura_CreateLSQFitContextPolynomialFloat(
					type_ext, order, num_data, &context);
		}
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), create_status);
		EXPECT_EQ(nullptr, context);
		Destroy(context, LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		std::cout << "ok." << std::flush << std::endl;
	}
}

/*
 * Test sakura_GetBaselineModelPolynomial
 * successful case (for normal polynomial model)
 */
TEST_F(Baseline, GetBaselineModelPolynomial) {
	uint16_t const order(20);
	size_t const num_chan(4096);
	size_t const num_model(order + 1);
	double out[num_chan * num_model];
	double answer[ELEMENTSOF(out)];

	// Setup answer----------------------------
	// Note that the model values are normalized to be 1 at the right end
	size_t idx = 0;
	for (size_t i = 0; i < num_chan; ++i) {
		double factor = static_cast<double>(i)/static_cast<double>(num_chan-1);
		for (size_t j = 0; j < num_model; ++j) {
			double value = 1.0;
			for (size_t k = 0; k < j; ++k) {
				value *= factor;
			}

			answer[idx] = value;
			idx++;
		}
	}
	//---------------------------------

	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_chan,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	for (size_t i = 0; i < ELEMENTSOF(out); ++i) {
		CheckAlmostEqual(answer[i], context->basis_data[i], 1e-10);
	}
	LSQFitTypeInternal type = context->lsqfit_type;
	EXPECT_EQ(LSQFitTypeInternal_kPolynomial, type);
	size_t num_bases = context->num_bases;
	ASSERT_EQ(num_model, num_bases);
	size_t num_basis_data = context->num_basis_data;
	ASSERT_EQ(num_chan, num_basis_data);

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetBaselineModelChebyshev
 * successful case (for Chebyshev polynomial model)
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

	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_chan,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	for (size_t i = 0; i < ELEMENTSOF(out); ++i) {
		CheckAlmostEqual(answer[i], context->basis_data[i], 1e-10);
	}

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_DestroyLSQFitContextFloat
 * successful case
 */
TEST_F(Baseline, DestroyLSQFitContextFloat) {
	uint16_t const order(20);
	size_t const num_chan(4096);

	double start, end;
	double elapsed_time = 0.0;
	size_t const num_repeat(NUM_REPEAT);
	for (size_t i = 0; i < num_repeat; ++i) {
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(
						LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order,
						num_chan, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		start = GetCurrentTime();
		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyLSQFitContextFloat(context);
		end = GetCurrentTime();
		elapsed_time += (end - start);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	cout << "Elapsed Time: " << elapsed_time << " sec." << endl;
}

/*
 * Test sakura_DestroyLSQFitContextFloat
 * failure case : context is a null pointer
 * returned value must be Status_kInvalidArgument
 */
TEST_F(Baseline, DestroyLSQFitContextFloatWithContextNullPointer) {
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;

	Destroy(context, LIBSAKURA_SYMBOL(Status_kInvalidArgument));
}

/*
 * Following the test code of sakura_SubtractBaseline
 * Test sakura_GetBestFitBaselineCoeffFromNormalDataWithoutClipping
 * successful case
 * the input data have smooth shape and no spiky feature, and
 * sakura_GetBestFitBaselineCoeff is executed without doing recursive
 * clipping.
 */
TEST_F(Baseline, GetBestFitBaselineCoeffFromSmoothDataWithoutClipping) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	double coeff_answer[num_coeff] = { 4.0, 3.0, 2.0, 1.0 };
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data, coeff_answer);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double coeff[num_coeff];

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, num_coeff,
			coeff, nullptr, nullptr, final_mask, &rms, &subbl_blstatus);

	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);

	for (size_t i = 0; i < num_coeff; ++i) {
		cout.precision();
		cout << "diff between coeff_answer and coeff "
				<< coeff_answer[i] - coeff[i] << endl;
		ASSERT_EQ((float )coeff_answer[i], (float )coeff[i]);
	}

	if (verbose) {
		Print1DArray("fmask ", num_data, final_mask);
	}

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Following the test code of sakura_SubtractBaseline
 * Test sakura_GetBestFitBaselineCoefficientsFloatPerformanceTest
 * successful case that will take a few seconds.
 * the input data have smooth shape and no spiky feature, and
 * sakura_GetBestFitBaselineCoeff is executed without doing recursive
 * clipping.
 */
TEST_F(Baseline, GetBestFitBaselineCoefficientsFloatPerformanceTest) {
	size_t const num_data(100000);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	double coeff_answer[num_coeff] = { pow(10.0, -10), pow(10.0, -10), pow(10.0,
			-10), pow(10.0, -10) };
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data, coeff_answer);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double coeff[num_coeff];

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;

	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);

	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	double elapsed_time = 0.0;
	size_t const num_repeat(1000);
	elapsed_time = 0.0;
	for (size_t i = 0; i < num_repeat; ++i) {
		double start = GetCurrentTime();
		LIBSAKURA_SYMBOL (Status) subbl_status =
		LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data,
				in_data, in_mask, clipping_threshold_sigma, num_fitting_max,
				num_coeff, coeff, nullptr, nullptr, final_mask, &rms,
				&subbl_blstatus);
		double end = GetCurrentTime();
		elapsed_time += (end - start);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	}
	std::cout << std::setprecision(5)
			<< "#x# benchmark Baseline_LSQFitPolynomialFloatPerformanceTest"
			<< " " << elapsed_time << std::endl;

	for (size_t i = 0; i < num_coeff; ++i) {
		cout.precision();
		cout << "diff between coeff_answer and coeff "
				<< coeff_answer[i] - coeff[i] << endl;
		EXPECT_LE(coeff_answer[i] - coeff[i], pow(10, -5));
	}

	if (verbose) {
		Print1DArray("fmask ", num_data, final_mask);
	}

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetBestFitBaselineCoeffFromNormalDataWithoutClippingWithDataNotAligned
 * failure case: data is not aligned
 * the input data have smooth shape and no spiky feature, and
 * sakura_GetBestFitBaselineCoeff is executed without doing recursive
 * clipping.
 */
TEST_F(Baseline, GetBestFitBaselineCoeffFromSmoothDataWithoutClippingWithDataNotAligned) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	double coeff_answer[num_coeff] = { 4.0, 3.0, 2.0, 1.0 };
	SIMD_ALIGN
	float in_data[num_data + 1];
	float *in_data_unaligned = in_data + 1;
	assert(!LIBSAKURA_SYMBOL(IsAligned)(in_data_unaligned));
	SetFloatPolynomial(num_data, in_data_unaligned, coeff_answer);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data_unaligned)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data_unaligned), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data_unaligned)];
	SIMD_ALIGN
	double coeff[num_coeff];
	float answer[ELEMENTSOF(in_data_unaligned)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data_unaligned), answer);
	if (verbose) {
		Print1DArray("in_data", num_data, in_data_unaligned);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;
	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) get_coeff_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data,
			in_data_unaligned, in_mask, clipping_threshold_sigma,
			num_fitting_max, num_coeff, coeff, nullptr, nullptr, final_mask,
			&rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), get_coeff_status);

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetBestFitBaselineCoeffFromSmoothDataWithoutClippingWithMaskNotAligned
 * failure case
 * the input data have smooth shape and no spiky feature, and
 * sakura_GetBestFitBaselineCoeff is executed without doing recursive
 * clipping.
 */
TEST_F(Baseline, GetBestFitBaselineCoeffFromSmoothDataWithoutClippingWithMaskNotAligned) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	double coeff_answer[num_coeff] = { 4.0, 3.0, 2.0, 1.0 };
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data, coeff_answer);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data) + 1];
	bool *in_mask_unaligned = in_mask + 1;
	assert(!LIBSAKURA_SYMBOL(IsAligned)(in_mask_unaligned));
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask_unaligned);
	in_mask_unaligned[3] = false;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double coeff[num_coeff];
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask_unaligned, clipping_threshold_sigma, num_fitting_max,
			num_coeff, coeff, nullptr, nullptr, final_mask, &rms,
			&subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetBestFitBaselineCoeffFromNormalDataWithoutClippingWithFinalMaskNotAligned
 * failure case
 * the input data have smooth shape and no spiky feature, and
 * sakura_GetBestFitBaselineCoeff is executed without doing recursive
 * clipping.
 */
TEST_F(Baseline, GetBestFitBaselineCoeffFromSmoothDataWithoutClippingWithFinalMaskNotAligned) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	double coeff_answer[num_coeff] = { 4.0, 3.0, 2.0, 1.0 };
	SIMD_ALIGN
	float in_data[num_data + 1];
	SetFloatPolynomial(num_data, in_data, coeff_answer);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data) + 1];
	bool *final_mask_unaligned = final_mask + 1;
	assert(!LIBSAKURA_SYMBOL(IsAligned)(final_mask_unaligned));
	SIMD_ALIGN
	double coeff[num_coeff];
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) get_final_mask_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, num_coeff,
			coeff, nullptr, nullptr, final_mask_unaligned, &rms,
			&subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), get_final_mask_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetBestFitBaselineCoeffFromNormalDataWithoutClippingWithCoeffNotAligned
 * failure case: coeff is not aligned
 * the input data have smooth shape and no spiky feature, and
 * sakura_GetBestFitBaselineCoeff is executed without doing recursive
 * clipping.
 */
TEST_F(Baseline, GetBestFitBaselineCoeffFromSmoothDataWithoutClippingWithCoeffNotAligned) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	double coeff_answer[num_coeff] = { 4.0, 3.0, 2.0, 1.0 };
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data, coeff_answer);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double coeff[num_coeff + 1];
	double *coeff_unaligned = coeff + 1;
	assert(!LIBSAKURA_SYMBOL(IsAligned)(coeff_unaligned));
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) get_coeff_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, num_coeff,
			coeff_unaligned, nullptr, nullptr, final_mask, &rms,
			&subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), get_coeff_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetBestFitBaselineCoeffFromNormalDataWithoutClippingWithDataNullPointer
 * failure case: data is a null pointer
 * the input data have smooth shape and no spiky feature, and
 * sakura_GetBestFitBaselineCoeff is executed without doing recursive
 * clipping.
 */
TEST_F(Baseline, GetBestFitBaselineCoeffFromSmoothDataWithoutClippingWithDataNullPointer) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	float *in_data = nullptr;
	SIMD_ALIGN
	bool in_mask[num_data];
	Set_XXX_Constant(true, ELEMENTSOF(in_mask), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_mask)];
	SIMD_ALIGN
	double coeff[num_coeff];
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) get_data_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, num_coeff,
			coeff, nullptr, nullptr, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), get_data_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetBestFitBaselineCoeffFromNormalDataWithoutClippingWithMaskNullPointer
 * failure case: mask is a null pointer
 * the input data have smooth shape and no spiky feature, and
 * sakura_GetBestFitBaselineCoeff is executed without doing recursive
 * clipping.
 */
TEST_F(Baseline, GetBestFitBaselineCoeffFromSmoothDataWithoutClippingWithMaskNullPointer) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	double coeff_answer[num_coeff] = { 4.0, 3.0, 2.0, 1.0 };
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data, coeff_answer);
	SIMD_ALIGN
	bool *in_mask = nullptr;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double coeff[num_coeff];
	float answer[ELEMENTSOF(in_data)];

	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) get_coeff_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, num_coeff,
			coeff, nullptr, nullptr, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), get_coeff_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetBestFitBaselineCoeffFromNormalDataWithoutClippingWithFinalMaskNullPointer
 * failure case: final_mask is a null pointer
 * the input data have smooth shape and no spiky feature, and
 * sakura_GetBestFitBaselineCoeff is executed without doing recursive
 * clipping.
 */
TEST_F(Baseline, GetBestFitBaselineCoeffFromSmoothDataWithoutClippingWithFinalMaskNullPointer) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	double coeff_answer[num_coeff] = { 4.0, 3.0, 2.0, 1.0 };
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data, coeff_answer);
	SIMD_ALIGN
	bool in_mask[num_data];
	Set_XXX_Constant(true, ELEMENTSOF(in_mask), in_mask);
	SIMD_ALIGN
	bool *final_mask = nullptr;
	SIMD_ALIGN
	double coeff[num_coeff];
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) get_coeff_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, num_coeff,
			coeff, nullptr, nullptr, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), get_coeff_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetBestFitBaselineCoeffFromNormalDataWithoutClippingWithCoeffNullPointer
 * failure case: coeff is a null pointer
 * the input data have smooth shape and no spiky feature, and
 * sakura_GetBestFitBaselineCoeff is executed without doing recursive
 * clipping.
 */
TEST_F(Baseline, GetBestFitBaselineCoeffFromSmoothDataWithoutClippingWithCoeffNullPointer) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	double coeff_answer[num_coeff] = { 4.0, 3.0, 2.0, 1.0 };
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data, coeff_answer);
	SIMD_ALIGN
	bool in_mask[num_data];
	Set_XXX_Constant(true, ELEMENTSOF(in_mask), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	double *coeff = nullptr;
	SIMD_ALIGN
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) get_coeff_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, num_coeff,
			coeff, nullptr, nullptr, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), get_coeff_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetBestFitBaselineCoeffWithZeroClipThreshold
 * failure case: clip_threshold_sigma == 0.0
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, GetBestFitBaselineCoeffWithZeroClipThreshold) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	double coeff_answer[num_coeff] = { 4.0, 3.0, 2.0, 1.0 };
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data, coeff_answer);
	SIMD_ALIGN
	bool in_mask[num_data];
	Set_XXX_Constant(true, ELEMENTSOF(in_mask), in_mask);
	SIMD_ALIGN
	double coeff[num_coeff];
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 0.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) get_coeff_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, num_coeff,
			coeff, nullptr, nullptr, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), get_coeff_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetBestFitBaselineCoeffWithNegativeClipThreshold
 * failure case: clip_threshold_sigma < 0.0
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, GetBestFitBaselineCoeffWithNegativeClipThreshold) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	double coeff_answer[num_coeff] = { 4.0, 3.0, 2.0, 1.0 };
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data, coeff_answer);
	SIMD_ALIGN
	bool in_mask[num_data];
	Set_XXX_Constant(true, ELEMENTSOF(in_mask), in_mask);
	SIMD_ALIGN
	double coeff[num_coeff];
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = -3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) get_coeff_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, num_coeff,
			coeff, nullptr, nullptr, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), get_coeff_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetBestFitBaselineCoeffWithZeroNumFittingMax
 * successful case: num_fitting_max == 0
 * in case num_fitting_max==0 is given, baseline subtraction
 * is not executed and values of coeff should not be modified.
 */
TEST_F(Baseline, GetBestFitBaselineCoeffWithZeroNumFittingMax) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	double coeff_answer[num_coeff] = { 4.0, 3.0, 2.0, 1.0 };
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data, coeff_answer);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	double coeff_orig[num_coeff] = { 1.1, 2.2, 3.3, 4.4 };
	SIMD_ALIGN
	double coeff[num_coeff];
	for (size_t i = 0; i < num_coeff; ++i) {
		coeff[i] = coeff_orig[i];
	}

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 0;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, num_coeff,
			coeff, nullptr, nullptr, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	for (size_t i = 0; i < num_coeff; ++i) {
		ASSERT_EQ(coeff_orig[i], coeff[i]);
	}
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineFromNormalDataWithoutClipping
 * successful case
 * the input data have smooth shape and no spiky feature, and
 * sakura_SubtractBaseline is executed without doing recursive
 * clipping.
 * the baseline-subtracted data should be zero throughout.
 */
TEST_F(Baseline, SubtractBaselineFromSmoothDataWithoutClipping) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		Print1DArray("fmask ", num_data, final_mask);
		Print1DArray("out   ", num_data, out);
		Print1DArray("answer", num_data, answer);
	}

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineFloatPerformanceTest
 * successful case
 * the input data have smooth shape and no spiky feature, and
 * sakura_SubtractBaseline is executed without doing recursive
 * clipping.
 * the baseline-subtracted data should be zero throughout.
 */
TEST_F(Baseline, SubtractBaselineFloatPerformanceTest) {
	size_t const num_data(NUM_DATA2 * 250);
	size_t const num_model(NUM_MODEL * 250);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	double start, end;
	double elapsed_time = 0.0;
	size_t const num_repeat(1);
	elapsed_time = 0.0;
	for (size_t i = 0; i < num_repeat; ++i) {
		start = GetCurrentTime();
		LIBSAKURA_SYMBOL (Status) subbl_status =
		LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data,
				in_data, in_mask, clipping_threshold_sigma, num_fitting_max,
				order + 1, nullptr, nullptr, out, final_mask, &rms,
				&subbl_blstatus);
		end = GetCurrentTime();
		elapsed_time += (end - start);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	}
	std::cout << std::setprecision(5)
			<< "#x# benchmark Baseline_LSQFitPolynomialFloatPerformanceTest"
			<< " " << elapsed_time << std::endl;

	if (verbose) {
		Print1DArray("fmask ", num_data, final_mask);
		Print1DArray("out   ", num_data, out);
		Print1DArray("answer", num_data, answer);
	}

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineFromSmoothDataWithClipping
 * successful case
 * the input data have smooth shape and no spiky feature, and
 * execute sakura_SubtractBaseline where recursive baseline fitting,
 * 10 times at maximum, with clipping data outside 3 sigma level.
 * since the input data haven't any outliers the baseline-subtracted
 * data should be identical with those in case without clipping,
 * namely, zero throughout.
 */
TEST_F(Baseline, SubtractBaselineFromSmoothDataWithClipping) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 10;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		Print1DArray("fmask ", num_data, final_mask);
		Print1DArray("out   ", num_data, out);
		Print1DArray("answer", num_data, answer);
	}

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineFromSpikyDataWithClipping
 * successful case
 * the input data have three outliers.
 * execute sakura_SubtractBaseline where recursive
 * baseline fitting, 5 times at maximum, with clipping data
 * outside 3 sigma level.
 * in the process of recursive clipping, the outliers will be
 * removed from baseline fitting procedure and the final result
 * should be zero throughout except for the outliers.
 */
TEST_F(Baseline, SubtractBaselineFromSpikyDataWithClipping) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	for (size_t i = 0; i < num_data; ++i) {
		if (i == 3)
			in_data[i] += 1000.0f;
		if (i == 6)
			in_data[i] -= 200.0f;
		if (i == 8)
			in_data[i] += 100.0f;
	}
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);
	for (size_t i = 0; i < num_data; ++i) {
		if (i == 3)
			answer[i] += 1000.0f;
		if (i == 6)
			answer[i] -= 200.0f;
		if (i == 8)
			answer[i] += 100.0f;
	}

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 5;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	for (size_t i = 0; i < num_data; ++i) {
		ASSERT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		Print1DArray("fmask ", num_data, final_mask);
		Print1DArray("out   ", num_data, out);
		Print1DArray("answer", num_data, answer);
	}

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithDataNullPointer
 * failure case : data is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithDataNullPointer) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	float *in_data = nullptr;
	SIMD_ALIGN
	bool in_mask[num_data];
	Set_XXX_Constant(true, ELEMENTSOF(in_mask), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_mask)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_mask)];

	if (verbose) {
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 5;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithDataNotAligned
 * failure case : data is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithDataNotAligned) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data + 1];
	float *in_data_unaligned = in_data + 1;
	assert(!LIBSAKURA_SYMBOL(IsAligned)(in_data_unaligned));
	SetFloatPolynomial(num_data, in_data_unaligned);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 5;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data,
			in_data_unaligned, in_mask, clipping_threshold_sigma,
			num_fitting_max, order + 1, nullptr, nullptr, out, final_mask, &rms,
			&subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithMaskNullPointer
 * failure case : mask is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithMaskNullPointer) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	bool *in_mask = nullptr;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 5;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithMaskNotAligned
 * failure case : mask is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithMaskNotAligned) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data) + 1];
	bool *in_mask_unaligned = in_mask + 1;
	assert(!LIBSAKURA_SYMBOL(IsAligned)(in_mask_unaligned));
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask_unaligned);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask_unaligned);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 5;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask_unaligned, clipping_threshold_sigma, num_fitting_max,
			order + 1, nullptr, nullptr, out, final_mask, &rms,
			&subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithLSQFitContextFloatNullPointer
 * failure case : context is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithLSQFitContextFloatNullPointer) {
	size_t const num_data(NUM_DATA2);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	uint16_t order(1);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 5;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
}

/*
 * Test sakura_SubtractBaselineWithNumDataNumBasisDataNotEqual
 * failure case : num_data != num_basis_data of baseline context
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithNumDataNumBasisDataNotEqual) {
	size_t const num_data_for_context(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data_for_context];
	SetFloatPolynomial(num_data_for_context, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		Print1DArray("in_data", num_data_for_context, in_data);
		Print1DArray("in_mask", num_data_for_context, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order,
					num_data_for_context, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	size_t const num_data = context->num_basis_data + 1;
	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 5;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithNumDataLessThanNumBases
 * failure case : num_data < num_bases of baseline context
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithNumDataLessThanNumBases) {
	size_t const num_data_for_context(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data_for_context];
	SetFloatPolynomial(num_data_for_context, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		Print1DArray("in_data", num_data_for_context, in_data);
		Print1DArray("in_mask", num_data_for_context, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order,
					num_data_for_context, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	size_t const num_data = context->num_bases - 1;
	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 5;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithFinalMaskNullPointer
 * failure case : final_mask is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithFinalMaskNullPointer) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	bool *final_mask = nullptr;
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithFinalMaskNotAligned
 * failure case : final_mask is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithFinalMaskNotAligned) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data) + 1];
	bool *final_mask_unaligned = final_mask + 1;
	assert(!LIBSAKURA_SYMBOL(IsAligned)(final_mask_unaligned));
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask_unaligned, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithOutNullPointer
 * failure case : out is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithOutNullPointer) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	float *out = nullptr;

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithOutNotAligned
 * failure case : out is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithOutNotAligned) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data) + 1];
	float *out_unaligned = out + 1;
	assert(!LIBSAKURA_SYMBOL(IsAligned)(out_unaligned));

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out_unaligned, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithLSQFitStatusNullPointer
 * failure case : baseline_status is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithLSQFitStatusNullPointer) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	LIBSAKURA_SYMBOL(LSQFitStatus) * subbl_blstatus_ptr = nullptr;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, subbl_blstatus_ptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithZeroClipThreshold
 * failure case : clip_threshold_sigma == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithZeroClipThreshold) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 0.0;
	uint16_t num_fitting_max = -5;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineWithNegativeClipThreshold
 * failure case : clip_threshold_sigma < 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineWithNegativeClipThreshold) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = -3.0;
	uint16_t num_fitting_max = 5;
	float rms;

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaseline
 * failure case : too many masked input data
 * returned value : Status_kNG
 */
TEST_F(Baseline, SubtractBaselineWithTooManyMaskedData) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(false, ELEMENTSOF(in_data), in_mask);
	in_mask[0] = true;
	in_mask[ELEMENTSOF(in_data) - 1] = true;

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 3;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float rms;

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kNG), subbl_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaseline
 * failure case : too many data clipped in the process of recursive
 * baseline fitting.
 * for this case, unmasked data should be less than the minimum needed
 * for fitting with the given baseline model just before the sixth
 * recursive baseline fitting.
 * returned value : Status_kNG
 */
TEST_F(Baseline, SubtractBaselineTooManyDataClipped) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data] = { 1.0, 2.5, 8.0, 1013.0, 21.5, 31.5, 42.5, 57.5,
			173.5, 90.5, 111.5, 132.5, 157.0, 182.5, 211.5 };
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	size_t order = 3;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 1.0;
	uint16_t num_fitting_max = 10;
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float rms;

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;
	LIBSAKURA_SYMBOL (Status) status = LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(
			context, order, num_data, in_data, in_mask,
			clipping_threshold_sigma, num_fitting_max, order + 1, nullptr,
			nullptr, out, final_mask, &rms, &subbl_blstatus);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kNG), status);
	EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kNotEnoughData), subbl_blstatus);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineFloatChebyshevPerformanceTest
 * successful case
 * subtract baseline from a realistic spectrum data with 3840 channels.
 * the data used here are taken from the first spectrum of the test
 * data used in the e2e test script for Sakura. clipping conditions
 * also are identical with those used in e2e test. the baseline model
 * here is a big set of Chebyshev polynomials with order up to 400.
 */
TEST_F(Baseline, SubtractBaselineFloatChebyshevPerformanceTest) {
	size_t const num_data(13840);
	size_t const num_model(400);

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<float> rnd(-0.01, 0.01);
	SIMD_ALIGN
	float in_data[num_data];
	float start_d = num_data / 2.0;
	for (size_t i = 0; i < num_data; ++i) {
		float x = i;
		in_data[i] = 0.001 * pow(x - start_d, 2.0) + rnd(mt);
	}

	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	Set_XXX_Constant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 2;
	float rms;

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
		Print1DArray("in_mask", num_data, in_mask);
	}

	LIBSAKURA_SYMBOL (LSQFitStatus) subbl_blstatus;

	double start, end;
	start = GetCurrentTime();
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, in_data,
			in_mask, clipping_threshold_sigma, num_fitting_max, order + 1,
			nullptr, nullptr, out, final_mask, &rms, &subbl_blstatus);
	end = GetCurrentTime();
	cout << "Elapsed Time: " << (end - start) << " sec." << endl;
	std::cout << std::setprecision(5)
			<< "#x# benchmark Baseline_SubtractBaselineFloatChebyshevPerformanceTest"
			<< " " << (end - start) << std::endl;
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	float limit_residual = 0.04f;
	for (size_t i = 0; i < num_data; ++i) {
		EXPECT_TRUE((out[i] > -limit_residual) && (out[i] < limit_residual));
	}

	if (verbose) {
		Print1DArray("fmask ", num_data, final_mask);
		Print1DArray("out   ", num_data, out);
		Print1DArray("answer", num_data, answer);
	}

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineUsingCoefficientsFloat
 * successful case :  : out[NUM_DATA2] == 0 with order is 2
 * subtract best fit model from input data using input coeff
 */
TEST_F(Baseline, SubtractBaselineUsingCoefficientsFloat) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);
	size_t const num_coeff(NUM_MODEL);

	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		in_data[i] = i * i + 2 * i + 3;
	}

	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	SIMD_ALIGN
	double coeff[num_coeff] = { 3.0f, 2.0f, 1.0f };
	LIBSAKURA_SYMBOL (Status) subbl_status =
	LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, num_data, in_data,
			num_coeff, coeff, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	for (size_t i = 0; i < num_data; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		Print1DArray("out   ", num_data, out);
		Print1DArray("answer", num_data, answer);
	}

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineUsingCoefficientsFloatWithNotAligned
 * failure case: data/coeff/out is not aligned
 * subtract best fit model from input data using input coeff
 */

TEST_F(Baseline, SubtractBaselineUsingCoefficientsFloatWithNotAligned) {
	{ // in_data is not aligned
		size_t const num_data(NUM_DATA2);
		size_t const num_model(NUM_MODEL);

		SIMD_ALIGN
		float in_data[num_data + 1];
		float *in_data_unaligned = in_data + 1;
		assert(!LIBSAKURA_SYMBOL(IsAligned)(in_data_unaligned));

		for (size_t i = 0; i < num_data; ++i) {
			in_data_unaligned[i] = 1.0 + i;
		}

		SIMD_ALIGN
		float out[ELEMENTSOF(in_data_unaligned)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data_unaligned)];
		Set_XXX_Constant(0.0f, ELEMENTSOF(in_data_unaligned), answer);

		if (verbose) {
			Print1DArray("in_data_unaligned", num_data, in_data_unaligned);
		}

		size_t order = num_model - 2;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(
						LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order,
						num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		size_t num_coeff = context->num_bases;
		SIMD_ALIGN
		double coeff[num_coeff];
		for (size_t i = 0; i < num_coeff; ++i) {
			coeff[i] = 1.0f;
		}
		LIBSAKURA_SYMBOL (Status) subbl_status =
		LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, num_data,
				in_data_unaligned, num_coeff, coeff, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
		Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
	}
	{ // coeff is not aligned
		size_t const num_data(NUM_DATA2);
		size_t const num_model(NUM_MODEL);

		SIMD_ALIGN
		float in_data[num_data];
		for (size_t i = 0; i < num_data; ++i) {
			in_data[i] = 1.0 + i;
		}

		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

		if (verbose) {
			Print1DArray("in_data", num_data, in_data);
		}

		size_t order = num_model - 2;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(
						LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order,
						num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		size_t num_coeff = context->num_bases;
		SIMD_ALIGN
		double coeff[num_coeff + 1];
		double *coeff_unaligned = coeff + 1;
		assert(!LIBSAKURA_SYMBOL(IsAligned)(coeff_unaligned));
		for (size_t i = 0; i < num_coeff; ++i) {
			coeff_unaligned[i] = 1.0f;
		}
		LIBSAKURA_SYMBOL (Status) subbl_status =
		LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, num_data, in_data,
				num_coeff, coeff_unaligned, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyLSQFitContextFloat(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	{ // out is not aligned
		size_t const num_data(NUM_DATA2);
		size_t const num_model(NUM_MODEL);

		SIMD_ALIGN
		float in_data[num_data];
		for (size_t i = 0; i < num_data; ++i) {
			in_data[i] = 1.0 + i;
		}
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data) + 1];
		float *out_unaligned = out + 1;
		assert(!LIBSAKURA_SYMBOL(IsAligned)(out_unaligned));
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

		if (verbose) {
			Print1DArray("in_data", num_data, in_data);
		}

		size_t order = num_model - 2;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(
						LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order,
						num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		size_t num_coeff = context->num_bases;
		SIMD_ALIGN
		double coeff[num_coeff];
		for (size_t i = 0; i < num_coeff; ++i) {
			coeff[i] = 1.0f;
		}
		LIBSAKURA_SYMBOL (Status) subbl_status =
		LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, num_data, in_data,
				num_coeff, coeff, out_unaligned);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
		Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
	}
}

/*
 * Test sakura_SubtractBaselineUsingCoefficientsFloatWithNullPointer
 * failure case : data/coeff/out is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineUsingCoefficientsFloatWithNullPointer) {
	{ // data is nullpointer
		size_t const num_data(NUM_DATA2);
		size_t const num_model(NUM_MODEL);
		float *in_data = nullptr;
		SIMD_ALIGN
		float out[num_data];
		size_t order = num_model - 2;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(
						LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order,
						num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		size_t num_coeff = context->num_bases;
		SIMD_ALIGN
		double coeff[num_coeff];
		for (size_t i = 0; i < num_coeff; ++i) {
			coeff[i] = 1.0f;
		}
		LIBSAKURA_SYMBOL (Status) subbl_status =
		LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, num_data, in_data,
				num_coeff, coeff, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
		Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
	}
	{ // coeff is nullpointer
		size_t const num_data(NUM_DATA2);
		size_t const num_model(NUM_MODEL);

		SIMD_ALIGN
		float in_data[num_data];
		for (size_t i = 0; i < num_data; ++i) {
			in_data[i] = 1.0 + i;
		}
		SIMD_ALIGN
		float out[num_data];
		size_t order = num_model - 2;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(
						LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order,
						num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		size_t num_coeff = context->num_bases;
		double *coeff = nullptr;
		LIBSAKURA_SYMBOL (Status) subbl_status =
		LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, num_data, in_data,
				num_coeff, coeff, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
		Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
	}
	{ // out is nullpointer
		size_t const num_data(NUM_DATA2);
		size_t const num_model(NUM_MODEL);
		SIMD_ALIGN
		float in_data[num_data];
		for (size_t i = 0; i < num_data; ++i) {
			in_data[i] = 1.0 + i;
		}
		float *out = nullptr;
		size_t order = num_model - 2;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(
						LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order,
						num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		size_t num_coeff = context->num_bases;
		SIMD_ALIGN
		double coeff[num_coeff];
		for (size_t i = 0; i < num_coeff; ++i) {
			coeff[i] = 1.0f;
		}
		LIBSAKURA_SYMBOL (Status) subbl_status =
		LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, num_data, in_data,
				num_coeff, coeff, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
		Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
	}
	{ // context is nullpointer
		size_t const num_data(NUM_DATA2);
		SIMD_ALIGN
		float in_data[num_data];
		for (size_t i = 0; i < num_data; ++i) {
			in_data[i] = 1.0 + i;
		}
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		size_t num_coeff = 2;
		SIMD_ALIGN
		double coeff[num_coeff];
		for (size_t i = 0; i < num_coeff; ++i) {
			coeff[i] = 1.0f;
		}
		LIBSAKURA_SYMBOL (Status) subbl_status =
		LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, num_data, in_data,
				num_coeff, coeff, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
		Destroy(context, LIBSAKURA_SYMBOL(Status_kInvalidArgument));
	}
}

/*
 * Test sakura_SubtractBaselineUsingCoefficientsFloatInvalidArguments
 * failure case :  num_data != context->num_basis_data
 * failure case :  num_data < context->num_bases
 * failure case :  num_coeff > context->num_bases
 * returned value : Status_kInvalidArgument
 */
TEST_F(Baseline, SubtractBaselineUsingCoefficientsFloatInvalidArguments) {
	{ // num_data != context->num_basis_data
		size_t const num_data(NUM_DATA2);
		size_t const num_model(NUM_MODEL);
		size_t const bad_num_data(NUM_DATA2 + 1);
		SIMD_ALIGN
		float in_data[num_data];
		for (size_t i = 0; i < num_data; ++i) {
			in_data[i] = 1.0 + i;
		}
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);
		size_t order = num_model - 2;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(
						LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order,
						num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		size_t num_coeff = context->num_bases;
		SIMD_ALIGN
		double coeff[num_coeff];
		for (size_t i = 0; i < num_coeff; ++i) {
			coeff[i] = 1.0f;
		}
		LIBSAKURA_SYMBOL (Status) subbl_status =
		LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, bad_num_data,
				in_data, num_coeff, coeff, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
		Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
	}
	{ // num_data < context->num_bases
		size_t const num_data(NUM_DATA2);
		size_t const num_model(NUM_MODEL);
		SIMD_ALIGN
		float in_data[num_data];
		for (size_t i = 0; i < num_data; ++i) {
			in_data[i] = 1.0 + i;
		}
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);
		size_t order = num_model - 2;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(
						LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order,
						num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		size_t num_coeff = context->num_bases;
		SIMD_ALIGN
		double coeff[num_coeff];
		for (size_t i = 0; i < num_coeff; ++i) {
			coeff[i] = 1.0f;
		}
		size_t bad_num_data = context->num_bases - 1;
		LIBSAKURA_SYMBOL (Status) subbl_status =
		LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, bad_num_data,
				in_data, num_coeff, coeff, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
		Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
	}
	{ // num_coeff > context->num_bases
		size_t const num_data(NUM_DATA2);
		size_t const num_model(NUM_MODEL);
		SIMD_ALIGN
		float in_data[num_data];
		for (size_t i = 0; i < num_data; ++i) {
			in_data[i] = 1.0 + i;
		}
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float answer[ELEMENTSOF(in_data)];
		Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);
		size_t order = num_model - 2;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(
						LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order,
						num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		size_t num_coeff = context->num_bases;
		SIMD_ALIGN
		double coeff[num_coeff];
		for (size_t i = 0; i < num_coeff; ++i) {
			coeff[i] = 1.0f;
		}
		size_t bad_num_coeff = context->num_bases + 1;
		LIBSAKURA_SYMBOL (Status) subbl_status =
		LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, num_data, in_data,
				bad_num_coeff, coeff, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
		Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
	}
}

/*
 * Test sakura_SubtractBaselineUsingCoefficientsFloatWithCoeffZeroPadding
 * successful case
 * subtract best fit model from input data using input coeff
 */
TEST_F(Baseline, SubtractBaselineUsingCoefficientsFloatWithCoeffZeroPadding) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	float in_data[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		in_data[i] = 1.0 + i;
	}
	SIMD_ALIGN
	float answer[ELEMENTSOF(in_data)];
	Set_XXX_Constant(0.0f, ELEMENTSOF(in_data), answer);

	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	size_t num_coeff = context->num_bases;
	SIMD_ALIGN
	double given_coeff[2] = { 1.0f, 1.0f };
	SIMD_ALIGN
	double coeff[num_coeff];
	coeff[0] = given_coeff[0];
	coeff[1] = given_coeff[1];
	coeff[2] = 0.0f;

	LIBSAKURA_SYMBOL (Status) subbl_status;
	subbl_status = LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, num_data,
			in_data, num_coeff, coeff, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);

	if (verbose) {
		Print1DArray("out", num_data, out);
	}
	for (size_t i = 0; i < num_data; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_SubtractBaselineUsingCoefficientsFloatPerformanceTest
 * successful case : out[NUM_DATA2] == 0
 * subtract best fit model from input data using input coeff
 */
TEST_F(Baseline, SubtractBaselineUsingCoefficientsFloatPerformanceTest) {
	size_t const num_data(500000);
	size_t const num_model(NUM_MODEL);
	size_t const num_coeff(NUM_MODEL);

	float *in_data = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_in_data(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*in_data) * num_data, &in_data));
	if (in_data == nullptr) {
		throw bad_alloc();
	}
	for (size_t i = 0; i < num_data; ++i) {
		in_data[i] = 2 * i + 3;
	}
	float *out = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_out(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*out) * num_data, &out));
	if (out == nullptr) {
		throw bad_alloc();
	}
	float *answer = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_answer(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*answer) * num_data, &answer));
	if (answer == nullptr) {
		throw bad_alloc();
	}
	Set_XXX_Constant(0.0f, num_data, answer);

	if (verbose) {
		Print1DArray("in_data", num_data, in_data);
	}

	size_t order = num_model - 1;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	SIMD_ALIGN
	double coeff[num_coeff] = { 3.0, 2.0 };

	size_t loop_max = 1000;
	double start_time, end_time;
	LIBSAKURA_SYMBOL (Status) subbl_status;
	start_time = GetCurrentTime();
	for (size_t i = 0; i < loop_max; ++i) {
		subbl_status = LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context,
				num_data, in_data, num_coeff, coeff, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
	}
	end_time = GetCurrentTime();
	std::cout << std::setprecision(5)
			<< "#x# benchmark Baseline_SubtractPolynomialFloatPerformanceTest"
			<< " " << (end_time - start_time) << std::endl;
	for (size_t i = 0; i < num_data; ++i) {
		EXPECT_EQ(answer[i], out[i]);
	}

	if (verbose) {
		Print1DArray("out   ", num_data, out);
		Print1DArray("answer", num_data, answer);
	}

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_GetNumberOfCoefficientsFloats
 * It can get number of Coefficients from context
 * successful case
 */
TEST_F(Baseline, GetNumberOfCoefficientsFloat) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL);
	size_t order = num_model - 2;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextPolynomialFloat(
					LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	size_t num_coeff = 0;
	LIBSAKURA_SYMBOL (Status) num_status = sakura_GetNumberOfCoefficientsFloat(
			context, order, &num_coeff);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), num_status);
	EXPECT_EQ(num_coeff, context->num_bases);

	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
}

/*
 * Test sakura_sakura_GetNumberOfCoefficientsFloatsWithNullPointer
 * It can get number of Coefficients from context
 * failed case : context is nullpointer
 */
TEST_F(Baseline, GetNumberOfCoefficientsFloatWithNullPointer) {
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	uint16_t order(1);
	size_t num_coeff = 0;
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument),
			LIBSAKURA_SYMBOL(GetNumberOfCoefficientsFloat)(context, order, &num_coeff));
}

