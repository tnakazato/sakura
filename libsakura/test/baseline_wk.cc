/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2015
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
 * Test cases to be implemented.
 * functions to be tested include:
 * - sakura_SubtractBaselineCubicSplineFloat
 * - sakura_GetBestFitBaselineCoefficientsCubicSplineFloat
 * test cases are as follows:
 * - for cubic spline functions:
 *     (1) simple successful case
 *     	 for cubic spline ones: (num_pieces=1,2,3, num_data=4*num_pieces+a where a=0,1,2,3,10)
 *     	 for sinusoidal ones: ()
 *     (2) time-consuming successful case for performance measurement
 *     (3) error cases
 *         (3-1) null pointer cases
 *         (3-2) non-aligned array cases
 *         (3-3) bad parameter value cases
 */

#include <cmath>
#include <iostream>
#include <random>
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
class BaselineWK: public ::testing::Test {
protected:

	BaselineWK() :
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

	//Set sinusoidal values of float into an array
	void SetFloatSinusoidal(size_t num_nwave, size_t const *nwave,
			double const *coeff, size_t num_data, float *data) {
		double factor = 2.0 * M_PI / (double) (num_data - 1);
		for (size_t i = 0; i < num_data; ++i) {
			data[i] = 0.0f;
			double x = factor * (double) i;
			size_t coeff_idx = 0;
			for (size_t j = 0; j < num_nwave; ++j) {
				//amplitude of each sinusoid is unity
				if (nwave[j] == 0) {
					data[i] += coeff[coeff_idx++];
				} else {
					double theta = nwave[j] * x;
					data[i] += (float) (coeff[coeff_idx++] * sin(theta));
					data[i] += (float) (coeff[coeff_idx++] * cos(theta));
				}
			}
		}
	}

	// Set constant float values into an array
	void SetFloatConstant(float value, size_t const num_data, float *data) {
		for (size_t i = 0; i < num_data; ++i) {
			data[i] = value;
		}
	}

	// Set constant float plus Gaussian noise values into an array
	void SetFloatConstantWithGaussianNoise(float value, float sigma,
			size_t const num_data, float *data) {
		std::random_device rd;
		std::mt19937 mt(rd());
		std::normal_distribution<> out(value, sigma);
		for (size_t i = 0; i < num_data; ++i) {
			data[i] = out(mt);
		}
	}
	// Set constant double values into an array
	void SetDoubleConstant(double value, size_t const num_data, double *data) {
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

	// Get number of sinusoidal coefficients from wave numbers to use
	size_t GetNumberOfSinusoidalCoefficients(size_t num_nwave,
			size_t const *nwave) {
		return (nwave[0] == 0) ? (2 * num_nwave - 1) : (2 * num_nwave);
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
	uint16_t const dummy_nwave = 0;

};

/*
 * Test GetBestFitBaselineCoefficientsCubicSpline
 * successful case
 * compute the best-fit baseline coefficients for cubic spline
 * for cases of combination of num_pieces=(1,2,3) and num_data=(4*num_pieces+a) where a=(0,1,2,3)
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplineSuccessfulCase) {
	uint16_t const dummy = 0;

	for (size_t num_pieces = 1; num_pieces <= 3; ++num_pieces) {
		cout << "    Testing for num_pieces = " << num_pieces
				<< " cases: num_data = ";
		size_t num_extra_max = 3;
		SIMD_ALIGN
		double boundary[num_pieces];
		SIMD_ALIGN
		double answer[4 * num_pieces];
		SetDoubleConstant(1.0, ELEMENTSOF(answer), answer);
		SIMD_ALIGN
		double out[ELEMENTSOF(answer)];
		float rms;
		LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

		for (size_t num_extra = 0; num_extra <= num_extra_max; ++num_extra) {
			size_t const num_data = ELEMENTSOF(answer) + num_extra;
			SIMD_ALIGN
			float data[num_data];
			SetFloatPolynomial(num_data, data, answer);
			SIMD_ALIGN
			bool mask[ELEMENTSOF(data)];
			SetBoolConstant(true, ELEMENTSOF(data), mask);
			cout << num_data << ((num_extra < num_extra_max) ? ", " : "");
			if (verbose) {
				PrintArray("data", num_data, data);
			}
			LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
			LIBSAKURA_SYMBOL (Status) create_status =
					sakura_CreateBaselineContext(
							LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), dummy,
							num_pieces, dummy, num_data, &context);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

			LIBSAKURA_SYMBOL (Status) coeff_status = LIBSAKURA_SYMBOL(
					GetBestFitBaselineCoefficientsCubicSplineFloat)(context,
					num_data, data, mask, 5.0f, 1, num_pieces, out, mask, &rms,
					boundary, &baseline_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

			for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
				CheckAlmostEqual(answer[i], out[i], 1.0e-6);
			}

			if (verbose) {
				PrintArray("data  ", ELEMENTSOF(data), data);
				PrintArray("out   ", ELEMENTSOF(out), out);
				PrintArray("answer", ELEMENTSOF(answer), answer);
			}

			LIBSAKURA_SYMBOL (Status) destroy_status =
					sakura_DestroyBaselineContext(context);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
		}
		cout << endl;
	}
}

/*
 * Test GetBestFitBaselineCoefficientsCubicSpline
 * successful cases with masked data
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplineSuccessfulCaseWithMaskedData) {
	enum MAItems {
		MA_kMiddle, MA_kLeft, MA_kRight, MA_kLeftRight, MA_kNumElems
	};
	vector<string> ma_param_names = { "middle", "left", "right", "leftright" };
	cout << "    Testing for ";

	uint16_t const dummy = 0;
	size_t num_pieces = 3;
	SIMD_ALIGN
	double boundary[num_pieces];
	SIMD_ALIGN
	double answer[4 * num_pieces];
	SetDoubleConstant(1.0, ELEMENTSOF(answer), answer);
	SIMD_ALIGN
	double out[ELEMENTSOF(answer)];
	float rms;
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;
	LIBSAKURA_SYMBOL (Status) coeff_status;

	size_t const num_data = 100;
	SIMD_ALIGN
	float data[num_data];
	SetFloatPolynomial(num_data, data, answer);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	if (verbose) {
		PrintArray("data", num_data, data);
	}
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), dummy, num_pieces,
			dummy, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	for (MAItems item = static_cast<MAItems>(0); item < MA_kNumElems; item =
			static_cast<MAItems>(item + 1)) {
		cout << ma_param_names[item] << ((item < MA_kNumElems - 1) ? ", " : "");
		size_t num_masked_idx = 10;

		num_masked_idx = (item == MA_kLeftRight) ? 20 : 10;
		size_t masked_idx[num_masked_idx];
		size_t ichan;
		switch (item) {
		case MA_kMiddle:
			for (size_t i = 0; i < num_masked_idx; ++i) {
				masked_idx[i] = (ELEMENTSOF(data) - num_masked_idx) / 2 + i;
			}
			break;
		case MA_kLeft:
			for (size_t i = 0; i < num_masked_idx; ++i) {
				masked_idx[i] = i;
			}
			break;
		case MA_kRight:
			for (size_t i = 0; i < num_masked_idx; ++i) {
				masked_idx[i] = ELEMENTSOF(data) - num_masked_idx + i;
			}
			break;
		case MA_kLeftRight:
			ichan = 0;
			for (; ichan < num_masked_idx / 2; ++ichan) {
				masked_idx[ichan] = ichan;
			}
			for (; ichan < num_masked_idx; ++ichan) {
				masked_idx[ichan] = ELEMENTSOF(data) - num_masked_idx + ichan;
			}
			break;
		default:
			assert(false);
			break;
		}
		for (size_t i = 0; i < num_masked_idx; ++i) {
			mask[masked_idx[i]] = false;
			data[masked_idx[i]] = -100000.0; // spike, which should not affect the fitting result
		}
		coeff_status = LIBSAKURA_SYMBOL(
				GetBestFitBaselineCoefficientsCubicSplineFloat)(context,
				num_data, data, mask, 5.0f, 1, num_pieces, out, mask, &rms,
				boundary, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

		if (verbose) {
			PrintArray("data  ", ELEMENTSOF(data), data);
			PrintArray("out   ", ELEMENTSOF(out), out);
			PrintArray("answer", ELEMENTSOF(answer), answer);
		}
		for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
			CheckAlmostEqual(answer[i], out[i], 1.0e-5);
		}
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test GetBestFitBaselineCoefficientsCubicSpline
 * time-consuming successful case for performance measurement
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplinePerformanceTest) {
	size_t const order = 0;
	size_t const num_repeat = 300;
	size_t const num_pieces = 2;
	SIMD_ALIGN
	double answer[4 * num_pieces];
	for (size_t i = 0; i < ELEMENTSOF(answer); i += 4) {
		answer[i] = 1.0;
		answer[i + 1] = 1e-4;
		answer[i + 2] = 1e-8;
		answer[i + 3] = 1e-12;
	}
	size_t const num_data = 70000;
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data, answer);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(in_data)];
	SetBoolConstant(true, ELEMENTSOF(in_data), mask);
	if (verbose) {
		PrintArray("in_data", num_data, in_data);
	}
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_pieces,
			dummy_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	double out[ELEMENTSOF(answer)];
	float rms;
	SIMD_ALIGN
	double boundary[num_pieces];
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;
	double start_time = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		LIBSAKURA_SYMBOL (Status) coeff_status = LIBSAKURA_SYMBOL(
				GetBestFitBaselineCoefficientsCubicSplineFloat)(context,
				num_data, in_data, mask, 5.0f, 1, num_pieces, out, mask, &rms,
				boundary, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
	}
	double end_time = sakura_GetCurrentTime();
	std::cout << std::setprecision(5)
			<< "#x# benchmark BaselineWK_GetBestFitBaselineCoefficientsCubicSplinePerformanceTest"
			<< " " << (end_time - start_time) << std::endl;

	EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

	if (verbose) {
		PrintArray("data  ", num_data, in_data);
		PrintArray("out   ", ELEMENTSOF(answer), out);
		PrintArray("answer", ELEMENTSOF(answer), answer);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test GetBestFitBaselineCoefficientsCubicSpline
 * erroneous cases: null pointer cases
 * parameters to be tested include context, data, mask, coeff, final_mask and baseline_status.
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplineErroneousCasesNullPointer) {
	enum NPItems {
		NP_kContext,
		NP_kData,
		NP_kMask,
		NP_kCoeff,
		NP_kFinalMask,
		NP_kBoundary,
		NP_kBaselineStatus,
		NP_kNumElems
	};
	vector<string> np_param_names = { "context", "data", "mask", "coeff",
			"final_mask", "boundary", "baseline_status" };
	cout << "    Testing for ";

	uint16_t const dummy = 0;
	size_t const num_pieces = 2;
	size_t const num_data = 10;
	double coeff[4];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	SIMD_ALIGN
	float data[num_data];
	SetFloatPolynomial(num_data, data, coeff);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	double out[ELEMENTSOF(coeff)];
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), dummy, num_pieces,
			dummy, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	float rms;
	SIMD_ALIGN
	double boundary[num_pieces];
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		float *data_ptr = data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		double *out_ptr = out;
		double *boundary_ptr = boundary;
		LIBSAKURA_SYMBOL(BaselineContext) *context_ptr = context;
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status_ptr = &baseline_status;

		switch (item) {
		case NP_kData:
			data_ptr = nullptr;
			break;
		case NP_kMask:
			mask_ptr = nullptr;
			break;
		case NP_kCoeff:
			out_ptr = nullptr;
			break;
		case NP_kFinalMask:
			final_mask_ptr = nullptr;
			break;
		case NP_kContext:
			context_ptr = nullptr;
			break;
		case NP_kBoundary:
			boundary_ptr = nullptr;
			break;
		case NP_kBaselineStatus:
			baseline_status_ptr = nullptr;
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL (Status) coeff_status = LIBSAKURA_SYMBOL(
				GetBestFitBaselineCoefficientsCubicSplineFloat)(context_ptr,
				num_data, data_ptr, mask_ptr, 5.0f, 1, num_pieces, out_ptr,
				final_mask_ptr, &rms, boundary_ptr, baseline_status_ptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test GetBestFitBaselineCoefficientsCubicSpline
 * erroneous cases: unaligned cases
 * parameters to be tested include data, mask, coeff and final_mask.
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplineErroneousCasesUnaligned) {
	enum UAItems {
		UA_kData, UA_kMask, UA_kCoeff, UA_kFinalMask, UA_kBoundary, UA_kNumElems
	};
	vector<string> ua_param_names = { "data", "mask", "coeff", "final_mask",
			"boundary" };
	cout << "    Testing for ";

	uint16_t const dummy = 0;
	size_t const num_pieces = 2;
	size_t const num_data = 10;
	double coeff[4];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	SIMD_ALIGN
	float data[num_data + 1];
	SetFloatPolynomial(num_data, data, coeff);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	double out[ELEMENTSOF(coeff) + 1];
	float rms;
	SIMD_ALIGN
	double boundary[num_pieces + 1];
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), dummy, num_pieces,
			dummy, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (UAItems item = static_cast<UAItems>(0); item < UA_kNumElems; item =
			static_cast<UAItems>(item + 1)) {
		cout << ua_param_names[item] << ((item < UA_kNumElems - 1) ? ", " : "");

		float *data_ptr = data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		double *out_ptr = out;
		double *boundary_ptr = boundary;

		switch (item) {
		case UA_kData:
			++data_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(data_ptr));
			break;
		case UA_kMask:
			++mask_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(mask_ptr));
			break;
		case UA_kCoeff:
			++out_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(out_ptr));
			break;
		case UA_kFinalMask:
			++final_mask_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(final_mask_ptr));
			break;
		case UA_kBoundary:
			++boundary_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(boundary_ptr));
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL (Status) coeff_status =
		LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsCubicSplineFloat)(
				context, num_data, data_ptr, mask_ptr, 5.0f, 1, num_pieces,
				out_ptr, final_mask_ptr, &rms, boundary_ptr, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test GetBestFitBaselineCoefficientsCubicSpline
 * erroneous cases: bad parameter value cases as follows:
 *     (1) num_data < context->num_bases
 *     (2) num_data < (!=) context->num_basis_data
 *     (3) num_data > (!=) context->num_basis_data
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplineErroneousCasesBadParameterValue) {
	enum BVItems {
		BV_kDataLTNumBases,
		BV_kDataLTNumBasisData,
		BV_kDataGTNumBasisData,
		BV_kNumElems
	};
	vector<string> bv_param_names = { "(num_data < context->num_bases)",
			"(num_data < context->num_basis_data)",
			"(num_data > context->num_basis_data)" };
	cout << "    Testing for cases " << endl;

	uint16_t const dummy = 0;
	size_t const num_basis_data = 10;
	size_t const num_boundary = 1;
	SIMD_ALIGN
	double boundary[num_boundary];
	double coeff[4 * num_boundary];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	SIMD_ALIGN
	double out[ELEMENTSOF(coeff)];
	float rms;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), dummy, num_boundary,
			dummy, num_basis_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (BVItems item = static_cast<BVItems>(0); item < BV_kNumElems; item =
			static_cast<BVItems>(item + 1)) {
		cout << "        " << bv_param_names[item]
				<< ((item < BV_kNumElems - 1) ? ", " : "") << endl;
		size_t num_data = 0;
		switch (item) {
		case BV_kDataLTNumBases:
			num_data = 2;
			break;
		case BV_kDataLTNumBasisData:
			num_data = 5;
			break;
		case BV_kDataGTNumBasisData:
			num_data = 15;
			break;
		default:
			assert(false);
		}

		SIMD_ALIGN
		float data[num_data];
		SetFloatPolynomial(num_data, data, coeff);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(data)];
		SetBoolConstant(true, ELEMENTSOF(data), mask);
		SIMD_ALIGN
		bool final_mask[ELEMENTSOF(data)];

		LIBSAKURA_SYMBOL (Status) coeff_status = LIBSAKURA_SYMBOL(
				GetBestFitBaselineCoefficientsCubicSplineFloat)(context,
				num_data, data, mask, 5.0f, 1, num_boundary, out, final_mask,
				&rms, boundary, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test sakura_SubtractBaselineCubicSplineFloat
 * successful case
 * subtract best fit model from input data
 */
TEST_F(BaselineWK, SubtractBaselineCubicSplineSuccessfulCase) {
	uint16_t const dummy = 0;
	double coeff[4];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	float rms;
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (size_t num_pieces = 1; num_pieces <= 3; ++num_pieces) {
		cout << "    Testing for num_pieces = " << num_pieces
				<< " cases: num_data = ";
		size_t num_extra_max = 3;
		SIMD_ALIGN
		double boundary[num_pieces];

		for (size_t num_extra = 0; num_extra <= num_extra_max; ++num_extra) {
			size_t const num_data = 4 * num_pieces + num_extra;
			cout << num_data << ((num_extra < num_extra_max) ? ", " : "");
			SIMD_ALIGN
			float data[num_data];
			SetFloatPolynomial(num_data, data, coeff);
			SIMD_ALIGN
			bool mask[ELEMENTSOF(data)];
			SetBoolConstant(true, ELEMENTSOF(data), mask);
			if (verbose) {
				PrintArray("data", num_data, data);
			}
			LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
			LIBSAKURA_SYMBOL (Status) create_status =
					sakura_CreateBaselineContext(
							LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), dummy,
							num_pieces, dummy, num_data, &context);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
			SIMD_ALIGN
			float out[ELEMENTSOF(data)];
			SIMD_ALIGN
			float answer[ELEMENTSOF(data)];
			SetFloatConstant(0.0, ELEMENTSOF(data), answer);

			LIBSAKURA_SYMBOL(Status) sub_status =
			LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(context,
					num_pieces, num_data, data, mask, 5.0f, 1, true, mask, out,
					&rms, boundary, &baseline_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

			for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
				CheckAlmostEqual(answer[i], out[i], 1.0e-6);
			}
			if (verbose) {
				PrintArray("data  ", ELEMENTSOF(data), data);
				PrintArray("out   ", ELEMENTSOF(out), out);
				PrintArray("answer", ELEMENTSOF(answer), answer);
			}

			LIBSAKURA_SYMBOL (Status) destroy_status =
					sakura_DestroyBaselineContext(context);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
		}
		cout << endl;
	}
}

/*
 * Test SubtractBaselineCubicSplineFloat
 * time-consuming successful case for performance measurement
 */
TEST_F(BaselineWK, SubtractBaselineCubicSplinePerformanceTest) {
	size_t const order = 0;
	size_t const num_repeat = 300;
	size_t const num_pieces = 2;
	size_t const num_data = 70000;
	SIMD_ALIGN
	float data[num_data];
	double coeff[4] = { 1.0, 1e-4, 1e-8, 1e-12 };
	SetFloatPolynomial(num_data, data, coeff);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	if (verbose) {
		PrintArray("data", num_data, data);
	}
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_pieces,
			dummy_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;
	SIMD_ALIGN
	double answer[ELEMENTSOF(data)];
	SetDoubleConstant(0.0, ELEMENTSOF(data), answer);
	float rms;
	SIMD_ALIGN
	double boundary[num_pieces];
	double start_time = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		LIBSAKURA_SYMBOL(Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(context, num_pieces,
				num_data, data, mask, 5.0f, 1, true, mask, out, &rms, boundary,
				&baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);
	}
	double end_time = sakura_GetCurrentTime();
	std::cout << std::setprecision(5)
			<< "#x# benchmark BaselineWK_SubtractBaselineCubicSplinePerformanceTest"
			<< " " << (end_time - start_time) << std::endl;
	EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

	if (verbose) {
		PrintArray("data  ", ELEMENTSOF(data), data);
		PrintArray("out   ", ELEMENTSOF(out), out);
		PrintArray("answer", ELEMENTSOF(answer), answer);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test SubtractBaselineCubicSplineFloat
 * erroneous cases: null pointer cases
 * parameters to be tested include context, data, mask, final_mask, out and baseline_status.
 */
TEST_F(BaselineWK, SubtractBaselineCubicSplineErroneousCasesNullPointer) {
	size_t const order = 0;
	enum NPItems {
		NP_kContext,
		NP_kData,
		NP_kMask,
		NP_kFinalMask,
		NP_kOut,
		NP_kBoundary,
		NP_kBaselineStatus,
		NP_kNumElems
	};
	vector<string> np_param_names = { "context", "data", "mask", "final_mask",
			"out", "boundary", "baseline_status" };
	cout << "    Testing for ";

	size_t const num_pieces = 2;
	size_t const num_data = 10;
	double coeff[4];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	SIMD_ALIGN
	float data[num_data];
	SetFloatPolynomial(num_data, data, coeff);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_pieces,
			dummy_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;
	float rms;
	SIMD_ALIGN
	double boundary[num_pieces];

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		LIBSAKURA_SYMBOL(BaselineContext) *context_ptr = context;
		float *data_ptr = data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		float *out_ptr = out;
		double *boundary_ptr = boundary;
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status_ptr = &baseline_status;

		switch (item) {
		case NP_kContext:
			context_ptr = nullptr;
			break;
		case NP_kData:
			data_ptr = nullptr;
			break;
		case NP_kMask:
			mask_ptr = nullptr;
			break;
		case NP_kFinalMask:
			final_mask_ptr = nullptr;
			break;
		case NP_kOut:
			out_ptr = nullptr;
			break;
		case NP_kBoundary:
			boundary_ptr = nullptr;
			break;
		case NP_kBaselineStatus:
			baseline_status_ptr = nullptr;
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL(Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(context_ptr,
				num_pieces, num_data, data_ptr, mask_ptr, 5.0f, 1, true,
				final_mask_ptr, out_ptr, &rms, boundary_ptr,
				baseline_status_ptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test SubtractBaselineCubicSplineFloat
 * erroneous cases: unaligned cases
 * parameters to be tested include data, mask, final_mask and out
 */
TEST_F(BaselineWK, SubtractBaselineCubicSplineErroneousCasesUnaligned) {
	size_t const order = 0;
	enum UAItems {
		UA_kData, UA_kMask, UA_kFinalMask, UA_kOut, UA_kBoundary, UA_kNumElems
	};
	vector<string> ua_param_names = { "data", "mask", "final_mask", "out",
			"boundary" };
	cout << "    Testing for ";

	size_t const num_pieces = 2;
	size_t const num_data = 10;
	double coeff[4];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	SIMD_ALIGN
	float data[num_data + 1];
	SetFloatPolynomial(num_data, data, coeff);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_pieces,
			dummy_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	float rms;
	SIMD_ALIGN
	double boundary[num_pieces + 1];
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (UAItems item = static_cast<UAItems>(0); item < UA_kNumElems; item =
			static_cast<UAItems>(item + 1)) {
		cout << ua_param_names[item] << ((item < UA_kNumElems - 1) ? ", " : "");

		float *data_ptr = data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		float *out_ptr = out;
		double *boundary_ptr = boundary;

		switch (item) {
		case UA_kData:
			++data_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(data_ptr));
			break;
		case UA_kMask:
			++mask_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(mask_ptr));
			break;
		case UA_kFinalMask:
			++final_mask_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(final_mask_ptr));
			break;
		case UA_kOut:
			++out_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(out_ptr));
			break;
		case UA_kBoundary:
			++boundary_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(boundary_ptr));
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL(Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(context, num_pieces,
				num_data, data_ptr, mask_ptr, 5.0f, 1, true, final_mask_ptr,
				out_ptr, &rms, boundary_ptr, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test SubtractBaselineCubicSplineFloat
 * erroneous cases: bad parameter value cases as follows:
 *     (1) num_data < context->num_bases
 *     (2) num_data < (!=) context->num_basis_data
 *     (3) num_data > (!=) context->num_basis_data
 */
TEST_F(BaselineWK, SubtractBaselineCubicSplineErroneousCasesBadParameterValue) {
	enum BVItems {
		BV_kDataLTNumBases,
		BV_kDataLTNumBasisData,
		BV_kDataGTNumBasisData,
		BV_kNumElems
	};
	vector<string> bv_param_names = { "(num_data < context->num_bases)",
			"(num_data < context->num_basis_data)",
			"(num_data > context->num_basis_data)" };
	cout << "    Testing for cases " << endl;

	size_t const num_basis_data = 10;
	uint16_t const dummy = 0;
	double coeff[4];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	float rms;
	size_t const num_pieces = 1;
	SIMD_ALIGN
	double boundary[num_pieces];
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), dummy, num_pieces,
			dummy, num_basis_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (BVItems item = static_cast<BVItems>(0); item < BV_kNumElems; item =
			static_cast<BVItems>(item + 1)) {
		cout << "        " << bv_param_names[item]
				<< ((item < BV_kNumElems - 1) ? ", " : "") << endl;
		size_t num_data = 0;
		switch (item) {
		case BV_kDataLTNumBases:
			num_data = 2;
			break;
		case BV_kDataLTNumBasisData:
			num_data = 5;
			break;
		case BV_kDataGTNumBasisData:
			num_data = 15;
			break;
		default:
			assert(false);
		}

		SIMD_ALIGN
		float data[num_data];
		SetFloatPolynomial(num_data, data, coeff);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(data)];
		SetBoolConstant(true, ELEMENTSOF(data), mask);
		SIMD_ALIGN
		bool final_mask[ELEMENTSOF(data)];
		SIMD_ALIGN
		float out[ELEMENTSOF(data)];

		LIBSAKURA_SYMBOL(Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(context, num_pieces,
				num_data, data, mask, 5.0f, 1, true, final_mask, out, &rms,
				boundary, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test GetBestFitBaselineCoefficientsSinusoid
 * successful cases
 * compute the best-fit baseline coefficients for sinusoid
 * context is created with maximum nwave of 2 in all cases, that is,
 * context contains bases with nwave of 0(constant), 1(sin/cos)
 * and 2(sin/cos).
 * test cases: all combinations of nwave between 0 and 2.
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsSinusoidSuccessfulCase) {
	uint16_t const dummy = 1;
	size_t const num_data = NUM_DATA2;
	size_t const context_nwave = 3;
	size_t const num_nwave_set = 4;
	size_t const nwave_set_length[num_nwave_set] = { 4, 6, 4, 1 };
	size_t nwave_set1[4][1] = { { 0 }, { 1 }, { 2 }, { 3 } };
	size_t nwave_set2[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 },
			{ 1, 3 }, { 2, 3 } };
	size_t nwave_set3[4][3] = { { 0, 1, 2 }, { 0, 1, 3 }, { 0, 2, 3 },
			{ 1, 2, 3 } };
	size_t nwave_set4[1][4] = { { 0, 1, 2, 3 } };
	cout << "    Baseline context has sinusoidal bases with nwave up to 3..."
			<< endl;

	SIMD_ALIGN
	bool mask[num_data];
	SetBoolConstant(true, ELEMENTSOF(mask), mask);
	float rms;
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (size_t num_nwave = 1; num_nwave <= num_nwave_set; ++num_nwave) {
		cout << "    Testing for num_nwave = " << num_nwave << ": nwave = ";

		size_t const num_cases = nwave_set_length[num_nwave - 1];
		size_t nwave[num_nwave];
		for (size_t j = 0; j < num_cases; ++j) {
			cout << "{";
			for (size_t k = 0; k < num_nwave; ++k) {
				if (num_nwave == 1) {
					nwave[k] = nwave_set1[j][k];
				} else if (num_nwave == 2) {
					nwave[k] = nwave_set2[j][k];
				} else if (num_nwave == 3) {
					nwave[k] = nwave_set3[j][k];
				} else if (num_nwave == 4) {
					nwave[k] = nwave_set4[j][k];
				}
				cout << nwave[k];
				if (k < num_nwave - 1) {
					cout << ", ";
				}
			}
			cout << "}";
			if (j < num_cases - 1) {
				cout << ", ";
			}
			size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave,
					nwave);
			double answer[num_coeff];
			SetDoubleConstant(1.0, num_coeff, answer);
			SIMD_ALIGN
			float data[ELEMENTSOF(mask)];
			SetFloatSinusoidal(num_nwave, nwave, answer, num_data, data);
			SIMD_ALIGN
			double out[ELEMENTSOF(answer)];

			LIBSAKURA_SYMBOL (Status) coeff_status = LIBSAKURA_SYMBOL(
					GetBestFitBaselineCoefficientsSinusoidFloat)(context,
					num_data, data, mask, 5.0f, 1, num_nwave, nwave, num_coeff,
					out, mask, &rms, &baseline_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

			for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
				CheckAlmostEqual(answer[i], out[i], 1.0e-6);
			}
			if (verbose) {
				PrintArray("data  ", ELEMENTSOF(data), data);
				PrintArray("out   ", ELEMENTSOF(out), out);
				PrintArray("answer", ELEMENTSOF(answer), answer);
			}
		}
		cout << endl;
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test GetBestFitBaselineCoefficientsSinusoid
 * successful cases with masked data
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsSinusoidSuccessfulCaseWithMaskedData) {
	enum MAItems {
		MA_kMiddle, MA_kLeft, MA_kRight, MA_kLeftRight, MA_kNumElems
	};
	vector<string> ma_param_names = { "middle", "left", "right", "leftright" };
	cout << "    Testing for ";

	size_t const num_data = NUM_DATA2;
	uint16_t const dummy = 1;
	size_t const context_nwave = 2;
	size_t const num_nwave = 3;
	size_t nwave[num_nwave] = { 0, 1, 2 };
	size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave, nwave);
	double answer[num_coeff];
	SetDoubleConstant(1.0, num_coeff, answer);
	SIMD_ALIGN
	float data[num_data];
	SetFloatSinusoidal(num_nwave, nwave, answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	double out[ELEMENTSOF(answer)];
	float rms;
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;
	LIBSAKURA_SYMBOL (Status) coeff_status;

	for (MAItems item = static_cast<MAItems>(0); item < MA_kNumElems; item =
			static_cast<MAItems>(item + 1)) {
		cout << ma_param_names[item] << ((item < MA_kNumElems - 1) ? ", " : "");
		size_t num_masked_idx = 10;

		num_masked_idx = (item == MA_kLeftRight) ? 2 : 1;
		size_t masked_idx[num_masked_idx];
		size_t ichan;
		switch (item) {
		case MA_kMiddle:
			for (size_t i = 0; i < num_masked_idx; ++i) {
				masked_idx[i] = (ELEMENTSOF(data) - num_masked_idx) / 2 + i;
			}
			break;
		case MA_kLeft:
			for (size_t i = 0; i < num_masked_idx; ++i) {
				masked_idx[i] = i;
			}
			break;
		case MA_kRight:
			for (size_t i = 0; i < num_masked_idx; ++i) {
				masked_idx[i] = ELEMENTSOF(data) - num_masked_idx + i;
			}
			break;
		case MA_kLeftRight:
			ichan = 0;
			for (; ichan < num_masked_idx / 2; ++ichan) {
				masked_idx[ichan] = ichan;
			}
			for (; ichan < num_masked_idx; ++ichan) {
				masked_idx[ichan] = ELEMENTSOF(data) - num_masked_idx + ichan;
			}
			break;
		default:
			assert(false);
			break;
		}
		for (size_t i = 0; i < num_masked_idx; ++i) {
			mask[masked_idx[i]] = false;
			data[masked_idx[i]] = 10.0; // spike, which should not affect the fitting result
		}
		coeff_status = LIBSAKURA_SYMBOL(
				GetBestFitBaselineCoefficientsSinusoidFloat)(context, num_data,
				data, mask, 5.0f, 1, num_nwave, nwave, num_coeff, out, mask,
				&rms, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

		if (verbose) {
			PrintArray("data  ", ELEMENTSOF(data), data);
			PrintArray("out   ", ELEMENTSOF(out), out);
			PrintArray("answer", ELEMENTSOF(answer), answer);
		}
		for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
			CheckAlmostEqual(answer[i], out[i], 1.0e-5);
		}
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test GetBestFitBaselineCoefficientsSinusoid
 * time-consuming successful case for performance measurement
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsSinusoidPerformanceTest) {
	size_t const num_repeat = 100;
	size_t const num_data = 50000;
	uint16_t const dummy = 1;
	size_t const context_nwave = 10;
	size_t const num_nwave = 10;
	size_t nwave[num_nwave] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave, nwave);
	double answer[num_coeff];
	SetDoubleConstant(1.0, num_coeff, answer);
	SIMD_ALIGN
	float data[num_data];
	SetFloatSinusoidal(num_nwave, nwave, answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	double out[ELEMENTSOF(answer)];
	float rms;
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;
	LIBSAKURA_SYMBOL (Status) coeff_status;
	double start_time = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		coeff_status = LIBSAKURA_SYMBOL(
				GetBestFitBaselineCoefficientsSinusoidFloat)(context, num_data,
				data, mask, 5.0f, 1, num_nwave, nwave, num_coeff, out, mask,
				&rms, &baseline_status);
	}
	double end_time = sakura_GetCurrentTime();
	std::cout << std::setprecision(5)
			<< "#x# benchmark BaselineWK_GetBestFitBaselineCoefficientsSinusoidPerformanceTest"
			<< " " << (end_time - start_time) << std::endl;
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
	EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

	for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
		CheckAlmostEqual(answer[i], out[i], 1.0e-6);
	}
	if (verbose) {
		PrintArray("data  ", ELEMENTSOF(data), data);
		PrintArray("out   ", ELEMENTSOF(out), out);
		PrintArray("answer", ELEMENTSOF(answer), answer);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test GetBestFitBaselineCoefficientsSinusoid
 * erroneous cases: null pointer cases
 * parameters to be tested include context, data, mask, nwave, coeff, final_mask and baseline_status.
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsSinusoidErroneousCasesNullPointer) {
	enum NPItems {
		NP_kContext,
		NP_kData,
		NP_kMask,
		NP_kNwave,
		NP_kCoeff,
		NP_kFinalMask,
		NP_kBaselineStatus,
		NP_kNumElems
	};
	vector<string> np_param_names = { "context", "data", "mask", "nwave",
			"coeff", "final_mask", "baseline_status" };
	cout << "    Testing for ";

	size_t const num_data = 10;
	uint16_t const dummy = 1;
	size_t const context_nwave = 2;
	size_t const num_nwave = 3;
	size_t nwave[num_nwave] = { 0, 1, 2 };
	size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave, nwave);
	double coeff_answer[num_coeff];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff_answer), coeff_answer);
	SIMD_ALIGN
	float data[num_data];
	SetFloatSinusoidal(num_nwave, nwave, coeff_answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	double coeff[ELEMENTSOF(coeff_answer)];
	float rms;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		LIBSAKURA_SYMBOL(BaselineContext) *context_ptr = context;
		float *data_ptr = data;
		bool *mask_ptr = mask;
		size_t *nwave_ptr = nwave;
		bool *final_mask_ptr = final_mask;
		double *coeff_ptr = coeff;
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status_ptr = &baseline_status;

		switch (item) {
		case NP_kContext:
			context_ptr = nullptr;
			break;
		case NP_kData:
			data_ptr = nullptr;
			break;
		case NP_kMask:
			mask_ptr = nullptr;
			break;
		case NP_kNwave:
			nwave_ptr = nullptr;
			break;
		case NP_kCoeff:
			coeff_ptr = nullptr;
			break;
		case NP_kFinalMask:
			final_mask_ptr = nullptr;
			break;
		case NP_kBaselineStatus:
			baseline_status_ptr = nullptr;
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL (Status) coeff_status =
		LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsSinusoidFloat)(
				context_ptr, num_data, data_ptr, mask_ptr, 5.0f, 1, num_nwave,
				nwave_ptr, num_coeff, coeff_ptr, final_mask_ptr, &rms,
				baseline_status_ptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test GetBestFitBaselineCoefficientsSinusoid
 * erroneous cases: unaligned cases
 * parameters to be tested include data, mask, coeff and final_mask.
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsSinusoidErroneousCasesUnaligned) {
	enum UAItems {
		UA_kData, UA_kMask, UA_kCoeff, UA_kFinalMask, UA_kNumElems
	};
	vector<string> ua_param_names = { "data", "mask", "coeff", "final_mask" };
	cout << "    Testing for ";

	size_t const num_data = 10;
	uint16_t const dummy = 1;
	size_t const context_nwave = 2;
	size_t const num_nwave = 3;
	size_t nwave[num_nwave] = { 0, 1, 2 };
	size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave, nwave);
	double coeff_answer[num_coeff];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff_answer), coeff_answer);
	SIMD_ALIGN
	float data[num_data + 1];
	SetFloatSinusoidal(num_nwave, nwave, coeff_answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	double coeff[ELEMENTSOF(coeff_answer) + 1];
	float rms;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (UAItems item = static_cast<UAItems>(0); item < UA_kNumElems; item =
			static_cast<UAItems>(item + 1)) {
		cout << ua_param_names[item] << ((item < UA_kNumElems - 1) ? ", " : "");

		float *data_ptr = data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		double *coeff_ptr = coeff;

		switch (item) {
		case UA_kData:
			++data_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(data_ptr));
			break;
		case UA_kMask:
			++mask_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(mask_ptr));
			break;
		case UA_kCoeff:
			++coeff_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(coeff_ptr));
			break;
		case UA_kFinalMask:
			++final_mask_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(final_mask_ptr));
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL (Status) coeff_status =
		LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsSinusoidFloat)(context,
				num_data, data_ptr, mask_ptr, 5.0f, 1, num_nwave, nwave,
				num_coeff, coeff_ptr, final_mask_ptr, &rms, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test GetBestFitBaselineCoefficientsSinusoid
 * erroneous cases: bad parameter value cases as follows:
 *     (1) num_data < context->num_bases
 *     (2) num_data < (!=) context->num_basis_data
 *     (3) num_data > (!=) context->num_basis_data
 *     (4) clip_threshold_sigma == 0
 *     (5) clip_threshold_sigma < 0
 *     (6) num_nwave == 0
 *     (7) nwave with duplicated elements
 *     (8) nwave with not in ascending order
 *     (9) maximum nwave to create context < maximum given nwave
 *     (10) num_coeff > context->num_bases
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsSinusoidErroneousCasesBadParameterValue) {
	enum BVItems {
		BV_kDataLTNumBases,
		BV_kDataLTNumBasisData,
		BV_kDataGTNumBasisData,
		BV_kClipThresholdZero,
		BV_kClipThresholdNegative,
		BV_kNumNWaveZero,
		BV_kNWaveDuplicated,
		BV_kNWaveNotAscending,
		BV_kGivenTooLargeMaximumNWave,
		BV_kNumCoeffGTNumBases,
		BV_kNumElems
	};
	vector<string> bv_param_names = { "(num_data < context->num_bases)",
			"(num_data < context->num_basis_data)",
			"(num_data > context->num_basis_data)",
			"(clip_threshold_sigma == 0)", "(clip_threshold_sigma < 0)",
			"(num_nwave == 0)", "(nwave duplicated)", "(nwave not ascending)",
			"(context_max_nwave < given_max_nwave)",
			"(num_coeff > context->num_bases)" };
	cout << "    Testing for cases " << endl;

	size_t const num_basis_data = 10;
	uint16_t const dummy = 1;
	size_t const context_nwave = 2;
	float rms;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_basis_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (BVItems item = static_cast<BVItems>(0); item < BV_kNumElems; item =
			static_cast<BVItems>(item + 1)) {
		cout << "        " << bv_param_names[item]
				<< ((item < BV_kNumElems - 1) ? ", " : "") << endl;
		size_t num_data = num_basis_data;
		float clip_threshold_sigma = 5.0f;
		size_t num_nwave = 3;
		size_t nwave[num_nwave];
		for (size_t i = 0; i < num_nwave; ++i) {
			nwave[i] = i;
		}
		size_t tmp_nwave = nwave[0];
		size_t num_coeff = 5;
		size_t num_coeff_by_nwave = GetNumberOfSinusoidalCoefficients(num_nwave,
				nwave);

		switch (item) {
		case BV_kDataLTNumBases:
			num_data = 2;
			break;
		case BV_kDataLTNumBasisData:
			num_data = 5;
			break;
		case BV_kDataGTNumBasisData:
			num_data = 15;
			break;
		case BV_kClipThresholdZero:
			clip_threshold_sigma = 0.0f;
			break;
		case BV_kClipThresholdNegative:
			clip_threshold_sigma = -5.0f;
			break;
		case BV_kNumNWaveZero:
			num_nwave = 0;
			break;
		case BV_kNWaveDuplicated:
			nwave[1] = nwave[0];
			break;
		case BV_kNWaveNotAscending:
			tmp_nwave = nwave[1];
			nwave[1] = nwave[0];
			nwave[0] = tmp_nwave;
			break;
		case BV_kGivenTooLargeMaximumNWave:
			nwave[num_nwave - 1] = context_nwave + 10;
			break;
		case BV_kNumCoeffGTNumBases:
			num_coeff = num_coeff_by_nwave + 1;
			break;
		default:
			assert(false);
		}

		double coeff[num_coeff];
		SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
		SIMD_ALIGN
		float data[num_data];
		SetFloatSinusoidal(num_nwave, nwave, coeff, num_data, data);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(data)];
		SetBoolConstant(true, ELEMENTSOF(data), mask);
		SIMD_ALIGN
		bool final_mask[ELEMENTSOF(data)];
		SIMD_ALIGN
		double out[ELEMENTSOF(coeff)];

		LIBSAKURA_SYMBOL (Status) coeff_status =
		LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsSinusoidFloat)(context,
				num_data, data, mask, clip_threshold_sigma, 1, num_nwave, nwave,
				num_coeff, out, final_mask, &rms, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test SubtractBaselineSinusoid
 * successful cases
 * compute the best-fit sinusoidal baseline and subtract it
 * context is created with maximum nwave of 2 in all cases, that is,
 * context contains bases with nwave of 0(constant), 1(sin/cos)
 * and 2(sin/cos).
 * test cases: all combinations of nwave between 0 and 2.
 */
TEST_F(BaselineWK, SubtractBaselineSinusoidSuccessfulCase) {
	uint16_t const dummy = 1;
	size_t const num_data = NUM_DATA2;
	size_t const context_nwave = 3;
	size_t const num_nwave_set = 4;
	size_t const nwave_set_length[num_nwave_set] = { 4, 6, 4, 1 };
	size_t nwave_set1[4][1] = { { 0 }, { 1 }, { 2 }, { 3 } };
	size_t nwave_set2[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 },
			{ 1, 3 }, { 2, 3 } };
	size_t nwave_set3[4][3] = { { 0, 1, 2 }, { 0, 1, 3 }, { 0, 2, 3 },
			{ 1, 2, 3 } };
	size_t nwave_set4[1][4] = { { 0, 1, 2, 3 } };
	cout << "    Baseline context has sinusoidal bases with nwave up to 3..."
			<< endl;

	SIMD_ALIGN
	bool mask[num_data];
	SetBoolConstant(true, ELEMENTSOF(mask), mask);
	SIMD_ALIGN
	float out[ELEMENTSOF(mask)];
	float answer[ELEMENTSOF(mask)];
	SetFloatConstant(0.0f, ELEMENTSOF(answer), answer);
	float rms;
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (size_t num_nwave = 1; num_nwave <= num_nwave_set; ++num_nwave) {
		cout << "    Testing for num_nwave = " << num_nwave << ": nwave = ";

		size_t const num_cases = nwave_set_length[num_nwave - 1];
		size_t nwave[num_nwave];
		for (size_t j = 0; j < num_cases; ++j) {
			cout << "{";
			for (size_t k = 0; k < num_nwave; ++k) {
				if (num_nwave == 1) {
					nwave[k] = nwave_set1[j][k];
				} else if (num_nwave == 2) {
					nwave[k] = nwave_set2[j][k];
				} else if (num_nwave == 3) {
					nwave[k] = nwave_set3[j][k];
				} else if (num_nwave == 4) {
					nwave[k] = nwave_set4[j][k];
				}
				cout << nwave[k];
				if (k < num_nwave - 1) {
					cout << ", ";
				}
			}
			cout << "}";
			if (j < num_cases - 1) {
				cout << ", ";
			}
			size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave,
					nwave);
			double coeff[num_coeff];
			SetDoubleConstant(1.0, num_coeff, coeff);
			SIMD_ALIGN
			float data[ELEMENTSOF(mask)];
			SetFloatSinusoidal(num_nwave, nwave, coeff, num_data, data);

			LIBSAKURA_SYMBOL (Status) sub_status =
			LIBSAKURA_SYMBOL(SubtractBaselineSinusoidFloat)(context, num_nwave,
					nwave, num_data, data, mask, 5.0f, 1,
					true, mask, out, &rms, &baseline_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

			for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
				CheckAlmostEqual(answer[i], out[i], 1.0e-6);
			}
			if (verbose) {
				PrintArray("data  ", ELEMENTSOF(data), data);
				PrintArray("out   ", ELEMENTSOF(out), out);
				PrintArray("answer", ELEMENTSOF(answer), answer);
			}
		}
		cout << endl;
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test SubtractBaselineSinusoid
 * successful cases with masked data
 */
TEST_F(BaselineWK, SubtractBaselineSinusoidSuccessfulCaseWithMaskedData) {
	enum MAItems {
		MA_kMiddle, MA_kLeft, MA_kRight, MA_kLeftRight, MA_kNumElems
	};
	vector<string> ma_param_names = { "middle", "left", "right", "leftright" };
	cout << "    Testing for ";

	size_t const num_data = NUM_DATA2;
	uint16_t const dummy = 1;
	size_t const context_nwave = 2;
	size_t const num_nwave = 3;
	size_t nwave[num_nwave] = { 0, 1, 2 };
	size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave, nwave);
	double coeff[num_coeff];
	SetDoubleConstant(1.0, num_coeff, coeff);
	SIMD_ALIGN
	float data[num_data];
	SetFloatSinusoidal(num_nwave, nwave, coeff, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	float answer[ELEMENTSOF(data)];
	for (size_t i = 0; i < ELEMENTSOF(coeff); ++i) {
		SetFloatConstant(0.0f, ELEMENTSOF(answer), answer);
	}
	float rms;
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;
	LIBSAKURA_SYMBOL (Status) sub_status;

	for (MAItems item = static_cast<MAItems>(0); item < MA_kNumElems; item =
			static_cast<MAItems>(item + 1)) {
		cout << ma_param_names[item] << ((item < MA_kNumElems - 1) ? ", " : "");
		size_t num_masked_idx = 10;

		num_masked_idx = (item == MA_kLeftRight) ? 2 : 1;
		size_t masked_idx[num_masked_idx];
		size_t ichan;
		switch (item) {
		case MA_kMiddle:
			for (size_t i = 0; i < num_masked_idx; ++i) {
				masked_idx[i] = (ELEMENTSOF(data) - num_masked_idx) / 2 + i;
			}
			break;
		case MA_kLeft:
			for (size_t i = 0; i < num_masked_idx; ++i) {
				masked_idx[i] = i;
			}
			break;
		case MA_kRight:
			for (size_t i = 0; i < num_masked_idx; ++i) {
				masked_idx[i] = ELEMENTSOF(data) - num_masked_idx + i;
			}
			break;
		case MA_kLeftRight:
			ichan = 0;
			for (; ichan < num_masked_idx / 2; ++ichan) {
				masked_idx[ichan] = ichan;
			}
			for (; ichan < num_masked_idx; ++ichan) {
				masked_idx[ichan] = ELEMENTSOF(data) - num_masked_idx + ichan;
			}
			break;
		default:
			assert(false);
			break;
		}
		for (size_t i = 0; i < num_masked_idx; ++i) {
			mask[masked_idx[i]] = false;
			data[masked_idx[i]] = 10.0; // spike, which should not affect the fitting result
		}
		sub_status = LIBSAKURA_SYMBOL(SubtractBaselineSinusoidFloat)(context,
				num_nwave, nwave, num_data, data, mask, 5.0f, 1, true, mask,
				out, &rms, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

		if (verbose) {
			PrintArray("data  ", ELEMENTSOF(data), data);
			PrintArray("out   ", ELEMENTSOF(out), out);
			PrintArray("answer", ELEMENTSOF(answer), answer);
		}
		for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
			if (mask[i]) {
				CheckAlmostEqual(answer[i], out[i], 1.0e-6);
			}
		}
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test SubtractBaselineSinusoid
 * time-consuming successful case for performance measurement
 */
TEST_F(BaselineWK, SubtractBaselineSinusoidPerformanceTest) {
	size_t const num_repeat = 100;
	size_t const num_data = 50000;
	uint16_t const dummy = 1;
	size_t const context_nwave = 10;
	size_t const num_nwave = 10;
	size_t nwave[num_nwave] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave, nwave);
	double coeff[num_coeff];
	SetDoubleConstant(1.0, num_coeff, coeff);
	SIMD_ALIGN
	float data[num_data];
	SetFloatSinusoidal(num_nwave, nwave, coeff, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	float answer[ELEMENTSOF(data)];
	for (size_t i = 0; i < ELEMENTSOF(coeff); ++i) {
		SetFloatConstant(0.0f, ELEMENTSOF(answer), answer);
	}
	float rms;
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;
	LIBSAKURA_SYMBOL (Status) sub_status;
	double start_time = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		sub_status = LIBSAKURA_SYMBOL(SubtractBaselineSinusoidFloat)(context,
				num_nwave, nwave, num_data, data, mask, 5.0f, 1, true, mask,
				out, &rms, &baseline_status);
	}
	double end_time = sakura_GetCurrentTime();
	std::cout << std::setprecision(5)
			<< "#x# benchmark BaselineWK_SubtractBaselineSinusoidPerformanceTest"
			<< " " << (end_time - start_time) << std::endl;
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);
	EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

	for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
		CheckAlmostEqual(answer[i], out[i], 5.0e-6);
	}
	if (verbose) {
		PrintArray("data  ", ELEMENTSOF(data), data);
		PrintArray("out   ", ELEMENTSOF(out), out);
		PrintArray("answer", ELEMENTSOF(answer), answer);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test SubtractBaselineSinusoid
 * erroneous cases: null pointer cases
 * parameters to be tested include context, data, mask, nwave, coeff, final_mask and baseline_status.
 */
TEST_F(BaselineWK, SubtractBaselineSinusoidErroneousCasesNullPointer) {
	enum NPItems {
		NP_kContext,
		NP_kNwave,
		NP_kData,
		NP_kMask,
		NP_kFinalMask,
		NP_kOut,
		NP_kBaselineStatus,
		NP_kNumElems
	};
	vector<string> np_param_names = { "context", "nwave", "data", "mask",
			"final_mask", "out", "baseline_status" };
	cout << "    Testing for ";

	size_t const num_data = 10;
	uint16_t const dummy = 1;
	size_t const context_nwave = 2;
	size_t const num_nwave = 3;
	size_t nwave[num_nwave] = { 0, 1, 2 };
	size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave, nwave);
	double coeff_answer[num_coeff];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff_answer), coeff_answer);
	SIMD_ALIGN
	float data[num_data];
	SetFloatSinusoidal(num_nwave, nwave, coeff_answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	float rms;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		LIBSAKURA_SYMBOL(BaselineContext) *context_ptr = context;
		float *data_ptr = data;
		bool *mask_ptr = mask;
		size_t *nwave_ptr = nwave;
		bool *final_mask_ptr = final_mask;
		float *out_ptr = out;
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status_ptr = &baseline_status;

		switch (item) {
		case NP_kContext:
			context_ptr = nullptr;
			break;
		case NP_kNwave:
			nwave_ptr = nullptr;
			break;
		case NP_kData:
			data_ptr = nullptr;
			break;
		case NP_kMask:
			mask_ptr = nullptr;
			break;
		case NP_kFinalMask:
			final_mask_ptr = nullptr;
			break;
		case NP_kOut:
			out_ptr = nullptr;
			break;
		case NP_kBaselineStatus:
			baseline_status_ptr = nullptr;
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL (Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineSinusoidFloat)(context_ptr, num_nwave,
				nwave_ptr, num_data, data_ptr, mask_ptr, 5.0f, 1, true,
				final_mask_ptr, out_ptr, &rms, baseline_status_ptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test SubtractBaselineSinusoidFloat
 * erroneous cases: unaligned cases
 * parameters to be tested include data, mask, final_mask and out
 */
TEST_F(BaselineWK, SubtractBaselineSinusoidErroneousCasesUnaligned) {
	enum UAItems {
		UA_kData, UA_kMask, UA_kFinalMask, UA_kOut, UA_kNumElems
	};
	vector<string> ua_param_names = { "data", "mask", "final_mask", "out" };
	cout << "    Testing for ";

	size_t const num_data = 10;
	uint16_t const dummy = 1;
	size_t const context_nwave = 2;
	size_t const num_nwave = 3;
	size_t nwave[num_nwave] = { 0, 1, 2 };
	size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave, nwave);
	double coeff_answer[num_coeff];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff_answer), coeff_answer);
	SIMD_ALIGN
	float data[num_data + 1];
	SetFloatSinusoidal(num_nwave, nwave, coeff_answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	float rms;
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (UAItems item = static_cast<UAItems>(0); item < UA_kNumElems; item =
			static_cast<UAItems>(item + 1)) {
		cout << ua_param_names[item] << ((item < UA_kNumElems - 1) ? ", " : "");

		float *data_ptr = data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		float *out_ptr = out;

		switch (item) {
		case UA_kData:
			++data_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(data_ptr));
			break;
		case UA_kMask:
			++mask_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(mask_ptr));
			break;
		case UA_kFinalMask:
			++final_mask_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(final_mask_ptr));
			break;
		case UA_kOut:
			++out_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(out_ptr));
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL(Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineSinusoidFloat)(context, num_nwave,
				nwave, num_data, data_ptr, mask_ptr, 5.0f, 1, true,
				final_mask_ptr, out_ptr, &rms, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test SubtractBaselineSinusoidFloat
 * erroneous cases: bad parameter value cases as follows:
 *     (1) num_data < context->num_bases
 *     (2) num_data < (!=) context->num_basis_data
 *     (3) num_data > (!=) context->num_basis_data
 *     (4) clip_threshold_sigma == 0
 *     (5) clip_threshold_sigma < 0
 *     (6) num_nwave == 0
 *     (7) nwave with duplicated elements
 *     (8) nwave with not in ascending order
 *     (9) max_nwave in context < given max_nwave
 *     (10) num_bases computed for given nwave > context->num_bases
 */
TEST_F(BaselineWK, SubtractBaselineSinusoidErroneousCasesBadParameterValue) {
	enum BVItems {
		BV_kDataLTNumBases,
		BV_kDataLTNumBasisData,
		BV_kDataGTNumBasisData,
		BV_kClipThresholdZero,
		BV_kClipThresholdNegative,
		BV_kNumNWaveZero,
		BV_kNWaveDuplicated,
		BV_kNWaveNotAscending,
		BV_kGivenTooLargeMaximumNWave,
		BV_kNumBasesByNWaveGTNumBases,
		BV_kNumElems
	};
	vector<string> bv_param_names = { "(num_data < context->num_bases)",
			"(num_data < context->num_basis_data)",
			"(num_data > context->num_basis_data)",
			"(clip_threshold_sigma == 0)", "(clip_threshold_sigma < 0)",
			"(num_nwave == 0)", "(nwave duplicated)", "(nwave not ascending)",
			"(context_max_nwave < given_max_nwave)",
			"(num_bases by given nwave > context->num_bases)" };
	cout << "    Testing for cases " << endl;

	size_t const num_basis_data = 10;
	uint16_t const dummy = 1;
	size_t const num_coeff = 5;
	double coeff[num_coeff];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);

	float rms;
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	for (BVItems item = static_cast<BVItems>(0); item < BV_kNumElems; item =
			static_cast<BVItems>(item + 1)) {
		cout << "        " << bv_param_names[item]
				<< ((item < BV_kNumElems - 1) ? ", " : "") << endl;

		size_t context_nwave = 2;
		size_t num_data = num_basis_data;
		float clip_threshold_sigma = 5.0f;
		size_t num_nwave = 3;
		size_t nwave[num_nwave];
		for (size_t i = 0; i < num_nwave; ++i) {
			nwave[i] = i;
		}
		size_t tmp_nwave;

		switch (item) {
		case BV_kDataLTNumBases:
			num_data = 2;
			break;
		case BV_kDataLTNumBasisData:
			num_data = 5;
			break;
		case BV_kDataGTNumBasisData:
			num_data = 15;
			break;
		case BV_kClipThresholdZero:
			clip_threshold_sigma = 0.0f;
			break;
		case BV_kClipThresholdNegative:
			clip_threshold_sigma = -5.0f;
			break;
		case BV_kNumNWaveZero:
			num_nwave = 0;
			break;
		case BV_kNWaveDuplicated:
			nwave[1] = nwave[0];
			break;
		case BV_kNWaveNotAscending:
			tmp_nwave = nwave[1];
			nwave[1] = nwave[0];
			nwave[0] = tmp_nwave;
			break;
		case BV_kGivenTooLargeMaximumNWave:
			nwave[num_nwave - 1] = context_nwave + 10;
			break;
		case BV_kNumBasesByNWaveGTNumBases:
			context_nwave = 1;
			break;
		default:
			assert(false);
		}

		SIMD_ALIGN
		float data[num_data];
		SetFloatSinusoidal(num_nwave, nwave, coeff, num_data, data);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(data)];
		SetBoolConstant(true, ELEMENTSOF(data), mask);
		SIMD_ALIGN
		bool final_mask[ELEMENTSOF(data)];
		SIMD_ALIGN
		float out[ELEMENTSOF(data)];
		LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
				context_nwave, num_basis_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		LIBSAKURA_SYMBOL(Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineSinusoidFloat)(context, num_nwave,
				nwave, num_data, data, mask, clip_threshold_sigma, 1, true,
				final_mask, out, &rms, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

/*
 * Test SubtractBaselineSinusoidUsingCoefficientsFloat
 * successful cases
 * subtract sinusoidal baseline with given coefficients
 * context is created with maximum nwave of 2 in all cases, that is,
 * context contains bases with nwave of 0(constant), 1(sin/cos)
 * and 2(sin/cos).
 * test cases: all combinations of nwave between 0 and 2.
 */
TEST_F(BaselineWK, SubtractBaselineSinusoidUsingCoefficientsFloatSuccessfulCase) {
	uint16_t const dummy = 1;
	size_t const num_data = NUM_DATA2;
	size_t const context_nwave = 3;
	size_t const num_nwave_set = 4;
	size_t const nwave_set_length[num_nwave_set] = { 4, 6, 4, 1 };
	size_t nwave_set1[4][1] = { { 0 }, { 1 }, { 2 }, { 3 } };
	size_t nwave_set2[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 },
			{ 1, 3 }, { 2, 3 } };
	size_t nwave_set3[4][3] = { { 0, 1, 2 }, { 0, 1, 3 }, { 0, 2, 3 },
			{ 1, 2, 3 } };
	size_t nwave_set4[1][4] = { { 0, 1, 2, 3 } };
	cout << "    Baseline context has sinusoidal bases with nwave up to 3..."
			<< endl;

	SIMD_ALIGN
	float out[num_data];
	float answer[ELEMENTSOF(out)];
	SetFloatConstant(0.0f, ELEMENTSOF(answer), answer);
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	for (size_t num_nwave = 1; num_nwave <= num_nwave_set; ++num_nwave) {
		cout << "    Testing for num_nwave = " << num_nwave << ": nwave = ";

		size_t const num_cases = nwave_set_length[num_nwave - 1];
		size_t nwave[num_nwave];
		for (size_t j = 0; j < num_cases; ++j) {
			cout << "{";
			for (size_t k = 0; k < num_nwave; ++k) {
				if (num_nwave == 1) {
					nwave[k] = nwave_set1[j][k];
				} else if (num_nwave == 2) {
					nwave[k] = nwave_set2[j][k];
				} else if (num_nwave == 3) {
					nwave[k] = nwave_set3[j][k];
				} else if (num_nwave == 4) {
					nwave[k] = nwave_set4[j][k];
				}
				cout << nwave[k];
				if (k < num_nwave - 1) {
					cout << ", ";
				}
			}
			cout << "}";
			if (j < num_cases - 1) {
				cout << ", ";
			}

			size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave,
					nwave);
			SIMD_ALIGN
			double coeff[num_coeff];
			SetDoubleConstant(1.0, num_coeff, coeff);
			SIMD_ALIGN
			float data[ELEMENTSOF(out)];
			SetFloatSinusoidal(num_nwave, nwave, coeff, num_data, data);

			LIBSAKURA_SYMBOL (Status) sub_status =
			LIBSAKURA_SYMBOL(SubtractBaselineSinusoidUsingCoefficientsFloat)(
					context, num_data, data, num_nwave, nwave, num_coeff, coeff,
					out);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);

			for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
				CheckAlmostEqual(answer[i], out[i], 1.0e-6);
			}
			if (verbose) {
				PrintArray("data  ", ELEMENTSOF(data), data);
				PrintArray("out   ", ELEMENTSOF(out), out);
				PrintArray("answer", ELEMENTSOF(answer), answer);
			}
		}
		cout << endl;
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test SubtractBaselineSinusoidUsingCoefficientsFloat
 * time-consuming successful case for performance measurement
 */
TEST_F(BaselineWK, SubtractBaselineSinusoidUsingCoefficientsPerformanceTest) {
	size_t const num_repeat = 2000;
	size_t const num_data = 50000;
	uint16_t const dummy = 1;
	size_t const context_nwave = 10;
	size_t const num_nwave = 10;
	size_t nwave[num_nwave] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave, nwave);
	SIMD_ALIGN
	double coeff[num_coeff];
	SetDoubleConstant(1.0, num_coeff, coeff);
	SIMD_ALIGN
	float data[num_data];
	SetFloatSinusoidal(num_nwave, nwave, coeff, num_data, data);
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	float answer[ELEMENTSOF(data)];
	for (size_t i = 0; i < ELEMENTSOF(coeff); ++i) {
		SetFloatConstant(0.0f, ELEMENTSOF(answer), answer);
	}
	LIBSAKURA_SYMBOL (Status) sub_status;
	double start_time = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineSinusoidUsingCoefficientsFloat)(
				context, num_data, data, num_nwave, nwave, num_coeff, coeff,
				out);
	}
	double end_time = sakura_GetCurrentTime();
	std::cout << std::setprecision(5)
			<< "#x# benchmark BaselineWK_SubtractBaselineSinusoidUsingCoefficientsPerformanceTest"
			<< " " << (end_time - start_time) << std::endl;
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);

	for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
		CheckAlmostEqual(answer[i], out[i], 5.0e-6);
	}
	if (verbose) {
		PrintArray("data  ", ELEMENTSOF(data), data);
		PrintArray("out   ", ELEMENTSOF(out), out);
		PrintArray("answer", ELEMENTSOF(answer), answer);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test SubtractBaselineSinusoidUsingCoefficients
 * erroneous cases: null pointer cases
 * parameters to be tested include context, data, nwave, coeff, and out.
 */
TEST_F(BaselineWK, SubtractBaselineSinusoidUsingCoefficientsErroneousCasesNullPointer) {
	enum NPItems {
		NP_kContext, NP_kData, NP_kNwave, NP_kCoeff, NP_kOut, NP_kNumElems
	};
	vector<string> np_param_names =
			{ "context", "data", "nwave", "coeff", "out" };
	cout << "    Testing for ";

	size_t const num_data = 10;
	uint16_t const dummy = 1;
	size_t const context_nwave = 2;
	size_t const num_nwave = 3;
	size_t nwave[num_nwave] = { 0, 1, 2 };
	size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave, nwave);
	SIMD_ALIGN
	double coeff[num_coeff];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	SIMD_ALIGN
	float data[num_data];
	SetFloatSinusoidal(num_nwave, nwave, coeff, num_data, data);
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		LIBSAKURA_SYMBOL(BaselineContext) *context_ptr = context;
		float *data_ptr = data;
		size_t *nwave_ptr = nwave;
		double *coeff_ptr = coeff;
		float *out_ptr = out;

		switch (item) {
		case NP_kContext:
			context_ptr = nullptr;
			break;
		case NP_kData:
			data_ptr = nullptr;
			break;
		case NP_kNwave:
			nwave_ptr = nullptr;
			break;
		case NP_kCoeff:
			coeff_ptr = nullptr;
			break;
		case NP_kOut:
			out_ptr = nullptr;
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL (Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineSinusoidUsingCoefficientsFloat)(
				context_ptr, num_data, data_ptr, num_nwave, nwave_ptr,
				num_coeff, coeff_ptr, out_ptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test SubtractBaselineSinusoidUsingCoefficients
 * erroneous cases: unaligned cases
 * parameters to be tested include data, coeff and out
 */
TEST_F(BaselineWK, SubtractBaselineSinusoidUsingCoefficientsErroneousCasesUnaligned) {
	enum UAItems {
		UA_kData, UA_kCoeff, UA_kOut, UA_kNumElems
	};
	vector<string> ua_param_names = { "data", "coeff", "out" };
	cout << "    Testing for ";

	size_t const num_data = 10;
	uint16_t const dummy = 1;
	size_t const context_nwave = 2;
	size_t const num_nwave = 3;
	size_t nwave[num_nwave] = { 0, 1, 2 };
	size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave, nwave);
	SIMD_ALIGN
	double coeff[num_coeff + 1];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	SIMD_ALIGN
	float data[num_data + 1];
	SetFloatSinusoidal(num_nwave, nwave, coeff, num_data, data);
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
			context_nwave, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	for (UAItems item = static_cast<UAItems>(0); item < UA_kNumElems; item =
			static_cast<UAItems>(item + 1)) {
		cout << ua_param_names[item] << ((item < UA_kNumElems - 1) ? ", " : "");

		float *data_ptr = data;
		double *coeff_ptr = coeff;
		float *out_ptr = out;

		switch (item) {
		case UA_kData:
			++data_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(data_ptr));
			break;
		case UA_kCoeff:
			++coeff_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(coeff_ptr));
			break;
		case UA_kOut:
			++out_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(out_ptr));
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL (Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineSinusoidUsingCoefficientsFloat)(
				context, num_data, data_ptr, num_nwave, nwave, num_coeff,
				coeff_ptr, out_ptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyBaselineContext(
			context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test SubtractBaselineSinusoidUsingCoefficientsFloat
 * erroneous cases: bad parameter value cases as follows:
 *     (1) num_data < context->num_bases
 *     (2) num_data < (!=) context->num_basis_data
 *     (3) num_data > (!=) context->num_basis_data
 *     (4) num_nwave == 0
 *     (5) nwave with duplicated elements
 *     (6) nwave with not in ascending order
 *     (7) max_nwave in context < given max_nwave
 *     (8) num_coeff > context->num_bases
 */
TEST_F(BaselineWK, SubtractBaselineSinusoidUsingCoefficientsErroneousCasesBadParameterValue) {
	enum BVItems {
		BV_kDataLTNumBases,
		BV_kDataLTNumBasisData,
		BV_kDataGTNumBasisData,
		BV_kNumNWaveZero,
		BV_kNWaveDuplicated,
		BV_kNWaveNotAscending,
		BV_kGivenTooLargeMaximumNWave,
		BV_kNumCoeffGTNumBases,
		BV_kNumElems
	};
	vector<string> bv_param_names = { "(num_data < context->num_bases)",
			"(num_data < context->num_basis_data)",
			"(num_data > context->num_basis_data)",
			"(clip_threshold_sigma == 0)", "(clip_threshold_sigma < 0)",
			"(num_nwave == 0)", "(nwave duplicated)", "(nwave not ascending)",
			"(context_max_nwave < given_max_nwave)",
			"(num_coeff > context->num_bases)" };
	cout << "    Testing for cases " << endl;

	size_t const num_basis_data = 10;
	uint16_t const dummy = 1;

	for (BVItems item = static_cast<BVItems>(0); item < BV_kNumElems; item =
			static_cast<BVItems>(item + 1)) {
		cout << "        " << bv_param_names[item]
				<< ((item < BV_kNumElems - 1) ? ", " : "") << endl;
		size_t num_data = num_basis_data;
		size_t context_nwave = 2;
		size_t num_nwave = 3;
		size_t nwave[num_nwave];
		for (size_t i = 0; i < num_nwave; ++i) {
			nwave[i] = i;
		}
		size_t num_coeff = 5;
		size_t num_coeff_by_nwave = GetNumberOfSinusoidalCoefficients(num_nwave,
				nwave);
		size_t tmp_nwave;

		switch (item) {
		case BV_kDataLTNumBases:
			num_data = 2;
			break;
		case BV_kDataLTNumBasisData:
			num_data = 5;
			break;
		case BV_kDataGTNumBasisData:
			num_data = 15;
			break;
		case BV_kNumNWaveZero:
			num_nwave = 0;
			break;
		case BV_kNWaveDuplicated:
			nwave[1] = nwave[0];
			break;
		case BV_kNWaveNotAscending:
			tmp_nwave = nwave[1];
			nwave[1] = nwave[0];
			nwave[0] = tmp_nwave;
			break;
		case BV_kGivenTooLargeMaximumNWave:
			nwave[num_nwave - 1] = context_nwave + 10;
			break;
		case BV_kNumCoeffGTNumBases:
			num_coeff = num_coeff_by_nwave + 1;
			break;
		default:
			assert(false);
		}

		SIMD_ALIGN
		double coeff[num_coeff];
		SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
		SIMD_ALIGN
		float data[num_data];
		SetFloatSinusoidal(num_nwave, nwave, coeff, num_data, data);
		SIMD_ALIGN
		float out[ELEMENTSOF(data)];
		LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kSinusoid), dummy, dummy,
				context_nwave, num_basis_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		LIBSAKURA_SYMBOL(Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineSinusoidUsingCoefficientsFloat)(
				context, num_data, data, num_nwave, nwave, num_coeff, coeff,
				out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}
