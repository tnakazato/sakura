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
 * Test cases to be implemented.
 * functions to be tested include:
 * - sakura_SubtractBaselineCubicSplineFloat
 * - sakura_GetBestFitBaselineCoefficientsCubicSplineFloat
 * - sakura_SubtractBaselineSinusoidFloat
 * - sakura_GetBestFitBaselineCoefficientsSinusoidFloat
 * - sakura_SubtractBaselineSinusoidUsingCoefficientsFloat
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
class BaselineWK: public ::testing::Test {
protected:

	BaselineWK() :
			verbose(false) {
	}

	virtual void SetUp() {
		// Initialize sakura
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Initialize)(nullptr,
				nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}

	virtual void TearDown() {
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	//Set (coeff[0]+coeff[1]*x+coeff[2]*x*x+...) float values into an array
	void SetFloatPolynomial(size_t num_coeff, double const *coeff,
			size_t num_data, float *data) {
		for (size_t i = 0; i < num_data; ++i) {
			double val = 0.0;
			double x = (double) i;
			for (size_t j = 0; j < num_coeff; ++j) {
				val *= x;
				val += coeff[num_coeff - 1 - j];
			}
			data[i] = static_cast<float>(val);
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
 * successful cases with masked data
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplineSuccessfulCaseWithMaskedData) {
	enum MAItems {
		MA_kMiddle, MA_kLeft, MA_kRight, MA_kLeftRight, MA_kNumElems
	};
	vector<string> ma_param_names = { "middle", "left", "right", "leftright" };
	cout << "    Testing for ";

	size_t num_pieces = 3;
	SIMD_ALIGN
	size_t boundary[num_pieces + 1];
	SIMD_ALIGN
	double answer[4 * num_pieces];
	SetDoubleConstant(1.0, ELEMENTSOF(answer), answer);
	SIMD_ALIGN
	double out[num_pieces][4];
	float rms;
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL (Status) coeff_status;

	size_t const num_data = 20;
	SIMD_ALIGN
	float data[num_data];
	SetFloatPolynomial(4, answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	if (verbose) {
		PrintArray("data", num_data, data);
	}
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextCubicSplineFloat(num_pieces, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	for (MAItems item = static_cast<MAItems>(0); item < MA_kNumElems; item =
			static_cast<MAItems>(item + 1)) {
		cout << ma_param_names[item] << ((item < MA_kNumElems - 1) ? ", " : "");
		size_t num_masked_idx = 2;//10;

		num_masked_idx = (item == MA_kLeftRight) ? 4 : 2;//20 : 10;
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
		coeff_status = LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(context,
				num_pieces, num_data, data, mask, 5.0f, 1, out, nullptr,
				nullptr, mask, &rms, boundary, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);
		for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
			CheckAlmostEqual(answer[i], out[i / 4][i % 4], 1.0e-5);
		}
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test GetBestFitBaselineCoefficientsCubicSpline
 * successful case with num_fitting_max == 0
 * coeff values should not be modified.
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplineZeroNumClippingMax) {
	size_t const num_pieces = 2;
	double coeff_orig[num_pieces][4] = { { 1.1, 2.2, 3.3, 4.4 }, { 5.5, 6.6,
			7.7, 8.8 } };
	SIMD_ALIGN
	double answer[4 * num_pieces];
	for (size_t i = 0; i < ELEMENTSOF(answer); i += 4) {
		answer[i] = 1.0;
		answer[i + 1] = 1e-4;
		answer[i + 2] = 1e-8;
		answer[i + 3] = 1e-12;
	}
	size_t const num_data = NUM_DATA2;
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(ELEMENTSOF(answer), answer, num_data, in_data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(in_data)];
	SetBoolConstant(true, ELEMENTSOF(in_data), mask);
	if (verbose) {
		PrintArray("in_data", num_data, in_data);
	}
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextCubicSplineFloat(num_pieces, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	double out[num_pieces][4];
	for (size_t i = 0; i < num_pieces; ++i) {
		for (size_t j = 0; j < 4; ++j) {
			out[i][j] = coeff_orig[i][j];
		}
	}
	float rms;
	SIMD_ALIGN
	size_t boundary[num_pieces + 1];
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) coeff_status;
	coeff_status = LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(context, num_pieces,
			num_data, in_data, mask, 5.0f, 0, out, nullptr, nullptr, mask,
			&rms, boundary, &baseline_status);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
	EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);
	for (size_t i = 0; i < num_pieces; ++i) {
		for (size_t j = 0; j < 4; ++j) {
			EXPECT_EQ(coeff_orig[i][j], out[i][j]);
		}
	}
	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test GetBestFitBaselineCoefficientsCubicSpline
 * time-consuming successful case for performance measurement
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplinePerformanceTest) {
	size_t const num_repeat = 300;
	size_t const num_pieces = 2;
	SIMD_ALIGN
	double answer[4 * num_pieces];
	for (size_t i = 0; i < ELEMENTSOF(answer); i += 4) {
		answer[i] = 1.0;
		answer[i + 1] = 0.0;
		answer[i + 2] = 0.0;
		answer[i + 3] = 0.0;
	}
	size_t const num_data = 70000;
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(4, answer, num_data, in_data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(in_data)];
	SetBoolConstant(true, ELEMENTSOF(in_data), mask);
	if (verbose) {
		PrintArray("in_data", num_data, in_data);
	}
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextCubicSplineFloat(num_pieces, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	double out[num_pieces][4];
	float rms;
	SIMD_ALIGN
	size_t boundary[num_pieces + 1];
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) coeff_status;
	double start_time = GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		coeff_status = LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(context,
				num_pieces, num_data, in_data, mask, 5.0f, 1, out, nullptr,
				nullptr, mask, &rms, boundary, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
	}
	double end_time = GetCurrentTime();
	std::cout << std::setprecision(5)
			<< "#x# benchmark BaselineWK_LSQFitCubicSplinePerformanceTest"
			<< " " << (end_time - start_time) << std::endl;

	EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
		NP_kLSQFitStatus,
		NP_kNumElems
	};
	vector<string> np_param_names = { "context", "data", "mask", "coeff",
			"final_mask", "boundary", "baseline_status" };
	cout << "    Testing for ";

	size_t const num_pieces = 2;
	size_t const num_data = 10;
	double coeff[4];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	SIMD_ALIGN
	float data[num_data];
	SetFloatPolynomial(ELEMENTSOF(coeff), coeff, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	double out[1][4];
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextCubicSplineFloat(num_pieces, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	float rms;
	SIMD_ALIGN
	size_t boundary[num_pieces + 1];
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		float *data_ptr = data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		double (*out_ptr)[4] = out;
		size_t *boundary_ptr = boundary;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_ptr = context;
		LIBSAKURA_SYMBOL(LSQFitStatus) *baseline_status_ptr = &baseline_status;
		LIBSAKURA_SYMBOL(Status) coeff_status;

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
		case NP_kLSQFitStatus:
			baseline_status_ptr = nullptr;
			break;
		default:
			assert(false);
		}

		coeff_status = LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(context_ptr,
				num_pieces, num_data, data_ptr, mask_ptr, 5.0f, 1, out_ptr,
				nullptr, nullptr, final_mask_ptr, &rms, boundary_ptr,
				baseline_status_ptr);
		LIBSAKURA_SYMBOL(Status) coeff_status_ref = LIBSAKURA_SYMBOL(
				Status_kInvalidArgument);
		if (item == NP_kCoeff) {
			coeff_status_ref = LIBSAKURA_SYMBOL(Status_kOK);
		}
		EXPECT_EQ(coeff_status_ref, coeff_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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

	size_t const num_pieces = 2;
	size_t const num_data = 10;
	double coeff[4];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	SIMD_ALIGN
	float data[num_data + 1];
	SetFloatPolynomial(ELEMENTSOF(coeff), coeff, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	double out[num_pieces + 1][4];
	float rms;
	SIMD_ALIGN
	size_t boundary[num_pieces + 1];
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextCubicSplineFloat(num_pieces, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;

	for (UAItems item = static_cast<UAItems>(0); item < UA_kNumElems; item =
			static_cast<UAItems>(item + 1)) {
		cout << ua_param_names[item] << ((item < UA_kNumElems - 1) ? ", " : "");

		float *data_ptr = data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		double *out_unaligned = reinterpret_cast<double *>(out) + 1;
		double (*out_ptr)[4] = out;
		size_t *boundary_ptr = boundary;

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
			out_ptr = reinterpret_cast<double (*)[4]>(out_unaligned);
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

		LIBSAKURA_SYMBOL (Status) coeff_status;
		coeff_status = LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(context,
				num_pieces, num_data, data_ptr, mask_ptr, 5.0f, 1, out_ptr,
				nullptr, nullptr, final_mask_ptr, &rms, boundary_ptr,
				&baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test GetBestFitBaselineCoefficientsCubicSpline
 * erroneous cases: bad parameter value cases as follows:
 *     (1) num_data < context->num_bases
 *     (2) num_data < (!=) context->num_basis_data
 *     (3) num_data > (!=) context->num_basis_data
 *     (4) num_pieces == 0
 *     (5) clip_threshold_sigma == 0
 *     (6) clip_threshold_sigma < 0
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplineErroneousCasesBadParameterValue) {
	enum BVItems {
		BV_kDataLTNumBases,
		BV_kDataLTNumBasisData,
		BV_kDataGTNumBasisData,
		BV_kNumPiecesZero,
		BV_kClipThresholdZero,
		BV_kClipThresholdNegative,
		BV_kNumElems
	};
	vector<string> bv_param_names = { "(num_data < context->num_bases)",
			"(num_data < context->num_basis_data)",
			"(num_data > context->num_basis_data)", "(num_pieces == 0)",
			"(clip_threshold_sigma == 0.0)", "(clip_threshold_sigma < 0.0)" };
	cout << "    Testing for cases " << endl;

	size_t const num_basis_data = 10;
	size_t const num_pieces_orig = 1;
	SIMD_ALIGN
	size_t boundary[num_pieces_orig + 1];
	double coeff[4 * num_pieces_orig];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	SIMD_ALIGN
	double out[num_pieces_orig][4];
	float rms;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextCubicSplineFloat(num_pieces_orig,
					num_basis_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) coeff_status;

	for (BVItems item = static_cast<BVItems>(0); item < BV_kNumElems; item =
			static_cast<BVItems>(item + 1)) {
		cout << "        " << bv_param_names[item]
				<< ((item < BV_kNumElems - 1) ? ", " : "") << endl;
		size_t num_data = num_basis_data;
		size_t num_pieces = num_pieces_orig;
		float clip_threshold = 5.0f;
		uint16_t num_fitting_max = 1;
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
		case BV_kNumPiecesZero:
			num_pieces = 0;
			break;
		case BV_kClipThresholdZero:
			clip_threshold = 0.0f;
			break;
		case BV_kClipThresholdNegative:
			clip_threshold = -5.0f;
			break;
		default:
			assert(false);
		}

		SIMD_ALIGN
		float data[num_data];
		SetFloatPolynomial(ELEMENTSOF(coeff), coeff, num_data, data);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(data)];
		SetBoolConstant(true, ELEMENTSOF(data), mask);
		SIMD_ALIGN
		bool final_mask[ELEMENTSOF(data)];

		coeff_status = LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(context,
				num_pieces, num_data, data, mask, clip_threshold,
				num_fitting_max, out, nullptr, nullptr, final_mask, &rms,
				boundary, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test SubtractBaselineCubicSplineFloat
 * time-consuming successful case for performance measurement
 */
TEST_F(BaselineWK, SubtractBaselineCubicSplinePerformanceTest) {
	size_t const num_repeat = 300;
	size_t const num_pieces = 2;
	size_t const num_data = 70000;
	SIMD_ALIGN
	float data[num_data];
	double coeff[4] = { 1.0, 1e-4, 1e-8, 1e-12 };
	SetFloatPolynomial(ELEMENTSOF(coeff), coeff, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	if (verbose) {
		PrintArray("data", num_data, data);
	}
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextCubicSplineFloat(num_pieces, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	SIMD_ALIGN
	double answer[ELEMENTSOF(data)];
	SetDoubleConstant(0.0, ELEMENTSOF(data), answer);
	float rms;
	SIMD_ALIGN
	size_t boundary[num_pieces + 1];
	LIBSAKURA_SYMBOL(Status) sub_status;

	double start_time = GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		sub_status = LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(context,
				num_pieces, num_data, data, mask, 5.0f, 1, nullptr, nullptr,
				out, mask, &rms, boundary, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);
	}
	double end_time = GetCurrentTime();
	std::cout << std::setprecision(5)
			<< "#x# benchmark BaselineWK_LSQFitCubicSplinePerformanceTest"
			<< " " << (end_time - start_time) << std::endl;
	EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);
	if (verbose) {
		PrintArray("data  ", ELEMENTSOF(data), data);
		PrintArray("out   ", ELEMENTSOF(out), out);
		PrintArray("answer", ELEMENTSOF(answer), answer);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test SubtractBaselineCubicSplineFloat
 * erroneous cases: null pointer cases
 * parameters to be tested include context, data, mask, final_mask, out and baseline_status.
 */
TEST_F(BaselineWK, SubtractBaselineCubicSplineErroneousCasesNullPointer) {
	enum NPItems {
		NP_kContext,
		NP_kData,
		NP_kMask,
		NP_kFinalMask,
		NP_kOut,
		NP_kBoundary,
		NP_kLSQFitStatus,
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
	SetFloatPolynomial(ELEMENTSOF(coeff), coeff, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextCubicSplineFloat(num_pieces, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	float rms;
	SIMD_ALIGN
	size_t boundary[num_pieces + 1];
	LIBSAKURA_SYMBOL(Status) sub_status;

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_ptr = context;
		float *data_ptr = data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		float *out_ptr = out;
		size_t *boundary_ptr = boundary;
		LIBSAKURA_SYMBOL(LSQFitStatus) *baseline_status_ptr = &baseline_status;

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
		case NP_kLSQFitStatus:
			baseline_status_ptr = nullptr;
			break;
		default:
			assert(false);
		}

		sub_status =
		LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(context_ptr, num_pieces,
				num_data, data_ptr, mask_ptr, 5.0f, 1, nullptr, nullptr,
				out_ptr, final_mask_ptr, &rms, boundary_ptr,
				baseline_status_ptr);
		LIBSAKURA_SYMBOL(Status) sub_status_ref = LIBSAKURA_SYMBOL(
				Status_kInvalidArgument);
		if (item == NP_kOut) {
			sub_status_ref = LIBSAKURA_SYMBOL(Status_kOK);
		}
		EXPECT_EQ(sub_status_ref, sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test SubtractBaselineCubicSplineFloat
 * erroneous cases: unaligned cases
 * parameters to be tested include data, mask, final_mask and out
 */
TEST_F(BaselineWK, SubtractBaselineCubicSplineErroneousCasesUnaligned) {
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
	SetFloatPolynomial(ELEMENTSOF(coeff), coeff, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextCubicSplineFloat(num_pieces, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	float rms;
	SIMD_ALIGN
	size_t boundary[num_pieces + 1];
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) sub_status;

	for (UAItems item = static_cast<UAItems>(0); item < UA_kNumElems; item =
			static_cast<UAItems>(item + 1)) {
		cout << ua_param_names[item] << ((item < UA_kNumElems - 1) ? ", " : "");

		float *data_ptr = data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		float *out_ptr = out;
		size_t *boundary_ptr = boundary;

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

		sub_status = LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(context,
				num_pieces, num_data, data_ptr, mask_ptr, 5.0f, 1, nullptr,
				nullptr, out_ptr, final_mask_ptr, &rms, boundary_ptr,
				&baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test SubtractBaselineCubicSplineFloat
 * erroneous cases: bad parameter value cases as follows:
 *     (1) num_data < context->num_bases
 *     (2) num_data < (!=) context->num_basis_data
 *     (3) num_data > (!=) context->num_basis_data
 *     (4) num_pieces == 0
 *     (5) clip_threshold_sigma == 0
 *     (6) clip_threshold_sigma < 0
 */
TEST_F(BaselineWK, SubtractBaselineCubicSplineErroneousCasesBadParameterValue) {
	enum BVItems {
		BV_kDataLTNumBases,
		BV_kDataLTNumBasisData,
		BV_kDataGTNumBasisData,
		BV_kNumPiecesZero,
		BV_kClipThresholdZero,
		BV_kClipThresholdNegative,
		BV_kNumElems
	};
	vector<string> bv_param_names = { "(num_data < context->num_bases)",
			"(num_data < context->num_basis_data)",
			"(num_data > context->num_basis_data)", "(num_pieces == 0)",
			"(clip_threshold_sigma == 0.0)", "(clip_threshold_sigma < 0.0)" };
	cout << "    Testing for cases " << endl;

	size_t const num_basis_data = 10;
	double coeff[4];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
	float rms;
	size_t const num_pieces_context = 1;
	SIMD_ALIGN
	size_t boundary[num_pieces_context + 1];
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextCubicSplineFloat(num_pieces_context,
					num_basis_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) sub_status;

	for (BVItems item = static_cast<BVItems>(0); item < BV_kNumElems; item =
			static_cast<BVItems>(item + 1)) {
		cout << "        " << bv_param_names[item]
				<< ((item < BV_kNumElems - 1) ? ", " : "") << endl;
		size_t num_data = num_basis_data;
		size_t num_pieces = num_pieces_context;
		float clip_threshold = 5.0f;
		uint16_t num_fitting_max = 1;
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
		case BV_kNumPiecesZero:
			num_pieces = 0;
			break;
		case BV_kClipThresholdZero:
			clip_threshold = 0.0f;
			break;
		case BV_kClipThresholdNegative:
			clip_threshold = -5.0f;
			break;
		default:
			assert(false);
		}

		SIMD_ALIGN
		float data[num_data];
		SetFloatPolynomial(ELEMENTSOF(coeff), coeff, num_data, data);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(data)];
		SetBoolConstant(true, ELEMENTSOF(data), mask);
		SIMD_ALIGN
		bool final_mask[ELEMENTSOF(data)];
		SIMD_ALIGN
		float out[ELEMENTSOF(data)];

		sub_status = LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(context,
				num_pieces, num_data, data, mask, clip_threshold,
				num_fitting_max, nullptr, nullptr, out, final_mask, &rms,
				boundary, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) coeff_status;

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

			coeff_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context,
					num_nwave, nwave, num_data, data, mask, 5.0f, 1, num_coeff,
					out, nullptr, nullptr, mask, &rms, &baseline_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);
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

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	double out[ELEMENTSOF(answer)];
	float rms;
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) coeff_status;

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

		coeff_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context, num_nwave,
				nwave, num_data, data, mask, 5.0f, 1, num_coeff, out, nullptr,
				nullptr, mask, &rms, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);
		if (verbose) {
			PrintArray("data  ", ELEMENTSOF(data), data);
			PrintArray("out   ", ELEMENTSOF(out), out);
			PrintArray("answer", ELEMENTSOF(answer), answer);
		}
		for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
			CheckAlmostEqual(answer[i], out[i], 1.0e-5);
		}
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test GetBestFitBaselineCoefficientsSinusoid
 * successful case with num_fitting_max == 0
 * coeff values should not be modified.
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsSinusoidWithZeroNumFittingMax) {
	size_t const num_data = NUM_DATA2;
	size_t const context_nwave = 2;
	size_t const num_nwave = 2;
	size_t nwave[num_nwave] = { 0, 1 };
	size_t num_coeff = GetNumberOfSinusoidalCoefficients(num_nwave, nwave);
	double answer[num_coeff];
	SetDoubleConstant(1.0, num_coeff, answer);
	SIMD_ALIGN
	float data[num_data];
	SetFloatSinusoidal(num_nwave, nwave, answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	double out_orig[3] = { 1.1, 2.2, 3.3 };
	SIMD_ALIGN
	double out[ELEMENTSOF(answer)];
	for (size_t i = 0; i < ELEMENTSOF(out); ++i) {
		out[i] = out_orig[i];
	}
	float rms;
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL (Status) coeff_status;
	coeff_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context, num_nwave,
			nwave, num_data, data, mask, 5.0f, 0, num_coeff, out, nullptr,
			nullptr, mask, &rms, &baseline_status);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
	EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);
	for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
		CheckAlmostEqual(out_orig[i], out[i], 1.0e-6);
	}
	if (verbose) {
		PrintArray("data  ", ELEMENTSOF(data), data);
		PrintArray("out   ", ELEMENTSOF(out), out);
		PrintArray("answer", ELEMENTSOF(answer), answer);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test GetBestFitBaselineCoefficientsSinusoid
 * time-consuming successful case for performance measurement
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsSinusoidPerformanceTest) {
	size_t const num_repeat = 100;
	size_t const num_data = 50000;
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
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	double out[ELEMENTSOF(answer)];
	float rms;
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL (Status) coeff_status;

	double start_time = GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		coeff_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context, num_nwave,
				nwave, num_data, data, mask, 5.0f, 1, num_coeff, out, nullptr,
				nullptr, mask, &rms, &baseline_status);
	}
	double end_time = GetCurrentTime();
	std::cout << std::setprecision(5)
			<< "#x# benchmark BaselineWK_LSQFitSinusoidPerformanceTest" << " "
			<< (end_time - start_time) << std::endl;
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
	EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);
	for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
		CheckAlmostEqual(answer[i], out[i], 1.0e-6);
	}
	if (verbose) {
		PrintArray("data  ", ELEMENTSOF(data), data);
		PrintArray("out   ", ELEMENTSOF(out), out);
		PrintArray("answer", ELEMENTSOF(answer), answer);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
		NP_kLSQFitStatus,
		NP_kNumElems
	};
	vector<string> np_param_names = { "context", "data", "mask", "nwave",
			"coeff", "final_mask", "baseline_status" };
	cout << "    Testing for ";

	size_t const num_data = 10;
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
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) coeff_status;

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_ptr = context;
		float *data_ptr = data;
		bool *mask_ptr = mask;
		size_t *nwave_ptr = nwave;
		bool *final_mask_ptr = final_mask;
		double *coeff_ptr = coeff;
		LIBSAKURA_SYMBOL(LSQFitStatus) *baseline_status_ptr = &baseline_status;

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
		case NP_kLSQFitStatus:
			baseline_status_ptr = nullptr;
			break;
		default:
			assert(false);
		}

		coeff_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context_ptr,
				num_nwave, nwave_ptr, num_data, data_ptr, mask_ptr, 5.0f, 1,
				num_coeff, coeff_ptr, nullptr, nullptr, final_mask_ptr, &rms,
				baseline_status_ptr);
		LIBSAKURA_SYMBOL(Status) status_ref = LIBSAKURA_SYMBOL(
				Status_kInvalidArgument);
		if (item == NP_kCoeff) {
			status_ref = LIBSAKURA_SYMBOL(Status_kOK);
		}
		EXPECT_EQ(status_ref, coeff_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL (Status) coeff_status;

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

		coeff_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context, num_nwave,
				nwave, num_data, data_ptr, mask_ptr, 5.0f, 1, num_coeff,
				coeff_ptr, nullptr, nullptr, final_mask_ptr, &rms,
				&baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
	size_t const context_nwave = 2;
	float rms;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave,
					num_basis_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) coeff_status;

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

		coeff_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context, num_nwave,
				nwave, num_data, data, mask, clip_threshold_sigma, 1, num_coeff,
				out, nullptr, nullptr, final_mask, &rms, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) sub_status;

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

			sub_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context,
					num_nwave, nwave, num_data, data, mask, 5.0f, 1, num_coeff,
					nullptr, nullptr, out, mask, &rms, &baseline_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);
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

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	float answer[ELEMENTSOF(data)];
	for (size_t i = 0; i < ELEMENTSOF(coeff); ++i) {
		SetFloatConstant(0.0f, ELEMENTSOF(answer), answer);
	}
	float rms;
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) sub_status;

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

		sub_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context, num_nwave,
				nwave, num_data, data, mask, 5.0f, 1, num_coeff, nullptr,
				nullptr, out, mask, &rms, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);
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

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	float answer[ELEMENTSOF(data)];
	for (size_t i = 0; i < ELEMENTSOF(coeff); ++i) {
		SetFloatConstant(0.0f, ELEMENTSOF(answer), answer);
	}
	float rms;
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL (Status) sub_status;

	double start_time = GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		sub_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context, num_nwave,
				nwave, num_data, data, mask, 5.0f, 1, num_coeff, nullptr,
				nullptr, out, mask, &rms, &baseline_status);
	}
	double end_time = GetCurrentTime();
	std::cout << std::setprecision(5)
			<< "#x# benchmark BaselineWK_LSQFitSinusoidPerformanceTest" << " "
			<< (end_time - start_time) << std::endl;
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);
	EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);
	for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
		CheckAlmostEqual(answer[i], out[i], 5.0e-6);
	}
	if (verbose) {
		PrintArray("data  ", ELEMENTSOF(data), data);
		PrintArray("out   ", ELEMENTSOF(out), out);
		PrintArray("answer", ELEMENTSOF(answer), answer);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
		NP_kLSQFitStatus,
		NP_kNumElems
	};
	vector<string> np_param_names = { "context", "nwave", "data", "mask",
			"final_mask", "out", "baseline_status" };
	cout << "    Testing for ";

	size_t const num_data = 10;
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
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) sub_status;

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_ptr = context;
		float *data_ptr = data;
		bool *mask_ptr = mask;
		size_t *nwave_ptr = nwave;
		bool *final_mask_ptr = final_mask;
		float *out_ptr = out;
		LIBSAKURA_SYMBOL(LSQFitStatus) *baseline_status_ptr = &baseline_status;

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
		case NP_kLSQFitStatus:
			baseline_status_ptr = nullptr;
			break;
		default:
			assert(false);
		}

		sub_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context_ptr,
				num_nwave, nwave_ptr, num_data, data_ptr, mask_ptr, 5.0f, 1,
				num_coeff, nullptr, nullptr, out_ptr, final_mask_ptr, &rms,
				baseline_status_ptr);
		LIBSAKURA_SYMBOL(Status) sub_status_ref = LIBSAKURA_SYMBOL(
				Status_kInvalidArgument);
		if (item == NP_kOut) {
			sub_status_ref = LIBSAKURA_SYMBOL(Status_kOK);
		}
		EXPECT_EQ(sub_status_ref, sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) sub_status;

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

		sub_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context, num_nwave,
				nwave, num_data, data_ptr, mask_ptr, 5.0f, 1, num_coeff,
				nullptr, nullptr, out_ptr, final_mask_ptr, &rms,
				&baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
	size_t const num_coeff = 5;
	double coeff[num_coeff];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);

	float rms;
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;
	LIBSAKURA_SYMBOL(Status) sub_status;

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
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextSinusoidFloat(context_nwave,
						num_basis_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		sub_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context, num_nwave,
				nwave, num_data, data, mask, clip_threshold_sigma, 1, num_coeff,
				nullptr, nullptr, out, final_mask, &rms, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyLSQFitContextFloat(context);
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
	LIBSAKURA_SYMBOL(Status) sub_status;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL(Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
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

			sub_status = LIBSAKURA_SYMBOL(SubtractSinusoidFloat)(context,
					num_data, data, num_nwave, nwave, num_coeff, coeff, out);
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

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
}

/*
 * Test SubtractBaselineSinusoidUsingCoefficientsFloat
 * time-consuming successful case for performance measurement
 */
TEST_F(BaselineWK, SubtractBaselineSinusoidUsingCoefficientsPerformanceTest) {
	size_t const num_repeat = 2000;
	size_t const num_data = 50000;
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
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	float out[ELEMENTSOF(data)];
	float answer[ELEMENTSOF(data)];
	for (size_t i = 0; i < ELEMENTSOF(coeff); ++i) {
		SetFloatConstant(0.0f, ELEMENTSOF(answer), answer);
	}
	LIBSAKURA_SYMBOL (Status) sub_status;

	double start_time = GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		sub_status = LIBSAKURA_SYMBOL(SubtractSinusoidFloat)(context, num_data,
				data, num_nwave, nwave, num_coeff, coeff, out);
	}
	double end_time = GetCurrentTime();
	std::cout << std::setprecision(5)
			<< "#x# benchmark BaselineWK_SubtractSinusoidPerformanceTest" << " "
			<< (end_time - start_time) << std::endl;
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);
	for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
		CheckAlmostEqual(answer[i], out[i], 5.0e-6);
	}
	if (verbose) {
		PrintArray("data  ", ELEMENTSOF(data), data);
		PrintArray("out   ", ELEMENTSOF(out), out);
		PrintArray("answer", ELEMENTSOF(answer), answer);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
	LIBSAKURA_SYMBOL (Status) sub_status;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_ptr = context;
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

		sub_status = LIBSAKURA_SYMBOL(SubtractSinusoidFloat)(context_ptr,
				num_data, data_ptr, num_nwave, nwave_ptr, num_coeff, coeff_ptr,
				out_ptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
	LIBSAKURA_SYMBOL (Status) sub_status;
	LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(context_nwave, num_data,
					&context);
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

		sub_status = LIBSAKURA_SYMBOL(SubtractSinusoidFloat)(context, num_data,
				data_ptr, num_nwave, nwave, num_coeff, coeff_ptr, out_ptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
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
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextSinusoidFloat(context_nwave,
						num_basis_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		LIBSAKURA_SYMBOL(Status) sub_status;
		sub_status = LIBSAKURA_SYMBOL(SubtractSinusoidFloat)(context, num_data,
				data, num_nwave, nwave, num_coeff, coeff, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyLSQFitContextFloat(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

TEST_F(BaselineWK, CreateLSQFitContextFloat) {
	LSQFitTypeInternal bltypes[] = { LSQFitTypeInternal_kPolynomial,
			LSQFitTypeInternal_kChebyshev, LSQFitTypeInternal_kCubicSpline,
			LSQFitTypeInternal_kSinusoid };
	uint16_t order_valid_min = 0;
	size_t num_data_valid_min = 0;
	for (size_t itype = 0; itype < ELEMENTSOF(bltypes); ++itype) {
		order_valid_min =
				(bltypes[itype] == LSQFitTypeInternal_kCubicSpline) ? 1 : 0;
		for (uint16_t order = 0; order < order_valid_min + 3; ++order) {
			if (bltypes[itype] == LSQFitTypeInternal_kPolynomial) {
				num_data_valid_min = order + 1;
			} else if (bltypes[itype] == LSQFitTypeInternal_kChebyshev) {
				num_data_valid_min = order + 1;
			} else if (bltypes[itype] == LSQFitTypeInternal_kCubicSpline) {
				num_data_valid_min = order + 3;
			} else if (bltypes[itype] == LSQFitTypeInternal_kSinusoid) {
				num_data_valid_min = 2 * order + 2;
			}
			for (size_t num_data = 0; num_data < num_data_valid_min + 3;
					++num_data) {
				if (verbose) {
					cout << "   testing [bltype=" << itype << ", order="
							<< order << ", num_data=" << num_data << "] ... ";
				}
				LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
				LIBSAKURA_SYMBOL (Status) status;
				if (bltypes[itype] == LSQFitTypeInternal_kSinusoid) {
					status = sakura_CreateLSQFitContextSinusoidFloat(order,
							num_data, &context);
				} else if (bltypes[itype]
						== LSQFitTypeInternal_kCubicSpline) {
					status = sakura_CreateLSQFitContextCubicSplineFloat(order,
							num_data, &context);
				} else {
					LIBSAKURA_SYMBOL(LSQFitType) type_ext = LIBSAKURA_SYMBOL(
							LSQFitType_kPolynomial);
					if (bltypes[itype] == LSQFitTypeInternal_kChebyshev) {
						type_ext = LIBSAKURA_SYMBOL(LSQFitType_kChebyshev);
					}
					status = sakura_CreateLSQFitContextPolynomialFloat(type_ext, order,
							num_data, &context);
				}
				if ((num_data < num_data_valid_min)
						|| (order < order_valid_min)) {
					EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument),
							status);
				} else {
					EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
				}
				if (verbose) {
					if (status == LIBSAKURA_SYMBOL(Status_kInvalidArgument)) {
						cout << "(InvalidArgument)";
					}
					cout << "OK" << endl;
				}
			}
		}
	}
}

/*
 * Test LSQFitPolynomial
 * successful case
 */
TEST_F(BaselineWK, LSQFitPolynomialSuccessfulCases) {
	enum NPCases {
		NP_kNo, NP_kCoeff, NP_kBestFit, NP_kResidual, NP_kAll, NP_kNumElems
	};
	vector<string> np_cases_names = { "no nullptr", "coeff=nullptr",
			"best_fit=nullptr", "residual=nullptr", "all nullptr" };
	cout << "    Testing for ";

	size_t const order = 3;
	SIMD_ALIGN
	double coeff_answer[order + 1];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff_answer), coeff_answer);
	SIMD_ALIGN
	double coeff[ELEMENTSOF(coeff_answer)];
	float rms;
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;

	size_t const num_data = ELEMENTSOF(coeff_answer);
	SIMD_ALIGN
	float data[num_data];
	SetFloatPolynomial(ELEMENTSOF(coeff_answer), coeff_answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	if (verbose) {
		PrintArray("data", num_data, data);
	}
	SIMD_ALIGN
	float best_fit[num_data];
	SIMD_ALIGN
	float residual[num_data];
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateLSQFitContextPolynomialFloat(
			sakura_LSQFitType_kPolynomial, order, num_data, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	for (NPCases item = static_cast<NPCases>(0); item < NP_kNumElems; item =
			static_cast<NPCases>(item + 1)) {
		cout << np_cases_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		double *coeff_ptr = coeff;
		float *best_fit_ptr = best_fit;
		float *residual_ptr = residual;

		switch (item) {
		case NP_kNo:
			break;
		case NP_kCoeff:
			coeff_ptr = nullptr;
			break;
		case NP_kBestFit:
			best_fit_ptr = nullptr;
			break;
		case NP_kResidual:
			residual_ptr = nullptr;
			break;
		case NP_kAll:
			coeff_ptr = nullptr;
			best_fit_ptr = nullptr;
			residual_ptr = nullptr;
			break;
		default:
			assert(false);
		}
		LIBSAKURA_SYMBOL (Status) fit_status =
		LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, order, num_data, data,
				mask, 5.0f, 1, ELEMENTSOF(coeff), coeff_ptr, best_fit_ptr,
				residual_ptr, mask, &rms, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), fit_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);

		bool check_coeff = true;
		bool check_best_fit = true;
		bool check_residual = true;
		if (item == NP_kCoeff || item == NP_kAll) {
			check_coeff = false;
		}
		if (item == NP_kBestFit || item == NP_kAll) {
			check_best_fit = false;
		}
		if (item == NP_kResidual || item == NP_kAll) {
			check_residual = false;
		}
		if (check_coeff) {
			for (size_t i = 0; i < ELEMENTSOF(coeff_answer); ++i) {
				CheckAlmostEqual(coeff_answer[i], coeff[i], 1.0e-6);
			}
		}
		if (check_best_fit) {
			for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
				CheckAlmostEqual(data[i], best_fit[i], 1.0e-6);
			}
		}
		if (check_residual) {
			for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
				CheckAlmostEqual(0.0, residual[i], 1.0e-6);
			}
		}
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test LSQFitCubicSpline
 * successful case
 */
TEST_F(BaselineWK, LSQFitCubicSplineSuccessfulCases) {
	enum NPCases {
		NP_kNo, NP_kCoeff, NP_kBestFit, NP_kResidual, NP_kAll, NP_kNumElems
	};
	vector<string> np_cases_names = { "no nullptr", "coeff=nullptr",
			"best_fit=nullptr", "residual=nullptr", "all nullptr" };
	cout << "    Testing for ";

	size_t const num_pieces = 3;
	SIMD_ALIGN
	size_t boundary[num_pieces + 1];
	SIMD_ALIGN
	double coeff_answer[4 * num_pieces];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff_answer), coeff_answer);
	SIMD_ALIGN
	double coeff[num_pieces][4];
	float rms;
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;

	size_t const num_data = ELEMENTSOF(coeff_answer);
	SIMD_ALIGN
	float data[num_data];
	SetFloatPolynomial(4, coeff_answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	if (verbose) {
		PrintArray("data", num_data, data);
	}
	SIMD_ALIGN
	float best_fit[num_data];
	SIMD_ALIGN
	float residual[num_data];
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextCubicSplineFloat(num_pieces, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	for (NPCases item = static_cast<NPCases>(0); item < NP_kNumElems; item =
			static_cast<NPCases>(item + 1)) {
		cout << np_cases_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		double (*coeff_ptr)[4] = coeff;
		float *best_fit_ptr = best_fit;
		float *residual_ptr = residual;

		switch (item) {
		case NP_kNo:
			break;
		case NP_kCoeff:
			coeff_ptr = nullptr;
			break;
		case NP_kBestFit:
			best_fit_ptr = nullptr;
			break;
		case NP_kResidual:
			residual_ptr = nullptr;
			break;
		case NP_kAll:
			coeff_ptr = nullptr;
			best_fit_ptr = nullptr;
			residual_ptr = nullptr;
			break;
		default:
			assert(false);
		}
		LIBSAKURA_SYMBOL (Status) fit_status =
		LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(context, num_pieces, num_data,
				data, mask, 5.0f, 1, coeff_ptr, best_fit_ptr, residual_ptr,
				mask, &rms, boundary, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), fit_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);

		bool check_coeff = true;
		bool check_best_fit = true;
		bool check_residual = true;
		if (item == NP_kCoeff || item == NP_kAll) {
			check_coeff = false;
		}
		if (item == NP_kBestFit || item == NP_kAll) {
			check_best_fit = false;
		}
		if (item == NP_kResidual || item == NP_kAll) {
			check_residual = false;
		}
		if (check_coeff) {
			for (size_t i = 0; i < ELEMENTSOF(coeff_answer); ++i) {
				CheckAlmostEqual(coeff_answer[i], coeff[i / 4][i % 4], 1.0e-5);
			}
		}
		if (check_best_fit) {
			for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
				CheckAlmostEqual(data[i], best_fit[i], 1.0e-5);
			}
		}
		if (check_residual) {
			for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
				CheckAlmostEqual(0.0, residual[i], 1.0e-5);
			}
		}
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}

/*
 * Test LSQFitSinusoid
 * successful case
 */
TEST_F(BaselineWK, LSQFitSinusoidSuccessfulCases) {
	enum NPCases {
		NP_kNo, NP_kCoeff, NP_kBestFit, NP_kResidual, NP_kAll, NP_kNumElems
	};
	vector<string> np_cases_names = { "no nullptr", "coeff=nullptr",
			"best_fit=nullptr", "residual=nullptr", "all nullptr" };
	cout << "    Testing for ";

	size_t const nwave_max = 3;
	size_t const num_nwave = nwave_max + 1;
	size_t const nwave[num_nwave] = { 0, 1, 2, 3 };
	SIMD_ALIGN
	double coeff_answer[nwave_max * 2 + 1];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff_answer), coeff_answer);
	SIMD_ALIGN
	double coeff[ELEMENTSOF(coeff_answer)];
	float rms;
	LIBSAKURA_SYMBOL(LSQFitStatus) baseline_status;

	size_t const num_data = ELEMENTSOF(coeff_answer) + 5;
	SIMD_ALIGN
	float data[num_data];
	SetFloatSinusoidal(num_nwave, nwave, coeff_answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	if (verbose) {
		PrintArray("data", num_data, data);
	}
	SIMD_ALIGN
	float best_fit[num_data];
	SIMD_ALIGN
	float residual[num_data];
	LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
			sakura_CreateLSQFitContextSinusoidFloat(nwave_max, num_data,
					&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	for (NPCases item = static_cast<NPCases>(0); item < NP_kNumElems; item =
			static_cast<NPCases>(item + 1)) {
		cout << np_cases_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		double *coeff_ptr = coeff;
		float *best_fit_ptr = best_fit;
		float *residual_ptr = residual;

		switch (item) {
		case NP_kNo:
			break;
		case NP_kCoeff:
			coeff_ptr = nullptr;
			break;
		case NP_kBestFit:
			best_fit_ptr = nullptr;
			break;
		case NP_kResidual:
			residual_ptr = nullptr;
			break;
		case NP_kAll:
			coeff_ptr = nullptr;
			best_fit_ptr = nullptr;
			residual_ptr = nullptr;
			break;
		default:
			assert(false);
		}
		LIBSAKURA_SYMBOL (Status) fit_status =
		LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context, num_nwave, nwave,
				num_data, data, mask, 5.0f, 1, ELEMENTSOF(coeff), coeff_ptr,
				best_fit_ptr, residual_ptr, mask, &rms, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), fit_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(LSQFitStatus_kOK), baseline_status);

		bool check_coeff = true;
		bool check_best_fit = true;
		bool check_residual = true;
		if (item == NP_kCoeff || item == NP_kAll) {
			check_coeff = false;
		}
		if (item == NP_kBestFit || item == NP_kAll) {
			check_best_fit = false;
		}
		if (item == NP_kResidual || item == NP_kAll) {
			check_residual = false;
		}
		if (check_coeff) {
			for (size_t i = 0; i < ELEMENTSOF(coeff_answer); ++i) {
				CheckAlmostEqual(coeff_answer[i], coeff[i], 1.0e-6);
			}
		}
		if (check_best_fit) {
			for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
				CheckAlmostEqual(data[i], best_fit[i], 1.0e-6);
			}
		}
		if (check_residual) {
			for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
				CheckAlmostEqual(0.0, residual[i], 1.0e-6);
			}
		}
	}

	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyLSQFitContextFloat(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;
}
