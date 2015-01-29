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
 * Test cases to be implemented.
 * functions to be tested include:
 * - sakura_SubtractBaselineCubicSplineFloat
 * - sakura_GetBestFitBaselineCoefficientsCubicSplineFloat
 * test cases are as follows:
 *     (1) simple successful case (num_pieces=1,2,3, num_data=4*num_pieces+a where a=0,1,2,3,10
 *     (2) time-consuming successful case for performance measurement
 *     (3) error cases
 *         (3-1) null pointer cases
 *         (3-2) non-aligned array cases
 *         (3-3) bad parameter value cases
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

	// Set constant float values into an array
	void SetFloatConstant(float value, size_t const num_data, float *data) {
		for (size_t i = 0; i < num_data; ++i) {
			data[i] = value;
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
 * Test GetBestFitBaselineCoefficientsCubicSpline
 * successful case
 * compute the best-fit baseline coefficients for cubic spline
 * for cases of combination of num_pieces=(1,2,3) and num_data=(4*num_pieces+a) where a=(0,1,2,3)
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplineSuccessfulCase) {
	for (size_t num_pieces = 1; num_pieces <= 3; ++num_pieces) {
		cout << "    Testing for num_pieces = " << num_pieces
				<< " cases: num_data = ";
		size_t num_extra_max = 3;
		for (size_t num_extra = 0; num_extra <= num_extra_max; ++num_extra) {
			SIMD_ALIGN
			double coeff[4 * num_pieces];
			SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
			size_t const num_data(ELEMENTSOF(coeff) + num_extra);
			SIMD_ALIGN
			float in_data[num_data];
			SetFloatPolynomial(num_data, in_data, coeff);
			SIMD_ALIGN
			bool mask[ELEMENTSOF(in_data)];
			SetBoolConstant(true, ELEMENTSOF(in_data), mask);
			cout << num_data << ((num_extra < num_extra_max) ? ", " : "");
			if (verbose) {
				PrintArray("in_data", num_data, in_data);
			}
			LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
			LIBSAKURA_SYMBOL (Status) create_status =
					sakura_CreateBaselineContext(
							LIBSAKURA_SYMBOL(BaselineType_kCubicSpline),
							num_pieces, num_data, &context);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
			SIMD_ALIGN
			double out[ELEMENTSOF(coeff)];
			double answer[ELEMENTSOF(coeff)];
			for (size_t i = 0; i < ELEMENTSOF(coeff); ++i) {
				answer[i] = coeff[i];
			}

			LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;
			LIBSAKURA_SYMBOL (Status) coeff_status = LIBSAKURA_SYMBOL(
					GetBestFitBaselineCoefficientsCubicSplineFloat)(context,
					num_data, in_data, mask, 5.0f, 1, num_pieces, out, mask,
					&baseline_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

			for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
				CheckAlmostEqual(answer[i], out[i], 1.0e-6);
			}

			if (verbose) {
				PrintArray("data  ", num_data, in_data);
				PrintArray("out   ", ELEMENTSOF(answer), out);
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
 * time-consuming successful case for performance measurement
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplinePerformanceTest) {
	size_t const num_repeat(300);
	size_t const num_pieces(2);
	SIMD_ALIGN
	double answer[4 * num_pieces];
	for (size_t i = 0; i < ELEMENTSOF(answer); i += 4) {
		answer[i] = 1.0;
		answer[i + 1] = 1e-4;
		answer[i + 2] = 1e-8;
		answer[i + 3] = 1e-12;
	}
	size_t const num_data(70000);
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
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), num_pieces, num_data,
			&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	double out[ELEMENTSOF(answer)];
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;
	for (size_t i = 0; i < num_repeat; ++i) {
		LIBSAKURA_SYMBOL (Status) coeff_status = LIBSAKURA_SYMBOL(
				GetBestFitBaselineCoefficientsCubicSplineFloat)(context,
				num_data, in_data, mask, 5.0f, 1, num_pieces, out, mask,
				&baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
	}
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
		NP_kBaselineStatus,
		NP_kNumElems
	};
	vector<string> np_param_names = { "context", "data", "mask", "coeff",
			"final_mask", "baseline_status" };
	cout << "    Testing for ";

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");
		size_t const num_pieces(2);
		size_t const num_data(10);
		double coeff[4];
		SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
		SIMD_ALIGN
		float in_data[num_data];
		SetFloatPolynomial(num_data, in_data, coeff);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(in_data)];
		SetBoolConstant(true, ELEMENTSOF(in_data), mask);
		SIMD_ALIGN
		bool final_mask[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		double out[ELEMENTSOF(coeff)];
		LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), num_pieces,
				num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

		float *in_data_ptr = in_data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		double *out_ptr = out;
		LIBSAKURA_SYMBOL(BaselineContext) *context_ptr = context;
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status_ptr = &baseline_status;

		switch (item) {
		case NP_kData:
			in_data_ptr = nullptr;
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
		case NP_kBaselineStatus:
			baseline_status_ptr = nullptr;
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL (Status) coeff_status = LIBSAKURA_SYMBOL(
				GetBestFitBaselineCoefficientsCubicSplineFloat)(context_ptr,
				num_data, in_data_ptr, mask_ptr, 5.0f, 1, num_pieces, out_ptr,
				final_mask_ptr, baseline_status_ptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	cout << endl;
}

/*
 * Test GetBestFitBaselineCoefficientsCubicSpline
 * erroneous cases: unaligned cases
 * parameters to be tested include data, mask, coeff and final_mask.
 */
TEST_F(BaselineWK, GetBestFitBaselineCoefficientsCubicSplineErroneousCasesUnaligned) {
	enum UAItems {
		UA_kData, UA_kMask, UA_kCoeff, UA_kFinalMask, UA_kNumElems
	};
	vector<string> ua_param_names = { "data", "mask", "coeff", "final_mask" };
	cout << "    Testing for ";

	for (UAItems item = static_cast<UAItems>(0); item < UA_kNumElems; item =
			static_cast<UAItems>(item + 1)) {
		cout << ua_param_names[item] << ((item < UA_kNumElems - 1) ? ", " : "");
		size_t const num_pieces(2);
		size_t const num_data(10);
		double coeff[4];
		SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
		SIMD_ALIGN
		float in_data[num_data + 1];
		SetFloatPolynomial(num_data, in_data, coeff);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(in_data)];
		SetBoolConstant(true, ELEMENTSOF(in_data), mask);
		SIMD_ALIGN
		bool final_mask[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		double out[ELEMENTSOF(coeff) + 1];
		LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), num_pieces,
				num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

		float *in_data_ptr = in_data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		double *out_ptr = out;

		switch (item) {
		case UA_kData:
			++in_data_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(in_data_ptr));
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
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL (Status) coeff_status =
		LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsCubicSplineFloat)(
				context, num_data, in_data_ptr, mask_ptr, 5.0f, 1, num_pieces,
				out_ptr, final_mask_ptr, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
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

	for (BVItems item = static_cast<BVItems>(0); item < BV_kNumElems; item =
			static_cast<BVItems>(item + 1)) {
		cout << "        " << bv_param_names[item]
				<< ((item < BV_kNumElems - 1) ? ", " : "") << endl;
		size_t const num_pieces(1);
		size_t num_data;
		size_t num_basis_data(10);
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

		double coeff[4];
		SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
		SIMD_ALIGN
		float in_data[num_data];
		SetFloatPolynomial(num_data, in_data, coeff);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(in_data)];
		SetBoolConstant(true, ELEMENTSOF(in_data), mask);
		SIMD_ALIGN
		bool final_mask[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		double out[ELEMENTSOF(coeff)];
		LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), num_pieces,
				num_basis_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

		LIBSAKURA_SYMBOL (Status) coeff_status = LIBSAKURA_SYMBOL(
				GetBestFitBaselineCoefficientsCubicSplineFloat)(context,
				num_data, in_data, mask, 5.0f, 1, num_pieces, out, final_mask,
				&baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

/*
 * Test sakura_SubtractBaselineCubicSplineFloat
 * successful case
 * subtract best fit model from input data
 */
TEST_F(BaselineWK, SubtractBaselineCubicSplineSuccessfulCase) {
	for (size_t num_pieces = 1; num_pieces <= 3; ++num_pieces) {
		cout << "    Testing for num_pieces = " << num_pieces
				<< " cases: num_data = ";
		size_t num_extra_max = 3;
		for (size_t num_extra = 0; num_extra <= num_extra_max; ++num_extra) {
			size_t const num_data(4 * num_pieces + num_extra);
			cout << num_data << ((num_extra < num_extra_max) ? ", " : "");
			double coeff[4];
			SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
			SIMD_ALIGN
			float in_data[num_data];
			SetFloatPolynomial(num_data, in_data, coeff);
			SIMD_ALIGN
			bool mask[ELEMENTSOF(in_data)];
			SetBoolConstant(true, ELEMENTSOF(in_data), mask);
			if (verbose) {
				PrintArray("in_data", num_data, in_data);
			}
			LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
			LIBSAKURA_SYMBOL (Status) create_status =
					sakura_CreateBaselineContext(
							LIBSAKURA_SYMBOL(BaselineType_kCubicSpline),
							num_pieces, num_data, &context);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
			SIMD_ALIGN
			float out[ELEMENTSOF(in_data)];
			SIMD_ALIGN
			float answer[ELEMENTSOF(in_data)];
			SetFloatConstant(0.0, ELEMENTSOF(in_data), answer);

			LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;
			LIBSAKURA_SYMBOL(Status) sub_status =
			LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(context,
					num_pieces, num_data, in_data, mask, 5.0f, 1, true, mask,
					out, &baseline_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

			for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
				CheckAlmostEqual(answer[i], out[i], 1.0e-6);
			}

			if (verbose) {
				PrintArray("data  ", num_data, in_data);
				PrintArray("out   ", ELEMENTSOF(answer), out);
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
	size_t const num_repeat(300);
	size_t const num_pieces(2);
	size_t const num_data(70000);
	SIMD_ALIGN
	float in_data[num_data];
	double coeff[4] = { 1.0, 1e-4, 1e-8, 1e-12 };
	SetFloatPolynomial(num_data, in_data, coeff);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(in_data)];
	SetBoolConstant(true, ELEMENTSOF(in_data), mask);
	if (verbose) {
		PrintArray("in_data", num_data, in_data);
	}
	LIBSAKURA_SYMBOL(BaselineContext) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
			LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), num_pieces, num_data,
			&context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;
	SIMD_ALIGN
	double answer[ELEMENTSOF(in_data)];
	SetDoubleConstant(0.0, ELEMENTSOF(in_data), answer);
	for (size_t i = 0; i < num_repeat; ++i) {
		LIBSAKURA_SYMBOL(Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(context, num_pieces,
				num_data, in_data, mask, 5.0f, 1, true, mask, out,
				&baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), sub_status);
	}
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
		NP_kBaselineStatus,
		NP_kNumElems
	};
	vector<string> np_param_names = { "context", "data", "mask", "final_mask",
			"out", "baseline_status" };
	cout << "    Testing for ";

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");
		size_t const num_pieces(2);
		size_t const num_data(10);
		double coeff[4];
		SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
		SIMD_ALIGN
		float in_data[num_data];
		SetFloatPolynomial(num_data, in_data, coeff);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(in_data)];
		SetBoolConstant(true, ELEMENTSOF(in_data), mask);
		SIMD_ALIGN
		bool final_mask[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), num_pieces,
				num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

		float *in_data_ptr = in_data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		float *out_ptr = out;
		LIBSAKURA_SYMBOL(BaselineContext) *context_ptr = context;
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status_ptr = &baseline_status;

		switch (item) {
		case NP_kData:
			in_data_ptr = nullptr;
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
		case NP_kContext:
			context_ptr = nullptr;
			break;
		case NP_kBaselineStatus:
			baseline_status_ptr = nullptr;
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL(Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(context_ptr,
				num_pieces, num_data, in_data_ptr, mask_ptr, 5.0f, 1, true,
				final_mask_ptr, out_ptr, baseline_status_ptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
	cout << endl;
}

/*
 * Test SubtractBaselineCubicSplineFloat
 * erroneous cases: unaligned cases
 * parameters to be tested include data, mask, final_mask and out
 */
TEST_F(BaselineWK, SubtractBaselineCubicSplineErroneousCasesUnaligned) {
	enum UAItems {
		UA_kData, UA_kMask, UA_kFinalMask, UA_kOut, UA_kNumElems
	};
	vector<string> ua_param_names = { "data", "mask", "final_mask", "out" };
	cout << "    Testing for ";

	for (UAItems item = static_cast<UAItems>(0); item < UA_kNumElems; item =
			static_cast<UAItems>(item + 1)) {
		cout << ua_param_names[item] << ((item < UA_kNumElems - 1) ? ", " : "");
		size_t const num_pieces(2);
		size_t const num_data(10);
		double coeff[4];
		SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
		SIMD_ALIGN
		float in_data[num_data + 1];
		SetFloatPolynomial(num_data, in_data, coeff);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(in_data)];
		SetBoolConstant(true, ELEMENTSOF(in_data), mask);
		SIMD_ALIGN
		bool final_mask[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), num_pieces,
				num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

		float *in_data_ptr = in_data;
		bool *mask_ptr = mask;
		bool *final_mask_ptr = final_mask;
		float *out_ptr = out;

		switch (item) {
		case UA_kData:
			++in_data_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(in_data_ptr));
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
		LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(context, num_pieces,
				num_data, in_data_ptr, mask_ptr, 5.0f, 1, true, final_mask_ptr,
				out_ptr, &baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
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

	for (BVItems item = static_cast<BVItems>(0); item < BV_kNumElems; item =
			static_cast<BVItems>(item + 1)) {
		cout << "        " << bv_param_names[item]
				<< ((item < BV_kNumElems - 1) ? ", " : "") << endl;
		size_t const num_pieces(1);
		size_t num_data;
		size_t num_basis_data(10);
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

		double coeff[4];
		SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
		SIMD_ALIGN
		float in_data[num_data];
		SetFloatPolynomial(num_data, in_data, coeff);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(in_data)];
		SetBoolConstant(true, ELEMENTSOF(in_data), mask);
		SIMD_ALIGN
		bool final_mask[ELEMENTSOF(in_data)];
		SIMD_ALIGN
		float out[ELEMENTSOF(in_data)];
		LIBSAKURA_SYMBOL(BaselineContext) *context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContext(
				LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), num_pieces,
				num_basis_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

		LIBSAKURA_SYMBOL(Status) sub_status =
		LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(context, num_pieces,
				num_data, in_data, mask, 5.0f, 1, true, final_mask, out,
				&baseline_status);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), sub_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyBaselineContext(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}
