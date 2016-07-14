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
 * test order parameters in
 * - sakura_GetNumberOfCoefficientsFloat
 * - sakura_SubtractBaselineFloat
 * the test cases necessary are,
 * - successful cases, order is smaller than max available number by context (polynomial, chebychev, cspline=only for GetNumberOfCoefficientsFloat)
 * - failure cases where order is larger than max available number by context (polynomial, chebychev)
 */

#include <cmath>
#include <iostream>
#include <string>
#include <map>
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
class BaselineKS: public ::testing::Test {
protected:

	BaselineKS() :
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

	void GenerateFromContext(
	LIBSAKURA_SYMBOL(LSQFitContextFloat) const *context, size_t const num_data,
			float *in_data, size_t const num_coeff, double *coeff) {
		assert(num_data == context->num_basis_data);
		assert(num_coeff <= context->num_bases);
		SetFloatConstant(0, num_data, in_data);
		size_t const num_bases(context->num_bases);
		SIMD_ALIGN
		double coeff_apply[num_coeff];
		std::copy(coeff, coeff + num_coeff, coeff_apply);
		if (context->lsqfit_type == LSQFitTypeInternal_kPolynomial) {
			double const max_data_x =
					static_cast<double>(context->num_basis_data - 1);
			double factor = 1.0;
			for (size_t i = 0; i < num_coeff; ++i) {
				coeff_apply[i] *= factor;
				factor *= max_data_x;
			}
		}
		for (size_t i = 0; i < num_data; ++i) {
			for (size_t k = 0; k < num_coeff; ++k) {
				in_data[i] += coeff_apply[k]
						* context->basis_data[i * num_bases + k];
			}
		}
	}

	bool verbose;
	uint16_t const nwave = 0;
};

/*
 * Test sakura_SubtractBaselineFloat with
 * {input order=3} < {the other for context generation=5}.
 * Fitting functions: polynomial, Chebyshev
 */
TEST_F(BaselineKS, SubtractBaselineOrder) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data];
	double coeff[] = { 4., 1., -1., 0.25 };
	SetFloatPolynomial(ELEMENTSOF(coeff), coeff, num_data, in_data);
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	SetBoolConstant(true, ELEMENTSOF(in_data), in_mask);
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(in_data)];
	SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);
	size_t const gen_order(5);
	size_t const in_order(3);
	LIBSAKURA_SYMBOL(LSQFitType) bltypes[] = { LIBSAKURA_SYMBOL(
			LSQFitType_kPolynomial), LIBSAKURA_SYMBOL(LSQFitType_kChebyshev) };

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
		cout << "order (context) = " << gen_order << endl;
		cout << "order (fitting) = " << in_order << endl;
	}

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	for (size_t i = 0; i < ELEMENTSOF(bltypes); ++i) {
		LIBSAKURA_SYMBOL(LSQFitType) type(bltypes[i]);
		cout << "Testing baseline type = " << type << endl;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(type, gen_order,
						num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		LIBSAKURA_SYMBOL (LSQFitStatus) op_blstatus;
		LIBSAKURA_SYMBOL (Status) op_status;
		op_status = LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, in_order,
				num_data, in_data, in_mask, clipping_threshold_sigma,
				num_fitting_max, in_order + 1, nullptr, nullptr, out,
				final_mask, &rms, &op_blstatus);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), op_status);
		for (size_t i = 0; i < num_data; ++i) {
			ASSERT_EQ(answer[i], out[i]);
		}
		if (verbose) {
			PrintArray("fmask ", num_data, final_mask);
			PrintArray("out   ", num_data, out);
			PrintArray("answer", num_data, answer);
		}

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyLSQFitContextFloat(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

/*
 * Test sakura_SubtractBaselineFloat with
 * invalid input order (> order for context generation).
 * Fitting functions: polynomial, Chebyshev
 */
TEST_F(BaselineKS, SubtractBaselineBadOrder) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	bool final_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	size_t const gen_order(5);
	size_t const in_order(6);
	LIBSAKURA_SYMBOL(LSQFitType) bltypes[] = { LIBSAKURA_SYMBOL(
			LSQFitType_kPolynomial), LIBSAKURA_SYMBOL(LSQFitType_kChebyshev) };

	if (verbose) {
		cout << "order (context) = " << gen_order << endl;
		cout << "order (fitting) = " << in_order << endl;
	}

	float clipping_threshold_sigma = 3.0;
	uint16_t num_fitting_max = 1;
	float rms;

	for (size_t i = 0; i < ELEMENTSOF(bltypes); ++i) {
		LIBSAKURA_SYMBOL(LSQFitType) type(bltypes[i]);
		cout << "Testing baseline type = " << type << endl;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(type, gen_order,
						num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		LIBSAKURA_SYMBOL (LSQFitStatus) op_blstatus;
		LIBSAKURA_SYMBOL (Status) op_status;
		op_status = LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context, in_order,
				num_data, in_data, in_mask, clipping_threshold_sigma,
				num_fitting_max, in_order + 1, nullptr, nullptr, out,
				final_mask, &rms, &op_blstatus);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), op_status);

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyLSQFitContextFloat(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

/*
 * Test sakura_GetNumberOfCoefficientsFloat with
 * {input order=2} < {the other for context generation=5}.
 * Fitting function: polynomial, Chebyshev, CubicSpline
 */
TEST_F(BaselineKS, GetNumberOfCoefficientsFloatOrder) {
	size_t const num_data(NUM_DATA2);
	size_t gen_order = 5; // the order to generate a context
	size_t test_order = 2; // the order to test
	//uint16_t const test_nwave = 2;
	map<LSQFitTypeInternal, size_t> answers;
	answers[LSQFitTypeInternal_kPolynomial] = test_order + 1;
	answers[LSQFitTypeInternal_kChebyshev] = test_order + 1;
	answers[LSQFitTypeInternal_kCubicSpline] = 4 * test_order;
	LSQFitTypeInternal bltypes[] = { LSQFitTypeInternal_kPolynomial,
			LSQFitTypeInternal_kChebyshev, LSQFitTypeInternal_kCubicSpline };
	for (size_t i = 0; i < ELEMENTSOF(bltypes); ++i) {
		LSQFitTypeInternal type(bltypes[i]);
		cout << "Testing baseline type = " << type << endl;
		if (answers.find(type) == answers.end())
			continue;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status;
		if (type == LSQFitTypeInternal_kCubicSpline) {
			create_status = sakura_CreateLSQFitContextCubicSplineFloat(
					gen_order, num_data, &context);
		} else {
			LIBSAKURA_SYMBOL(LSQFitType) type_external = LIBSAKURA_SYMBOL(
					LSQFitType_kPolynomial);
			if (type == LSQFitTypeInternal_kChebyshev) {
				type_external = LIBSAKURA_SYMBOL(LSQFitType_kChebyshev);
			}
			create_status = sakura_CreateLSQFitContextPolynomialFloat(
					type_external, gen_order, num_data, &context);
		}

		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		size_t reference(answers[type]);
		size_t num_coeff = 0;
		LIBSAKURA_SYMBOL (Status) num_status;
		num_status = sakura_GetNumberOfCoefficientsFloat(context, test_order,
				&num_coeff);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), num_status);
		EXPECT_EQ(num_coeff, reference);
		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyLSQFitContextFloat(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

/*
 * Test sakura_GetNumberOfCoefficientsFloat with
 * invalid input order (> order for context generation).
 * Fitting functions: polynomial, Chebyshev
 */
TEST_F(BaselineKS, GetNumberOfCoefficientsFloatBadOrder) {
	size_t const num_data(NUM_DATA2);
	size_t gen_order = 5; // the order to generate a context
	size_t test_order = 6; // the order to test
	LIBSAKURA_SYMBOL(LSQFitType) bltypes[] = { LIBSAKURA_SYMBOL(
			LSQFitType_kPolynomial), LIBSAKURA_SYMBOL(LSQFitType_kChebyshev) //,
//					LIBSAKURA_SYMBOL(LSQFitType_kCubicSpline), // this never fails
			};
	for (size_t i = 0; i < ELEMENTSOF(bltypes); ++i) {
		LIBSAKURA_SYMBOL(LSQFitType) type(bltypes[i]);
		cout << "Testing baseline type = " << type << endl;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateLSQFitContextPolynomialFloat(type, gen_order,
						num_data, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		size_t num_coeff = 0;
		LIBSAKURA_SYMBOL (Status) num_status;
		num_status = sakura_GetNumberOfCoefficientsFloat(context, test_order,
				&num_coeff);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), num_status);
		EXPECT_EQ(num_coeff, 0);
		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyLSQFitContextFloat(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

/*
 * Test sakura_SubtractBaselineUsingCoefficientsFloat
 * successful case with
 * {input num_coeff=4} < {the num_bases in context=6}.
 * Fitting function: polynomial, Chebyshev
 */
TEST_F(BaselineKS, SubtractBaselineUsingCoefficientsFloatNumCoeff) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	double coeff[] = { 4., 1., -1., 0.25 };
	size_t const num_coeff(ELEMENTSOF(coeff));
	size_t const order(num_coeff + 1); // make it just bigger than num_coeff will do.

	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	float answer[ELEMENTSOF(in_data)];
	SetFloatConstant(0.0f, ELEMENTSOF(in_data), answer);

	LIBSAKURA_SYMBOL(LSQFitType) bltypes[] = { LIBSAKURA_SYMBOL(
			LSQFitType_kPolynomial), LIBSAKURA_SYMBOL(LSQFitType_kChebyshev) };
	for (size_t i = 0; i < ELEMENTSOF(bltypes); ++i) {
		LIBSAKURA_SYMBOL(LSQFitType) type(bltypes[i]);
		cout << "Testing baseline type = " << type << endl;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status;
		create_status = sakura_CreateLSQFitContextPolynomialFloat(type, order,
				num_data, &context);

		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		// generate input data from basis in context
		GenerateFromContext(context, num_data, in_data, num_coeff, coeff);
		if (verbose) {
			PrintArray("in_data", num_data, in_data);
			PrintArray("coeff", num_coeff, coeff);
		}
		cout << "Fitting with num_coeff = " << num_coeff << " (num_bases = "
				<< context->num_bases << ")" << endl;
		LIBSAKURA_SYMBOL (Status) subbl_status;
		subbl_status = LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context,
				num_data, in_data, num_coeff, coeff, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), subbl_status);
		if (verbose) {
			PrintArray("out", num_data, out);
		}
		for (size_t i = 0; i < num_data; ++i) {
			CheckAlmostEqual(answer[i], out[i], 1.0e-6);
		}
		if (verbose) {
			PrintArray("out   ", num_data, out);
			PrintArray("answer", num_data, answer);
		}

		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyLSQFitContextFloat(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

/*
 * Test sakura_SubtractBaselineUsingCoefficientsFloat
 * failure case with
 * num_coeff > num_model
 * num_coeff = 0
 * Fitting function: Polynomial, Chebyshev, Sinusoid
 */
TEST_F(BaselineKS, SubtractBaselineUsingCoefficientsFloatBadNumCoeff) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data];
	size_t const order(4);
	size_t const bad_coeffs[] = { 2 * order + 5, 0 }; // Set large enough coeff
	SIMD_ALIGN
	double coeff[2 * order + 5]; // allocate big enough array
	SIMD_ALIGN
	float out[ELEMENTSOF(in_data)];

	LSQFitTypeInternal bltypes[] = { LSQFitTypeInternal_kPolynomial,
			LSQFitTypeInternal_kChebyshev, LSQFitTypeInternal_kSinusoid };

	for (size_t i = 0; i < ELEMENTSOF(bltypes); ++i) {
		LSQFitTypeInternal type(bltypes[i]);
		cout << "Testing baseline type = " << type << endl;
		LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status;
		if (type == LSQFitTypeInternal_kSinusoid) {
			create_status = sakura_CreateLSQFitContextSinusoidFloat(nwave,
					num_data, &context);
		} else {
			LIBSAKURA_SYMBOL(LSQFitType) type_external = LIBSAKURA_SYMBOL(
					LSQFitType_kPolynomial);
			if (type == LSQFitTypeInternal_kChebyshev) {
				type_external = LIBSAKURA_SYMBOL(LSQFitType_kChebyshev);
			}
			create_status = sakura_CreateLSQFitContextPolynomialFloat(
					type_external, order, num_data, &context);
		}
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		for (size_t j = 0; j < ELEMENTSOF(bad_coeffs); ++j) {
			size_t num_coeff = bad_coeffs[j];
			cout << "Fitting with num_coeff = " << num_coeff << " (num_bases = "
					<< context->num_bases << ")" << endl;
			LIBSAKURA_SYMBOL (Status) subbl_status;
			if (type == LSQFitTypeInternal_kSinusoid) {
				size_t nwave_list[nwave];
				subbl_status =
				LIBSAKURA_SYMBOL(SubtractSinusoidFloat)(context, num_data,
						in_data, nwave, nwave_list, num_coeff, coeff, out);
			} else {
				subbl_status =
				LIBSAKURA_SYMBOL(SubtractPolynomialFloat)(context, num_data,
						in_data, num_coeff, coeff, out);
			}
			ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), subbl_status);
		}
		LIBSAKURA_SYMBOL (Status) destroy_status =
				sakura_DestroyLSQFitContextFloat(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);
	}
}

