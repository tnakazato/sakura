/*
 * numeric_operation.cc
 *
 *  Created on: 2013/11/11
 *      Author: wataru
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

/* the number of elements in input/output array to test */
#define NUM_DATA 4096
#define NUM_MODEL 20
#define NUM_REPEAT 3000
#define NUM_REPEAT2 300
#define NUM_REPEAT3 1500000
#define NUM_EXCLUDE 5

using namespace std;

/*
 * A super class to test numeric operations of array(s)
 */
class NumericOperation: public ::testing::Test {
protected:

	NumericOperation() :
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

	// Set (1+x+x*x+x*x*x+x*x*x*x) float values into an array
	void SetFloatPolynomial(size_t const num_data, float *data) {
		for (size_t i = 0; i < num_data; ++i) {
			double x = (double) i;
			data[i] = (float) (1.0 + x + x * x + x * x * x + x * x * x * x);
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

	// Set input data using a polynomial
	void SetInputData(size_t const num_data, float *data) {
		for (size_t i = 0; i < num_data; ++i) {
			double x = (double) i;
			double val = 4.0 + 0.000056 * x - 0.000037 * x * x
					+ 0.0000012 * x * x * x + 0.0000009 * x * x * x * x
					+ 0.0000006 * x * x * x * x * x;
			data[i] = (float) val;
		}
	}

	// Set polynomial model
	void SetPolynomialModel(size_t const num_data, size_t const num_model,
			double *model) {
		size_t idx = 0;
		for (size_t j = 0; j < num_data; ++j) {
			double value = 1.0;
			for (size_t i = 0; i < num_model; ++i) {
				model[idx] = value;
				value *= (double) j;
				idx++;
			}
		}
	}

	// Set Chebyshev polynomial model
	void SetChebyshevModel(size_t const num_data, size_t const num_model,
			double *model) {
		size_t idx = 0;
		for (size_t j = 0; j < num_data; ++j) {
			for (size_t i = 0; i < num_model; ++i) {
				double value = 1.0;
				double x = 2.0 * (double) j / (double) (num_data - 1) - 1.0;
				if (i == 0) {
					value = 1.0;
				} else if (i == 1) {
					value = x;
				} else {
					value = 2.0 * x * model[idx - 1] - model[idx - 2];
				}
				model[idx] = value;
				idx++;
			}
		}
	}

	// Set reference values
	void SetAnswers(size_t const num_data, float const *data,
			size_t const num_model, double const *model, double *answer_matrix,
			double *answer_vector) {
		size_t idx = 0;
		for (size_t i = 0; i < num_model; ++i) {
			for (size_t j = 0; j < num_model; ++j) {
				double val = 0.0;
				for (size_t k = 0; k < num_data; ++k) {
					val += model[num_model * k + i] * model[num_model * k + j];
				}
				answer_matrix[idx] = val;
				idx++;
			}
		}
		for (size_t i = 0; i < num_model; ++i) {
			double val = 0.0;
			for (size_t j = 0; j < num_data; ++j) {
				val += model[num_model * j + i] * data[j];
			}
			answer_vector[i] = val;
		}
	}

	// Set reference values for testing Update
	void SetAnswers(size_t const num_data, float const *data,
	bool const *mask, size_t const num_model, double const *model,
			double *answer_matrix, double *answer_vector) {
		size_t idx = 0;
		for (size_t i = 0; i < num_model; ++i) {
			for (size_t j = 0; j < num_model; ++j) {
				double val = 0.0;
				for (size_t k = 0; k < num_data; ++k) {
					if (mask[k]) {
						val += model[num_model * k + i]
								* model[num_model * k + j];
					}
				}
				answer_matrix[idx] = val;
				idx++;
			}
		}
		for (size_t i = 0; i < num_model; ++i) {
			double val = 0.0;
			for (size_t j = 0; j < num_data; ++j) {
				if (mask[j]) {
					val += model[num_model * j + i] * data[j];
				}
			}
			answer_vector[i] = val;
		}
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
				cout << ", ";
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
				cout << ", ";
		}
		cout << " ]" << endl;
	}

	bool verbose;

};

/*
 * Test sakura_GetLeastSquareFittingCoefficients
 * successful case
 * note: repeating NUM_REPEAT times for performance measurement
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficients) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double answer[num_model * num_model];
	SIMD_ALIGN
	double answer_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, answer, answer_vector);

	if (verbose) {
		PrintArray("in_mask", num_data, in_mask);
		PrintArray("model  ", num_data, num_model, model);
	}

	size_t const num_repeat(NUM_REPEAT);
	double start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		LIBSAKURA_SYMBOL(Status) status =
				sakura_GetLeastSquareFittingCoefficients(num_data, in_data,
						in_mask, num_model, model, out, out_vector);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}
	double end = sakura_GetCurrentTime();

	for (size_t i = 0; i < num_model * num_model; ++i) {
		double deviation;
		if (answer[i] != 0.0) {
			deviation = fabs((out[i] - answer[i]) / answer[i]);
		} else if (out[i] != 0.0) {
			deviation = fabs((out[i] - answer[i]) / out[i]);
		} else {
			deviation = fabs(out[i] - answer[i]);
		}
		ASSERT_LE(deviation, 1e-7);
	}
	for (size_t i = 0; i < num_model; ++i) {
		double deviation;
		if (answer_vector[i] != 0.0) {
			deviation = fabs(
					(out_vector[i] - answer_vector[i]) / answer_vector[i]);
		} else if (out_vector[i] != 0.0) {
			deviation = fabs(
					(out_vector[i] - answer_vector[i]) / out_vector[i]);
		} else {
			deviation = fabs(out_vector[i] - answer_vector[i]);
		}
		ASSERT_LE(deviation, 1e-7);
	}

	if (verbose) {
		cout << "Elapse time of " << num_repeat << " repetition: "
				<< end - start << " sec." << endl;
		PrintArray("out   ", num_model, num_model, out);
		PrintArray("answer", num_model, num_model, answer);
	}
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsNumDataZero
 * failure case : num_data == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsNumDataZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask, num_model, model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsWithDataNullPointer
 * failure case : data is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsWithDataNullPointer) {
	size_t const num_data(NUM_DATA);
	float *in_data = nullptr;
	SIMD_ALIGN
	bool in_mask[num_data];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_mask)];

	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask, num_model, model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsWithDataNotAligned
 * failure case : data is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsWithDataNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data + 1];
	float *in_data_unaligned = in_data + 1;
	SIMD_ALIGN
	bool in_mask[num_data];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_mask)];

	SetInputData(num_data, in_data_unaligned);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data_unaligned, in_mask, num_model, model, out,
			out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsWithMaskNullPointer
 * failure case : mask is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsWithMaskNullPointer) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	bool *in_mask = nullptr;
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask, num_model, model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsWithMaskNotAligned
 * failure case : mask is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsWithMaskNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data) + 1];
	bool *in_mask_unaligned = in_mask + 1;
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask_unaligned);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask_unaligned, num_model, model, out,
			out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsWithNumModelBasesZero
 * failure case : num_model_bases == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsWithNumModelBasesZero) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(0);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask, num_model, model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsWithNumDataLessThanNumModelBases
 * failure case : num_model_bases > num_data
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsWithNumDataLessThanNumModelBases) {
	size_t const num_data(NUM_MODEL - 1);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask, num_model, model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsWithBasisDataNullPointer
 * failure case : basis_data is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsWithBasisDataNullPointer) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	double *model = nullptr;

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask, num_model, model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsWithBasisDataNotAligned
 * failure case : basis_data is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsWithBasisDataNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data) + 1];
	double *model_unaligned = model + 1;

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model_unaligned);

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask, num_model, model_unaligned, out,
			out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsWithLsqMatrixNullPointer
 * failure case : lsq_matrix is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsWithLsqMatrixNullPointer) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	double *out = nullptr;
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask, num_model, model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsWithLsqMatrixNotAligned
 * failure case : lsq_matrix is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsWithLsqMatrixNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model + 1];
	double *out_unaligned = out + 1;
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask, num_model, model, out_unaligned,
			out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsWithLsqVectorNullPointer
 * failure case : lsq_vector is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsWithLsqVectorNullPointer) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double *out_vector = nullptr;
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask, num_model, model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLeastSquareFittingCoefficientsWithLsqVectorNotAligned
 * failure case : lsq_vector is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsWithLsqVectorNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model + 1];
	double *out_vector_unaligned = out_vector + 1;
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask, num_model, model, out,
			out_vector_unaligned);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test GetLeastSquareFittingCoefficientsTooManyMaskedData:
 * failure case of too many masked data
 * returned value : Status_kNG
 */
TEST_F(NumericOperation, GetLeastSquareFittingCoefficientsTooManyMaskedData) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(false, num_data, in_mask);
	for (size_t i = 0; i < num_model / 2; ++i) {
		in_mask[i] = true;
	}
	SetChebyshevModel(num_data, num_model, model);

	if (verbose) {
		PrintArray("in_mask", num_data, in_mask);
		PrintArray("model  ", num_data, num_model, model);
	}

	LIBSAKURA_SYMBOL(Status) status = sakura_GetLeastSquareFittingCoefficients(
			num_data, in_data, in_mask, num_model, model, out, out_vector);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kNG), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficients
 * successful case
 * note : repeating NUM_REPEAT2 times for performance measurement
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficients) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];
	SIMD_ALIGN
	double answer[num_model * num_model];
	SIMD_ALIGN
	double answer_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	LIBSAKURA_SYMBOL(Status) status_getlsq =
			sakura_GetLeastSquareFittingCoefficients(num_data, in_data, in_mask,
					num_model, model, answer, answer_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_getlsq);

	if (verbose) {
		PrintArray("in_mask", num_data, in_mask);
		PrintArray("model  ", num_data, num_model, model);
	}

	size_t const num_repeat(NUM_REPEAT2);
	double elapse_time = 0.0;
	for (size_t i = 0; i < num_repeat; ++i) {
		SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
				in_lsq_vector);
		double start = sakura_GetCurrentTime();
		LIBSAKURA_SYMBOL(Status) status =
				sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
						num_exclude_indices, exclude_indices, num_model, model,
						in_lsq_matrix, in_lsq_vector);
		double end = sakura_GetCurrentTime();
		elapse_time += (end - start);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}

	for (size_t i = 0; i < num_model * num_model; ++i) {
		double deviation;
		if (answer[i] != 0.0) {
			deviation = fabs((in_lsq_matrix[i] - answer[i]) / answer[i]);
		} else if (in_lsq_matrix[i] != 0.0) {
			deviation = fabs((in_lsq_matrix[i] - answer[i]) / in_lsq_matrix[i]);
		} else {
			deviation = fabs(in_lsq_matrix[i] - answer[i]);
		}
		ASSERT_LE(deviation, 1e-7);
	}
	for (size_t i = 0; i < num_model; ++i) {
		double deviation;
		if (answer_vector[i] != 0.0) {
			deviation = fabs(
					(in_lsq_vector[i] - answer_vector[i]) / answer_vector[i]);
		} else if (in_lsq_vector[i] != 0.0) {
			deviation = fabs(
					(in_lsq_vector[i] - answer_vector[i]) / in_lsq_vector[i]);
		} else {
			deviation = fabs(in_lsq_vector[i] - answer_vector[i]);
		}
		ASSERT_LE(deviation, 1e-7);
	}

	if (verbose) {
		cout << "Elapse time of " << num_repeat << " repetition: "
				<< elapse_time << " sec." << endl;
		PrintArray("out   ", num_model, num_model, in_lsq_matrix);
		PrintArray("answer", num_model, num_model, answer);
	}
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsNumDataZero
 * failure case : num_data == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsNumDataZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithDataNullPointer
 * failure case : data is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithDataNullPointer) {
	size_t const num_data(NUM_DATA);
	float *data_np = nullptr;
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[num_data];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, data_np,
					num_exclude_indices, exclude_indices, num_model, model,
					in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithDataNotAligned
 * failure case : data is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithDataNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data + 1];
	float *in_data_unaligned = in_data + 1;
	SIMD_ALIGN
	bool in_mask[num_data];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_mask)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data_unaligned);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data_unaligned, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data,
					in_data_unaligned, num_exclude_indices, exclude_indices,
					num_model, model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithExcludeIndicesNullPointer
 * failure case : exclude_indices is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithExcludeIndicesNullPointer) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	size_t *exclude_indices_np = nullptr;
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices_np, num_model, model,
					in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithExcludeIndicesNotAligned
 * failure case : exclude_indices is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithExcludeIndicesNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE + 1];
	size_t *exclude_indices_unaligned = exclude_indices + 1;
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices_unaligned[i] = i * i;
		in_mask[exclude_indices_unaligned[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices_unaligned, num_model,
					model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithNumExcludeIndicesGreaterThanNumData
 * failure case : num_exclude_indices > num_data
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithNumExcludeIndicesGreaterThanNumData) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	size_t const num_exclude_indices_toolarge = num_data + 1;
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE + 1];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices_toolarge, exclude_indices, num_model,
					model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithExcludeIndicesHasValueEqualToNumData
 * failure case : exclude_indices[i] == num_data (where 0 <= i < num_exclude_indices)
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithExcludeIndicesHasValueEqualToNumData) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE + 1];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	exclude_indices[0] = num_data;
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithExcludeIndicesHasValueGreaterThanNumData
 * failure case : exclude_indices[i] > num_data (where 0 <= i < num_exclude_indices)
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithExcludeIndicesHasValueGreaterThanNumData) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE + 1];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	exclude_indices[0] = num_data + 1;
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithNumModelBasisZero
 * failure case : num_model == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithNumModelBasisZero) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE + 1];
	size_t const num_model(NUM_MODEL);
	size_t const num_model_zero = 0;
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model_zero, model,
					in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithNumModelBasisGreaterThanNumData
 * failure case : num_model > num_data
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithNumModelBasisGreaterThanNumData) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE + 1];
	size_t const num_model(NUM_MODEL);
	size_t const num_model_toolarge = num_data + 1;
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model_toolarge,
					model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithBasisDataNullPointer
 * failure case : model is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithBasisDataNullPointer) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	double *model_np = nullptr;
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model_np,
					in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithBasisDataNotAligned
 * failure case : model is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithBasisDataNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data) + 1];
	double *model_unaligned = model + 1;
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model_unaligned);
	SetAnswers(num_data, in_data, num_model, model_unaligned, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model,
					model_unaligned, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithLsqMatrixNullPointer
 * failure case : in_lsq_matrix is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithLsqMatrixNullPointer) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	double *in_lsq_matrix_np = nullptr;
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					in_lsq_matrix_np, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithLsqMatrixNotAligned
 * failure case : in_lsq_matrix is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithLsqMatrixNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model + 1];
	double *in_lsq_matrix_unaligned = in_lsq_matrix + 1;
	SIMD_ALIGN
	double in_lsq_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix_unaligned,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					in_lsq_matrix_unaligned, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithLsqVectorNullPointer
 * failure case : in_lsq_vector is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithLsqVectorNullPointer) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];
	double *in_lsq_vector_np = nullptr;

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					in_lsq_matrix, in_lsq_vector_np);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLeastSquareFittingCoefficientsWithLsqVectorNotAligned
 * failure case : in_lsq_vector is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLeastSquareFittingCoefficientsWithLsqVectorNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_matrix[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model + 1];
	double *in_lsq_vector_unaligned = in_lsq_vector + 1;

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, in_lsq_matrix,
			in_lsq_vector_unaligned);

	LIBSAKURA_SYMBOL(Status) status =
			sakura_UpdateLeastSquareFittingCoefficients(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					in_lsq_matrix, in_lsq_vector_unaligned);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

#define NUM_DATA_TESTLU 50
#define NUM_MODEL_TESTLU 5
/*
 * Test sakura_SolveSimultaneousEquationsByLU
 * successful case
 * note : repeating NUM_REPEAT3 times for performance measurement
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLU) {
	size_t const num_data(NUM_DATA_TESTLU);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL_TESTLU);
	SIMD_ALIGN
	double lsq_vector[num_model];
	SIMD_ALIGN
	double model[ELEMENTSOF(lsq_vector) * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double lsq_matrix[ELEMENTSOF(lsq_vector) * ELEMENTSOF(lsq_vector)];
	SIMD_ALIGN
	double out[ELEMENTSOF(lsq_vector)];
	SIMD_ALIGN
	float answer[ELEMENTSOF(lsq_vector)];

	SetFloatPolynomial(num_data, in_data);
	SetBoolConstant(true, ELEMENTSOF(in_data), in_mask);
	SetPolynomialModel(ELEMENTSOF(in_data), ELEMENTSOF(lsq_vector), model);
	SetFloatConstant(1.0, ELEMENTSOF(lsq_vector), answer);

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
		PrintArray("model", num_data, num_model, model);
	}

	LIBSAKURA_SYMBOL(Status) coeff_status =
			sakura_GetLeastSquareFittingCoefficients(num_data, in_data, in_mask,
					num_model, model, lsq_matrix, lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);

	if (verbose) {
		PrintArray("lsq_matrix", num_model, num_model, lsq_matrix);
		PrintArray("lsq_vector", num_model, lsq_vector);
	}

	size_t const num_repeat = NUM_REPEAT3;
	for (size_t i = 0; i < num_repeat; ++i) {
		LIBSAKURA_SYMBOL(Status) solve_status =
				sakura_SolveSimultaneousEquationsByLU(num_model, lsq_matrix,
						lsq_vector, out);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), solve_status);
	}

	for (size_t i = 0; i < num_model; ++i) {
		float deviation = (out[i] - answer[i]) / answer[i];
		ASSERT_LE(deviation, 1e-7);
	}

	if (verbose) {
		PrintArray("out   ", num_model, out);
		PrintArray("answer", num_model, answer);
	}
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUWithInMatrixNullPointer
 * failure case : in_matrix is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUWithInMatrixNullPointer) {
	size_t const num_data(NUM_DATA_TESTLU);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL_TESTLU);
	SIMD_ALIGN
	double lsq_vector[num_model];
	SIMD_ALIGN
	double model[ELEMENTSOF(lsq_vector) * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double lsq_matrix[ELEMENTSOF(lsq_vector) * ELEMENTSOF(lsq_vector)];
	double *in_matrix_np = nullptr;
	SIMD_ALIGN
	double out[ELEMENTSOF(lsq_vector)];

	SetFloatPolynomial(num_data, in_data);
	SetBoolConstant(true, ELEMENTSOF(in_data), in_mask);
	SetPolynomialModel(ELEMENTSOF(in_data), ELEMENTSOF(lsq_vector), model);
	SetAnswers(num_data, in_data, in_mask, num_model, model, lsq_matrix,
			lsq_vector);

	LIBSAKURA_SYMBOL(Status) solve_status =
			sakura_SolveSimultaneousEquationsByLU(num_model, in_matrix_np,
					lsq_vector, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUWithInMatrixNotAligned
 * failure case : in_matrix is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUWithInMatrixNotAligned) {
	size_t const num_data(NUM_DATA_TESTLU);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL_TESTLU);
	SIMD_ALIGN
	double lsq_vector[num_model];
	SIMD_ALIGN
	double model[ELEMENTSOF(lsq_vector) * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double lsq_matrix[ELEMENTSOF(lsq_vector) * ELEMENTSOF(lsq_vector) + 1];
	double *lsq_matrix_unaligned = lsq_matrix + 1;
	SIMD_ALIGN
	double out[ELEMENTSOF(lsq_vector)];

	SetFloatPolynomial(num_data, in_data);
	SetBoolConstant(true, ELEMENTSOF(in_data), in_mask);
	SetPolynomialModel(ELEMENTSOF(in_data), ELEMENTSOF(lsq_vector), model);
	SetAnswers(num_data, in_data, in_mask, num_model, model,
			lsq_matrix_unaligned, lsq_vector);

	LIBSAKURA_SYMBOL(Status) solve_status =
			sakura_SolveSimultaneousEquationsByLU(num_model,
					lsq_matrix_unaligned, lsq_vector, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUWithInVectorNullPointer
 * failure case : in_vector is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUWithInVectorNullPointer) {
	size_t const num_data(NUM_DATA_TESTLU);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL_TESTLU);
	SIMD_ALIGN
	double lsq_vector[num_model];
	double *in_vector_np = nullptr;
	SIMD_ALIGN
	double model[ELEMENTSOF(lsq_vector) * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double lsq_matrix[ELEMENTSOF(lsq_vector) * ELEMENTSOF(lsq_vector)];
	SIMD_ALIGN
	double out[ELEMENTSOF(lsq_vector)];

	SetFloatPolynomial(num_data, in_data);
	SetBoolConstant(true, ELEMENTSOF(in_data), in_mask);
	SetPolynomialModel(ELEMENTSOF(in_data), ELEMENTSOF(lsq_vector), model);
	SetAnswers(num_data, in_data, in_mask, num_model, model, lsq_matrix,
			lsq_vector);

	LIBSAKURA_SYMBOL(Status) solve_status =
			sakura_SolveSimultaneousEquationsByLU(num_model, lsq_matrix,
					in_vector_np, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUWithInVectorNotAligned
 * failure case : in_vector is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUWithInVectorNotAligned) {
	size_t const num_data(NUM_DATA_TESTLU);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL_TESTLU);
	SIMD_ALIGN
	double lsq_vector[num_model + 1];
	double *in_vector_unaligned = lsq_vector + 1;
	SIMD_ALIGN
	double model[ELEMENTSOF(lsq_vector) * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double lsq_matrix[ELEMENTSOF(lsq_vector) * ELEMENTSOF(lsq_vector)];
	SIMD_ALIGN
	double out[ELEMENTSOF(lsq_vector)];

	SetFloatPolynomial(num_data, in_data);
	SetBoolConstant(true, ELEMENTSOF(in_data), in_mask);
	SetPolynomialModel(ELEMENTSOF(in_data), ELEMENTSOF(lsq_vector), model);
	SetAnswers(num_data, in_data, in_mask, num_model, model, lsq_matrix,
			in_vector_unaligned);

	LIBSAKURA_SYMBOL(Status) solve_status =
			sakura_SolveSimultaneousEquationsByLU(num_model, lsq_matrix,
					in_vector_unaligned, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUWithOutNullPointer
 * failure case : out is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUWithOutNullPointer) {
	size_t const num_data(NUM_DATA_TESTLU);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL_TESTLU);
	SIMD_ALIGN
	double lsq_vector[num_model];
	SIMD_ALIGN
	double model[ELEMENTSOF(lsq_vector) * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double lsq_matrix[ELEMENTSOF(lsq_vector) * ELEMENTSOF(lsq_vector)];
	double *out = nullptr;

	SetFloatPolynomial(num_data, in_data);
	SetBoolConstant(true, ELEMENTSOF(in_data), in_mask);
	SetPolynomialModel(ELEMENTSOF(in_data), ELEMENTSOF(lsq_vector), model);
	SetAnswers(num_data, in_data, in_mask, num_model, model, lsq_matrix,
			lsq_vector);

	LIBSAKURA_SYMBOL(Status) solve_status =
			sakura_SolveSimultaneousEquationsByLU(num_model, lsq_matrix,
					lsq_vector, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUWithOutNotAligned
 * failure case : out is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUWithOutNotAligned) {
	size_t const num_data(NUM_DATA_TESTLU);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL_TESTLU);
	SIMD_ALIGN
	double lsq_vector[num_model];
	SIMD_ALIGN
	double model[ELEMENTSOF(lsq_vector) * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double lsq_matrix[ELEMENTSOF(lsq_vector) * ELEMENTSOF(lsq_vector)];
	SIMD_ALIGN
	double out[ELEMENTSOF(lsq_vector) + 1];
	double *out_unaligned = out + 1;

	SetFloatPolynomial(num_data, in_data);
	SetBoolConstant(true, ELEMENTSOF(in_data), in_mask);
	SetPolynomialModel(ELEMENTSOF(in_data), ELEMENTSOF(lsq_vector), model);
	SetAnswers(num_data, in_data, in_mask, num_model, model, lsq_matrix,
			lsq_vector);

	LIBSAKURA_SYMBOL(Status) solve_status =
			sakura_SolveSimultaneousEquationsByLU(num_model, lsq_matrix,
					lsq_vector, out_unaligned);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}
