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
 * numeric_operation.cc
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

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"

/* the number of elements in input/output array to test */
#define NUM_DATA0 65536
#define NUM_DATA 4096
#define NUM_DATA2 50
#define NUM_DATA3 500
#define NUM_MODEL 20
#define NUM_MODEL2 5
#define NUM_MODEL3 499
#define NUM_REPEAT0 200
#define NUM_REPEAT 3000
#define NUM_REPEAT2 3000000
#define NUM_REPEAT3 15
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
		LIBSAKURA_SYMBOL (Status) status = LIBSAKURA_SYMBOL(Initialize)(nullptr,
				nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}

	virtual void TearDown() {
		LIBSAKURA_SYMBOL (CleanUp)();
	}

	// Set (1+x+x*x+x*x*x+x*x*x*x) float values into an array
	void SetFloatPolynomial(size_t const num_data, float *data) {
		for (size_t i = 0; i < num_data; ++i) {
			double x = (double) i;
			data[i] = (float) (1.0 + x + x * x + x * x * x + x * x * x * x);
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
			size_t const num_model, double const *model, size_t const num_coeff,
			double *answer_matrix, double *answer_vector) {
		size_t idx = 0;
		for (size_t i = 0; i < num_coeff; ++i) {
			for (size_t j = 0; j < num_coeff; ++j) {
				double val = 0.0;
				for (size_t k = 0; k < num_data; ++k) {
					val += model[num_model * k + i] * model[num_model * k + j];
				}
				answer_matrix[idx] = val;
				idx++;
			}
		}
		for (size_t i = 0; i < num_coeff; ++i) {
			double val = 0.0;
			for (size_t j = 0; j < num_data; ++j) {
				val += model[num_model * j + i] * data[j];
			}
			answer_vector[i] = val;
		}
	}

	// Set reference values for testing Update
	void SetAnswers(size_t const num_data, float const *data, bool const *mask,
			size_t const num_model, double const *model, size_t const num_coeff,
			double *answer_matrix, double *answer_vector) {
		size_t idx = 0;
		for (size_t i = 0; i < num_coeff; ++i) {
			for (size_t j = 0; j < num_coeff; ++j) {
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
		for (size_t i = 0; i < num_coeff; ++i) {
			double val = 0.0;
			for (size_t j = 0; j < num_data; ++j) {
				if (mask[j]) {
					val += model[num_model * j + i] * data[j];
				}
			}
			answer_vector[i] = val;
		}
	}

	// Set reference values for Cubic Spline
	void SetAnswersCubicSpline(size_t const num_data, float const *data,
	bool const *mask, size_t const num_boundary, double const *boundary,
			double *answer_matrix, double *answer_vector) {
		size_t num_vector = 3 + num_boundary;
		double *basis_data = nullptr;
		unique_ptr<void, DefaultAlignedMemory> storage_for_basisdata(
				DefaultAlignedMemory::AlignedAllocateOrException(
						sizeof(*basis_data) * num_data * 4, &basis_data));
		if (basis_data == nullptr) {
			throw bad_alloc();
		}
		for (size_t i = 0; i < num_data; ++i) {
			double x = (double) i;
			basis_data[i] = 1.0;
			basis_data[i + num_data] = x;
			basis_data[i + num_data * 2] = x * x;
			basis_data[i + num_data * 3] = x * x * x;
		}
		double *aux_data = nullptr;
		unique_ptr<void, DefaultAlignedMemory> storage_for_auxdata(
				DefaultAlignedMemory::AlignedAllocateOrException(
						sizeof(*aux_data) * num_data * num_boundary,
						&aux_data));
		if (aux_data == nullptr) {
			throw bad_alloc();
		}
		for (size_t i = 0; i < num_data; ++i) {
			for (size_t j = 0; j < num_boundary; ++j) {
				double x = (double) i - boundary[j];
				aux_data[i + num_data * j] = max(x * x * x, 0.0);
			}
		}

		for (size_t i = 0; i < num_vector; ++i) {
			// computing vector components
			double vec = 0.0;
			if (i < 3) {
				double *basis_data_i = &basis_data[num_data * i];
				for (size_t j = 0; j < num_data; ++j) {
					if (mask[j]) {
						vec += basis_data_i[j] * data[j];
					}
				}
			} else if (i == num_vector - 1) {
				size_t bidx = i - 3;
				double *aux_data_bidx = &aux_data[num_data * bidx];
				for (size_t j = 0; j < num_data; ++j) {
					if (mask[j]) {
						vec += aux_data_bidx[j] * data[j];
					}
				}
			} else if (num_boundary > 1) {
				size_t bidx = i - 3;
				double *aux_data_bidx = &aux_data[num_data * bidx];
				double *aux_data_bidxr = &aux_data[num_data * (bidx + 1)];
				for (size_t j = 0; j < num_data; ++j) {
					if (mask[j]) {
						vec += (aux_data_bidx[j] - aux_data_bidxr[j]) * data[j];
					}
				}
			}
			answer_vector[i] = vec;
			// computing matrix components
			double *answer_matrix_i = &answer_matrix[num_vector * i];
			for (size_t icol = 0; icol < num_vector; ++icol) {
				double mtx = 0.0;
				if (i < 3) {
					double *basis_data_i = &basis_data[num_data * i];
					double *basis_data_icol = &basis_data[num_data * icol];
					if (icol < 3) {
						for (size_t j = 0; j < num_data; ++j) {
							if (mask[j])
								mtx += basis_data_i[j] * basis_data_icol[j];
						}
					} else if (icol == num_vector - 1) {
						size_t bidx2 = icol - 3;
						double *aux_data_bidx2 = &aux_data[num_data * bidx2];
						for (size_t j = 0; j < num_data; ++j) {
							if (mask[j])
								mtx += basis_data_i[j] * aux_data_bidx2[j];
						}
					} else if (num_boundary > 1) {
						size_t bidx2 = icol - 3;
						double *aux_data_bidx2 = &aux_data[num_data * bidx2];
						double *aux_data_bidx2r = &aux_data[num_data * (bidx2 + 1)];
						for (size_t j = 0; j < num_data; ++j) {
							if (mask[j])
								mtx += basis_data_i[j] * (aux_data_bidx2[j] - aux_data_bidx2r[j]);
						}
					}
				} else if (i == num_vector - 1) {
					size_t bidx = i - 3;
					double *aux_data_bidx = &aux_data[num_data * bidx];
					double *basis_data_icol = &basis_data[num_data * icol];
					if (icol < 3) {
						for (size_t j = 0; j < num_data; ++j) {
							if (mask[j])
								mtx += aux_data_bidx[j] * basis_data_icol[j];
						}
					} else if (icol == num_vector - 1) {
						size_t bidx2 = icol - 3;
						double *aux_data_bidx2 = &aux_data[num_data * bidx2];
						for (size_t j = 0; j < num_data; ++j) {
							if (mask[j])
								mtx += aux_data_bidx[j] * aux_data_bidx2[j];
						}
					} else if (num_boundary > 1) {
						size_t bidx2 = icol - 3;
						double *aux_data_bidx2 = &aux_data[num_data * bidx2];
						double *aux_data_bidx2r = &aux_data[num_data * (bidx2 + 1)];
						for (size_t j = 0; j < num_data; ++j) {
							if (mask[j])
								mtx += aux_data_bidx[j] * (aux_data_bidx2[j] - aux_data_bidx2r[j]);
						}
					}
				} else if (num_boundary > 1) {
					size_t bidx = i - 3;
					double *aux_data_bidx = &aux_data[num_data * bidx];
					double *aux_data_bidxr = &aux_data[num_data * (bidx + 1)];
					double *basis_data_icol = &basis_data[num_data * icol];
					if (icol < 3) {
						for (size_t j = 0; j < num_data; ++j) {
							if (mask[j])
								mtx += (aux_data_bidx[j] - aux_data_bidxr[j])
										* basis_data_icol[j];
						}
					} else if (icol == num_vector - 1) {
						size_t bidx2 = icol - 3;
						double *aux_data_bidx2 = &aux_data[num_data * bidx2];
						for (size_t j = 0; j < num_data; ++j) {
							if (mask[j])
								mtx += (aux_data_bidx[j] - aux_data_bidxr[j])
										* aux_data_bidx2[j];
						}
					} else if (num_boundary > 1) {
						size_t bidx2 = icol - 3;
						double *aux_data_bidx2 = &aux_data[num_data * bidx2];
						double *aux_data_bidx2r = &aux_data[num_data * (bidx2 + 1)];
						for (size_t j = 0; j < num_data; ++j) {
							if (mask[j])
								mtx += (aux_data_bidx[j] - aux_data_bidxr[j])
										* (aux_data_bidx2[j] - aux_data_bidx2r[j]);
						}
					}
				}
				answer_matrix_i[icol] = mtx;
			}
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
 * Test sakura_GetLSQFittingCoefficients
 * successful case
 * note: repeating NUM_REPEAT times for performance measurement
 */
TEST_F(NumericOperation, GetLSQFittingCoefficients) {
	size_t const num_data(NUM_DATA0);
	size_t const num_model(NUM_MODEL);
	size_t const num_model_num_data(num_model * num_data);
	size_t const num_model_num_model(num_model * num_model);

	float *in_data = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_indata(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*in_data) * num_data, &in_data));
	if (in_data == nullptr) {
		throw bad_alloc();
	}

	bool *in_mask = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_inmask(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*in_mask) * num_data, &in_mask));
	if (in_mask == nullptr) {
		throw bad_alloc();
	}

	double *out = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_out(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*out) * num_model_num_model, &out));
	if (out == nullptr) {
		throw bad_alloc();
	}

	double *out_vector = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_outvector(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*out_vector) * num_model, &out_vector));
	if (out_vector == nullptr) {
		throw bad_alloc();
	}

	double *model = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_model(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*model) * num_model_num_data, &model));
	if (model == nullptr) {
		throw bad_alloc();
	}

	double *answer = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_answer(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*answer) * num_model_num_model, &answer));
	if (answer == nullptr) {
		throw bad_alloc();
	}

	double *answer_vector = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_answervector(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*answer_vector) * num_model, &answer_vector));
	if (answer_vector == nullptr) {
		throw bad_alloc();
	}

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, num_model, answer,
			answer_vector);

	if (verbose) {
		PrintArray("in_mask", num_data, in_mask);
		PrintArray("model  ", num_data, num_model, model);
	}

	size_t const num_repeat(NUM_REPEAT0);
	double start, end;
	double elapsed_time = 0.0;
	for (size_t i = 0; i < num_repeat; ++i) {
		start = sakura_GetCurrentTime();
		LIBSAKURA_SYMBOL (Status) status =
				sakura_GetLSQFittingCoefficientsDouble(num_data,
						in_data, in_mask, num_model, model, num_model, out,
						out_vector);
		end = sakura_GetCurrentTime();
		elapsed_time += (end - start);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}
	cout << "Elapsed Time: " << elapsed_time << " sec." << endl;

	for (size_t i = 0; i < num_model * num_model; ++i) {
		CheckAlmostEqual(answer[i], out[i], 1e-10);
	}
	for (size_t i = 0; i < num_model; ++i) {
		CheckAlmostEqual(answer_vector[i], out_vector[i], 1e-10);
	}

	if (verbose) {
		PrintArray("out   ", num_model, num_model, out);
		PrintArray("answer", num_model, num_model, answer);
	}
}

/*
 * Test sakura_GetLSQFittingCoefficientsNumDataZero
 * failure case : num_data == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsNumDataZero) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQFittingCoefficientsWithDataNullPointer
 * failure case : data is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsWithDataNullPointer) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQFittingCoefficientsWithDataNotAligned
 * failure case : data is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsWithDataNotAligned) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data,
					in_data_unaligned, in_mask, num_model, model, num_model,
					out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQFittingCoefficientsWithMaskNullPointer
 * failure case : mask is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsWithMaskNullPointer) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQFittingCoefficientsWithMaskNotAligned
 * failure case : mask is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsWithMaskNotAligned) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask_unaligned, num_model, model, num_model, out,
					out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQFittingCoefficientsWithNumModelBasesZero
 * failure case : num_model_bases == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsWithNumModelBasesZero) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQFittingCoefficientsWithNumDataLessThanNumModelBases
 * failure case : num_model_bases > num_data
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsWithNumDataLessThanNumModelBases) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQFittingCoefficientsWithBasisDataNullPointer
 * failure case : basis_data is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsWithBasisDataNullPointer) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQFittingCoefficientsWithBasisDataNotAligned
 * failure case : basis_data is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsWithBasisDataNotAligned) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model_unaligned, num_model, out,
					out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQFittingCoefficientsWithLsqMatrixNullPointer
 * failure case : lsq_matrix is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsWithLsqMatrixNullPointer) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQFittingCoefficientsWithLsqMatrixNotAligned
 * failure case : lsq_matrix is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsWithLsqMatrixNotAligned) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, out_unaligned,
					out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQFittingCoefficientsWithLsqVectorNullPointer
 * failure case : lsq_vector is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsWithLsqVectorNullPointer) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQFittingCoefficientsWithLsqVectorNotAligned
 * failure case : lsq_vector is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsWithLsqVectorNotAligned) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, out,
					out_vector_unaligned);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test GetLSQFittingCoefficientsTooManyMaskedData:
 * failure case of too many masked data
 * returned value : Status_kNG
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsTooManyMaskedData) {
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

	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, out, out_vector);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kNG), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficients
 * successful case
 * note : repeating NUM_REPEAT2 times for performance measurement
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficients) {
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
	double in_lsq_matrix_orig[num_model * num_model];
	SIMD_ALIGN
	double in_lsq_vector[num_model];
	SIMD_ALIGN
	double in_lsq_vector_orig[num_model];
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
	LIBSAKURA_SYMBOL (Status) status_getlsq =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, answer,
					answer_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_getlsq);

	if (verbose) {
		PrintArray("in_mask", num_data, in_mask);
		PrintArray("model  ", num_data, num_model, model);
	}

	SetAnswers(num_data, in_data, num_model, model, num_model,
			in_lsq_matrix_orig, in_lsq_vector_orig);

	size_t const num_repeat(NUM_REPEAT2);
	double elapsed_time = 0.0;
	for (size_t i = 0; i < num_repeat; ++i) {
		for (size_t j = 0; j < ELEMENTSOF(in_lsq_matrix); ++j) {
			in_lsq_matrix[j] = in_lsq_matrix_orig[j];
		}
		for (size_t j = 0; j < ELEMENTSOF(in_lsq_vector); ++j) {
			in_lsq_vector[j] = in_lsq_vector_orig[j];
		}
		double start = sakura_GetCurrentTime();
		LIBSAKURA_SYMBOL (Status) status =
				sakura_UpdateLSQFittingCoefficientsDouble(num_data,
						in_data, num_exclude_indices, exclude_indices,
						num_model, model, num_model, in_lsq_matrix,
						in_lsq_vector);
		double end = sakura_GetCurrentTime();
		elapsed_time += (end - start);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}
	cout << "Elapsed Time: " << elapsed_time << " sec." << endl;

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
		PrintArray("out   ", num_model, num_model, in_lsq_matrix);
		PrintArray("answer", num_model, num_model, answer);
	}
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsNumDataZero
 * failure case : num_data == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsNumDataZero) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					num_model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithDataNullPointer
 * failure case : data is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithDataNullPointer) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, data_np,
					num_exclude_indices, exclude_indices, num_model, model,
					num_model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithDataNotAligned
 * failure case : data is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithDataNotAligned) {
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
	SetAnswers(num_data, in_data_unaligned, num_model, model, num_model,
			in_lsq_matrix, in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data,
					in_data_unaligned, num_exclude_indices, exclude_indices,
					num_model, model, num_model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithExcludeIndicesNullPointer
 * failure case : exclude_indices is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithExcludeIndicesNullPointer) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices_np, num_model, model,
					num_model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithExcludeIndicesNotAligned
 * failure case : exclude_indices is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithExcludeIndicesNotAligned) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices_unaligned, num_model,
					model, num_model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithNumExcludeIndicesGreaterThanNumData
 * failure case : num_exclude_indices > num_data
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithNumExcludeIndicesGreaterThanNumData) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices_toolarge, exclude_indices, num_model,
					model, num_model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithExcludeIndicesHasValueEqualToNumData
 * failure case : exclude_indices[i] == num_data (where 0 <= i < num_exclude_indices)
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithExcludeIndicesHasValueEqualToNumData) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					num_model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithExcludeIndicesHasValueGreaterThanNumData
 * failure case : exclude_indices[i] > num_data (where 0 <= i < num_exclude_indices)
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithExcludeIndicesHasValueGreaterThanNumData) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					num_model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithNumModelBasisZero
 * failure case : num_model == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithNumModelBasisZero) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model_zero, model,
					num_model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithNumModelBasisGreaterThanNumData
 * failure case : num_model > num_data
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithNumModelBasisGreaterThanNumData) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model_toolarge,
					model, num_model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithBasisDataNullPointer
 * failure case : model is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithBasisDataNullPointer) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model_np,
					num_model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithBasisDataNotAligned
 * failure case : model is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithBasisDataNotAligned) {
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
	SetAnswers(num_data, in_data, num_model, model_unaligned, num_model,
			in_lsq_matrix, in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model,
					model_unaligned, num_model, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithLsqMatrixNullPointer
 * failure case : in_lsq_matrix is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithLsqMatrixNullPointer) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					num_model, in_lsq_matrix_np, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithLsqMatrixNotAligned
 * failure case : in_lsq_matrix is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithLsqMatrixNotAligned) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model,
			in_lsq_matrix_unaligned, in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					num_model, in_lsq_matrix_unaligned, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithLsqVectorNullPointer
 * failure case : in_lsq_vector is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithLsqVectorNullPointer) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					num_model, in_lsq_matrix, in_lsq_vector_np);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQFittingCoefficientsWithLsqVectorNotAligned
 * failure case : in_lsq_vector is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsWithLsqVectorNotAligned) {
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
	SetAnswers(num_data, in_data, num_model, model, num_model, in_lsq_matrix,
			in_lsq_vector_unaligned);

	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					num_model, in_lsq_matrix, in_lsq_vector_unaligned);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_SolveSimultaneousEquationsByLU
 * successful case
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLU) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL2);
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

	LIBSAKURA_SYMBOL (Status) coeff_status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, lsq_matrix,
					lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);

	if (verbose) {
		PrintArray("lsq_matrix", num_model, num_model, lsq_matrix);
		PrintArray("lsq_vector", num_model, lsq_vector);
	}

	LIBSAKURA_SYMBOL (Status) solve_status =
			sakura_SolveSimultaneousEquationsByLUDouble(num_model, lsq_matrix,
					lsq_vector, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), solve_status);

	for (size_t i = 0; i < num_model; ++i) {
		float deviation = (out[i] - answer[i]) / out[i];
		ASSERT_LE(deviation, 1e-7);
	}

	if (verbose) {
		PrintArray("out   ", num_model, out);
		PrintArray("answer", num_model, answer);
	}
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUBigOrderModel
 * successful case
 * note : repeating NUM_REPEAT3 times for performance measurement
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUBigOrderModel) {
	size_t const num_data(NUM_DATA3);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL3);
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

	LIBSAKURA_SYMBOL (Status) coeff_status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_model, lsq_matrix,
					lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);

	if (verbose) {
		PrintArray("lsq_matrix", num_model, num_model, lsq_matrix);
		PrintArray("lsq_vector", num_model, lsq_vector);
	}

	double start, end;
	double elapsed_time = 0.0;
	size_t const num_repeat = NUM_REPEAT3;
	for (size_t i = 0; i < num_repeat; ++i) {
		start = LIBSAKURA_SYMBOL(GetCurrentTime)();
		LIBSAKURA_SYMBOL (Status) solve_status =
				sakura_SolveSimultaneousEquationsByLUDouble(num_model,
						lsq_matrix, lsq_vector, out);
		end = LIBSAKURA_SYMBOL(GetCurrentTime)();
		elapsed_time += (end - start);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), solve_status);
	}
	cout << "Elapsed Time: " << elapsed_time << " sec." << endl;

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
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL2);
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
	SetAnswers(num_data, in_data, in_mask, num_model, model, num_model,
			lsq_matrix, lsq_vector);

	LIBSAKURA_SYMBOL (Status) solve_status =
			sakura_SolveSimultaneousEquationsByLUDouble(num_model, in_matrix_np,
					lsq_vector, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUWithInMatrixNotAligned
 * failure case : in_matrix is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUWithInMatrixNotAligned) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL2);
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
	SetAnswers(num_data, in_data, in_mask, num_model, model, num_model,
			lsq_matrix_unaligned, lsq_vector);

	LIBSAKURA_SYMBOL (Status) solve_status =
			sakura_SolveSimultaneousEquationsByLUDouble(num_model,
					lsq_matrix_unaligned, lsq_vector, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUWithInVectorNullPointer
 * failure case : in_vector is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUWithInVectorNullPointer) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL2);
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
	SetAnswers(num_data, in_data, in_mask, num_model, model, num_model,
			lsq_matrix, lsq_vector);

	LIBSAKURA_SYMBOL (Status) solve_status =
			sakura_SolveSimultaneousEquationsByLUDouble(num_model, lsq_matrix,
					in_vector_np, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUWithInVectorNotAligned
 * failure case : in_vector is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUWithInVectorNotAligned) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL2);
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
	SetAnswers(num_data, in_data, in_mask, num_model, model, num_model,
			lsq_matrix, in_vector_unaligned);

	LIBSAKURA_SYMBOL (Status) solve_status =
			sakura_SolveSimultaneousEquationsByLUDouble(num_model, lsq_matrix,
					in_vector_unaligned, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUWithOutNullPointer
 * failure case : out is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUWithOutNullPointer) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL2);
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
	SetAnswers(num_data, in_data, in_mask, num_model, model, num_model,
			lsq_matrix, lsq_vector);

	LIBSAKURA_SYMBOL (Status) solve_status =
			sakura_SolveSimultaneousEquationsByLUDouble(num_model, lsq_matrix,
					lsq_vector, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUWithOutNotAligned
 * failure case : out is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUWithOutNotAligned) {
	size_t const num_data(NUM_DATA2);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL2);
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
	SetAnswers(num_data, in_data, in_mask, num_model, model, num_model,
			lsq_matrix, lsq_vector);

	LIBSAKURA_SYMBOL (Status) solve_status =
			sakura_SolveSimultaneousEquationsByLUDouble(num_model, lsq_matrix,
					lsq_vector, out_unaligned);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_GetLSQFittingCoefficients
 * successful case with num_lsq_bases < num_model
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsOrder) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2 + 1);
	size_t const num_coeff(NUM_MODEL2);

	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double out_vector[num_coeff];
	SIMD_ALIGN
	double out[ELEMENTSOF(out_vector) * ELEMENTSOF(out_vector)];
	SIMD_ALIGN
	double answer_vector[ELEMENTSOF(out_vector)];
	SIMD_ALIGN
	double answer[ELEMENTSOF(out)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);
	SetAnswers(num_data, in_data, num_model, model, num_coeff, answer,
			answer_vector);

	if (verbose) {
		PrintArray("in_mask", num_data, in_mask);
		PrintArray("model  ", num_data, num_model, model);
	}

	cout << "number of fitting coefficients = " << num_coeff << " (model = "
			<< num_model << ")" << endl;
	LIBSAKURA_SYMBOL (Status) status =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_coeff, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);

	for (size_t i = 0; i < ELEMENTSOF(out); ++i) {
		CheckAlmostEqual(answer[i], out[i], 1e-10);
	}
	for (size_t i = 0; i < ELEMENTSOF(out_vector); ++i) {
		CheckAlmostEqual(answer_vector[i], out_vector[i], 1e-10);
	}

	if (verbose) {
		PrintArray("out   ", num_coeff, num_coeff, out);
		PrintArray("answer", num_coeff, num_coeff, answer);
	}
}

/*
 * Test sakura_GetLSQFittingCoefficients
 * failure case with
 * num_lsq_bases > num_model
 * num_lsq_bases = 0
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsBadOrder) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const bad_coeffs[] = { NUM_MODEL2 + 1, 0 };

	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	for (size_t i = 0; i < ELEMENTSOF(bad_coeffs); ++i) {
		size_t num_coeff = bad_coeffs[i];
		SIMD_ALIGN
		double out[num_coeff * num_coeff];
		SIMD_ALIGN
		double out_vector[num_coeff];
		cout << "number of fitting coefficients = " << num_coeff << " (model = "
				<< num_model << ")" << endl;
		LIBSAKURA_SYMBOL (Status) status =
				sakura_GetLSQFittingCoefficientsDouble(num_data,
						in_data, in_mask, num_model, model, num_coeff, out,
						out_vector);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
	}
}

/*
 * Test sakura_UpdateLSQFittingCoefficients
 * successful case with num_coeff < num_model
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsOrder) {
	size_t const num_data(NUM_DATA);
	size_t const num_model(NUM_MODEL + 1);
	size_t const num_coeff(NUM_MODEL);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_vector[num_coeff];
	SIMD_ALIGN
	double in_lsq_vector_orig[ELEMENTSOF(in_lsq_vector)];
	SIMD_ALIGN
	double in_lsq_matrix[ELEMENTSOF(in_lsq_vector) * ELEMENTSOF(in_lsq_vector)];
	SIMD_ALIGN
	double in_lsq_matrix_orig[ELEMENTSOF(in_lsq_matrix)];
	SIMD_ALIGN
	double answer[ELEMENTSOF(in_lsq_matrix)];
	SIMD_ALIGN
	double answer_vector[ELEMENTSOF(in_lsq_vector)];

	SetInputData(ELEMENTSOF(in_data), in_data);
	SetBoolConstant(true, ELEMENTSOF(in_data), in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	LIBSAKURA_SYMBOL (Status) status_getlsq =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, num_coeff, answer,
					answer_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_getlsq);

	if (verbose) {
		PrintArray("in_mask", ELEMENTSOF(in_mask), in_mask);
		PrintArray("model  ", num_data, num_model, model);
	}
	cout << "number of fitting coefficients = " << num_coeff << " (model = "
			<< num_model << ")" << endl;

	SetAnswers(num_data, in_data, num_model, model, num_coeff,
			in_lsq_matrix_orig, in_lsq_vector_orig);

	for (size_t j = 0; j < ELEMENTSOF(in_lsq_matrix); ++j) {
		in_lsq_matrix[j] = in_lsq_matrix_orig[j];
	}
	for (size_t j = 0; j < ELEMENTSOF(in_lsq_vector); ++j) {
		in_lsq_vector[j] = in_lsq_vector_orig[j];
	}
	LIBSAKURA_SYMBOL (Status) status =
			sakura_UpdateLSQFittingCoefficientsDouble(num_data, in_data,
					num_exclude_indices, exclude_indices, num_model, model,
					num_coeff, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < ELEMENTSOF(answer); ++i) {
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
	for (size_t i = 0; i < ELEMENTSOF(answer_vector); ++i) {
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
		PrintArray("out   ", num_coeff, num_coeff, in_lsq_matrix);
		PrintArray("answer", num_coeff, num_coeff, answer);
	}
}

/*
 * Test sakura_UpdateLSQFittingCoefficients
 * failure case with
 * num_coeff > num_model
 * num_coeff = 0
 */
TEST_F(NumericOperation, UpdateLSQFittingCoefficientsBadOrder) {
	size_t const num_data(NUM_DATA);
	size_t const num_model(NUM_MODEL);
	size_t const bad_coeffs[] = { num_model + 1, 0 };
	size_t const good_coeff(num_model);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double in_lsq_vector[good_coeff];
	SIMD_ALIGN
	double in_lsq_matrix[ELEMENTSOF(in_lsq_vector) * ELEMENTSOF(in_lsq_vector)];

	SetInputData(ELEMENTSOF(in_data), in_data);
	SetBoolConstant(true, ELEMENTSOF(in_data), in_mask);
	SetChebyshevModel(num_data, num_model, model);
	LIBSAKURA_SYMBOL (Status) status_getlsq =
			sakura_GetLSQFittingCoefficientsDouble(num_data, in_data,
					in_mask, num_model, model, good_coeff, in_lsq_matrix,
					in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_getlsq);

	if (verbose) {
		PrintArray("in_mask", ELEMENTSOF(in_mask), in_mask);
		PrintArray("model  ", num_data, num_model, model);
	}
	for (size_t i = 0; i < ELEMENTSOF(bad_coeffs); ++i) {
		size_t num_coeff(bad_coeffs[i]);
		cout << "number of fitting coefficients = " << num_coeff << " (model = "
				<< num_model << ")" << endl;

		LIBSAKURA_SYMBOL (Status) status =
				sakura_UpdateLSQFittingCoefficientsDouble(num_data,
						in_data, num_exclude_indices, exclude_indices,
						num_model, model, num_coeff, in_lsq_matrix,
						in_lsq_vector);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
	}
}

/*
 * Test GetLSQFittingCoefficientsCubicSpline
 * successful case
 * compute the best-fit baseline coefficients for cubic spline
 * for cases of combination of two conditions: one is num_boundary
 * variation (1,2,3) and the other is boundary position (just on
 * or between data positions)
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsCubicSplineSuccessfulCase) {
	for (size_t num_boundary = 1; num_boundary <= 3; ++num_boundary) {
		cout << "    Testing for num_boundary = " << num_boundary
				<< " cases: num_data = ";
		size_t ipos_max = 0;
		if (num_boundary > 1)
			ipos_max = 1;
		for (size_t ipos = 0; ipos <= ipos_max; ++ipos) {
			// -------------------------------------------
			// ipos==0 : boundary on data points
			// ipos==1 : boundary between data points
			// -------------------------------------------
			SIMD_ALIGN
			double coeff[4 * num_boundary];
			SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
			size_t const num_data(ELEMENTSOF(coeff));
			SIMD_ALIGN
			float in_data[num_data];
			SetFloatPolynomial(num_data, in_data, coeff);
			SIMD_ALIGN
			bool mask[ELEMENTSOF(in_data)];
			SetBoolConstant(true, ELEMENTSOF(in_data), mask);
			if (ipos == 0)
				cout << num_data;
			if ((ipos == 0) || (ipos < ipos_max))
				cout << ", ";
			SIMD_ALIGN
			double boundary[num_boundary];
			for (size_t i = 0; i < num_boundary; ++i) {
				boundary[i] = floor(
						(double) i * (double) ELEMENTSOF(in_data)
								/ (double) num_boundary);
				if ((ipos > 0) && (i > 0))
					boundary[i] += 0.5;
			}
			cout << ((ipos == 0) ? "boundary =  " : "") << "[";
			for (size_t i = 0; i < num_boundary; ++i) {
				cout << boundary[i] << ((i < num_boundary - 1) ? ", " : "");
			}
			cout << "]" << ((ipos < ipos_max) ? ", " : "");
			SIMD_ALIGN
			double basis_data[4 * ELEMENTSOF(in_data)];
			size_t num_vector(3 + num_boundary);
			SIMD_ALIGN
			double out_vector[num_vector];
			SIMD_ALIGN
			double out_matrix[num_vector * num_vector];
			double answer_vector[ELEMENTSOF(out_vector)];
			double answer_matrix[ELEMENTSOF(out_matrix)];
			SetAnswersCubicSpline(num_data, in_data, mask, num_boundary,
					boundary, answer_matrix, answer_vector);
			SetPolynomialModel(num_data, 4, basis_data);
			LIBSAKURA_SYMBOL (Status) coeff_status =
					LIBSAKURA_SYMBOL(GetLSQFittingCoefficientsCubicSplineDouble)(
							num_data, in_data, mask, num_boundary, boundary,
							basis_data, out_matrix, out_vector);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);

			for (size_t i = 0; i < ELEMENTSOF(answer_vector); ++i) {
				CheckAlmostEqual(answer_vector[i], out_vector[i], 1.0e-6);
			}
			for (size_t i = 0; i < ELEMENTSOF(answer_matrix); ++i) {
				CheckAlmostEqual(answer_matrix[i], out_matrix[i], 1.0e-6);
			}

			if (verbose) {
				PrintArray("data       ", num_data, in_data);
				PrintArray("vector:out ", ELEMENTSOF(answer_vector),
						out_vector);
				PrintArray("vector:ans ", ELEMENTSOF(answer_vector),
						answer_vector);
				PrintArray("matrix:out ", ELEMENTSOF(answer_matrix),
						out_matrix);
				PrintArray("matrix:ans ", ELEMENTSOF(answer_matrix),
						answer_matrix);
			}

		}
		cout << endl;
	}
}

/*
 * Test GetLSQFittingCoefficientsCubicSpline
 * time-consuming successful case for performance measurement
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsCubicSplinePerformanceTest) {
	size_t const num_repeat(300);
	size_t const num_boundary(2);
	SIMD_ALIGN
	double coeff[4 * num_boundary];
	for (size_t i = 0; i < ELEMENTSOF(coeff); i += 4) {
		coeff[i] = 1.0;
		coeff[i + 1] = 1e-4;
		coeff[i + 2] = 1e-8;
		coeff[i + 3] = 1e-12;
	}
	size_t const num_data(70000);
	SIMD_ALIGN
	float in_data[num_data];
	SetFloatPolynomial(num_data, in_data, coeff);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(in_data)];
	SetBoolConstant(true, ELEMENTSOF(in_data), mask);
	SIMD_ALIGN
	double boundary[num_boundary];
	for (size_t i = 0; i < num_boundary; ++i) {
		boundary[i] = floor(
				(double) i * (double) ELEMENTSOF(in_data)
						/ (double) num_boundary);
	}
	SIMD_ALIGN
	double basis_data[4 * ELEMENTSOF(in_data)];
	SetPolynomialModel(num_data, 4, basis_data);
	size_t num_vector(3 + num_boundary);
	SIMD_ALIGN
	double out_vector[num_vector];
	SIMD_ALIGN
	double out_matrix[num_vector * num_vector];
	LIBSAKURA_SYMBOL (Status) coeff_status;
	for (size_t i = 0; i < num_repeat; ++i) {
		coeff_status =
		LIBSAKURA_SYMBOL(GetLSQFittingCoefficientsCubicSplineDouble)(
				num_data, in_data, mask, num_boundary, boundary, basis_data,
				out_matrix, out_vector);
	}
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);
}

/*
 * Test GetLSQFittingCoefficientsCubicSpline
 * erroneous cases: null pointer cases
 * parameters to be tested include data, mask, boundary, basis_data, out_matrix and out_vector.
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsCubicSplineErroneousCasesNullPointer) {
	enum NPItems {
		NP_kData,
		NP_kMask,
		NP_kBoundary,
		NP_kBasisData,
		NP_kOutMatrix,
		NP_kOutVector,
		NP_kNumElems
	};
	vector<string> np_param_names = { "data", "mask", "boundary", "basis_data",
			"out_matrix", "out_vector" };
	cout << "    Testing for ";

	for (NPItems item = static_cast<NPItems>(0); item < NP_kNumElems; item =
			static_cast<NPItems>(item + 1)) {
		cout << np_param_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");
		size_t const num_boundary(2);
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
		double boundary[num_boundary];
		for (size_t i = 0; i < num_boundary; ++i) {
			boundary[i] = floor(
					(double) i * (double) ELEMENTSOF(in_data)
							/ (double) num_boundary);
		}
		SIMD_ALIGN
		double basis_data[4 * ELEMENTSOF(in_data)];
		SetPolynomialModel(num_data, 4, basis_data);
		size_t num_vector(3 + num_boundary);
		SIMD_ALIGN
		double out_vector[num_vector];
		SIMD_ALIGN
		double out_matrix[num_vector * num_vector];

		float *in_data_ptr = in_data;
		bool *mask_ptr = mask;
		double *boundary_ptr = boundary;
		double *basis_data_ptr = basis_data;
		double *out_matrix_ptr = out_matrix;
		double *out_vector_ptr = out_vector;

		switch (item) {
		case NP_kData:
			in_data_ptr = nullptr;
			break;
		case NP_kMask:
			mask_ptr = nullptr;
			break;
		case NP_kBoundary:
			boundary_ptr = nullptr;
			break;
		case NP_kBasisData:
			basis_data_ptr = nullptr;
			break;
		case NP_kOutMatrix:
			out_matrix_ptr = nullptr;
			break;
		case NP_kOutVector:
			out_vector_ptr = nullptr;
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL (Status) coeff_status =
		LIBSAKURA_SYMBOL(GetLSQFittingCoefficientsCubicSplineDouble)(
				num_data, in_data_ptr, mask_ptr, num_boundary, boundary_ptr,
				basis_data_ptr, out_matrix_ptr, out_vector_ptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);
	}
	cout << endl;
}

/*
 * Test GetLSQFittingCoefficientsCubicSpline
 * erroneous cases: unaligned cases
 * parameters to be tested include data, mask, boundary, basis_data, out_matrix and out_vector.
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsCubicSplineErroneousCasesUnaligned) {
	enum UAItems {
		UA_kData,
		UA_kMask,
		UA_kBoundary,
		UA_kBasisData,
		UA_kOutMatrix,
		UA_kOutVector,
		UA_kNumElems
	};
	vector<string> ua_param_names = { "data", "mask", "boundary", "basis_data",
			"out_matrix", "out_vector" };
	cout << "    Testing for ";

	for (UAItems item = static_cast<UAItems>(0); item < UA_kNumElems; item =
			static_cast<UAItems>(item + 1)) {
		cout << ua_param_names[item] << ((item < UA_kNumElems - 1) ? ", " : "");
		size_t const num_boundary(2);
		size_t const num_data(10);
		double coeff[4];
		SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
		SIMD_ALIGN
		float in_data[num_data+1];
		SetFloatPolynomial(num_data, in_data, coeff);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(in_data)];
		SetBoolConstant(true, ELEMENTSOF(in_data), mask);
		SIMD_ALIGN
		double boundary[num_boundary+1];
		for (size_t i = 0; i < num_boundary; ++i) {
			boundary[i] = floor(
					(double) i * (double) ELEMENTSOF(in_data)
							/ (double) num_boundary);
		}
		SIMD_ALIGN
		double basis_data[4 * ELEMENTSOF(in_data)+1];
		SetPolynomialModel(num_data, 4, basis_data);
		size_t num_vector(3 + num_boundary);
		SIMD_ALIGN
		double out_vector[num_vector+1];
		SIMD_ALIGN
		double out_matrix[num_vector * num_vector + 1];

		float *in_data_ptr = in_data;
		bool *mask_ptr = mask;
		double *boundary_ptr = boundary;
		double *basis_data_ptr = basis_data;
		double *out_matrix_ptr = out_matrix;
		double *out_vector_ptr = out_vector;

		switch (item) {
		case UA_kData:
			++in_data_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(in_data_ptr));
			break;
		case UA_kMask:
			++mask_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(mask_ptr));
			break;
		case UA_kBoundary:
			++boundary_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(boundary_ptr));
			break;
		case UA_kBasisData:
			++basis_data_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(basis_data_ptr));
			break;
		case UA_kOutMatrix:
			++out_matrix_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(out_matrix_ptr));
			break;
		case UA_kOutVector:
			++out_vector_ptr;
			assert(!LIBSAKURA_SYMBOL(IsAligned)(out_vector_ptr));
			break;
		default:
			assert(false);
		}

		LIBSAKURA_SYMBOL (Status) coeff_status =
		LIBSAKURA_SYMBOL(GetLSQFittingCoefficientsCubicSplineDouble)(
				num_data, in_data_ptr, mask_ptr, num_boundary, boundary_ptr,
				basis_data_ptr, out_matrix_ptr, out_vector_ptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);
	}
	cout << endl;
}

/*
 * Test GetLSQFittingCoefficientsCubicSpline
 * erroneous cases: bad parameter value cases as follows:
 *     (1) num_data == 0
 *     (2) num_boundary == 0
 */
TEST_F(NumericOperation, GetLSQFittingCoefficientsCubicSplineErroneousCasesBadParameterValue) {
	enum BVItems {
		BV_kNumDataZero,
		BV_kNumBoundaryZero,
		BV_kNumElems
	};
	vector<string> bv_param_names = { "(num_data == 0)", "(num_boundary == 0)" };
	cout << "    Testing for cases " << endl;

	for (BVItems item = static_cast<BVItems>(0); item < BV_kNumElems; item =
			static_cast<BVItems>(item + 1)) {
		cout << "        " << bv_param_names[item]
				<< ((item < BV_kNumElems - 1) ? ", " : "") << endl;
		size_t num_boundary = 1;
		size_t num_data = 10;
		switch (item) {
		case BV_kNumDataZero:
			num_data = 0;
			break;
		case BV_kNumBoundaryZero:
			num_boundary = 0;
			break;
		default:
			assert(false);
		}
		SIMD_ALIGN
		double coeff[4];
		SetDoubleConstant(1.0, ELEMENTSOF(coeff), coeff);
		SIMD_ALIGN
		float in_data[num_data];
		SetFloatPolynomial(num_data, in_data, coeff);
		SIMD_ALIGN
		bool mask[ELEMENTSOF(in_data)];
		SetBoolConstant(true, ELEMENTSOF(in_data), mask);
		SIMD_ALIGN
		double boundary[num_boundary];
		for (size_t i = 0; i < num_boundary; ++i) {
			boundary[i] = floor(
					(double) i * (double) ELEMENTSOF(in_data)
							/ (double) num_boundary);
		}
		SIMD_ALIGN
		double basis_data[4 * ELEMENTSOF(in_data)];
		SetPolynomialModel(num_data, 4, basis_data);
		size_t num_vector(3 + num_boundary);
		SIMD_ALIGN
		double out_vector[num_vector];
		SIMD_ALIGN
		double out_matrix[num_vector * num_vector];

		LIBSAKURA_SYMBOL (Status) coeff_status =
		LIBSAKURA_SYMBOL(GetLSQFittingCoefficientsCubicSplineDouble)(
				num_data, in_data, mask, num_boundary, boundary,
				basis_data, out_matrix, out_vector);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), coeff_status);

	}
}

