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
#include <random>
#include <string>
#include <sys/time.h>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"
#include "testutil.h"

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
		LIBSAKURA_SYMBOL (Status)
		status = LIBSAKURA_SYMBOL(Initialize)(nullptr, nullptr);
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

	//Set (A[0]+A[1]*x+A[2]*x*x+A[3]*x*x*x) float values into an array
	void SetFloatPolynomial(size_t num_data, float *data,
			double *coeff_answer) {
		for (size_t i = 0; i < num_data; ++i) {
			double x = (double) i;
			data[i] = (float) (coeff_answer[0] + coeff_answer[1] * x
					+ coeff_answer[2] * x * x + coeff_answer[3] * x * x * x);
		}
	}

	// Set constant float plus Gaussian noise values into an array
	void SetFloatConstantWithGaussianNoise(float value, float sigma,
			size_t const num_data, float *data) {
		// using fixed seed value to get identical datasets
		std::mt19937 mt(3148285944);
		//--------
		// using random seed value
		//std::random_device rd;
		//std::mt19937 mt(rd());
		//--------
		std::normal_distribution<> out(value, sigma);
		for (size_t i = 0; i < num_data; ++i) {
			data[i] = out(mt);
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
						double *aux_data_bidx2r = &aux_data[num_data
								* (bidx2 + 1)];
						for (size_t j = 0; j < num_data; ++j) {
							if (mask[j])
								mtx += basis_data_i[j]
										* (aux_data_bidx2[j]
												- aux_data_bidx2r[j]);
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
						double *aux_data_bidx2r = &aux_data[num_data
								* (bidx2 + 1)];
						for (size_t j = 0; j < num_data; ++j) {
							if (mask[j])
								mtx += aux_data_bidx[j]
										* (aux_data_bidx2[j]
												- aux_data_bidx2r[j]);
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
						double *aux_data_bidx2r = &aux_data[num_data
								* (bidx2 + 1)];
						for (size_t j = 0; j < num_data; ++j) {
							if (mask[j])
								mtx += (aux_data_bidx[j] - aux_data_bidxr[j])
										* (aux_data_bidx2[j]
												- aux_data_bidx2r[j]);
						}
					}
				}
				answer_matrix_i[icol] = mtx;
			}
		}
	}

	// Set Gaussian profile(s) + noise data
	void SetFloatGaussian(size_t num_data, float const sigma_noise,
			size_t const num_line, double const *ans_height,
			double const *ans_center, double const *ans_sigma, float *data) {
		SetFloatConstantWithGaussianNoise(0.0f, sigma_noise, num_data, data);
		// Add Gaussian profile to the input data
		for (size_t iline = 0; iline < num_line; ++iline) {
			for (size_t i = 0; i < num_data; ++i) {
				double value = ((double) i - ans_center[iline])
						/ ans_sigma[iline];
				data[i] += ans_height[iline] * exp(-0.5 * value * value);
			}
		}
	}

	// Check if the expected and actual values are enough close to each other
	void CheckAlmostEqual(double expected, double actual, double tolerance) {
		double deviation = fabs(actual - expected);
		double val = max(fabs(actual), fabs(expected)) * tolerance + tolerance;
		ASSERT_LE(deviation, val);
	}

	// Check if the actual values are in range [expected-error, expected+error]
	void CheckInRange(double expected, double actual, double error) {
		double const factor = 10.0; //to cope with modest accuracy in LM estimate
		ASSERT_GE(actual, expected - factor * error);
		ASSERT_LE(actual, expected + factor * error);
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
 * Test sakura_GetLSQCoefficients
 * successful case
 * note: repeating NUM_REPEAT times for performance measurement
 */
TEST_F(NumericOperation, GetLSQCoefficients) {
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

	size_t *use_idx = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_use_idx(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*use_idx) * num_model, &use_idx));
	if (use_idx == nullptr) {
		throw bad_alloc();
	}
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
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
		start = GetCurrentTime();
		LIBSAKURA_SYMBOL (Status)
		status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
				num_model, model, num_model, use_idx, out, out_vector);
		end = GetCurrentTime();
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
 * Test sakura_GetLSQCoefficientsNumDataZero
 * failure case : num_data == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsNumDataZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQCoefficientsWithDataNullPointer
 * failure case : data is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsWithDataNullPointer) {
	size_t const num_data(NUM_DATA);
	float *in_data = nullptr;
	SIMD_ALIGN
	bool in_mask[num_data];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_mask)];

	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQCoefficientsWithDataNotAligned
 * failure case : data is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsWithDataNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data + 1];
	float *in_data_unaligned = in_data + 1;
	SIMD_ALIGN
	bool in_mask[num_data];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_mask)];

	SetInputData(num_data, in_data_unaligned);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data_unaligned,
			in_mask, num_model, model, num_model, use_idx, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQCoefficientsWithMaskNullPointer
 * failure case : mask is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsWithMaskNullPointer) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	bool *in_mask = nullptr;
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQCoefficientsWithMaskNotAligned
 * failure case : mask is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsWithMaskNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data) + 1];
	bool *in_mask_unaligned = in_mask + 1;
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask_unaligned);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data,
			in_mask_unaligned, num_model, model, num_model, use_idx, out,
			out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQCoefficientsWithNumModelBasesZero
 * failure case : num_model_bases == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsWithNumModelBasesZero) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(0);
	SIMD_ALIGN
	size_t use_idx[num_model]; // not initialize since it has zero length
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQCoefficientsWithNumDataLessThanNumModelBases
 * failure case : num_model_bases > num_data
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsWithNumDataLessThanNumModelBases) {
	size_t const num_data(NUM_MODEL - 1);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQCoefficientsWithBasisDataNullPointer
 * failure case : basis_data is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsWithBasisDataNullPointer) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double out_vector[num_model];
	double *model = nullptr;

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQCoefficientsWithBasisDataNotAligned
 * failure case : basis_data is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsWithBasisDataNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model_unaligned, num_model, use_idx, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQCoefficientsWithLsqMatrixNullPointer
 * failure case : lsq_matrix is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsWithLsqMatrixNullPointer) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	double *out = nullptr;
	SIMD_ALIGN
	double out_vector[num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQCoefficientsWithLsqMatrixNotAligned
 * failure case : lsq_matrix is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsWithLsqMatrixNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, out_unaligned, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQCoefficientsWithLsqVectorNullPointer
 * failure case : lsq_vector is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsWithLsqVectorNullPointer) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double *out_vector = nullptr;
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	SetChebyshevModel(num_data, num_model, model);

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, out, out_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_GetLSQCoefficientsWithLsqVectorNotAligned
 * failure case : lsq_vector is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, GetLSQCoefficientsWithLsqVectorNotAligned) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, out, out_vector_unaligned);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test GetLSQCoefficientsTooManyMaskedData:
 * failure case of too many masked data
 * returned value : Status_kNG
 */
TEST_F(NumericOperation, GetLSQCoefficientsTooManyMaskedData) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, out, out_vector);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kNG), status);
}

/*
 * Test sakura_UpdateLSQCoefficients
 * successful case
 * note : repeating NUM_REPEAT2 times for performance measurement
 */
TEST_F(NumericOperation, UpdateLSQCoefficients) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[num_exclude_indices] = { };
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	size_t const num_lsq_matrix(num_model * num_model);
	SIMD_ALIGN
	double in_lsq_matrix[num_lsq_matrix];
	SIMD_ALIGN
	double in_lsq_matrix_orig[num_lsq_matrix];
	SIMD_ALIGN
	double in_lsq_vector[num_model];
	SIMD_ALIGN
	double in_lsq_vector_orig[num_model];
	SIMD_ALIGN
	size_t use_idx[num_model] = { };
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	SIMD_ALIGN
	double answer[num_lsq_matrix];
	SIMD_ALIGN
	double answer_vector[num_model];

	SetInputData(num_data, in_data);
	SetBoolConstant(true, num_data, in_mask);
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
		in_mask[exclude_indices[i]] = false;
	}
	SetChebyshevModel(num_data, num_model, model);
	LIBSAKURA_SYMBOL (Status)
	status_getlsq = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, answer, answer_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_getlsq);

	if (verbose) {
		PrintArray("in_mask", num_data, in_mask);
		PrintArray("model  ", num_data, num_model, model);
	}

	SetAnswers(num_data, in_data, num_model, model, num_model,
			in_lsq_matrix_orig, in_lsq_vector_orig);

	SetBoolConstant(true, num_data, in_mask);
	size_t const num_repeat(NUM_REPEAT2);
	double elapsed_time = 0.0;
	for (size_t i = 0; i < num_repeat; ++i) {
		for (size_t j = 0; j < ELEMENTSOF(in_lsq_matrix); ++j) {
			in_lsq_matrix[j] = in_lsq_matrix_orig[j];
		}
		for (size_t j = 0; j < ELEMENTSOF(in_lsq_vector); ++j) {
			in_lsq_vector[j] = in_lsq_vector_orig[j];
		}
		double start = GetCurrentTime();
		LIBSAKURA_SYMBOL (Status)
		status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
				num_exclude_indices, exclude_indices, num_model, model,
				num_model, use_idx, in_lsq_matrix, in_lsq_vector);
		double end = GetCurrentTime();
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
 * Test sakura_UpdateLSQCoefficientsNumDataZero
 * failure case : num_data == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsNumDataZero) {
	size_t const num_data(0);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[NUM_DATA2];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices, num_model, model, num_model,
			use_idx, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithDataNullPointer
 * failure case : data is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithDataNullPointer) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, data_np, in_mask,
			num_exclude_indices, exclude_indices, num_model, model, num_model,
			use_idx, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithDataNotAligned
 * failure case : data is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithDataNotAligned) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data_unaligned,
			in_mask, num_exclude_indices, exclude_indices, num_model, model,
			num_model, use_idx, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithExcludeIndicesNullPointer
 * failure case : exclude_indices is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithExcludeIndicesNullPointer) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices_np, num_model, model,
			num_model, use_idx, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithExcludeIndicesNotAligned
 * failure case : exclude_indices is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithExcludeIndicesNotAligned) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices_unaligned, num_model, model,
			num_model, use_idx, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithNumExcludeIndicesGreaterThanNumData
 * failure case : num_exclude_indices > num_data
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithNumExcludeIndicesGreaterThanNumData) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices_toolarge, exclude_indices, num_model, model,
			num_model, use_idx, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithExcludeIndicesHasValueEqualToNumData
 * failure case : exclude_indices[i] == num_data (where 0 <= i < num_exclude_indices)
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithExcludeIndicesHasValueEqualToNumData) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices, num_model, model, num_model,
			use_idx, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithExcludeIndicesHasValueGreaterThanNumData
 * failure case : exclude_indices[i] > num_data (where 0 <= i < num_exclude_indices)
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithExcludeIndicesHasValueGreaterThanNumData) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices, num_model, model, num_model,
			use_idx, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithNumModelBasisZero
 * failure case : num_model == 0
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithNumModelBasisZero) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices, num_model_zero, model,
			num_model, use_idx, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithNumModelBasisGreaterThanNumData
 * failure case : num_model > num_data
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithNumModelBasisGreaterThanNumData) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices, num_model_toolarge, model,
			num_model, use_idx, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithBasisDataNullPointer
 * failure case : model is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithBasisDataNullPointer) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices, num_model, model_np,
			num_model, use_idx, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithBasisDataNotAligned
 * failure case : model is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithBasisDataNotAligned) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices, num_model, model_unaligned,
			num_model, use_idx, in_lsq_matrix, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithLsqMatrixNullPointer
 * failure case : in_lsq_matrix is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithLsqMatrixNullPointer) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices, num_model, model, num_model,
			use_idx, in_lsq_matrix_np, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithLsqMatrixNotAligned
 * failure case : in_lsq_matrix is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithLsqMatrixNotAligned) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices, num_model, model, num_model,
			use_idx, in_lsq_matrix_unaligned, in_lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithLsqVectorNullPointer
 * failure case : in_lsq_vector is a null pointer
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithLsqVectorNullPointer) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices, num_model, model, num_model,
			use_idx, in_lsq_matrix, in_lsq_vector_np);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_UpdateLSQCoefficientsWithLsqVectorNotAligned
 * failure case : in_lsq_vector is not aligned
 * returned value : Status_kInvalidArgument
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsWithLsqVectorNotAligned) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices, num_model, model, num_model,
			use_idx, in_lsq_matrix, in_lsq_vector_unaligned);
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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

	LIBSAKURA_SYMBOL (Status)
	coeff_status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, lsq_matrix, lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);

	if (verbose) {
		PrintArray("lsq_matrix", num_model, num_model, lsq_matrix);
		PrintArray("lsq_vector", num_model, lsq_vector);
	}

	LIBSAKURA_SYMBOL (Status)
	solve_status = sakura_SolveSimultaneousEquationsByLUDouble(num_model,
			lsq_matrix, lsq_vector, out);
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
	double *model = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_model(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*model)
							* ELEMENTSOF(lsq_vector) * ELEMENTSOF(in_data),
					&model));
	if (model == nullptr) {
		throw bad_alloc();
	}
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	double *lsq_matrix = nullptr;
	unique_ptr<void, DefaultAlignedMemory> storage_for_lsq_matrix(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(*lsq_matrix)
							* ELEMENTSOF(lsq_vector) * ELEMENTSOF(lsq_vector),
					&lsq_matrix));
	if (lsq_matrix == nullptr) {
		throw bad_alloc();
	}
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

	LIBSAKURA_SYMBOL (Status)
	coeff_status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, lsq_matrix, lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);

	if (verbose) {
		PrintArray("lsq_matrix", num_model, num_model, lsq_matrix);
		PrintArray("lsq_vector", num_model, lsq_vector);
	}

	double start, end;
	double elapsed_time = 0.0;
	size_t const num_repeat = NUM_REPEAT3;
	for (size_t i = 0; i < num_repeat; ++i) {
		start = GetCurrentTime();
		LIBSAKURA_SYMBOL (Status)
		solve_status = sakura_SolveSimultaneousEquationsByLUDouble(num_model,
				lsq_matrix, lsq_vector, out);
		end = GetCurrentTime();
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

	LIBSAKURA_SYMBOL (Status)
	solve_status = sakura_SolveSimultaneousEquationsByLUDouble(num_model,
			in_matrix_np, lsq_vector, out);
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

	LIBSAKURA_SYMBOL (Status)
	solve_status = sakura_SolveSimultaneousEquationsByLUDouble(num_model,
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

	LIBSAKURA_SYMBOL (Status)
	solve_status = sakura_SolveSimultaneousEquationsByLUDouble(num_model,
			lsq_matrix, in_vector_np, out);
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

	LIBSAKURA_SYMBOL (Status)
	solve_status = sakura_SolveSimultaneousEquationsByLUDouble(num_model,
			lsq_matrix, in_vector_unaligned, out);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_SolveSimultaneousEquationsByLUWithInVectorAndOutOfSameAddress
 * failure case, with a parameter 'out' having the same address with 'in_vector'.
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLUWithInVectorAndOutOfSameAddress) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	SIMD_ALIGN
	double lsq_matrix[ELEMENTSOF(lsq_vector) * ELEMENTSOF(lsq_vector)];
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

	LIBSAKURA_SYMBOL (Status)
	coeff_status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_model, use_idx, lsq_matrix, lsq_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), coeff_status);

	double *out = lsq_vector;
	LIBSAKURA_SYMBOL (Status)
	solve_status = sakura_SolveSimultaneousEquationsByLUDouble(num_model,
			lsq_matrix, lsq_vector, out);
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

	LIBSAKURA_SYMBOL (Status)
	solve_status = sakura_SolveSimultaneousEquationsByLUDouble(num_model,
			lsq_matrix, lsq_vector, out);
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

	LIBSAKURA_SYMBOL (Status)
	solve_status = sakura_SolveSimultaneousEquationsByLUDouble(num_model,
			lsq_matrix, lsq_vector, out_unaligned);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), solve_status);
}

/*
 * Test sakura_GetLSQCoefficients
 * successful case with num_lsq_bases < num_model
 */
TEST_F(NumericOperation, GetLSQCoefficientsOrder) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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
	LIBSAKURA_SYMBOL (Status)
	status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_coeff, use_idx, out, out_vector);
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
 * Test sakura_GetLSQCoefficients
 * failure case with
 * num_lsq_bases > num_model
 * num_lsq_bases = 0
 */
TEST_F(NumericOperation, GetLSQCoefficientsBadOrder) {
	size_t const num_data(NUM_DATA2);
	size_t const num_model(NUM_MODEL2);
	size_t const bad_coeffs[] = { NUM_MODEL2 + 1, 0 };

	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}

	for (size_t i = 0; i < ELEMENTSOF(bad_coeffs); ++i) {
		size_t num_coeff = bad_coeffs[i];
		SIMD_ALIGN
		double out[num_coeff * num_coeff];
		SIMD_ALIGN
		double out_vector[num_coeff];
		cout << "number of fitting coefficients = " << num_coeff << " (model = "
				<< num_model << ")" << endl;
		LIBSAKURA_SYMBOL (Status)
		status = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
				num_model, model, num_coeff, use_idx, out, out_vector);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
	}
}

/*
 * Test sakura_UpdateLSQCoefficients
 * successful case with num_coeff < num_model
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsOrder) {
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
	size_t use_idx[num_model];
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
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
	LIBSAKURA_SYMBOL (Status)
	status_getlsq = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, num_coeff, use_idx, answer, answer_vector);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_getlsq);

	if (verbose) {
		PrintArray("in_mask", ELEMENTSOF(in_mask), in_mask);
		PrintArray("model  ", num_data, num_model, model);
	}
	cout << "number of fitting coefficients = " << num_coeff << " (model = "
			<< num_model << ")" << endl;

	SetAnswers(num_data, in_data, num_model, model, num_coeff,
			in_lsq_matrix_orig, in_lsq_vector_orig);

	SetBoolConstant(true, ELEMENTSOF(in_data), in_mask);
	for (size_t j = 0; j < ELEMENTSOF(in_lsq_matrix); ++j) {
		in_lsq_matrix[j] = in_lsq_matrix_orig[j];
	}
	for (size_t j = 0; j < ELEMENTSOF(in_lsq_vector); ++j) {
		in_lsq_vector[j] = in_lsq_vector_orig[j];
	}
	LIBSAKURA_SYMBOL (Status)
	status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_exclude_indices, exclude_indices, num_model, model, num_coeff,
			use_idx, in_lsq_matrix, in_lsq_vector);
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
 * Test sakura_UpdateLSQCoefficients
 * failure case with
 * num_coeff > num_model
 * num_coeff = 0
 */
TEST_F(NumericOperation, UpdateLSQCoefficientsBadOrder) {
	size_t const num_data(NUM_DATA);
	size_t const num_model(NUM_MODEL);
	size_t const bad_coeffs[2] = { num_model + 1, 0 };
	size_t const good_coeff(num_model);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_exclude_indices(NUM_EXCLUDE);
	SIMD_ALIGN
	size_t exclude_indices[NUM_EXCLUDE];
	for (size_t i = 0; i < num_exclude_indices; ++i) {
		exclude_indices[i] = i * i;
	}
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_data)];
	SIMD_ALIGN
	size_t use_idx[num_model] = { };
	for (size_t i = 0; i < num_model; ++i) {
		use_idx[i] = i;
	}
	SIMD_ALIGN
	double in_lsq_vector[good_coeff];
	SIMD_ALIGN
	double in_lsq_matrix[ELEMENTSOF(in_lsq_vector) * ELEMENTSOF(in_lsq_vector)];

	SetInputData(ELEMENTSOF(in_data), in_data);
	SetBoolConstant(true, ELEMENTSOF(in_data), in_mask);
	SetChebyshevModel(num_data, num_model, model);
	LIBSAKURA_SYMBOL (Status)
	status_getlsq = sakura_GetLSQCoefficientsDouble(num_data, in_data, in_mask,
			num_model, model, good_coeff, use_idx, in_lsq_matrix,
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

		LIBSAKURA_SYMBOL (Status)
		status = sakura_UpdateLSQCoefficientsDouble(num_data, in_data, in_mask,
				num_exclude_indices, exclude_indices, num_model, model,
				num_coeff, use_idx, in_lsq_matrix, in_lsq_vector);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
	}
}

/*
 * Test sakura_LMFitGaussianFloat
 * successful case with single Gaussian + noise data with some masked out regions
 */
TEST_F(NumericOperation, LMFitGaussianFloatSingle) {
	size_t const num_data = 100;
	float sigma_noise = 1.0;
	// Answer of Gaussian parameters
	size_t const num_lines = 1;
	double ans_height[num_lines] = { 7.0 };
	double ans_center[num_lines] = { 37.0 };
	double ans_sigma[num_lines] = { 3.5 };
	SIMD_ALIGN
	float data[num_data];
	SetFloatGaussian(num_data, sigma_noise, num_lines, ans_height, ans_center,
			ans_sigma, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	// Add bad data and mask them
	for (size_t i = 0; i < 10; ++i) {
		data[i] = 100.0;
		mask[i] = false;
	}
	for (size_t i = 50; i < 52; ++i) {
		data[i] = -300.0;
		mask[i] = false;
	}
	for (size_t i = 80; i < 100; ++i) {
		data[i] = 300.0;
		mask[i] = false;
	}

	// Initial guess of Gaussian parameters
	double out_height[1] = { 5.0 };
	double out_center[1] = { 40.0 };
	double out_sigma[1] = { 5.0 };

	double err_height[1];
	double err_center[1];
	double err_sigma[1];

	LIBSAKURA_SYMBOL (Status)
	status = sakura_LMFitGaussianFloat(num_data, data, mask, num_lines,
			out_height, out_center, out_sigma, err_height, err_center,
			err_sigma);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);

	std::cout << "(height) answer: " << ans_height[0] << ", result: "
			<< out_height[0] << ", error(1sigma): " << err_height[0]
			<< std::endl;
	std::cout << "(center) answer: " << ans_center[0] << ", result: "
			<< out_center[0] << ", error(1sigma): " << err_center[0]
			<< std::endl;
	std::cout << "(sigma)  answer: " << ans_sigma[0] << ", result: "
			<< out_sigma[0] << ", error(1sigma): " << err_sigma[0] << std::endl;

	CheckInRange(ans_height[0], out_height[0], err_height[0]);
	CheckInRange(ans_center[0], out_center[0], err_center[0]);
	CheckInRange(ans_sigma[0], out_sigma[0], err_sigma[0]);
}

/*
 * Test sakura_LMFitGaussianFloat
 * successful case with two Gaussians + noise data with some masked out regions
 * Note: the input data contain two Gaussian profiles, but fitting is done
 * for each peak, since it is quite difficult to fit double Gaussians at once
 * using Eigen's LM solver...
 */
TEST_F(NumericOperation, LMFitGaussianFloatDouble) {
	size_t const num_data = 100;
	float const sigma_noise = 0.10;
	// Answer of Gaussian parameters
	size_t const num_line = 2;
	double const ans_height[num_line] = { 7.0, 6.0 };
	double const ans_center[num_line] = { 27.0, 70.0 };
	double const ans_sigma[num_line] = { 1.0, 1.5 };
	SIMD_ALIGN
	float data[num_data];
	SetFloatGaussian(num_data, sigma_noise, num_line, ans_height, ans_center,
			ans_sigma, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	// Add bad data and mask them
	for (size_t i = 0; i < 10; ++i) {
		data[i] = 100.0;
		mask[i] = false;
	}
	for (size_t i = 50; i < 52; ++i) {
		data[i] = -300.0;
		mask[i] = false;
	}
	for (size_t i = 90; i < 100; ++i) {
		data[i] = 300.0;
		mask[i] = false;
	}

	// Initial guess of Gaussian parameters
	double out_height[num_line] = { 10.0, 10.0 };
	double out_center[num_line] = { 25.0, 69.0 };
	double out_sigma[num_line] = { 1.0, 1.0 };

	double err_height[num_line];
	double err_center[num_line];
	double err_sigma[num_line];

	for (size_t iline = 0; iline < num_line; ++iline) {
		for (size_t i = 0; i < num_data; ++i) {
			mask[i] = false;
		}
		size_t i_min = out_center[iline] - 5.0 * out_sigma[iline];
		size_t i_max = out_center[iline] + 5.0 * out_sigma[iline];
		for (size_t i = i_min; i < i_max; ++i) {
			mask[i] = true;
		}
		LIBSAKURA_SYMBOL (Status)
		status = sakura_LMFitGaussianFloat(num_data, data, mask, 1,
				&out_height[iline], &out_center[iline], &out_sigma[iline],
				&err_height[iline], &err_center[iline], &err_sigma[iline]);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);

		std::cout << "line profile #" << iline << ":" << std::endl;
		std::cout << "(height) answer: " << ans_height[iline] << ", result: "
				<< out_height[iline] << ", error(1sigma): " << err_height[iline]
				<< std::endl;
		std::cout << "(center) answer: " << ans_center[iline] << ", result: "
				<< out_center[iline] << ", error(1sigma): " << err_center[iline]
				<< std::endl;
		std::cout << "(sigma)  answer: " << ans_sigma[iline] << ", result: "
				<< out_sigma[iline] << ", error(1sigma): " << err_sigma[iline]
				<< std::endl;

		CheckInRange(ans_height[iline], out_height[iline], err_height[iline]);
		CheckInRange(ans_center[iline], out_center[iline], err_center[iline]);
		CheckInRange(ans_sigma[iline], out_sigma[iline], err_sigma[iline]);
	}
}

/*
 * Test sakura_LMFitGaussianFloat
 * failure case: num_data is less than needed value (= 3*num_lines)
 */
TEST_F(NumericOperation, LMFitGaussianFloatTooFewNumData) {
	size_t const num_data = 2;
	float sigma_noise = 1.0;
	// Answer of Gaussian parameters
	size_t const num_lines = 1;
	double ans_height[num_lines] = { 7.0 };
	double ans_center[num_lines] = { 37.0 };
	double ans_sigma[num_lines] = { 3.5 };
	SIMD_ALIGN
	float data[num_data];
	SetFloatGaussian(num_data, sigma_noise, num_lines, ans_height, ans_center,
			ans_sigma, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);

	// Initial guess of Gaussian parameters
	double out_height[1] = { 1.0 };
	double out_center[1] = { 50.0 };
	double out_sigma[1] = { 1.0 };

	double err_height[1];
	double err_center[1];
	double err_sigma[1];

	LIBSAKURA_SYMBOL (Status)
	status = sakura_LMFitGaussianFloat(num_data, data, mask, num_lines,
			out_height, out_center, out_sigma, err_height, err_center,
			err_sigma);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

/*
 * Test sakura_LMFitGaussianFloat
 * failure case: number of effective data is less than needed value (= 3*num_lines)
 */
TEST_F(NumericOperation, LMFitGaussianFloatTooFewNumDataDueToMasking) {
	size_t const num_data = 10;
	float sigma_noise = 1.0;
	// Answer of Gaussian parameters
	size_t const num_lines = 1;
	double ans_height[num_lines] = { 7.0 };
	double ans_center[num_lines] = { 37.0 };
	double ans_sigma[num_lines] = { 3.5 };
	SIMD_ALIGN
	float data[num_data];
	SetFloatGaussian(num_data, sigma_noise, num_lines, ans_height, ans_center,
			ans_sigma, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	// Give large mask (only 2 data effective at index of 0 and 4)
	for (size_t i = 1; i < 4; ++i) {
		mask[i] = false;
	}
	for (size_t i = 5; i < num_data; ++i) {
		mask[i] = false;
	}

	// Initial guess of Gaussian parameters
	double out_height[1] = { 1.0 };
	double out_center[1] = { 50.0 };
	double out_sigma[1] = { 1.0 };

	double err_height[1];
	double err_center[1];
	double err_sigma[1];

	LIBSAKURA_SYMBOL (Status)
	status = sakura_LMFitGaussianFloat(num_data, data, mask, num_lines,
			out_height, out_center, out_sigma, err_height, err_center,
			err_sigma);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kNG), status);
}
