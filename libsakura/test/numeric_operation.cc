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
 * Test sakura_GetMatrixCoefficientsForLeastSquareFitting
 * RESULT:
 * out = []
 */
TEST_F(NumericOperation, GetMatrixCoefficientsForLeastSquareFitting) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	bool in_mask[num_data];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model * num_model];
	SIMD_ALIGN
	double model[num_model * ELEMENTSOF(in_mask)];
	SIMD_ALIGN
	double answer[num_model * num_model];

	// Setup in_mask--------------------------
	for (size_t i = 0; i < num_data; ++i) {
		in_mask[i] = true;
	}
	// Setup model----------------------------
	size_t idx = 0;
	for (size_t j = 0; j < num_data; ++j) {
		for (size_t i = 0; i < num_model; ++i) {
			double value = 1.0;
			//poly----
			//for (size_t k = 0; k < i; ++k) { value *= (double)j; }
			//chebyshev----
			double x = 2.0 * (double) j / (double) (num_data - 1) - 1.0;
			if (i == 0) {
				value = 1.0;
			} else if (i == 1) {
				value = x;
			} else {
				value = 2.0 * x * model[idx - 1] - model[idx - 2];
			}
			//----------
			model[idx] = value;
			idx++;
		}
	}
	// Setup answer---------------------------
	idx = 0;
	for (size_t i = 0; i < num_model; ++i) {
		for (size_t j = 0; j < num_model; ++j) {
			double val = 0.0;
			for (size_t k = 0; k < num_data; ++k) {
				val += model[num_model * k + i] * model[num_model * k + j];
			}
			answer[idx] = val;
			idx++;
		}
	}
	//----------------------------------------------

	if (verbose) {
		PrintArray("in_mask", num_data, in_mask);
		PrintArray("model  ", num_data, num_model, model);
	}

	double start, end;
	size_t const num_repeat = NUM_REPEAT;
	LIBSAKURA_SYMBOL(Status) status;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_GetMatrixCoefficientsForLeastSquareFitting(num_data,
				in_mask, num_model, model, out);
	}
	end = sakura_GetCurrentTime();

	cout << "Elapse time of " << num_repeat << " repetition: " << end - start
			<< " sec." << endl;
	if (verbose) {
		PrintArray("out   ", num_model, num_model, out);
		PrintArray("answer", num_model, num_model, answer);
	}

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
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
}

/*
 * Test sakura_GetVectorCoefficientsForLeastSquareFitting
 * RESULT:
 * out = []
 */
TEST_F(NumericOperation, GetVectorCoefficientsForLeastSquareFitting) {
	size_t const num_data(NUM_DATA);
	SIMD_ALIGN
	float in_data[num_data];
	SIMD_ALIGN
	bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN
	double out[num_model];
	SIMD_ALIGN
	double model[ELEMENTSOF(out) * ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(out)];

	for (size_t i = 0; i < num_data; ++i) {
		double x = (double) i;
		double val = 4.0 + 0.000056 * x - 0.000037 * x * x
				+ 0.0000012 * x * x * x + 0.0000009 * x * x * x * x
				+ 0.0000006 * x * x * x * x * x;
		in_data[i] = (float) val;
		in_mask[i] = true;
	}
	size_t idx = 0;
	for (size_t j = 0; j < num_data; ++j) {
		for (size_t i = 0; i < num_model; ++i) {
			double value = 1.0;
			//poly----
			//for (size_t k = 0; k < i; ++k) { value *= (double)j; }
			//chebyshev----
			double x = 2.0 * (double) j / (double) (num_data - 1) - 1.0;
			if (i == 0) {
				value = 1.0;
			} else if (i == 1) {
				value = x;
			} else {
				value = 2.0 * x * model[idx - 1] - model[idx - 2];
			}
			//----------
			model[idx] = value;
			idx++;
		}
	}
	for (size_t i = 0; i < num_model; ++i) {
		double val = 0.0;
		for (size_t j = 0; j < num_data; ++j) {
			val += model[num_model * j + i] * in_data[j];
		}
		answer[i] = val;
	}

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
		PrintArray("model  ", num_data, num_model, model);
	}

	double start, end;
	size_t const num_repeat = NUM_REPEAT;
	LIBSAKURA_SYMBOL(Status) status;
	start = sakura_GetCurrentTime();
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_GetVectorCoefficientsForLeastSquareFitting(num_data,
				in_data, in_mask, num_model, model, out);
	}
	end = sakura_GetCurrentTime();

	cout << "Elapse time of " << num_repeat
			<< " repetition:                                     "
			<< (end - start) << " sec." << endl;
	if (verbose) {
		PrintArray("out   ", num_model, out);
		PrintArray("answer", num_model, answer);
	}

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_model; ++i) {
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
}

#define NUM_DATA_TESTLU 50
#define NUM_MODEL_TESTLU 5
/*
 * Test sakura_SolveSimultaneousEquationsByLU
 * RESULT:
 * out = []
 *//*
 TEST_F(NumericOperation, SolveSimultaneousEquationsByLU) {
 size_t const num_data(NUM_DATA_TESTLU);
 SIMD_ALIGN float in_data[num_data];
 SIMD_ALIGN bool in_mask[ELEMENTSOF(in_data)];
 size_t const num_model(NUM_MODEL_TESTLU);
 SIMD_ALIGN double lsq_vector[num_model];
 SIMD_ALIGN double model[ELEMENTSOF(lsq_vector)*ELEMENTSOF(in_data)];
 SIMD_ALIGN double lsq_matrix[ELEMENTSOF(lsq_vector)*ELEMENTSOF(lsq_vector)];
 SIMD_ALIGN double out[ELEMENTSOF(lsq_vector)];
 SIMD_ALIGN float answer[ELEMENTSOF(lsq_vector)] = {1.0, 1.0, 1.0, 1.0, 1.0};

 for (size_t i = 0; i < num_data; ++i) {
 double x = (double)i;

 double val = answer[0]
 + answer[1] * x
 + answer[2] * x * x
 + answer[3] * x * x * x
 + answer[4] * x * x * x * x;

 in_data[i] = (float)val;
 in_mask[i] = true;
 }
 size_t idx = 0;
 for (size_t i = 0; i < num_model; ++i) {
 for (size_t j = 0; j < num_data; ++j) {
 double value = 1.0;
 for (size_t k = 0; k < i; ++k) {
 value *= (double)j;
 }
 model[idx] = value;
 idx++;
 }
 }


 if (verbose) {
 PrintArray("in_data",    num_data,  in_data);
 PrintArray("in_mask",    num_data,  in_mask);
 PrintArray("model",      num_model, num_data, model);
 }

 sakura_GetMatrixCoefficientsForLeastSquareFitting(
 num_data, in_mask, num_model, model, lsq_matrix);
 sakura_GetVectorCoefficientsForLeastSquareFitting(
 num_data, in_data, in_mask, num_model, model, lsq_vector);

 if (verbose) {
 PrintArray("lsq_matrix", num_model, num_model, lsq_matrix);
 PrintArray("lsq_vector", num_model, lsq_vector);
 }

 double start, end;
 size_t const num_repeat = 1;//NUM_REPEAT;
 start = sakura_GetCurrentTime();
 LIBSAKURA_SYMBOL(Status) status;
 for (size_t i = 0; i < num_repeat; ++i) {
 status = sakura_SolveSimultaneousEquationsByLU(
 num_model, lsq_matrix, lsq_vector, out);
 }
 end = sakura_GetCurrentTime();

 if (verbose) {
 cout << "Elapse time of " << num_repeat << " repetition: " << end - start << " sec." << endl;
 PrintArray("out   ", num_model, out);
 PrintArray("answer", num_model, answer);
 }

 // Verification
 EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
 for (size_t i = 0 ; i < num_model; ++i){
 float deviation = (out[i] - answer[i]) / answer[i];
 ASSERT_LE(deviation, 1e-7);
 }
 }
 */
