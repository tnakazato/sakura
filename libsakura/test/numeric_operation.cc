#include <iostream>
#include <string>
#include <sys/time.h>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "aligned_memory.h"
#include "gtest/gtest.h"

/* the number of elements in input/output array to test */
#define NUM_DATA 4096
#define NUM_MODEL 6
#define NUM_REPEAT 10000

#define NUM_MODEL_SMALL 3

using namespace std;

/*
 * A super class to test numeric operations of array(s)
 */
class NumericOperation : public ::testing::Test {
protected:

	NumericOperation() : verbose(false)
	{}

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
		if (print_name) cout << name << " = ";
		cout << "[";
		for (size_t i = start_idx; i < start_idx+print_length-1; ++i)
			cout << data[i] << ", " ;
		cout << data[start_idx+print_length-1];
		cout << " ]";
		if (newline) cout << endl;
	}
	//1D double array
	void PrintArray(char const *name, size_t print_length, double const *data,
			size_t start_idx = 0, bool print_name = true, bool newline = true) {
		if (print_name) cout << name << " = ";
		cout << "[";
		for (size_t i = start_idx; i < start_idx+print_length-1; ++i)
			cout << data[i] << ", " ;
		cout << data[start_idx+print_length-1];
		cout << " ]";
		if (newline) cout << endl;
	}
	//1D bool array
	void PrintArray(char const *name, size_t print_length, bool const *data,
			size_t start_idx = 0, bool print_name = true, bool newline = true) {
		if (print_name) cout << name << " = ";
		cout << "[";
		for (size_t i = start_idx; i < start_idx+print_length-1; ++i)
			cout << (data[i] ? "T" : "F") << ", " ;
		cout << (data[start_idx+print_length-1] ? "T" : "F");
		cout << " ]";
		if (newline) cout << endl;
	}
	//given as 1D float array but actually stores (num_row * num_column) 2D data
	//for which column loop comes inside row loop.
	void PrintArray(char const *name, size_t num_row, size_t num_column, float const *data) {
		cout << name << " = [";
		for (size_t i = 0; i < num_row; ++i) {
			PrintArray(name, num_column, data, num_column*i, false, false);
			if (i < num_row-1) cout << ", ";
		}
		cout << " ]" << endl;
	}
	//given as 1D double array but actually stores (num_row * num_column) 2D data
	//for which column loop comes inside row loop.
	void PrintArray(char const *name, size_t num_row, size_t num_column, double const *data) {
		cout << name << " = [";
		for (size_t i = 0; i < num_row; ++i) {
			PrintArray(name, num_column, data, num_column*i, false, false);
			if (i < num_row-1) cout << ", ";
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
	SIMD_ALIGN bool in_mask[num_data];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN double out[num_model*num_model];
	SIMD_ALIGN double model[num_model*ELEMENTSOF(in_mask)];
	double answer[num_model*num_model];

	// Setup in_mask--------------------------
	for (size_t i = 0; i < num_data; ++i) {
		in_mask[i] = true;
	}
	// Setup model----------------------------
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
	// Setup answer---------------------------
	idx = 0;
	for (size_t i = 0; i < num_model; ++i) {
		for (size_t j = 0; j < num_model; ++j) {
			double val = 0.0;
			for (size_t k = 0; k < num_data; ++k) {
				val += model[num_data * i + k] * model[num_data * j + k];
			}
			answer[idx] = val;
			idx++;
		}
	}
	//----------------------------------------------

	if (verbose) {
		PrintArray("in_mask", num_data, in_mask);
		PrintArray("model  ", num_model, num_data, model);
	}

	double start, end;
	size_t const num_repeat = NUM_REPEAT;
	start = sakura_GetCurrentTime();
	LIBSAKURA_SYMBOL(Status) status;
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_GetMatrixCoefficientsForLeastSquareFitting(
				num_data, in_mask, num_model, model, out);
	}
	end = sakura_GetCurrentTime();

	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

	if (verbose) {
		PrintArray("out   ", num_model, num_model, out);
		PrintArray("answer", num_model, num_model, answer);
	}

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0 ; i < num_model*num_model; ++i){
		float deviation = (out[i] - answer[i]) / answer[i];
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
	SIMD_ALIGN float in_data[num_data];
	SIMD_ALIGN bool in_mask[ELEMENTSOF(in_data)];
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN double out[num_model];
	SIMD_ALIGN double model[ELEMENTSOF(out)*ELEMENTSOF(in_data)];
	float answer[ELEMENTSOF(out)];

	for (size_t i = 0; i < num_data; ++i) {
		double x = (double)i;
		double val = 4.0
				+ 0.000056  * x
				- 0.000037  * x * x
				+ 0.0000012 * x * x * x
				+ 0.0000009 * x * x * x * x
				+ 0.0000006 * x * x * x * x * x;
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
	for (size_t i = 0; i < num_model; ++i) {
		double val = 0.0;
		for (size_t j = 0; j < num_data; ++j) {
			val += model[num_data * i + j] * in_data[j];
		}
		answer[i] = val;
	}

/*
	cout << "------------------" << endl;
	cout.precision(40);
	cout << "in_data[4095] = " << in_data[4095] << endl;
	cout << "answer: ";
	for (size_t i = 0; i < num_model; ++i) {
		cout << answer[i];
		if (i < num_model-1) {
			cout << ", ";
		}
	}
	cout << endl;
*/
/*
	cout << "------------------" << endl;
	cout.precision(20);
	cout << "model[0]     = " << model[0]     << ", ..., [4095]  = " << model[4095]  << endl;
	cout << "model[4096]  = " << model[4096]  << ", ..., [8191]  = " << model[8191]  << endl;
	cout << "model[8192]  = " << model[8192]  << ", ..., [12287] = " << model[12287] << endl;
	cout << "model[12288] = " << model[12288] << ", ..., [16383] = " << model[16383] << endl;
	cout << "model[16384] = " << model[16384] << ", ..., [20479] = " << model[20479] << endl;
	cout << "------------------" << endl;
*/

	if (verbose) {
		PrintArray("in_data", num_data, in_data);
		PrintArray("in_mask", num_data, in_mask);
		PrintArray("model  ", num_model, num_data, model);
	}

	double start, end;
	size_t const num_repeat = NUM_REPEAT;
	start = sakura_GetCurrentTime();
	LIBSAKURA_SYMBOL(Status) status;
	for (size_t i = 0; i < num_repeat; ++i) {
		status = sakura_GetVectorCoefficientsForLeastSquareFitting(
				num_data, in_data, in_mask, num_model, model, out);
	}
	end = sakura_GetCurrentTime();
	cout << "Elapse time of actual operation: " << end - start << " sec"
			<< endl;

/*
	cout << "------------------" << endl;
	cout << "   out: ";
	for (size_t i = 0; i < num_model; ++i) {
		cout << out[i];
		if (i < num_model-1) {
			cout << ", ";
		}
	}
	cout << endl;
	cout.precision(4);
	cout << "  diff: ";
	for (size_t i = 0; i < num_model; ++i) {
		cout << ((out[i] - answer[i]) / answer[i]);
		if (i < num_model-1) {
			cout << ", ";
		}
	}
	cout << endl;
	cout << "------------------" << endl;
*/

	if (verbose) {
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

/*
 * Test sakura_SolveSimultaneousEquationsByLU
 * RESULT:
 * out = []
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLU) {
	size_t const num_model(NUM_MODEL_SMALL);
	SIMD_ALIGN double lsq_vector[num_model] = {32.0, 101.0, 367.0};
	SIMD_ALIGN double lsq_matrix[ELEMENTSOF(lsq_vector)*ELEMENTSOF(lsq_vector)] = {4.0, 7.0, 21.0, 7.0, 21.0, 73.0, 21.0, 73.0, 273.0};
	SIMD_ALIGN double out[ELEMENTSOF(lsq_vector)];
	float answer[ELEMENTSOF(lsq_vector)] = {1.0, 1.0, 1.0};

	if (verbose) {
		PrintArray("lsq_matrix", num_model, num_model, lsq_matrix);
		PrintArray("lsq_vector", num_model, lsq_vector);
	}

	LIBSAKURA_SYMBOL(Status) status =
				sakura_SolveSimultaneousEquationsByLU(num_model, lsq_matrix, lsq_vector, out);

	if (verbose) {
		PrintArray("out   ", num_model, out);
		PrintArray("answer", num_model, answer);
	}

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0 ; i < num_model; ++i){
		ASSERT_EQ((float)out[i], answer[i]);
	}
}
