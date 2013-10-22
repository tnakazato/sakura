#include <iostream>
#include <string>

#include <libsakura/sakura.h>
#include "aligned_memory.h"
#include "gtest/gtest.h"

/* the number of elements in input/output array to test */
#define NUM_IN 5
#define NUM_MODEL 3

using namespace std;

/*
 * A super class to test numeric operations of array(s)
 * INPUTS:
 * - in1 = [ 0.0, 1.0, 2.0, -3.0, -4.0 ]
 * - in2 = [ 5.0, 4.0, -1.0, 2.0, -3.0 ]
 */
class NumericOperation : public ::testing::Test
{
protected:

	NumericOperation() : verbose(false)
	{}

	virtual void SetUp()
	{
		in1_[0] = 0.0;
		in1_[1] = 1.0;
		in1_[2] = 2.0;
		in1_[3] = -3.0;
		in1_[4] = -4.0;

		in2_[0] = 5.0;
		in2_[1] = 4.0;
		in2_[2] = -1.0;
		in2_[3] = 2.0;
		in2_[4] = -3.0;

		// Initialize sakura
		LIBSAKURA_SYMBOL(Status) status =
				LIBSAKURA_SYMBOL(Initialize)(nullptr, nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}

	virtual void TearDown()
	{
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	void PrintInputs(){
		PrintArray("in1", NUM_IN, in1_);
		PrintArray("in2", NUM_IN, in2_);
	}

	void PrintArray(char const *name, size_t num_in, float *in){
		cout << name << " = [";
		for (size_t i = 0; i < num_in-1; i++)
			cout << in[i] << ", " ;
		cout << in[num_in-1];
		cout << " ]" << endl;
	}

	SIMD_ALIGN float in1_[NUM_IN];
	SIMD_ALIGN float in2_[NUM_IN];
	bool verbose;

};

/*
 * Test subtraction in1 - in2 by sakura_OperateFloatSubtraction
 * RESULT:
 * out = [ -5.0, -3.0, 3.0, -5.0, -1.0 ]
 */
TEST_F(NumericOperation, Subtraction) {
	SIMD_ALIGN float out[NUM_IN];
	float result[NUM_IN] = {-5.0, -3.0, 3.0, -5.0, -1.0};
	size_t const num_in(NUM_IN);

	if (verbose) PrintInputs();

	LIBSAKURA_SYMBOL(Status) status =
			sakura_OperateFloatSubtraction(num_in, in1_, in2_, out);

	if (verbose) PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0 ; i < num_in; ++i){
		ASSERT_EQ(out[i], result[i]);
	}
}

/*
 * Test sakura_GetBestFitModel
 * RESULT:
 * out = []
 */
TEST_F(NumericOperation, GetBestFitModel) {
	size_t const num_in(NUM_IN);
	SIMD_ALIGN float in_data[NUM_IN] = {1.0, 3.0, 7.0, 130.0, 21.0};
	SIMD_ALIGN bool in_mask[NUM_IN] = {true, true, true, false, true};
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN double model[NUM_MODEL*NUM_IN] = {1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 4.0, 9.0, 16.0};
	SIMD_ALIGN float out[NUM_IN];
	float result[NUM_IN] = {1.0, 3.0, 7.0, 13.0, 21.0};

	LIBSAKURA_SYMBOL(Status) status =
			sakura_GetBestFitModel(num_in, in_data, in_mask, num_model, model, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0 ; i < num_model; ++i){
		ASSERT_EQ(out[i], result[i]);
	}
}

/*
 * Test sakura_GetLeastSquareMatrix
 * RESULT:
 * out = []
 */
TEST_F(NumericOperation, GetLeastSquareMatrix) {
	size_t const num_in(NUM_IN);
	SIMD_ALIGN float in_data[NUM_IN] = {1.0, 3.0, 7.0, 130.0, 21.0};
	SIMD_ALIGN bool in_mask[NUM_IN] = {true, true, true, false, true};
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN double model[NUM_MODEL*NUM_IN] = {1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 1.0, 4.0, 9.0, 16.0};
	SIMD_ALIGN double out[NUM_MODEL*NUM_MODEL];
	SIMD_ALIGN double out_vector[NUM_MODEL];
	float result_matrix[NUM_MODEL*NUM_MODEL] = {4.0, 7.0, 21.0, 7.0, 21.0, 73.0, 21.0, 73.0, 273.0};
	float result_vector[NUM_MODEL] = {32.0, 101.0, 367.0};

	LIBSAKURA_SYMBOL(Status) status =
			sakura_GetLeastSquareMatrix(num_in, in_data, in_mask, num_model, model, out, out_vector);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0 ; i < num_model*num_model; ++i){
		ASSERT_EQ(out[i], result_matrix[i]);
	}
	for (size_t i = 0 ; i < num_model; ++i){
		ASSERT_EQ(out_vector[i], result_vector[i]);
	}
}

/*
 * Test sakura_SolveSimultaneousEquationsByLU
 * RESULT:
 * out = []
 */
TEST_F(NumericOperation, SolveSimultaneousEquationsByLU) {
	size_t const num_model(NUM_MODEL);
	SIMD_ALIGN double lsq_matrix[NUM_MODEL*NUM_MODEL] = {4.0, 7.0, 21.0, 7.0, 21.0, 73.0, 21.0, 73.0, 273.0};
	SIMD_ALIGN double lsq_vector[NUM_MODEL] = {32.0, 101.0, 367.0};
	SIMD_ALIGN double out[NUM_MODEL];
	float result[NUM_MODEL] = {1.0, 1.0, 1.0};

	LIBSAKURA_SYMBOL(Status) status =
				sakura_SolveSimultaneousEquationsByLU(num_model, lsq_matrix, lsq_vector, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0 ; i < num_model; ++i){
		ASSERT_EQ((float)out[i], result[i]);
	}
}
