#include <iostream>
#include <string>

#include <libsakura/sakura.h>
#include "gtest/gtest.h"

/* the number of elements in input/output array to test */
#define NUM_IN 5

using namespace std;

/*
 * A super class to test logical operations of array(s)
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
	}

	virtual void TearDown()
	{}

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

	float in1_[NUM_IN];
	float in2_[NUM_IN];
	bool verbose;

};

/*
 * Test subtraction in1 - in2 by sakura_OperateFloatSubtraction
 * RESULT:
 * out = [ -5.0, -3.0, 3.0, -5.0, -1.0 ]
 */
TEST_F(NumericOperation, Subtraction) {
	float out[NUM_IN];
	float result[NUM_IN] = {-5.0, -3.0, 3.0, -5.0, -1.0};
	size_t const num_in(NUM_IN);

	if (verbose) PrintInputs();

	sakura_OperateFloatSubtraction(num_in, in1_, in2_, out);

	if (verbose) PrintArray("out", num_in, out);

	// Verification
	for (size_t i = 0 ; i < num_in ; i++){
		ASSERT_EQ(out[i], result[i]);
	}
}
