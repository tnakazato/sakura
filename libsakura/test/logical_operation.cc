#include <iostream>
#include <string>
#include <math.h>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "aligned_memory.h"
#include "gtest/gtest.h"

/* the number of elements in input/output array to test */
#define NUM_IN 4

using namespace std;

/*
 * A super class to test logical operations of array(s)
 * INPUTS:
 * - in1 = [ F, F, T, T ]
 * - in2 = [ F, T, F, T ]
 */
class LogicalOperation : public ::testing::Test
{
protected:

	LogicalOperation() : verbose(false)
	{}

	virtual void SetUp()
	{
		size_t const a(2);
		for (size_t i = 0; i < NUM_IN; ++i){
			in1_[i] = (i / a > 0);
			in2_[i] = (i % a > 0);
		}

		// Initialize sakura
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Initialize)(nullptr,
				nullptr);
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

	void PrintArray(char const *name, size_t num_in, bool *in){
		cout << name << " = [";
		for (size_t i = 0; i < num_in-1; i++)
			cout << BToS(in[i]) << ", " ;
		cout << BToS(in[num_in-1]) ;
		cout << " ]" << endl;
	}

	/* Converts an input boolean value to a string.*/
	string BToS(bool in_value) {
		return in_value ? "T" : "F";
	}

	bool in1_[NUM_IN];
	bool in2_[NUM_IN];
	bool verbose;

};

/*
 * Test logical operation AND by sakura_OperateLogicalAnd
 * RESULT:
 * out = [false, false, false, true]
 */
TEST_F(LogicalOperation, And) {
	size_t const num_in(NUM_IN);
	SIMD_ALIGN bool in1[num_in];
	SIMD_ALIGN bool in2[ELEMENTSOF(in1)];
	SIMD_ALIGN bool out[ELEMENTSOF(in1)];
	bool answer[ELEMENTSOF(in1)] = {false, false, false, true};

	if (verbose) PrintInputs();

	for (size_t i = 0; i < num_in; ++i) {
		in1[i] = in1_[i];
		in2[i] = in2_[i];
	}
	LIBSAKURA_SYMBOL(Status) status = sakura_OperateLogicalAnd(num_in, in1, in2, out);

	if (verbose) PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0 ; i < num_in ; ++i){
		ASSERT_EQ(out[i], answer[i]);
	}
}
