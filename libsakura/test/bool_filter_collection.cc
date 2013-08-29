#include <iostream>
#include <string>
#include <sys/time.h>

#include <libsakura/sakura.h>
#include "aligned_memory.h"
#include "gtest/gtest.h"

/* the number of elements in input/output array to test */
#define NUM_IN 16
#define NUM_RANGE 2
#define NUM_IN_LONG 262144 //2**18

using namespace std;

/*
 * A super class to test various bit operation of an value and array
 * INPUTS:
 * - data = []
 * - upper_bound = []
 * - lower_bound = []
 */
template<typename DataType>
class BoolFilter : public ::testing::Test
{
protected:

	BoolFilter()
	: verbose(false)
	{
	}

	virtual void SetUp()
	{
//		for (size_t i = 0; i < NUM_IN; ++i){
//			data_[i] = i % ntype; /* repeat bit pattern of *00, *01, *10, *11,... */
//		}

		// Initialize sakura
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Initialize)();
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}

	virtual void TearDown()
	{
		// Clean-up sakura
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	SIMD_ALIGN DataType data_[NUM_IN];
	DataType upper_[NUM_RANGE];
	DataType lower_[NUM_RANGE];
	bool verbose;
};

/*
 * Tests various bool filter generation using float array
 * INPUTS:
 * - data = []
 * - upper_bound = []
 * - lower_bound = []
 */
class BoolFilterFloat : public BoolFilter<float>
{};

/*
 * Tests various bool filter generation using int array
 * INPUTS:
 * - data = []
 * - upper_bound = []
 * - lower_bound = []
 */
class BoolFilterInt : public BoolFilter<int>
{};

/*
 * Test bool filter generation sakura_SetTrueFloatInRangesInclusive
 * RESULT:
 * out = []
 */
TEST_F(BoolFilterFloat, RangesInclusive) {
	SIMD_ALIGN bool out[NUM_IN];
	bool answer[NUM_IN];// = {0, 1, 2, 3, 0, 0, 2, 2};
	//uint8_t result[NUM_IN];
	size_t const num_in(NUM_IN), num_range(NUM_RANGE);

//	if (verbose) PrintInputs();

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueFloatInRangesInclusive(num_in, data_, num_range, lower_, upper_, out);

//	if (verbose) PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
//	for (size_t i = 0 ; i < num_in ; ++i){
//		ASSERT_EQ(out[i], answer[i]);
//	}
}

/*
 * Test bool filter generation sakura_SetTrueIntInRangesInclusive
 * RESULT:
 * out = []
 */
TEST_F(BoolFilterInt, RangesInclusive) {
	SIMD_ALIGN bool out[NUM_IN];
	bool answer[NUM_IN];// = {0, 1, 2, 3, 0, 0, 2, 2};
	//uint8_t result[NUM_IN];
	size_t const num_in(NUM_IN), num_range(NUM_RANGE);

//	if (verbose) PrintInputs();

	LIBSAKURA_SYMBOL(Status) status = sakura_SetTrueIntInRangesInclusive(num_in, data_, num_range, lower_, upper_, out);

//	if (verbose) PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
//	for (size_t i = 0 ; i < num_in ; ++i){
//		ASSERT_EQ(out[i], answer[i]);
//	}
}


