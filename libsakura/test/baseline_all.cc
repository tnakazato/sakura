#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <sys/time.h>
#include <random>
#include <stdio.h>
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
class Baseline: public ::testing::Test {
protected:

	Baseline() :
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

	bool verbose;

};


void Destroy(LIBSAKURA_SYMBOL(BaselineContextFloat) *context, LIBSAKURA_SYMBOL (Status) status) {
	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyBaselineContextFloat(context);
	EXPECT_EQ(status, destroy_status);
}


////////////////////////////////////////////////////////////////////////////////
// TEST
////////////////////////////////////////////////////////////////////////////////

/*
 * Test sakura_DestroyBaselineContextFloat
 * successful case
 */
TEST_F(Baseline, DestroyBaselineContextFloat) {
	uint16_t const order(20);
	size_t const num_chan(4096);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	LIBSAKURA_SYMBOL (Status) create_status =
	sakura_CreateBaselineContextFloat(
						LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order,
						num_chan, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));

}

/*
 * Test sakura_DestroyBaselineContextFloat
 * failure case : context is a null pointer
 * returned value must be Status_kInvalidArgument
 */
TEST_F(Baseline, DestroyBaselineContextFloatWithContextNullPointer) {
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Destroy(context, LIBSAKURA_SYMBOL(Status_kInvalidArgument));
}


