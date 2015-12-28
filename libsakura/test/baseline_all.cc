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


void Destroy(LIBSAKURA_SYMBOL (Status) status, LIBSAKURA_SYMBOL(BaselineContextFloat) *context){
	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyBaselineContextFloat(context);
	EXPECT_EQ(status, destroy_status);
}


void Create(LIBSAKURA_SYMBOL (Status) status, sakura_BaselineType const baseline_type, int16_t const order,
		size_t const num_data, struct sakura_BaselineContextFloat** context){

	LIBSAKURA_SYMBOL (Status) create_status;

	if(baseline_type == LIBSAKURA_SYMBOL(BaselineType_kPolynomial) ||
			baseline_type == LIBSAKURA_SYMBOL(BaselineType_kChebyshev)){
		create_status = sakura_CreateBaselineContextFloat(baseline_type, order,
					num_data, context);
/*
	}else if(baseline_type == LIBSAKURA_SYMBOL(BaselineType_kChebyshev)){
		create_status = sakura_CreateBaselineContextFloat(baseline_type, order,
					num_data, context);
*/
	}else if(baseline_type == LIBSAKURA_SYMBOL(BaselineType_kCubicSpline)){
		create_status = sakura_CreateBaselineContextCubicSplineFloat(
					order, num_data, context);

	}else if(baseline_type == LIBSAKURA_SYMBOL(BaselineType_kSinusoid)){
		create_status = sakura_CreateBaselineContextSinusoidFloat(order,
					num_data, context);
	}

	EXPECT_EQ(status, create_status);

}

////////////////////////////////////////////////////////////////////////////////
// TEST (Create)
////////////////////////////////////////////////////////////////////////////////

/*
 * Test sakura_CreateBaselineContextFloatWithPolynomial/ChebyshevPolynomial
 * successful case (with polynomial/chebyshevpolynomial model using boundaries consisting of order=0 and num_data=1)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithPolynomial_ChebyshevPolynomial) {
	uint16_t const order(0);
	size_t const num_chan(1);

	size_t const num_func_types = 2;
	LIBSAKURA_SYMBOL(BaselineType) func_types[num_func_types]={LIBSAKURA_SYMBOL(BaselineType_kPolynomial),
			LIBSAKURA_SYMBOL(BaselineType_kChebyshev)
	};

	for(size_t i=0; i < num_func_types ; ++i){
		LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
		Create(LIBSAKURA_SYMBOL(Status_kOK), func_types[i], order, num_chan, &context);
		Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
	}

}


/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using boundaries consisting of npiece=10 and num_data=30)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece10_Numdata30) {
	uint16_t const npiece(10);
	size_t const num_data(30);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kOK),LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), npiece, num_data, &context);
	Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}


/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using boundaries consisting of npiece=1LB and num_data=30)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece1LB_Numdata30) {
	uint16_t const npiece(1);
	size_t const num_data(30);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kOK),LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), npiece, num_data, &context);
	Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}


/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using boundaries consisting of npiece=1LB and num_data=4LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece1LB_Numdata4LB) {
	uint16_t const npiece(1);
	size_t const num_data(4);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kOK),LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), npiece, num_data, &context);
	Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}


/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using npiece=65535(max))
 */
/*
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece65535UB_Num_data65538LB) {
	uint16_t const npiece(65535);
	size_t const num_data(65538);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kOK),LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), npiece, num_data, &context);
	Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}
*/

/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using npiece=10 and num_data=13LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece10IR_Numdata13LB) {
	uint16_t const npiece(10);
	size_t const num_data(13);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kOK),LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), npiece, num_data, &context);
	Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}


/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * failure case (with cubicspline model using npiece=10IR and num_data=12OR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSplineUsingNpiece10IR_Numdata12OR) {
	uint16_t const npiece(10);
	size_t const num_data(12);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kInvalidArgument),LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), npiece, num_data, &context);
	//Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}


/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * failure case (with cubicspline model using context=NULL)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_ContextNULL) {
	uint16_t const npiece(10);
	size_t const num_data(30);

	LIBSAKURA_SYMBOL(BaselineContextFloat) ** context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kInvalidArgument),LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), npiece, num_data, context);
	//Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}


/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of nwave=10IR and num_data=30IR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave10IR_Numdata30IR) {
	uint16_t const nwave(10);
	size_t const num_data(30);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kOK),LIBSAKURA_SYMBOL(BaselineType_kSinusoid), nwave, num_data, &context);
	Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of nwave=0LB and num_data=30IR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave0LB_Numdata30IR) {
	uint16_t const nwave(0);
	size_t const num_data(30);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kOK),LIBSAKURA_SYMBOL(BaselineType_kSinusoid), nwave, num_data, &context);
	Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}


/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of order=0LB and num_chan=1LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave0LB_Numdata1LB) {
	uint16_t const nwave(0);
	size_t const num_data(1);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kOK),LIBSAKURA_SYMBOL(BaselineType_kSinusoid), nwave, num_data, &context);
	Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of nwave=10IR and num_data=21LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave10IR_Numdata21LB) {
	uint16_t const nwave(10);
	size_t const num_data(21);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kOK),LIBSAKURA_SYMBOL(BaselineType_kSinusoid), nwave, num_data, &context);
	Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}


/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * failure case (with sinusoid model using values consisting of nwave=10IR and num_data=20OR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave10IR_Numdata20OR) {
	uint16_t const nwave(10);
	size_t const num_data(20);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kInvalidArgument),LIBSAKURA_SYMBOL(BaselineType_kSinusoid), nwave, num_data, &context);
	//Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}


/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * failure case (with sinusoid model using context=NULL)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_NULLcontext) {
	uint16_t const nwave(10);
	size_t const num_data(30);

	LIBSAKURA_SYMBOL(BaselineContextFloat) **context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kInvalidArgument),LIBSAKURA_SYMBOL(BaselineType_kSinusoid), nwave, num_data, context);
	//Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}



////////////////////////////////////////////////////////////////////////////////
// TEST (Destroy)
////////////////////////////////////////////////////////////////////////////////

/*
 * Test sakura_DestroyBaselineContextFloat
 * successful case
 */
TEST_F(Baseline, DestroyBaselineContextFloat) {
	uint16_t const order(20);
	size_t const num_chan(4096);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Create(LIBSAKURA_SYMBOL(Status_kOK),LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), order, num_chan, &context);
	Destroy(LIBSAKURA_SYMBOL(Status_kOK), context);
}

/*
 * Test sakura_DestroyBaselineContextFloat
 * failure case : context is a null pointer
 * returned value must be Status_kInvalidArgument
 */
TEST_F(Baseline, DestroyBaselineContextFloatWithContextNullPointer) {
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	Destroy(LIBSAKURA_SYMBOL(Status_kInvalidArgument), context);
}



