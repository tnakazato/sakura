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

void Destroy(LIBSAKURA_SYMBOL (Status) status,
LIBSAKURA_SYMBOL(BaselineContextFloat) *context) {
	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyBaselineContextFloat(context);
	EXPECT_EQ(status, destroy_status);
}

void Create(LIBSAKURA_SYMBOL (Status) status,
		BaselineTypeInternal const mybaseline_type, int16_t const order,
		size_t const num_data, struct sakura_BaselineContextFloat** context) {
	LIBSAKURA_SYMBOL (Status) create_status;

	//cout << "context " << context << endl;
	//cout << "*context "<< *context << endl;
	cout << "&contest " << &context << endl;

	if (mybaseline_type == BaselineTypeInternal_kPolynomial) {
		create_status = sakura_CreateBaselineContextFloat(
				LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order, num_data,
				context);
	} else if (mybaseline_type == BaselineTypeInternal_kChebyshev) {
		create_status = sakura_CreateBaselineContextFloat(
				LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data,
				context);
	} else if (mybaseline_type == BaselineTypeInternal_kCubicSpline) {
		create_status = sakura_CreateBaselineContextCubicSplineFloat(order,
				num_data, context);
	} else if (mybaseline_type == BaselineTypeInternal_kSinusoid) {
		create_status = sakura_CreateBaselineContextSinusoidFloat(order,
				num_data, context);
	}

	EXPECT_EQ(status, create_status);
}

void Create2(LIBSAKURA_SYMBOL (Status) status,
		BaselineTypeInternal const mybaseline_type, int16_t const order,
		size_t const num_data, struct sakura_BaselineContextFloat** context) {
	LIBSAKURA_SYMBOL (Status) create_status;


	if (mybaseline_type == BaselineTypeInternal_kPolynomial) {
		create_status = sakura_CreateBaselineContextFloat(
				LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order, num_data,
				context);
	} else if (mybaseline_type == BaselineTypeInternal_kChebyshev) {
		create_status = sakura_CreateBaselineContextFloat(
				LIBSAKURA_SYMBOL(BaselineType_kChebyshev), order, num_data,
				context);
	} else if (mybaseline_type == BaselineTypeInternal_kCubicSpline) {
		create_status = sakura_CreateBaselineContextCubicSplineFloat(order,
				num_data, context);
	} else if (mybaseline_type == BaselineTypeInternal_kSinusoid) {
		create_status = sakura_CreateBaselineContextSinusoidFloat(order,
				num_data, context);
	}

	EXPECT_EQ(status, create_status);
}


struct FitExecute{
	static void execute(){ cout << "FitExecute" << endl;}
};

struct NullExecute{
	static void execute(){cout << "NullExecute0"<< endl;}

	static void execute(LIBSAKURA_SYMBOL (Status) status, BaselineTypeInternal const mybaseline_type,
			int16_t const order, size_t const num_data, struct sakura_BaselineContextFloat** context){
			cout << "NullExecute1"<< endl;
		}

	static void execute(LIBSAKURA_SYMBOL (Status) status, struct sakura_BaselineContextFloat** context){
				cout << "NullExecute2"<< endl;
			}

	static void execute(LIBSAKURA_SYMBOL (Status) status, struct sakura_BaselineContextFloat* context){
			cout << "NullExecute3" << endl;
	}

};

struct CreatExecute{
	static void execute(LIBSAKURA_SYMBOL (Status) status, BaselineTypeInternal const mybaseline_type,
			int16_t const order, size_t const num_data, struct sakura_BaselineContextFloat** context){
		Create(status,mybaseline_type,order,num_data, context);
		cout << "CreatExecute" << endl;

	}
};

struct CreatExecute2{
	static void execute(LIBSAKURA_SYMBOL (Status) status, BaselineTypeInternal const mybaseline_type,
			int16_t const order, size_t const num_data, struct sakura_BaselineContextFloat** context){
		Create2(status,mybaseline_type,order,num_data, context);
		cout << "CreatExecute2" << endl;

	}
};



struct DestroyExecute{
	static void execute(LIBSAKURA_SYMBOL (Status) status, struct sakura_BaselineContextFloat* context){
		Destroy(status, context);
		cout << "DestroyExecute" << endl;
	}

	//static void execute(){}


};

template<class T_creator, class T_fitter, class T_destroyer>
void TestRun(LIBSAKURA_SYMBOL (Status) status,BaselineTypeInternal const mybaseline_type,
		int16_t const order,size_t const num_data,
		struct sakura_BaselineContextFloat** context){

	T_creator::execute(status,mybaseline_type,order,num_data, context);
	T_fitter::execute();
	T_destroyer::execute(status, *context);

/*
	if(context==nullptr){
		T_destroyer::execute(status, context);
	}else{
		T_destroyer::execute(status, *context);
	}
*/

}


template<class T_creator, class T_fitter, class T_destroyer>
void TestRun2(LIBSAKURA_SYMBOL (Status) status,BaselineTypeInternal const mybaseline_type,
		int16_t const order,size_t const num_data,
		struct sakura_BaselineContextFloat** context){

	T_creator::execute(status,mybaseline_type,order,num_data, context);
	T_fitter::execute();
	//T_destroyer::execute(status, *context);



	if(context==nullptr){
		T_destroyer::execute(status, context);
		//T_destroyer::execute();
	}else{
		T_destroyer::execute(status, *context);
	}



}






////////////////////////////////////////////////////////////////////////////////
// TEST: Create & Destroy
////////////////////////////////////////////////////////////////////////////////

/*
 * Test sakura_CreateBaselineContextFloatWithPolynomial/ChebyshevPolynomial
 * successful case (with polynomial/chebyshevpolynomial model using boundaries consisting of order=0 and num_data=1)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithPolynomial_ChebyshevPolynomial) {
	uint16_t const order(0);
	size_t const num_data(1);
	size_t const num_func_types = 2;
	BaselineTypeInternal func_types[num_func_types] = {
			BaselineTypeInternal_kPolynomial, BaselineTypeInternal_kChebyshev };

	for (size_t i = 0; i < num_func_types; ++i) {
		LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
		TestRun<CreatExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),func_types[i],order,num_data,&context);
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
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data,&context);
}

/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using boundaries consisting of npiece=1LB and num_data=30)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece1LB_Numdata30) {
	uint16_t const npiece(1);
	size_t const num_data(30);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data,&context);
}

/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using boundaries consisting of npiece=1LB and num_data=4LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece1LB_Numdata4LB) {
	uint16_t const npiece(1);
	size_t const num_data(4);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data,&context);
}

/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using npiece=10 and num_data=13LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece10IR_Numdata13LB) {
	uint16_t const npiece(10);
	size_t const num_data(13);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data,&context);
}


////////////////////////////////////////////////////////////////////////////////
// TEST: Create
////////////////////////////////////////////////////////////////////////////////


/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * failure case (with cubicspline model using npiece=10IR and num_data=12OR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSplineUsingNpiece10IR_Numdata12OR) {
	uint16_t const npiece(10);
	size_t const num_data(12);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, NullExecute>
	(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kCubicSpline,npiece,num_data,&context);
}

/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * failure case (with cubicspline model using context=NULL)

TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_ContextNULL) {
	uint16_t const npiece(10);
	size_t const num_data(30);

	LIBSAKURA_SYMBOL(BaselineContextFloat) ** context = nullptr;
	TestRun2<CreatExecute, NullExecute, NullExecute>
	(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kCubicSpline,npiece,num_data,context);




	//LIBSAKURA_SYMBOL(BaselineContextFloat) ** context = nullptr;
	//cout <<"####### context " << context << endl;
	//LIBSAKURA_SYMBOL (Status) create_status = sakura_CreateBaselineContextCubicSplineFloat(npiece,
	//				num_data, context);
	//EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), create_status);

}
*/



/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of nwave=10IR and num_data=30IR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave10IR_Numdata30IR) {
	uint16_t const nwave(10);
	size_t const num_data(30);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data,&context);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of nwave=0LB and num_data=30IR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave0LB_Numdata30IR) {
	uint16_t const nwave(0);
	size_t const num_data(30);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data,&context);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of order=0LB and num_chan=2LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave0LB_Numdata2LB) {
	uint16_t const nwave(0);
	size_t const num_data(2);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data,&context);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of nwave=10IR and num_data=22LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave10IR_Numdata22LB) {
	uint16_t const nwave(10);
	size_t const num_data(22);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data,&context);

}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * failure case (with sinusoid model using values consisting of nwave=10IR and num_data=21OR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave10IR_Numdata21OR) {
	uint16_t const nwave(10);
	size_t const num_data(21);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, NullExecute>
	(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kSinusoid,nwave,num_data,&context);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * failure case (with sinusoid model using context=NULL)

TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_NULLcontext) {
	uint16_t const nwave(10);
	size_t const num_data(30);

	LIBSAKURA_SYMBOL(BaselineContextFloat) **context = nullptr;
	TestRun2<CreatExecute, NullExecute, NullExecute>
	(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kSinusoid,nwave,num_data,context);

}
*/
////////////////////////////////////////////////////////////////////////////////
// TEST: Destroy
////////////////////////////////////////////////////////////////////////////////

/*
 * Test sakura_DestroyBaselineContextFloat
 * successful case
 */
TEST_F(Baseline, DestroyBaselineContextFloat) {
	uint16_t const order(20);
	size_t const num_data(4096);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,order,num_data,&context);
}

/*
 * Test sakura_DestroyBaselineContextFloat
 * failure case : context is a null pointer
 * returned value must be Status_kInvalidArgument
 */
TEST_F(Baseline, DestroyBaselineContextFloatWithContextNullPointer) {
	uint16_t const order(20);
	size_t const num_data(4096);

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<NullExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kCubicSpline,order,num_data,&context);
}

