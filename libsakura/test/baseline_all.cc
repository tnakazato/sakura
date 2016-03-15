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
		cout <<"BaselineTypeInternal_kSinusoid in Create" << endl;
	}


	EXPECT_EQ(status, create_status);

}


struct FitExecute{
	static void execute(){ cout << "FitExecute" << endl;}

	static void execute(BaselineTypeInternal const mybaseline_type, struct sakura_BaselineContextFloat const * context,
				size_t order, size_t const * nwave, size_t num_data, float const data[],
				bool const mask[], float clip_threshold_sigma, uint16_t num_fitting_max,
				size_t num_coeff, double coeff_ptr[], float best_fit_ptr[], float residual_ptr[],
				bool final_mask[], float * rms, size_t const * boundary, sakura_BaselineStatus  *baseline_status){

		if(mybaseline_type==BaselineTypeInternal_kPolynomial){
			LIBSAKURA_SYMBOL (Status) fit_status =
				LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context,
						order, num_data, data,
						mask, clip_threshold_sigma, num_fitting_max,
						num_coeff, coeff_ptr, best_fit_ptr, residual_ptr,
						final_mask,
						rms,
						baseline_status);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), fit_status);
			cout << "LSQFitPolynomialFloat" << endl;
		}

		if(mybaseline_type==BaselineTypeInternal_kSinusoid){
				LIBSAKURA_SYMBOL (Status) fit_status =
					LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context,
							order, nwave, num_data, data,
							mask, clip_threshold_sigma, num_fitting_max,
							num_coeff, coeff_ptr, best_fit_ptr, residual_ptr,
							final_mask,
							rms,
							baseline_status);
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), fit_status);
				cout << "LSQFitSinusoidFloat" << endl;
		}


		if(mybaseline_type==BaselineTypeInternal_kCubicSpline){


			cout << "LSQFitCubicSplineFloat" << endl;
		}

	}

};//Fit Execute (struct)

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


struct CreatExecuteNuLL{
	static void execute(LIBSAKURA_SYMBOL (Status) status, BaselineTypeInternal const mybaseline_type,
			int16_t const order, size_t const num_data, struct sakura_BaselineContextFloat** context){

		//if(context==nullptr){
		//	Create(status,mybaseline_type,order,num_data, nullptr);
		//}else{
			Create(status,mybaseline_type,order,num_data, nullptr);
		//cout << "CreatExecute" << endl;
		//}
		cout << "CreatExecuteNuLL" << endl;

	};
};



/*
struct CreatExecute{
	static void execute(LIBSAKURA_SYMBOL (Status) status, BaselineTypeInternal const mybaseline_type,
			int16_t const order, size_t const num_data, struct sakura_BaselineContextFloat* context){

		//if(context==nullptr){
		//	Create(status,mybaseline_type,order,num_data, nullptr);
		//}else{
			Create(status,mybaseline_type,order,num_data, context);
		//cout << "CreatExecute" << endl;
		//}
		cout << "CreatExecute" << endl;
	}
};
*/


struct DestroyExecute{
	static void execute(LIBSAKURA_SYMBOL (Status) status, struct sakura_BaselineContextFloat* context){
		Destroy(status, context);
		cout << "DestroyExecute" << endl;
	}
};

struct DestroyExecuteNuLL{
	static void execute(LIBSAKURA_SYMBOL (Status) status,struct sakura_BaselineContextFloat* context){
		Destroy(status, nullptr);
		cout << "DestroyExecute" << endl;
	}
};


struct ContextStatusNULL{
	static void execute(struct sakura_BaselineContextFloat* context){
		context=nullptr;
		cout << "ContextStatusNULL" << endl;

	}

};

struct ContextStatusNotNULL{
	static void execute(struct sakura_BaselineContextFloat* context){
	cout << "ContextStatusNotNULL" << endl;
	}
};



template<class T_creator, class T_fitter, class T_destroyer>
void TestRun(LIBSAKURA_SYMBOL (Status) status,BaselineTypeInternal const mybaseline_type,
		int16_t const order,size_t const num_data,
		struct sakura_BaselineContextFloat** context){


	//LIBSAKURA_SYMBOL(BaselineContextFloat) * context2 = nullptr;
	//LIBSAKURA_SYMBOL(BaselineContextFloat) ** context2 = nullptr;


	T_creator::execute(status,mybaseline_type,order,num_data, context);
	T_fitter::execute();
	if(context!=nullptr){
		T_destroyer::execute(status, *context);
		//T_destroyer::execute(status, context2);
	}

}

template<class T_creator, class T_fitter, class T_destroyer>
void TestRun0(LIBSAKURA_SYMBOL (Status) status,BaselineTypeInternal const mybaseline_type,
		int16_t const order,size_t const num_data){

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;

	T_creator::execute(status,mybaseline_type,order,num_data, &context);
	//T_creator::execute(status,mybaseline_type,order,num_data, context);


	T_fitter::execute();
	if(context!=nullptr){
		T_destroyer::execute(status, context);
		//T_destroyer::execute(status, context2);
	}

}




template<class T_creator, class T_fitter, class T_destroyer>
void TestRun2(LIBSAKURA_SYMBOL (Status) status,
		BaselineTypeInternal const mybaseline_type,
		int16_t const order,
		size_t const * nwave,//sinusoid
		size_t const num_data,
		float *data_ptr,
		bool *mask_ptr,
		float clip_threshold_sigma,
		uint16_t num_fitting_max,
		size_t num_coeff,
		double * coeff_ptr,
		float * best_fit_ptr,
		float * residual_ptr,
		bool * final_mask,
		float * rms_ptr,
		size_t const * boundary, //cubic
		LIBSAKURA_SYMBOL(BaselineStatus) * baseline_status
		){

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;

	T_creator::execute(status,mybaseline_type,order,num_data, &context);
	T_fitter::execute();
	if(context!=nullptr){
		T_fitter::execute(mybaseline_type, context,order, nwave, num_data, data_ptr, mask_ptr, clip_threshold_sigma,
					num_fitting_max,num_coeff, coeff_ptr, best_fit_ptr, residual_ptr,
					final_mask,rms_ptr,boundary, baseline_status);
		T_destroyer::execute(status, context);
	}
}


template<class T_creator, class T_fitter, class T_destroyer>
void TestRun3(LIBSAKURA_SYMBOL (Status) status,
		BaselineTypeInternal const mybaseline_type,
		int16_t const order,
		size_t const * nwave,//sinusoid
		size_t const num_data,
		float *data_ptr,
		bool *mask_ptr,
		float clip_threshold_sigma,
		uint16_t num_fitting_max,
		size_t num_coeff,
		double * coeff_ptr,
		float * best_fit_ptr,
		float * residual_ptr,
		bool * final_mask,
		float * rms_ptr,
		size_t const * boundary, //cubic
		LIBSAKURA_SYMBOL(BaselineStatus) * baseline_status
		){

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;


	T_creator::execute(status,mybaseline_type,order,num_data, &context);
	T_fitter::execute();
	if(context!=nullptr){
		T_fitter::execute(mybaseline_type, *context,order, nwave, num_data, data_ptr, mask_ptr, clip_threshold_sigma,
					num_fitting_max,num_coeff, coeff_ptr, best_fit_ptr, residual_ptr,
					final_mask,rms_ptr,boundary, baseline_status);
		T_destroyer::execute(status, context);
	}
}









//Set (coeff[0]+coeff[1]*x+coeff[2]*x*x+...) float values into an array
void SetFloatPolynomial(size_t num_coeff, double const *coeff,
		size_t num_data, float *data) {
	for (size_t i = 0; i < num_data; ++i) {
		double val = 0.0;
		double x = (double) i;
		for (size_t j = 0; j < num_coeff; ++j) {
			val *= x;
			val += coeff[num_coeff - 1 - j];
		}
		data[i] = static_cast<float>(val);
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


// Check if the expected and actual values are enough close to each other
void CheckAlmostEqual(double expected, double actual, double tolerance) {
	double deviation = fabs(actual - expected);
	double val = max(fabs(actual), fabs(expected)) * tolerance + tolerance;
	ASSERT_LE(deviation, val);
}

//Set sinusoidal values of float into an array
void SetFloatSinusoidal(size_t num_nwave, size_t const *nwave,
		double const *coeff, size_t num_data, float *data) {
	double factor = 2.0 * M_PI / (double) (num_data - 1);
	for (size_t i = 0; i < num_data; ++i) {
		data[i] = 0.0f;
		double x = factor * (double) i;
		size_t coeff_idx = 0;
		for (size_t j = 0; j < num_nwave; ++j) {
			//amplitude of each sinusoid is unity
			if (nwave[j] == 0) {
				data[i] += coeff[coeff_idx++];
			} else {
				double theta = nwave[j] * x;
				data[i] += (float) (coeff[coeff_idx++] * sin(theta));
				data[i] += (float) (coeff[coeff_idx++] * cos(theta));
			}
		}
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
		/*
		LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
		TestRun<CreatExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),func_types[i],order,num_data,&context);
		*/

		TestRun0<CreatExecute, NullExecute, DestroyExecute>
				(LIBSAKURA_SYMBOL(Status_kOK),func_types[i],order,num_data);
	}
}

/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using boundaries consisting of npiece=10 and num_data=30)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece10_Numdata30) {
	uint16_t const npiece(10);
	size_t const num_data(30);

	/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data,&context);
*/

	TestRun0<CreatExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data);


}

/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using boundaries consisting of npiece=1LB and num_data=30)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece1LB_Numdata30) {
	uint16_t const npiece(1);
	size_t const num_data(30);

	/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data,&context);
*/

	TestRun0<CreatExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data);

}

/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using boundaries consisting of npiece=1LB and num_data=4LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece1LB_Numdata4LB) {
	uint16_t const npiece(1);
	size_t const num_data(4);

	/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data,&context);
*/

	TestRun0<CreatExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data);

}

/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using npiece=10 and num_data=13LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece10IR_Numdata13LB) {
	uint16_t const npiece(10);
	size_t const num_data(13);

	/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data,&context);
*/

	TestRun0<CreatExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data);

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

	/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, NullExecute>
	(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kCubicSpline,npiece,num_data,&context);
	*/

	TestRun0<CreatExecute, NullExecute, NullExecute>
		(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kCubicSpline,npiece,num_data);


}


/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * failure case (with cubicspline model using context=NULL)
*/
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_ContextNULL) {
	uint16_t const npiece(10);
	size_t const num_data(30);

	/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) ** context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kCubicSpline,npiece,num_data,context);
	*/
	TestRun0<CreatExecuteNuLL, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kCubicSpline,npiece,num_data);


}




/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of nwave=10IR and num_data=30IR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave10IR_Numdata30IR) {
	uint16_t const nwave(10);
	size_t const num_data(30);

	/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data,&context);
	*/

	TestRun0<CreatExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of nwave=0LB and num_data=30IR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave0LB_Numdata30IR) {
	uint16_t const nwave(0);
	size_t const num_data(30);

	/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data,&context);
}
*/

TestRun0<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of order=0LB and num_chan=2LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave0LB_Numdata2LB) {
	uint16_t const nwave(0);
	size_t const num_data(2);
/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data,&context);
*/
	TestRun0<CreatExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of nwave=10IR and num_data=22LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave10IR_Numdata22LB) {
	uint16_t const nwave(10);
	size_t const num_data(22);
/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data,&context);
*/

	TestRun0<CreatExecute, NullExecute, DestroyExecute>
			(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * failure case (with sinusoid model using values consisting of nwave=10IR and num_data=21OR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave10IR_Numdata21OR) {
	uint16_t const nwave(10);
	size_t const num_data(21);

	/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, NullExecute>
	(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kSinusoid,nwave,num_data,&context);
*/

	TestRun0<CreatExecute, NullExecute, NullExecute>
		(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kSinusoid,nwave,num_data);

}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * failure case (with sinusoid model using context=NULL)
*/
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_NULLcontext) {
	uint16_t const nwave(10);
	size_t const num_data(30);
/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) **context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kSinusoid,nwave,num_data,context);
*/
	TestRun0<CreatExecuteNuLL, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kSinusoid,nwave,num_data);
}

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
/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<CreatExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,order,num_data,&context);
*/

	TestRun0<CreatExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,order,num_data);
}

/*
 * Test sakura_DestroyBaselineContextFloat
 * failure case : context is a null pointer
 * returned value must be Status_kInvalidArgument
 */
TEST_F(Baseline, DestroyBaselineContextFloatWithContextNullPointer) {
	uint16_t const order(20);
	size_t const num_data(4096);
/*
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	TestRun<NullExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kCubicSpline,order,num_data,&context);
*/

	TestRun0<NullExecute, NullExecute, DestroyExecuteNuLL>
		(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kCubicSpline,order,num_data);

}


////////////////////////////////////////////////////////////////////////////////
// TEST: Fit
////////////////////////////////////////////////////////////////////////////////
/*
 * Test LSQFitPolynomial
 * successful case
 */

TEST_F(Baseline, LSQFitPolynomialSuccessfulCases) {
	enum NPCases {
		NP_kNo, NP_kCoeff, NP_kBestFit, NP_kResidual, NP_kAll, NP_kNumElems
	};
	vector<string> np_cases_names = { "no nullptr", "coeff=nullptr",
			"best_fit=nullptr", "residual=nullptr", "all nullptr" };
	cout << "    Testing for ";

	size_t const order = 3;
	SIMD_ALIGN
	double coeff_answer[order + 1];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff_answer), coeff_answer);
	SIMD_ALIGN
	double coeff[ELEMENTSOF(coeff_answer)];
	float rms;
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	size_t const num_data = ELEMENTSOF(coeff_answer);
	SIMD_ALIGN
	float data[num_data];
	SetFloatPolynomial(ELEMENTSOF(coeff_answer), coeff_answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	//if (verbose) {
	//	PrintArray("data", num_data, data);
	//}
	SIMD_ALIGN
	float best_fit[num_data];
	SIMD_ALIGN
	float residual[num_data];

	size_t const * dummy=nullptr;

	//LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;

	for (NPCases item = static_cast<NPCases>(0); item < NP_kNumElems; item =
			static_cast<NPCases>(item + 1)) {
		cout << np_cases_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		double *coeff_ptr = coeff;
		float *best_fit_ptr = best_fit;
		float *residual_ptr = residual;

		switch (item) {
		case NP_kNo:
			break;
		case NP_kCoeff:
			coeff_ptr = nullptr;
			break;
		case NP_kBestFit:
			best_fit_ptr = nullptr;
			break;
		case NP_kResidual:
			residual_ptr = nullptr;
			break;
		case NP_kAll:
			coeff_ptr = nullptr;
			best_fit_ptr = nullptr;
			residual_ptr = nullptr;
			break;
		default:
			assert(false);
		}


		TestRun2<CreatExecute, FitExecute, DestroyExecute>
			(LIBSAKURA_SYMBOL(Status_kOK),
					BaselineTypeInternal_kPolynomial,
					order,
					dummy,
					num_data,
					data,
					mask,
					5.0f,
					1,
					ELEMENTSOF(coeff),
					coeff_ptr,
					best_fit_ptr,
					residual_ptr,
					mask,
					&rms,
					dummy,
					&baseline_status
					);

		bool check_coeff = true;
		bool check_best_fit = true;
		bool check_residual = true;
		if (item == NP_kCoeff || item == NP_kAll) {
			check_coeff = false;
		}
		if (item == NP_kBestFit || item == NP_kAll) {
			check_best_fit = false;
		}
		if (item == NP_kResidual || item == NP_kAll) {
			check_residual = false;
		}
		if (check_coeff) {
			for (size_t i = 0; i < ELEMENTSOF(coeff_answer); ++i) {
				CheckAlmostEqual(coeff_answer[i], coeff[i], 1.0e-6);
			}
		}
		if (check_best_fit) {
			for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
				CheckAlmostEqual(data[i], best_fit[i], 1.0e-6);
			}
		}
		if (check_residual) {
			for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
				CheckAlmostEqual(0.0, residual[i], 1.0e-6);
			}
		}
	}


	cout << endl;
}

/*
 * Test LSQFitSinusoid
 * successful case
 */

TEST_F(Baseline, LSQFitSinusoidSuccessfulCases) {
	enum NPCases {
		NP_kNo, NP_kCoeff, NP_kBestFit, NP_kResidual, NP_kAll, NP_kNumElems
	};
	vector<string> np_cases_names = { "no nullptr", "coeff=nullptr",
			"best_fit=nullptr", "residual=nullptr", "all nullptr" };
	cout << "    Testing for ";

	size_t const nwave_max = 3;
	size_t const num_nwave = nwave_max + 1;
	size_t const nwave[num_nwave] = { 0, 1, 2, 3 };
	SIMD_ALIGN
	double coeff_answer[nwave_max * 2 + 1];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff_answer), coeff_answer);
	SIMD_ALIGN
	double coeff[ELEMENTSOF(coeff_answer)];
	float rms;
	LIBSAKURA_SYMBOL(BaselineStatus) baseline_status;

	size_t const num_data = ELEMENTSOF(coeff_answer) + 5;
	SIMD_ALIGN
	float data[num_data];
	SetFloatSinusoidal(num_nwave, nwave, coeff_answer, num_data, data);
	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);
	//if (verbose) {
	//	PrintArray("data", num_data, data);
	//}
	SIMD_ALIGN
	float best_fit[num_data];
	SIMD_ALIGN
	float residual[num_data];
	//LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
	//LIBSAKURA_SYMBOL (Status) create_status =
	//		sakura_CreateBaselineContextSinusoidFloat(nwave_max, num_data,
	//				&context);
	//EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
	for (NPCases item = static_cast<NPCases>(0); item < NP_kNumElems; item =
			static_cast<NPCases>(item + 1)) {
		cout << np_cases_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

		double *coeff_ptr = coeff;
		float *best_fit_ptr = best_fit;
		float *residual_ptr = residual;

		switch (item) {
		case NP_kNo:
			break;
		case NP_kCoeff:
			coeff_ptr = nullptr;
			break;
		case NP_kBestFit:
			best_fit_ptr = nullptr;
			break;
		case NP_kResidual:
			residual_ptr = nullptr;
			break;
		case NP_kAll:
			coeff_ptr = nullptr;
			best_fit_ptr = nullptr;
			residual_ptr = nullptr;
			break;
		default:
			assert(false);
		}




		//LIBSAKURA_SYMBOL (Status) fit_status =
		//LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context, num_nwave, nwave,
		//		num_data, data, mask, 5.0f, 1, ELEMENTSOF(coeff), coeff_ptr,
		//		best_fit_ptr, residual_ptr, mask, &rms, &baseline_status);
		//EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), fit_status);
		//EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);


		size_t const * dummy=nullptr;

		TestRun2<CreatExecute, FitExecute, DestroyExecute>
			(LIBSAKURA_SYMBOL(Status_kOK),
					BaselineTypeInternal_kSinusoid,
					num_nwave,
					nwave,
					num_data,
					data,
					mask,
					5.0f,
					1,
					ELEMENTSOF(coeff),
					coeff_ptr,
					best_fit_ptr,
					residual_ptr,
					mask,
					&rms,
					dummy,
					&baseline_status
					);












		bool check_coeff = true;
		bool check_best_fit = true;
		bool check_residual = true;
		if (item == NP_kCoeff || item == NP_kAll) {
			check_coeff = false;
		}
		if (item == NP_kBestFit || item == NP_kAll) {
			check_best_fit = false;
		}
		if (item == NP_kResidual || item == NP_kAll) {
			check_residual = false;
		}
		if (check_coeff) {
			for (size_t i = 0; i < ELEMENTSOF(coeff_answer); ++i) {
				CheckAlmostEqual(coeff_answer[i], coeff[i], 1.0e-6);
			}
		}
		if (check_best_fit) {
			for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
				CheckAlmostEqual(data[i], best_fit[i], 1.0e-6);
			}
		}
		if (check_residual) {
			for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
				CheckAlmostEqual(0.0, residual[i], 1.0e-6);
			}
		}
	}

	//LIBSAKURA_SYMBOL (Status) destroy_status =
	//		sakura_DestroyBaselineContextFloat(context);
	//EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), destroy_status);

	cout << endl;

}



