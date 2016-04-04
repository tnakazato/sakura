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


// Check if the expected and actual values are enough close to each other
void CheckAlmostEqual(double expected, double actual, double tolerance) {
	double deviation = fabs(actual - expected);
	double val = max(fabs(actual), fabs(expected)) * tolerance + tolerance;
	ASSERT_LE(deviation, val);
}


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
	}


	EXPECT_EQ(status, create_status);

}


struct FitExecute{
	static void execute(){ cout << "FitExecute" << endl;}

	static void execute(BaselineTypeInternal const mybaseline_type, struct sakura_BaselineContextFloat const * context,
				size_t order, size_t const * nwave, size_t num_data, float const data[],
				bool const mask[], float clip_threshold_sigma, uint16_t num_fitting_max,
				size_t num_coeff, double coeff[], float best_fit[], float residual[],
				bool final_mask[],
				size_t const * boundary,
				double *coeff_answer_ptr){

		sakura_BaselineStatus  baseline_status;
		float rms;
		double *coeff_ptr = coeff;
		float *best_fit_ptr = best_fit;
		float *residual_ptr = residual;
		enum NPCases {NP_kNo, NP_kCoeff, NP_kBestFit, NP_kResidual, NP_kAll, NP_kNumElems};
		vector<string> np_cases_names = { "no nullptr", "coeff=nullptr",
								"best_fit=nullptr", "residual=nullptr", "all nullptr" };

		if(mybaseline_type==BaselineTypeInternal_kPolynomial){

			cout << "    Testing for ";

			for (NPCases item = static_cast<NPCases>(0); item < NP_kNumElems; item =
						static_cast<NPCases>(item + 1)) {
					cout << np_cases_names[item] << ((item < NP_kNumElems - 1) ? ", " : "");

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
				LIBSAKURA_SYMBOL (Status) fit_status =
								LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(context,
								order, num_data, data,
								mask, clip_threshold_sigma, num_fitting_max,
								num_coeff, coeff_ptr, best_fit_ptr, residual_ptr,
								final_mask,
								&rms,
								&baseline_status);

				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), fit_status);
				EXPECT_EQ(LIBSAKURA_SYMBOL(BaselineStatus_kOK), baseline_status);

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
					for (size_t i = 0; i < ELEMENTSOF(coeff_answer_ptr); ++i) {
						CheckAlmostEqual(coeff_answer_ptr[i], coeff[i], 1.0e-6);
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
			}//for

		}






		//todo
		if(mybaseline_type==BaselineTypeInternal_kSinusoid){
				/*
				cout << "#############order "<< order << endl;
				for(int i=0; i < 4; ++i){
						cout << "#############nwave " << *nwave << endl;
						++nwave;
				}
				cout << "#############num_data " << num_data << endl;
				*/

				LIBSAKURA_SYMBOL (Status) fit_status =
					LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(context,
							order+1, nwave, num_data, data,
							mask, clip_threshold_sigma, num_fitting_max,
							num_coeff, coeff, best_fit, residual,
							final_mask,
							&rms,
							&baseline_status);
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), fit_status);
				cout << "LSQFitSinusoidFloat" << endl;
		}

		//todo
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

struct CreateExecute{
	static void execute(LIBSAKURA_SYMBOL (Status) status, BaselineTypeInternal const mybaseline_type,
			int16_t const order, size_t const num_data, struct sakura_BaselineContextFloat** context){

			Create(status,mybaseline_type,order,num_data, context);
		cout << "CreatExecute" << endl;
	}
};


struct CreateExecuteNuLL{
	static void execute(LIBSAKURA_SYMBOL (Status) status, BaselineTypeInternal const mybaseline_type,
			int16_t const order, size_t const num_data, struct sakura_BaselineContextFloat** context){

			Create(status,mybaseline_type,order,num_data, nullptr);
			cout << "CreatExecuteNuLL" << endl;

	};
};





struct DestroyExecute{
	static void execute(LIBSAKURA_SYMBOL (Status) status, struct sakura_BaselineContextFloat* context){
		Destroy(status, context);
		cout << "DestroyExecute" << endl;
	}
};

struct DestroyExecuteNuLL{
	static void execute(LIBSAKURA_SYMBOL (Status) status,struct sakura_BaselineContextFloat* context){
		Destroy(status, nullptr);
		cout << "DestroyExecuteNuLL" << endl;
	}
};

/*
struct InitializeExecute{
		static void execute(LIBSAKURA_SYMBOL (Status) &status,
					BaselineTypeInternal &mybaseline_type,
					int16_t  *order,
					size_t  * nwave, //sinusoid
					float * clip_threshold_sigma,
					uint16_t * num_fitting_max,
					size_t * num_coeff,
					size_t  * boundary //cubic
					)

					{

					*order=3;
					*num_fitting_max=2;
					*num_coeff=*order+1;
					*clip_threshold_sigma =3;

			cout << "InitializerExecute" << endl;

		}

};
*/

struct InitializePoly{
		static void execute(LIBSAKURA_SYMBOL (Status) &status,
					BaselineTypeInternal &mybaseline_type,
					int16_t  *order,
					size_t  * nwave, //sinusoid
					float * clip_threshold_sigma,
					uint16_t * num_fitting_max,
					size_t * num_coeff,
					size_t  * boundary, //cubic
					size_t * num_coeff_answer
					){

					*order=3;
					*num_fitting_max=2;
					*num_coeff=*order+1;
					*clip_threshold_sigma =3;
					*num_coeff_answer=*order+1;

			cout << "InitializePolyExecute" << endl;

		}

};

/*
struct InitializeSinusoid{
		static void execute(LIBSAKURA_SYMBOL (Status) &status,
					BaselineTypeInternal &mybaseline_type,
					int16_t  *nwave_max,
					size_t  * nwave, //sinusoid
					float * clip_threshold_sigma,
					uint16_t * num_fitting_max,
					size_t * num_coeff,
					size_t  * boundary, //cubic
					size_t * num_coeff_answer
					){

					*nwave_max = 3;
					*num_nwave=*nwave_max + 1;

					for(size_t i=0; i < 4;++i){
						nwave[i]=i;
					}


					*num_fitting_max=2;
					*num_coeff=*nwave_max+1;
					*clip_threshold_sigma =3;
					*num_coeff_answer=*nwave_max * 2 + 1;

			cout << "InitializeSinusoid" << endl;

		}

};
*/





template<class T_creator, class T_fitter, class T_destroyer>
void TestRun0(LIBSAKURA_SYMBOL (Status) status,BaselineTypeInternal const mybaseline_type,
		int16_t const order,size_t const num_data){

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;

	T_creator::execute(status,mybaseline_type,order,num_data, &context);
	T_fitter::execute();
	if(context!=nullptr){
		T_destroyer::execute(status, context);
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


	//todo this is dummy (added this to be able to run T_fitter::execute(...)below
	double coeff_answer[7];


	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;

	T_creator::execute(status,mybaseline_type,order,num_data, &context);
	T_fitter::execute();
	if(context!=nullptr){
		T_fitter::execute(mybaseline_type, context,order, nwave, num_data, data_ptr, mask_ptr, clip_threshold_sigma,
					num_fitting_max,num_coeff, coeff_ptr, best_fit_ptr, residual_ptr,
					final_mask,
					//rms_ptr,
					boundary,
					//baseline_status,
					coeff_answer);
		T_destroyer::execute(status, context);
	}
}


//todo (TestRun2->TestRun3)
template<class T_initializer, class T_creator, class T_fitter, class T_destroyer>
void TestRun3(LIBSAKURA_SYMBOL(Status)  status, BaselineTypeInternal mybaseline_type
		){

	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;

	uint16_t num_fitting_max;
	float clip_threshold_sigma;
	//size_t  *nwave=nullptr; //sinusoid
	//bool * final_mask=nullptr;
	size_t * boundary=nullptr;//cubic
	int16_t order;
	size_t num_coeff;

	//sinusoid
	size_t nwave_max = 3; //=order
	size_t const num_nwave = nwave_max + 1; //=order+1

	//size_t const nwave[num_nwave] = { 0, 1, 2, 3 };
	//size_t  nwave[num_nwave] = { 0, 1, 2, 3 };
	size_t nwave[num_nwave];
	size_t num_coeff_answer;



	T_initializer::execute(status,
			mybaseline_type,
			&order,
			nwave,
			&clip_threshold_sigma,
			&num_fitting_max,
			&num_coeff,
			boundary,
			&num_coeff_answer
			);



	SIMD_ALIGN
	//double coeff_answer[order+1];
	double coeff_answer[num_coeff_answer];
	SetDoubleConstant(1.0, ELEMENTSOF(coeff_answer), coeff_answer);
	SIMD_ALIGN
	double coeff[ELEMENTSOF(coeff_answer)];

	size_t  num_data=ELEMENTSOF(coeff_answer);
	if (mybaseline_type==BaselineTypeInternal_kSinusoid){
		num_data = ELEMENTSOF(coeff_answer) + 5;
	}

	SIMD_ALIGN
	float best_fit[num_data];
	SIMD_ALIGN
	float residual[num_data];
	SIMD_ALIGN
	float data[num_data];

	if(mybaseline_type==BaselineTypeInternal_kPolynomial){
		SetFloatPolynomial(ELEMENTSOF(coeff_answer), coeff_answer, num_data, data);
	}else if(mybaseline_type==BaselineTypeInternal_kSinusoid){
		SetFloatSinusoidal(num_nwave, nwave, coeff_answer, num_data, data);


	}

	SIMD_ALIGN
	bool mask[ELEMENTSOF(data)];
	SetBoolConstant(true, ELEMENTSOF(data), mask);

	T_creator::execute(status,mybaseline_type,order,num_data, &context);
	T_fitter::execute();
	if(context!=nullptr){
		T_fitter::execute(mybaseline_type, context,order, nwave, num_data, data, mask, clip_threshold_sigma,
							num_fitting_max, num_coeff, coeff, best_fit, residual,
							mask, boundary, coeff_answer);
		T_destroyer::execute(status, context);
	}
}













/*
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
*/





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
		TestRun0<CreateExecute, NullExecute, DestroyExecute>
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

	TestRun0<CreateExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data);


}

/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using boundaries consisting of npiece=1LB and num_data=30)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece1LB_Numdata30) {
	uint16_t const npiece(1);
	size_t const num_data(30);

	TestRun0<CreateExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data);

}

/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using boundaries consisting of npiece=1LB and num_data=4LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece1LB_Numdata4LB) {
	uint16_t const npiece(1);
	size_t const num_data(4);

	TestRun0<CreateExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kCubicSpline,npiece,num_data);

}

/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * successful case (with cubicspline model using npiece=10 and num_data=13LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_Npiece10IR_Numdata13LB) {
	uint16_t const npiece(10);
	size_t const num_data(13);

	TestRun0<CreateExecute, NullExecute, DestroyExecute>
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

	TestRun0<CreateExecute, NullExecute, NullExecute>
		(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kCubicSpline,npiece,num_data);
}


/*
 * Test sakura_CreateBaselineContextFloatWithCubicSpline
 * failure case (with cubicspline model using context=NULL)
*/
TEST_F(Baseline, CreateBaselineContextFloatWithCubicSpline_ContextNULL) {
	uint16_t const npiece(10);
	size_t const num_data(30);

	TestRun0<CreateExecuteNuLL, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kCubicSpline,npiece,num_data);
}




/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of nwave=10IR and num_data=30IR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave10IR_Numdata30IR) {
	uint16_t const nwave(10);
	size_t const num_data(30);

	TestRun0<CreateExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of nwave=0LB and num_data=30IR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave0LB_Numdata30IR) {
	uint16_t const nwave(0);
	size_t const num_data(30);

TestRun0<CreateExecute, NullExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of order=0LB and num_chan=2LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave0LB_Numdata2LB) {
	uint16_t const nwave(0);
	size_t const num_data(2);

	TestRun0<CreateExecute, NullExecute, DestroyExecute>
		(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * successful case (with sinusoid model using boundaries consisting of nwave=10IR and num_data=22LB)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave10IR_Numdata22LB) {
	uint16_t const nwave(10);
	size_t const num_data(22);
	TestRun0<CreateExecute, NullExecute, DestroyExecute>
			(LIBSAKURA_SYMBOL(Status_kOK),BaselineTypeInternal_kSinusoid,nwave,num_data);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * failure case (with sinusoid model using values consisting of nwave=10IR and num_data=21OR)
 */
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_Nwave10IR_Numdata21OR) {
	uint16_t const nwave(10);
	size_t const num_data(21);

	TestRun0<CreateExecute, NullExecute, NullExecute>
		(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kSinusoid,nwave,num_data);
}

/*
 * Test sakura_CreateBaselineContextFloatWithSinusoid
 * failure case (with sinusoid model using context=NULL)
*/
TEST_F(Baseline, CreateBaselineContextFloatWithSinusoid_NULLcontext) {
	uint16_t const nwave(10);
	size_t const num_data(30);

	TestRun0<CreateExecuteNuLL, NullExecute, DestroyExecute>
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

	TestRun0<CreateExecute, NullExecute, DestroyExecute>
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

	TestRun0<NullExecute, NullExecute, DestroyExecuteNuLL>
		(LIBSAKURA_SYMBOL(Status_kInvalidArgument),BaselineTypeInternal_kCubicSpline,order,num_data);
}


////////////////////////////////////////////////////////////////////////////////
// TEST: Fit
////////////////////////////////////////////////////////////////////////////////
/*
 * Test LSQFitPolynomial
 * successful case using TestRun3
 */
TEST_F(Baseline, LSQFitPolynomialSuccessfulCases_TestRun3) {

	TestRun3<InitializePoly, CreateExecute, FitExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK), BaselineTypeInternal_kPolynomial);

}

/*
 * Test LSQFitSinusoid
 * successful case using TestRun3
 */
/*
TEST_F(Baseline, LSQFitSinusoidSuccessfulCases_TestRun3) {

	TestRun3<InitializeSinusoid, CreateExecute, FitExecute, DestroyExecute>
	(LIBSAKURA_SYMBOL(Status_kOK), BaselineTypeInternal_kSinusoid);
}
*/




/*
 * Test LSQFitPolynomial
 * successful case
 */










/*
 * Test LSQFitSinusoid
 * successful case
 */
/*
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

		size_t const * dummy=nullptr;

		TestRun2<CreateExecute, FitExecute, DestroyExecute>
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





		//TestRun3<InitializeSinusoid, CreateExecute, FitExecute,DestroyExecute>
		//	(LIBSAKURA_SYMBOL(Status_kOK), BaselineTypeInternal_kSinusoid);

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
*/


