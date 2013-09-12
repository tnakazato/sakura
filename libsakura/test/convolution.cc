#include <iostream>
#include <string>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fftw3.h>
#include <libsakura/sakura.h>
#include "gtest/gtest.h"

/* the number of elements in input/output array to test */
#define NUM_CHANNEL 128
#define NUM_WIDTH 5
#define NUM_IN 11
#define NUM_CENTER 1

extern "C" {
struct LIBSAKURA_SYMBOL(Convolve1DContext) {
	fftwf_plan plan_real_to_complex_float;
	fftwf_plan plan_complex_to_real_float;
	fftwf_complex *input_complex_spectrum;
	fftwf_complex *output_complex_spectrum;
	fftwf_complex *fft_applied_complex_kernel;
	size_t num_channel;
	float input_real_array[0];
};
}

using namespace std;

/*
 * A super class to test creating kernel of array(s) width=3, channels=128
 * INPUTS:
 * - in1 = [0.000141569,0.00226511, 0.0195716, 0.0913234,0.230121, 0.313146,
 *          0.230121,0.0913234,0.0195716,0.00226511,0.000141569]
 */
class CreateConvolve1DContext : public ::testing::Test
{
protected:

	CreateConvolve1DContext () : verbose(true)
	{}

	virtual void SetUp()
	{
		in1_[0]  = 0.000141569;
		in1_[1]  = 0.00226511;
		in1_[2]  = 0.0195716;
		in1_[3]  = 0.0913234;
		in1_[4]  = 0.230121;
		in1_[5]  = 0.313146;
		in1_[6]  = 0.230121;
		in1_[7]  = 0.0913234;
		in1_[8]  = 0.0195716;
		in1_[9]  = 0.00226511;
		in1_[10] = 0.000141569;

		// Initialize sakura
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Initialize)();
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}

	virtual void TearDown()
	{
		// Clean-up sakura
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	void PrintInputs(){
		PrintArray("in1", NUM_IN, in1_);
	}

	void PrintArray(char const *name, size_t num_in, float *in){
		cout << name << " = [";
		for (size_t i = 0; i < num_in-1; i++)
			cout << in[i] << ", " ;
		cout << in[num_in-1];
		cout << " ]" << endl;
	}

	float in1_[NUM_IN];
	bool verbose;

};

/*
 * Test creating gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * gaussian kernel is symmetric
 */
TEST_F(CreateConvolve1DContext , GaussianKernelShape) {
        LIBSAKURA_SYMBOL(Convolve1DContext) *context;
	    float out_left[NUM_CHANNEL]={};
        float out_right[NUM_CHANNEL]={};
        float center[NUM_CENTER]={};
	    size_t num_channel=128;
        size_t kernel_width=3;
        bool fftuse = true;
        size_t kernel_center;
        kernel_center = num_channel*0.5;

	if (verbose) PrintInputs();

	LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_channel,LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian),kernel_width,fftuse, &context);
	std::cout << "center = " << kernel_center << std::endl;
	center[0] = context->input_real_array[kernel_center];

	// Verification
	//EXPECT_EQ(in1_[5],center[0]) << "center verification" ;
	for ( uint i = 0 ; i < 5 ; ++i){
	       out_left[i]  = context->input_real_array[kernel_center -1 -i];
           out_right[i] = context->input_real_array[kernel_center +1 +i];
           //ASSERT_EQ(out_left[i],out_right[i]);
           EXPECT_EQ(out_left[i],out_right[i]);
   		   std::cout << "left[" << i << "]=" << out_left[i] << "   right[" << i << "]=" << out_right[i] << std::endl;
	}

	if(context != NULL){
		free(context);
		context = NULL;
	}
}
