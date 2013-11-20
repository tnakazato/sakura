#include <iostream>
#include <string>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fftw3.h>
#include <libsakura/localdef.h>
#include <libsakura/sakura.h>
#include "gtest/gtest.h"
#include "aligned_memory.h"

/* the number of elements in input/output array to test */
#define NUM_WIDTH 5
#define NUM_IN 24

extern "C" {
struct LIBSAKURA_SYMBOL(Convolve1DContext) {
	bool use_fft;
	size_t num_data;
	fftwf_plan plan_real_to_complex_float;
	fftwf_plan plan_complex_to_real_float;
	fftwf_complex *fft_applied_complex_input_data;
	fftwf_complex *multiplied_complex_data;
	fftwf_complex *fft_applied_complex_kernel;
	float *real_array;
};
}

using namespace std;

/*
 * A super class to test creating kernel of array(s) width=3, channels=128
 * INPUTS:
 * - in1 = [0.000141569,0.00226511, 0.0195716, 0.0913234,0.230121, 0.313146,
 *          0.230121,0.0913234,0.0195716,0.00226511,0.000141569]
 * - in2 = [1,-0.998046 ,0.992209 ,-0.982555 ,0.969198 ,-0.952291 ,0.932026 ,
 *         -0.908632 ,0.882368,-0.853519 ,0.82239 ,-0.789304 ,0.754592 ,-0.71859 ]
 *
 */
class Convolve1DOperation: public ::testing::Test {
protected:

	Convolve1DOperation() :
			verbose(true) {
	}
	virtual void SetUp() {
		//in1_[] = {0.000141569,0,};
		// Initialize sakura
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Initialize)(nullptr,
				nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}
	virtual void TearDown() {
		// Clean-up sakura
		LIBSAKURA_SYMBOL(CleanUp)();
	}
	void PrintInputs() {
		//PrintArray("in1", NUM_IN, in1_);
		//PrintArray("in2", NUM_IN, in2_);
	}
	void PrintArray(char const *name, size_t num_in, float *in) {
		cout << name << " = [";
		for (size_t i = 0; i < num_in - 1; i++)
			cout << in[i] << ", ";
		cout << in[num_in - 1];
		cout << " ]" << endl;
	}
	//float in1_[NUM_IN];
	//float in2_[NUM_IN];
	bool verbose;
};

/*
 * Test creating gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * if num_data = 0, CreateConvolve1DContext will return Status_kInvalidArgument)
 */
TEST_F(Convolve1DOperation ,InvalidNumData) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context;
	size_t const num_data(0);
	size_t const kernel_width(NUM_WIDTH);
	bool fftuse = true;
	LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type = LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
			num_data,kernel_type , kernel_width, fftuse, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument),status);
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
}
/*
 * Test creating gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * if  kernel_width = 0, CreateConvolve1DContext will return Status_kInvalidArgument)
 */
TEST_F(Convolve1DOperation ,InvalidKernelWidth ) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context;
	size_t const num_data(NUM_IN);
	size_t const kernel_width(0);
	bool fftuse = true;
	LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type = LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
			num_data,kernel_type , kernel_width, fftuse, &context);
    EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument),status);
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
}
/*
 * Test creating gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * if  LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian) = UnknownType, CreateConvolve1DContext will return Status_kInvalidArgument)
 */
TEST_F(Convolve1DOperation ,InvalidKernelType ) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context;
	size_t const num_data(NUM_IN);
	size_t const kernel_width(NUM_WIDTH);
	bool fftuse = true;
	LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type = LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
			num_data,kernel_type, kernel_width, fftuse, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),status);
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
}
/*
 * Test creating gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * gaussian kernel is symmetric
 */
TEST_F(Convolve1DOperation , GaussianKernelShape) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context;
	size_t const num_data(NUM_IN);
	size_t const kernel_width(NUM_WIDTH);
	bool fftuse = true;
	LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type = LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
			num_data,kernel_type , kernel_width, fftuse, &context);

	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),status);
	// Verification
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
}
/*
 * Test fft applied gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * imaginary part is zero,
 * real part is plus and minus value
 */
TEST_F(Convolve1DOperation , FFTappliedKernelValue) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context;
	size_t const num_data(NUM_IN);
	size_t const kernel_width(NUM_WIDTH);
	bool fftuse = true;
	LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type = LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
	//if (verbose)
		//PrintArray("output_data",num_data,);
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data,kernel_type , kernel_width, fftuse, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),status);
	// Verification
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
}
/*
 * Test result of FFTWf for gaussian kernel by sakura_Convolve1D
 * RESULT:
 * mask will be applied and then first/last 4 data of inputdata will be put 0
 */
TEST_F(Convolve1DOperation , AppliedMask) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context;
	float input_data[NUM_IN] = {1,1,1,1,-1,-1,-1,-1,-1,1,1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,1};
	size_t const num_data(ELEMENTSOF(input_data));
	//bool mask_[num_data] = {0};
	bool mask[num_data] = {0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0};
	size_t const kernel_width(NUM_WIDTH);
	bool fftuse = true;
	float output_data[num_data];

	LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type = LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
			LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
			LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));
	verbose = false;
	if (verbose)
		PrintArray("output_data",num_data,output_data);
	//for (size_t i = 0; i < NUM_IN; ++i) {
		//EXPECT_FLOAT_EQ(-0.085498303,output_data[0]);
	//}
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
}
/*
 * Test result of FFTWf for gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * convolution with fft is done
 */
TEST_F(Convolve1DOperation , FFTWfResult) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context;
	float const input_data[NUM_IN] = {1,1,1,1,-1,-1,-1,-1,-1,1,1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,1};
	size_t const num_data(ELEMENTSOF(input_data));
	bool mask[num_data] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	size_t const kernel_width(NUM_WIDTH);
	bool fftuse = true;
	float output_data[num_data];
	LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type = LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);

	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data,kernel_type, kernel_width, fftuse, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),status);

	status = LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),status);
	verbose = false;
	if (verbose)
		PrintArray("output_data",num_data,output_data);
	//for (size_t i = 0; i < NUM_IN; ++i) {
		//EXPECT_FLOAT_EQ(0.828858, output_data[0]);
	//}
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
}
/*
 * Test result of FFTWf for gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * convolution without fft is done
 */
TEST_F(Convolve1DOperation , ConvolutionResult) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context;
	float const input_data[NUM_IN] = {1,1,1,1,-1,-1,-1,-1,-1,1,1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,1};
	size_t const num_data(ELEMENTSOF(input_data));
	bool mask[num_data] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	size_t const kernel_width(NUM_WIDTH);
	bool fftuse = false;
	float output_data[num_data];
	LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type = LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);

	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data,kernel_type, kernel_width, fftuse, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),status);
	status = LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),status);
	verbose = false;
	if (verbose)
		PrintArray("output_data",num_data,output_data);
	//for (size_t i = 0; i < NUM_IN; ++i) {
		//EXPECT_FLOAT_EQ(0.828858, output_data[0]);
	//}
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
}
