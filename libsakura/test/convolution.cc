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
#define NUM_CHANNEL 24
#define NUM_WIDTH 5
#define NUM_IN 11
#define NUM_CENTER 1

extern "C" {
struct LIBSAKURA_SYMBOL(Convolve1DContext) {
	size_t num_data;
	fftwf_plan *plan_real_to_complex_float;
	fftwf_plan *plan_complex_to_real_float;
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
class CreateConvolve1DContext: public ::testing::Test {
protected:

	CreateConvolve1DContext() :
			verbose(true) {
	}

	virtual void SetUp() {
		in1_[0] = 0.000141569;
		in1_[1] = 0.00226511;
		in1_[2] = 0.0195716;
		in1_[3] = 0.0913234;
		in1_[4] = 0.230121;
		in1_[5] = 0.313146; // center
		in1_[6] = 0.230121;
		in1_[7] = 0.0913234;
		in1_[8] = 0.0195716;
		in1_[9] = 0.00226511;
		in1_[10] = 0.000141569;

		in2_[0] = 1;
		in2_[1] = -0.998046;
		in2_[2] = 0.992209;
		in2_[3] = -0.982555;
		in2_[4] = 0.969198;
		in2_[5] = -0.952291;
		in2_[6] = 0.932026;
		in2_[7] = -0.908632;
		in2_[8] = 0.882368;
		in2_[9] = -0.853519;
		in2_[10] = 0.82239;

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
		PrintArray("in1", NUM_IN, in1_);
		PrintArray("in2", NUM_IN, in2_);
	}

	void PrintArray(char const *name, size_t num_in, float *in) {
		cout << name << " = [";
		for (size_t i = 0; i < num_in - 1; i++)
			cout << in[i] << ", ";
		cout << in[num_in - 1];
		cout << " ]" << endl;
	}

	float in1_[NUM_IN];
	float in2_[NUM_IN];bool verbose;

};

/*
 * Test creating gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * gaussian kernel is symmetric
 */
TEST_F(CreateConvolve1DContext , GaussianKernelShape) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context;
	size_t num_data = 24;
	size_t kernel_width = 5;
	bool fftuse = true;

	//if (verbose)
	//	PrintInputs();

	LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data,
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian), kernel_width, fftuse,
			&context);
	// Verification
	//EXPECT_EQ(in1_[5],center_[0]) << "center verification" ;
	if (context != nullptr) {
	//	EXPECT_FLOAT_EQ(0.8495121, context->fft_applied_complex_kernel[1][0]);
	//	EXPECT_FLOAT_EQ(0.11184037, context->fft_applied_complex_kernel[1][1]);
	}
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
}

/*
 * Test fft applied gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * imaginary part is zero,
 * real part is plus and minus value
 */
TEST_F(CreateConvolve1DContext , FFTappliedKernelValue) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context;

	size_t num_data = 24;
	size_t kernel_width = 5;
	bool fftuse = true;

	//if (verbose)
	//	PrintInputs();

	LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data,
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian), kernel_width, fftuse,
			&context);
	// Verification
	for (size_t i = 0; i < 11; ++i) {
		//out1_[i] = context->fft_applied_complex_kernel[i][0];
		//EXPECT_EQ(in2_[i],out1_[i]);
		// std::cout << "in2_[" << i << "]=" << in2_[i] << "   fft_applied[" << i << "]=" << out1_[i] << std::endl;
	}
	//out2_[0] = context->fft_applied_complex_kernel[0][1];
	//EXPECT_EQ(0, out2_[0]);

	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);

}

/*
 * Test result of FFTWf for gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * smoothing is done
 */
TEST_F(CreateConvolve1DContext , FFTWfResult) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context;

	float outspec_[NUM_CHANNEL] = { 0, };
	float inspec_[NUM_CHANNEL] = { -0.0834709, -0.11101, -0.0766486, 0.0074527,
			0.10254, 0.137752, 0.100273, 0.0281713, -0.101949, -0.2157,
			-0.106586, 0.020014, -0.00773457, 0.014356, 0.0639472, 0.0489224,
			0.0463871, 0.0547668, -0.0214724, -0.0627455, -0.0823716, -0.149673,
			-0.138719, -0.100872 };
	bool mask_[NUM_CHANNEL] = { 0, };

	size_t num_data = 24;
	size_t kernel_width = 5;
	bool fftuse = true;

	//if (verbose)
	//	PrintInputs();

	LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data,
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian), kernel_width, fftuse,
			&context);

	LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, inspec_, mask_, outspec_);
	for (size_t i = 0; i < NUM_CHANNEL; ++i)
		std::cout << "outspec_[i] = " << outspec_[i] << endl;

	EXPECT_FLOAT_EQ(-0.074806437, outspec_[0]);
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
}
