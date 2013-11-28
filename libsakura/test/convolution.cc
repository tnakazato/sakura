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
#define NUM_WIDTH_THIN 2
#define NUM_IN 24
#define NUM_IN_ODD 23
#define NUM_IN_LARGE 128

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
			verbose(false) {
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
		for (size_t i = 0; i < num_in - 1; ++i)
			cout << in[i] << ", ";
		cout << in[num_in - 1];
		cout << " ]" << endl;
	}
	void PrintArray2(char const *name, size_t num_in, float *in) {
		cout << name << "";
		for (size_t i = 0; i < num_in; ++i) {
			cout << in[i] << "\n";
		}
		//cout << in[num_in - 1] << endl;
	}
	void PrintArrayBool(char const *name, size_t num_in, bool *in) {
		cout << name << "";
		for (size_t i = 0; i < num_in - 1; ++i)
			cout << in[i] << "\n";
		cout << in[num_in - 1] << endl;
	}
	//float in1_[NUM_IN];
	//float in2_[NUM_IN];
	bool verbose;
};

/*
 * Test creating gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * if num_data = 0,kernel_width = 0, LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian) = UnknownType,
 * CreateConvolve1DContext will return Status_kInvalidArgument)
 */
TEST_F(Convolve1DOperation ,InvalidArguments) {
	{ // num_data = 0
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		size_t const num_data(0);
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);

		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);

		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // kernel_width = 0
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		size_t const num_data(NUM_IN);
		size_t const kernel_width(0);
		bool fftuse = true;
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);

		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // KernelType = undefined
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		size_t const num_data(NUM_IN);
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		size_t invalid_kernel = 4; // undefined kernel number of enum
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
				static_cast<LIBSAKURA_SYMBOL(Convolve1DKernelType)>(invalid_kernel);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}

/*
 * Test different number of data between context->num_data and num_data
 * with running Convolve1D
 * RESULT:
 * kInvalidArgument status will be returned and then
 * message "num_data does't equal to context->num_data" will be shown.
 */
TEST_F(Convolve1DOperation , DifferentNumdata) {
	{ // num_data != context->num_data
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		size_t const bad_num_data(NUM_IN_ODD);
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			mask[i] = 1; // no mask data
			input_data[i] = 0.0;
		}
		input_data[NUM_IN / 2] = 1.0; // one spike
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Status) status_Convolve1d =
		LIBSAKURA_SYMBOL(Convolve1D)(context, bad_num_data, input_data, mask,
				output_data);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Convolve1d); // Status_kInvalidArgument
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}

/*
 * Test Odd Number of data
 * RESULT:
 * convolved data will be output without problem by Convolve1D
 */
TEST_F(Convolve1DOperation , OddNumdata) {
	{ // num_data is odd
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_ODD];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			mask[i] = 1;
			input_data[i] = 0.0;
		}
		input_data[num_data / 2] = 1.0; // one spike
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Status) status_Convolve1d =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask,
				output_data);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve1d); // Status_kOK
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}

/*
 * Test creating gaussian kernel by sakura_CreateConvolve1DContext
 * and checking kernel shape,value,even,odd.
 * RESULT:
 * intended gaussian kernel shape,value will be obtained.
 */
TEST_F(Convolve1DOperation , ValidateGaussianKernel) {
	{ // [even],withFFT, Gaussian Kernel Shape,input only 1 spike at center
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float gaussian_kernel_even[10] = { 0.0198866, 0.04829242, 0.09394373,
				0.1463953, 0.1827497, 0.1827497, 0.1463953, 0.09394372,
				0.04829243, 0.0198866 }; // calculated data beforehand
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
			mask[i] = 1; // no mask
		}
		input_data[num_data / 2] = 1.0; // center of kernel
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask,
				output_data);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + i]);
			EXPECT_FLOAT_EQ(gaussian_kernel_even[5 + (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [odd],withFFT, Gaussian Kernel Shape,input only 1 spike at center
		LIBSAKURA_SYMBOL(Convolve1DContext)*context = nullptr;
		float gaussian_kernel_even[11] = { 0.01174296, 0.03186111, 0.06924918,
				0.1205698, 0.168164, 0.1878875, 0.168164, 0.1205698, 0.06924917,
				0.03186112, 0.01174296 }; // calculated data beforehand
		float input_data[NUM_IN_ODD];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
			mask[i] = 1; // no mask
		}
		input_data[num_data / 2] = 1.0; // center of kernel
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		float output_data[num_data];
		LIBSAKURA_SYMBOL(
				Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(
				Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(
				Convolve1D)(context, num_data, input_data, mask, output_data);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
			EXPECT_FLOAT_EQ(gaussian_kernel_even[5 + (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
		}
		EXPECT_FLOAT_EQ(
				gaussian_kernel_even[ELEMENTSOF(gaussian_kernel_even)/2],
				output_data[(num_data / 2)]);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(
				DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [even],without FFT, Gaussian Kernel Shape,input only 1 spike at center
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float gaussian_kernel_even[10] = { 0.0198866, 0.04829242, 0.09394373,
				0.1463953, 0.1827497, 0.1827497, 0.1463953, 0.09394372,
				0.04829243, 0.0198866 }; // calculated data beforehand
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
			mask[i] = 1; // no mask
		}
		input_data[num_data / 2] = 1.0; // center of kernel
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = false; // without FFT
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask,
				output_data);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + i]);
			EXPECT_FLOAT_EQ(gaussian_kernel_even[5 + (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [odd],without FFT, Gaussian Kernel Shape,input only 1 spike at center
		LIBSAKURA_SYMBOL(Convolve1DContext)*context = nullptr;
		float gaussian_kernel_even[11] = { 0.01174296, 0.03186111, 0.06924918,
				0.1205698, 0.168164, 0.1878875, 0.168164, 0.1205698, 0.06924917,
				0.03186112, 0.01174296 }; // calculated data beforehand
		float input_data[NUM_IN_ODD];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
			mask[i] = 1; // no mask
		}
		input_data[num_data / 2] = 1.0; // center of kernel
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = false; // without FFT
		float output_data[num_data];
		LIBSAKURA_SYMBOL(
				Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(
				Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(
				Convolve1D)(context, num_data, input_data, mask, output_data);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
			EXPECT_FLOAT_EQ(gaussian_kernel_even[5 + (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
		}
		EXPECT_FLOAT_EQ(
				gaussian_kernel_even[ELEMENTSOF(gaussian_kernel_even)/2],
				output_data[(num_data / 2)]);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(
				DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}
/*
 * Test result of FFTWf for gaussian kernel by sakura_Convolve1D
 * RESULT:
 * mask will be applied and then first/last 4 data of inputdata will be put 0
 */
TEST_F(Convolve1DOperation , AppliedMask) {
	{ // with mask
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN] = { 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1,
				1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1 };
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data] = { 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
				1, 1, 1, 1, 0, 0, 0, 0 };
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		float output_data[num_data];

		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));
        //verbose = true;
		if (verbose)
			PrintArray("output_data", num_data, output_data);
//for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
//EXPECT_FLOAT_EQ(-0.085498303,output_data[0]);
//}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // without mask
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float const input_data[NUM_IN] = { 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1,
				1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1 };
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
				1, 1, 1, 1, 1, 1, 1, 1 };
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		status = LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				mask, output_data);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
//verbose = true ;
		if (verbose)
			PrintArray("output_data", num_data, output_data);
//for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
//EXPECT_FLOAT_EQ(0.828858, output_data[0]);
//}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}

/*
 * Test 128ch input_data with mask by sakura_CreateConvolve1DContext
 * RESULT:
 * mask applied
 */
TEST_F(Convolve1DOperation ,ConvolutionWithMaskOnOff) {
	{ // without mask 128ch
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			input_data[i] = 0.0;
		for (size_t i = 32; i < 98; ++i)
			input_data[i] = 1.0;
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			mask[i] = 1;
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));
//verbose = true;
		if (verbose)
			PrintArray2("\n", num_data, output_data);
//PrintArrayBool("\n",num_data,mask);
//EXPECT_FLOAT_EQ(-0.085498303,output_data[0]);

		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // with mask 128ch
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			input_data[i] = 0.0;
		for (size_t i = 32; i < 98; ++i)
			input_data[i] = 1.0;
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			mask[i] = 1;
		for (size_t i = 58; i < 72; ++i)
			mask[i] = 0;
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));
//verbose = true;
		if (verbose)
			PrintArray2("\n", num_data, output_data);
//PrintArrayBool("\n",num_data,mask);
//EXPECT_FLOAT_EQ(-0.085498303,output_data[0]);

		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // input_data all 1, mask all 0
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			input_data[i] = 1.0;
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			mask[i] = 0;
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));

		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			EXPECT_FLOAT_EQ(0, output_data[i]);

		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // input_data all 1, mask all 1 (with FFT)
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			input_data[i] = 1.0;
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			mask[i] = 1;
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));

		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			EXPECT_FLOAT_EQ(1, output_data[i]);

		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // input_data 1 spike, mask all 1 then output_data should be equal to kernel (without FFT)
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			input_data[i] = 0.0;
		input_data[63] = 1.0; // center of kernel

		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			mask[i] = 1;
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = false;
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));
        //verbose = true;
		if (verbose)
			PrintArray2("\n", num_data, output_data);
          //PrintArray2("\n",num_data,context->real_array);
		verbose = false;
		for (size_t i = 0; i < ELEMENTSOF(input_data) / 2 - 1; ++i) {
			EXPECT_FLOAT_EQ(
					context->real_array[ELEMENTSOF(input_data)/2 + 1 + i],
					output_data[i]);
			EXPECT_FLOAT_EQ(context->real_array[i],
					output_data[ELEMENTSOF(input_data)/2 -1 + i]);
		}

		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // input_data 32-92  delta, mask 44-56, 68-80
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			input_data[i] = 0.0;
		for (size_t i = 33; i < 92; ++i)
			input_data[i] = 1.0;
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			mask[i] = 1;
		for (size_t i = 44; i < 56; ++i)
			mask[i] = 0;
		for (size_t i = 68; i < 80; ++i)
			mask[i] = 0;
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			input_data[i] = 0.0;
		for (size_t i = 32; i < 44; ++i)
			input_data[i] = 1.0;
		for (size_t i = 56; i < 68; ++i)
			input_data[i] = 1.0;
		for (size_t i = 80; i < 92; ++i)
			input_data[i] = 1.0;

		size_t const kernel_width(20); //kernel_width(NUM_WIDTH);
		bool fftuse = true;
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));
        //verbose = true;
		if (verbose)
			PrintArray2("\n", num_data, output_data);
		verbose = false;
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // with mask 128ch all zero against spike --> all 0
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			input_data[i] = 0.0;
		for (size_t i = 58; i < 72; ++i)
			input_data[i] = 1.0;
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			mask[i] = 0;
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));
        //verbose = true;
		if (verbose)
			PrintArray2("\n", num_data, output_data);
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			EXPECT_FLOAT_EQ(0, output_data[i]);

		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // with mask 128ch , spike( 58 < ch < 72)
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			input_data[i] = 0.0;
		for (size_t i = 58; i < 72; ++i)
			input_data[i] = 1.0;
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			mask[i] = 1;
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));
//verbose = true;
		if (verbose)
			PrintArray2("\n", num_data, output_data);
//verbose = false;
		EXPECT_GT(output_data[64], 0.6);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // without mask 128ch , and then input this output again as input_data
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			input_data[i] = 0.0;
		for (size_t i = 32; i < 98; ++i)
			input_data[i] = 1.0;
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i)
			mask[i] = 1;
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));
//verbose = true;
		if (verbose)
			PrintArray2("\n", num_data, output_data);
//PrintArrayBool("\n",num_data,mask);
//EXPECT_FLOAT_EQ(-0.085498303,output_data[0]);
		float output_data_reuse[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			if ((i > 44 && i < 54))
				mask[i] = 0;
			else
				mask[i] = 1;
			if ((i > 44 && i < 54))
				output_data_reuse[i] = 0.0;
			else
				output_data_reuse[i] = output_data[i];
		}
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, output_data_reuse, mask, output_data));
//verbose = true;
		if (verbose)
			PrintArray2("\n", num_data, output_data);
//PrintArrayBool("\n",num_data,mask);
//EXPECT_FLOAT_EQ(-0.085498303,output_data[0]);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}

}
/*
 * Test result of FFTWf for gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * convolution without fft is done
 */
TEST_F(Convolve1DOperation , ConvolutionResultWithoutFFT) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
	float const input_data[NUM_IN] = { 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1,
			1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1 };
	size_t const num_data(ELEMENTSOF(input_data));
	bool mask[num_data] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1 };
	size_t const kernel_width(NUM_WIDTH);
	bool fftuse = false;
	float output_data[num_data];
	LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
			num_data, kernel_type, kernel_width, fftuse, &context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	status = LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask,
			output_data);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	//verbose = true ;
	if (verbose)
		PrintArray("output_data", num_data, output_data);
	//for (size_t i = 0; i < ELEMENTSOF(input_data) ; ++i) {
	//EXPECT_FLOAT_EQ(0.828858, output_data[0]);
	//}
	LIBSAKURA_SYMBOL(Status) status_Destroy =
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
}

/*
 * Test Destroy context by sakura_DestroyConvolve1DContext
 * RESULT:
 * pointer of context and its member should be nil.
 */
TEST_F(Convolve1DOperation , DestroyConvolve1DContext) {
	{
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float const input_data[NUM_IN] = { 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1,
				1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1 };
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
				1, 1, 1, 1, 1, 1, 1, 1 };
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		status = LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				mask, output_data);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		if (verbose)
			PrintArray("output_data", num_data, output_data);

		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}
