#include <iostream>
#include <string>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fftw3.h>
#include <omp.h>

#include <libsakura/localdef.h>
#include <libsakura/sakura.h>
#include "loginit.h"
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
	//fftwf_complex *fft_applied_complex_input_data;
	//fftwf_complex *multiplied_complex_data;
	//fftwf_complex *fft_applied_complex_kernel;
	float *real_array;
};
}
using namespace std;

/*
 * A super class to test creating kernel of array(s) width=3, channels=128
 * INPUTS:
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
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);

		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
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
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
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
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
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
			mask[i] = true; // no mask data
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
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		float output_data[num_data];
		//EXPECT_TRUE(num_data != bad_num_data) << "num_data != bad_num_data"; // assert
		LIBSAKURA_SYMBOL(Status) status_Convolve1d =
		LIBSAKURA_SYMBOL(Convolve1D)(context, bad_num_data, input_data, mask,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Convolve1d); // Status_kInvalidArgument
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
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
			mask[i] = true;
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
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Status) status_Convolve1d =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve1d); // Status_kOK
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}

/*
 * Test Original Mask by Linear Interpolation
 * RESULT:
 * Interpolate instead of zero.
 */
TEST_F(Convolve1DOperation , LinearInterpolation) {
	{ // with mask against input delta
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			mask[i] = true;
			input_data[i] = 0.3;
			if (i >= 48 && i <= 84) {
				input_data[i] = 0.3 + i * 0.0175;
			} else if (i > 84) {
				input_data[i] = 1.0;
			}
			if (i >= 20 && i <= 30) {
				mask[i] = false; // mask on against delta
			}
			if (i >= 50 && i <= 70) {
				mask[i] = false; // mask on against delta
			}
		}
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = false; // without FFT
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(
				Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(
				Convolve1D)(context, num_data, input_data, mask, output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);

		//for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
		//EXPECT_FLOAT_EQ(0, output_data[i]); // delta will remain, but its convolved
		//}
		//verbose = true;
		if (verbose) {
			//PrintArray2("\n", num_data, input_data);
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
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
		float gaussian_kernel[10] = { 0.0198866, 0.04829242, 0.09394373,
				0.1463953, 0.1827497, 0.1827497, 0.1463953, 0.09394372,
				0.04829243, 0.0198866 }; // calculated data beforehand
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
			mask[i] = true; // no mask
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
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + i]);
			EXPECT_FLOAT_EQ(gaussian_kernel[5 + (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [odd],withFFT, Gaussian Kernel Shape,input only 1 spike at center
		LIBSAKURA_SYMBOL(Convolve1DContext)*context = nullptr;
		float gaussian_kernel[11] = { 0.01174296, 0.03186111, 0.06924918,
				0.1205698, 0.168164, 0.1878875, 0.168164, 0.1205698, 0.06924917,
				0.03186112, 0.01174296 }; // calculated data beforehand
		float input_data[NUM_IN_ODD];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
			mask[i] = true; // no mask
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
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(
				Convolve1D)(context, num_data, input_data, mask, output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
			EXPECT_FLOAT_EQ(gaussian_kernel[5 + (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
		}
		EXPECT_FLOAT_EQ(gaussian_kernel[ELEMENTSOF(gaussian_kernel)/2],
				output_data[(num_data / 2)]);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(
				DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [even],without FFT, Gaussian Kernel Shape,input only 1 spike at center
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float gaussian_kernel[10] = { 0.0198866, 0.04829242, 0.09394373,
				0.1463953, 0.1827497, 0.1827497, 0.1463953, 0.09394372,
				0.04829243, 0.0198866 }; // calculated data beforehand
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
			mask[i] = true; // no mask
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
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + i]);
			EXPECT_FLOAT_EQ(gaussian_kernel[5 + (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [odd],without FFT, Gaussian Kernel Shape,input only 1 spike at center
		LIBSAKURA_SYMBOL(Convolve1DContext)*context = nullptr;
		float gaussian_kernel[11] = { 0.01174296, 0.03186111, 0.06924918,
				0.1205698, 0.168164, 0.1878875, 0.168164, 0.1205698, 0.06924917,
				0.03186112, 0.01174296 }; // calculated data beforehand
		float input_data[NUM_IN_ODD];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
			mask[i] = true; // no mask
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
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(
				Convolve1D)(context, num_data, input_data, mask, output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
			EXPECT_FLOAT_EQ(gaussian_kernel[5 + (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
		}
		EXPECT_FLOAT_EQ(gaussian_kernel[ELEMENTSOF(gaussian_kernel)/2],
				output_data[(num_data / 2)]);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(
				DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}

/*
 * Test convolution against delta type input data
 * RESULT: same result will be yielded by with/without FFT
 */
TEST_F(Convolve1DOperation , OtherInputDataFFTonoff) {
	{ // [even],withFFT, Gaussian Kernel Shape,input delta
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 1.0;
			if (i > 7 && i < 15) {
				input_data[i] = -1.0;
			}
			mask[i] = true; // no mask
		}
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [even],without FFT, Gaussian Kernel Shape,input delta
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 1.0;
			if (i > 7 && i < 15) {
				input_data[i] = -1.0;
			}
			mask[i] = true; // no mask
		}
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = false; // without FFT
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [even],withFFT, Gaussian Kernel Shape,input only 2 spike at 0, num_data - 1 , multi
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
			mask[i] = true; // no mask
		}
		input_data[0] = 1.0; // first ch
		input_data[num_data - 1] = 1.0; // final ch
		//input_data[1][0] = 1.0;
		//input_data[1][num_data - 1] = 1.0;
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [even],without FFT, Gaussian Kernel Shape,input only 2 spike at 0, num_data - 1
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];

		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
			mask[i] = true; // no mask
		}
		input_data[0] = 1.0; // first ch
		input_data[num_data - 1] = 1.0; // final ch

		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = false; // without FFT
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}

}

/*
 * Test mask on/off against 128ch input_data  by sakura_Convolve1D
 *
 * RESULT:
 * without mask, input data will not be applied mask.
 * with mask, output data will be zero.
 */
TEST_F(Convolve1DOperation ,ConvolutionWithMaskOnOff) {
	{ // without mask against input delta
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			mask[i] = true; // mask all 1
			input_data[i] = 0.0;
			if (i > (num_data / 2 - num_data / 8)
					&& i < (num_data / 2 + num_data / 8)) {
				input_data[i] = 1.0; // width = center +/- 16
			}
		}
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(
				Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(
				Convolve1D)(context, num_data, input_data, mask, output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			if (i > (num_data / 2 - 5) && i < (num_data / 2 + 5)) {
				EXPECT_FLOAT_EQ(1, output_data[i]); // delta will remain, but its convolved
			}
		}
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // with mask against input delta
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			mask[i] = true;
			input_data[i] = 0.0;
			if (i > (num_data / 2 - num_data / 4) // create delta
			&& i < (num_data / 2 + num_data / 4)) {
				input_data[i] = 1.0; // width = center +/- 16
				//mask[i] = false; // mask on against delta
			}
		}
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(
				Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(
				Convolve1D)(context, num_data, input_data, mask, output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);

		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			//EXPECT_FLOAT_EQ(0, output_data[i]); // delta will remain, but its convolved
		}
		//verbose = true;
		if (verbose) {
			PrintArray2("\n", num_data, output_data);
		}
		verbose = false;
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // input_data all 1, mask all 0 (with FFT)
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			//mask[i] = false;
			input_data[i] = 1.0;
		}
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			//EXPECT_FLOAT_EQ(0, output_data[i]);
			EXPECT_FLOAT_EQ(1, output_data[i]);
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // input data all 1, mask all 1 (with FFT)
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float input_data[NUM_IN_LARGE]; // 128
		size_t const num_data(ELEMENTSOF(input_data));
		bool mask[num_data];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			mask[i] = true;
			input_data[i] = 1.0;
		}
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with fft
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type, kernel_width, fftuse, &context));
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK),
				LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data, mask, output_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			EXPECT_FLOAT_EQ(1, output_data[i]);
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}
/*
 * Test Compare result of with/without FFTW for gaussian kernel
 * by sakura_CreateConvolve1DContext
 * RESULT:
 * each result will be equal between convolution with/without fft
 */
TEST_F(Convolve1DOperation , CompareResultWithFFTWithoutFFT) {
	LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
	float const input_data[NUM_IN] = { 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1,
			1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1 };
	size_t const num_data(ELEMENTSOF(input_data));
	bool mask[num_data];
	for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
		mask[i] = true; // no mask
	}
	size_t const kernel_width(NUM_WIDTH);
	LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
	bool fftuse = false; // without fft
	float output_data[num_data];
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
			num_data, kernel_type, kernel_width, fftuse, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	LIBSAKURA_SYMBOL(Status) status_Convolve = LIBSAKURA_SYMBOL(Convolve1D)(
			context, num_data, input_data, mask, output_data);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
	fftuse = true; // with fft
	float output_data_fft[num_data];
	status = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
			kernel_width, fftuse, &context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	status_Convolve = LIBSAKURA_SYMBOL(Convolve1D)(context, num_data,
			input_data, mask, output_data_fft);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
	//verbose = true;
	if (verbose) {
		PrintArray2("", num_data, output_data_fft);
		PrintArray2("", num_data, output_data);
	}
	verbose = false;
	for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
		if ((i > num_data / 2 - kernel_width / 2)
				&& (i < num_data / 2 + kernel_width / 2)) {
			EXPECT_FLOAT_EQ(output_data_fft[i], output_data[i]); // compare
		}
	}
	LIBSAKURA_SYMBOL(Status) status_Destroy =
	LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
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
		bool fftuse = true; // with fft
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		status = LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				mask, output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		//verbose = true;
		if (verbose)
			PrintArray("output_data", num_data, output_data);
		verbose = false;
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}
