#include <iostream>
#include <string>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fftw3.h>
#include <climits>
#include <memory>

#include <libsakura/localdef.h>
#include <libsakura/sakura.h>
#include "loginit.h"
#include "gtest/gtest.h"
#include "aligned_memory.h"

/* the number of elements in input/output array to test */
#define NUM_WIDTH 5
#define NUM_WIDTH_THIN 2
#define NUM_IN 24
#define NUM_IN_ODD 25
#define NUM_IN_LARGE 64
#define NUM_IN_MAX 8192
#define LOOP_MAX 1000

extern "C" {
struct LIBSAKURA_SYMBOL(Convolve1DContext) {
	bool use_fft;
	size_t num_data;
	size_t kernel_width;
	fftwf_plan plan_real_to_complex_float;
	fftwf_plan plan_complex_to_real_float;
	fftwf_complex *fft_applied_complex_kernel;
	float *real_array;
	void *real_array_work;
	float *real_array_kernel;
	void *real_kernel_array_work;
};
}
namespace {
void *DummyAlloc(size_t size) {
	return nullptr;
}
} /* namespace */

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
		//PrintArray("in2", NUM_IN, in2_);
	}
	void PrintArray(char const *name, size_t num_in, float *in) {
		cout << name << "";
		for (size_t i = num_in / 2 - 10; i < num_in / 2 + 10; ++i) {
			cout << setprecision(10) << in[i] << "\n";
		}
		//cout << in[num_in - 1] << endl;
	}

	bool verbose;
};
/*
 * Test Alignment Check
 * RESULT:
 * input/output_data and real_array/real_array_kernel should be aligned
 */
TEST_F(Convolve1DOperation , AlignmentCheck) {
	{
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		SIMD_ALIGN
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
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

		ASSERT_TRUE(sakura_IsAligned(context->real_array))<< "real_array is not aligned";
		ASSERT_TRUE(sakura_IsAligned(context->real_array_kernel))<< "real_array_kernel is not aligned";
		SIMD_ALIGN
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Status) status_Convolve1d =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve1d);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}

/*
 * Test creating gaussian kernel by sakura_CreateConvolve1DContext
 * RESULT:
 * if num_data = 0,kernel_width = 0, LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian) = UnknownType,
 * CreateConvolve1DContext will return Status_kInvalidArgument)
 */
TEST_F(Convolve1DOperation ,InvalidArguments) {
	{ // num_data > INT_MAX
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		size_t const num_data(INT_MAX + 1);
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
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Destroy);
	}
	{ // num_data == 0
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
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Destroy);
	}
	{ // kernel_width == 0
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
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Destroy);
	}
	{ // (FFT) kernel_width > num_data
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		SIMD_ALIGN
		float input_data[NUM_IN_LARGE];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		input_data[num_data / 2] = 1.0; // center ch
		size_t const kernel_width( NUM_IN_LARGE + 1);
		bool fftuse = true; // FFT
		SIMD_ALIGN
		float output_data[ELEMENTSOF(input_data)];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create;
		status_Create = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data,
				kernel_type, kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve;
		status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("\n", num_data, output_data);
		}
		//verbose = false;
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // (Without FFT) kernel_width > num_data
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		SIMD_ALIGN
		float input_data[NUM_IN_LARGE];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		input_data[num_data / 2] = 1.0; // center ch
		size_t const kernel_width( NUM_IN_LARGE + 1);
		bool fftuse = false; // without FFT
		SIMD_ALIGN
		float output_data[ELEMENTSOF(input_data)];
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			output_data[i] = 0.0;
		}
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create;
		status_Create = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data,
				kernel_type, kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve;
		status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("without FFT\n", num_data, output_data);
		}
		verbose = false;
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // KernelType == undefined
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		size_t const num_data(NUM_IN);
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		auto invalid_kernel_type =
				static_cast<LIBSAKURA_SYMBOL(Convolve1DKernelType)>(-1);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, invalid_kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Destroy);
	}
	{ // context == nullptr
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		size_t const num_data(NUM_IN);
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true;
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
				LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, nullptr);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Destroy);
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
		SIMD_ALIGN
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		size_t const bad_num_data(NUM_IN_ODD);
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
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
		SIMD_ALIGN
		float output_data[num_data];
		EXPECT_TRUE(num_data != bad_num_data)
				<< "In this test, num_data != bad_num_data is expected"; // assert
		LIBSAKURA_SYMBOL(Status) status_Convolve1d =
		LIBSAKURA_SYMBOL(Convolve1D)(context, bad_num_data, input_data,
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
		SIMD_ALIGN
		float input_data[NUM_IN_ODD];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
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
		SIMD_ALIGN
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Status) status_Convolve1d =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve1d); // Status_kOK
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}

/*
 * Test Failed malloc of context
 * RESULT:
 * Convolve1D will return Status_kInvalidArgument
 */
TEST(Convolve1DOperationFailed , FailedMallocContext) {
	{
		LIBSAKURA_SYMBOL(Status) status_Init = LIBSAKURA_SYMBOL(Initialize)(
				DummyAlloc, nullptr);
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		SIMD_ALIGN
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Init);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kNoMemory), status_Create);
		SIMD_ALIGN
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Status) status_Convolve1d =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				output_data);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Convolve1d);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Destroy);
		LIBSAKURA_SYMBOL(CleanUp)();
	}

}

/*
 * Test creating gaussian kernel by sakura_CreateConvolve1DContext
 * and checking kernel shape,value,even,odd.
 * RESULT:
 * intended gaussian kernel shape,value will be obtained.
 */
TEST_F(Convolve1DOperation , ValidateGaussianKernel) {
	{ // [even],FFT, Gaussian Kernel Shape,input only 1 spike at center
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float gaussian_kernel[11] = { 0.011742971, 0.031861119, 0.069249183,
				0.12056981, 0.16816399, 0.187887, 0.16816399, 0.12056981,
				0.069249183, 0.031861119, 0.011742971 }; // calculated data beforehand
		SIMD_ALIGN
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		input_data[num_data / 2] = 1.0; // center of kernel
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		SIMD_ALIGN
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			std::cout << " even with FFT" << std::endl;
			PrintArray("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
			//EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
			//		output_data[(num_data / 2) + (i + 1)]);
			//EXPECT_FLOAT_EQ(gaussian_kernel[kernel_width + (i + 1)],
			//		output_data[(num_data / 2) + (i + 1)]);
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [odd],FFT, Gaussian Kernel Shape,input only 1 spike at center
		LIBSAKURA_SYMBOL(Convolve1DContext)*context = nullptr;
		float gaussian_kernel[11] = { 0.011742957, 0.031861119, 0.069249183,
				0.12056981, 0.16816399, 0.18788746, 0.16816399, 0.12056981,
				0.069249183, 0.031861119, 0.011742957 }; // calculated data beforehand
		SIMD_ALIGN
		float input_data[NUM_IN_ODD];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		input_data[num_data / 2] = 1.0; // center of kernel
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		SIMD_ALIGN
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
				Convolve1D)(context, num_data, input_data, output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
			//EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
			//		output_data[(num_data / 2) + (i + 1)]);
			//EXPECT_FLOAT_EQ(gaussian_kernel[kernel_width + (i + 1)],
			//		output_data[(num_data / 2) + (i + 1)]);
		}
		//EXPECT_FLOAT_EQ(gaussian_kernel[ELEMENTSOF(gaussian_kernel)/2],
		//		output_data[(num_data / 2)]);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(
				DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [even],without FFT, Gaussian Kernel Shape,input only 1 spike at center
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		float gaussian_kernel[11] = { 0.011742963, 0.0318611, 0.0692492,
				0.12056981, 0.168164, 0.187887, 0.168164, 0.12056981, 0.0692492,
				0.0318611, 0.011742963 }; // calculated data beforehand
		SIMD_ALIGN
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		input_data[num_data / 2] = 1.0; // center of kernel
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = false; // without FFT
		SIMD_ALIGN
		float output_data[ELEMENTSOF(input_data)];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			std::cout << " even without FFT" << std::endl;
			PrintArray("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width; ++i) {
			//EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
			//		output_data[(num_data / 2) + (i + 1)]);
			//EXPECT_FLOAT_EQ(gaussian_kernel[kernel_width + (i + 1)],
			//		output_data[(num_data / 2) + (i + 1)]);
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [odd],without FFT, Gaussian Kernel Shape,input only 1 spike at center
		LIBSAKURA_SYMBOL(Convolve1DContext)*context = nullptr;
		float gaussian_kernel[11] = { 0.011742963, 0.0318611, 0.0692492,
				0.12056981, 0.168164, 0.18788746, 0.168164, 0.12056981,
				0.0692492, 0.0318611, 0.011742963 }; // calculated data beforehand
		SIMD_ALIGN
		float input_data[NUM_IN_ODD];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		input_data[num_data / 2] = 1.0; // center of kernel
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = false; // without FFT
		SIMD_ALIGN
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
				Convolve1D)(context, num_data, input_data, output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width; ++i) {
			//EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
			//		output_data[(num_data / 2) + (i + 1)]);
			//EXPECT_FLOAT_EQ(gaussian_kernel[kernel_width + (i + 1)],
			//		output_data[(num_data / 2) + (i + 1)]);
		}
		//EXPECT_FLOAT_EQ(gaussian_kernel[ELEMENTSOF(gaussian_kernel)/2],
		//		output_data[(num_data / 2)]);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(
				DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}

/*
 * Test convolution against input spike
 * RESULT: same result will be yielded by with/without FFT
 */
TEST_F(Convolve1DOperation , OtherInputDataFFTonoff) {
	{ // [even],FFT, Gaussian Kernel Shape,input delta
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		SIMD_ALIGN
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 1.0;
			if (i > 7 && i < 15) {
				input_data[i] = -1.0;
			}
		}
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		SIMD_ALIGN
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("\n", num_data, output_data);
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
		SIMD_ALIGN
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 1.0;
			if (i > 7 && i < 15) {
				input_data[i] = -1.0;
			}
		}
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = false; // without FFT
		SIMD_ALIGN
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
	{ // [even],FFT, Gaussian Kernel Shape,input only 2 spike at 0, num_data - 1 , multi
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		SIMD_ALIGN
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		input_data[0] = 1.0; // first ch
		input_data[num_data - 1] = 1.0; // final ch
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		SIMD_ALIGN
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("\n", num_data, output_data);
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
		SIMD_ALIGN
		float input_data[NUM_IN];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		input_data[0] = 1.0; // first ch
		input_data[num_data - 1] = 1.0; // final ch

		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = false; // without FFT
		SIMD_ALIGN
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		LIBSAKURA_SYMBOL(Status) status_Convolve =
		LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("\n", num_data, output_data);
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
 * Test Compare result of with/without FFTW for gaussian kernel
 * by sakura_CreateConvolve1DContext
 * RESULT:
 * each result will be equal between convolution with/without fft
 */
TEST_F(Convolve1DOperation , CompareResultWithFFTWithoutFFT) {
	{
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		SIMD_ALIGN
		float input_data[NUM_IN_LARGE];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		input_data[num_data / 2] = 1.0; // center
		size_t const kernel_width(NUM_WIDTH);
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		bool fftuse = false; // without fft
		SIMD_ALIGN
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		LIBSAKURA_SYMBOL(Status) status_Convolve = LIBSAKURA_SYMBOL(Convolve1D)(
				context, num_data, input_data, output_data);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("output_data\n", num_data, output_data);
		}
		fftuse = true; // with fft
		float output_data_fft[num_data];
		status = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data,
				kernel_type, kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		status_Convolve = LIBSAKURA_SYMBOL(Convolve1D)(context, num_data,
				input_data, output_data_fft);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("output_data_fft\n", num_data, output_data_fft);
		}
		verbose = false;
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			if ((i > num_data / 2 - kernel_width / 2)
					&& (i < num_data / 2 + kernel_width / 2)) {
				//EXPECT_FLOAT_EQ(output_data_fft[i], output_data[i]); // compare
			}
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);
	}
}

/*
 * Test Destroy context by sakura_DestroyConvolve1DContext
 * RESULT:
 * pointer of context and its member should be nil.
 */
TEST_F(Convolve1DOperation , DestroyConvolve1DContext) {
	{
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		SIMD_ALIGN
		float const input_data[NUM_IN] = { 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1,
				1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1 };
		size_t const num_data(ELEMENTSOF(input_data));
		SIMD_ALIGN
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with fft
		SIMD_ALIGN
		float output_data[num_data];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data, kernel_type,
				kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		status = LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
				output_data);
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

/*
 * Test [Without FFT] Performance of Convolve1D
 * RESULT: elapsed time will be output
 */
TEST_F(Convolve1DOperation , PerformanceTestWithoutFFT) {
	{ // [even],without FFT, Gaussian Kernel Shape,input delta
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		SIMD_ALIGN
		float input_data[NUM_IN_MAX];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		input_data[num_data / 2] = 1.0; // center ch
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = false; // without FFT
		SIMD_ALIGN
		float output_data[ELEMENTSOF(input_data)];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create;
		double start = sakura_GetCurrentTime();
		status_Create = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data,
				kernel_type, kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		double end = sakura_GetCurrentTime();
		LIBSAKURA_SYMBOL(Status) status_Convolve;
		double start_time = sakura_GetCurrentTime();
		for (size_t i = 0; i < LOOP_MAX; ++i) {
			status_Convolve =
			LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
					output_data);
		}
		double end_time = sakura_GetCurrentTime();
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("without FFT\n", num_data, output_data);
		}
		verbose = false;
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);

		std::cout << "without FFT create = " << num_data << ", elapsed time = "
				<< end - start << "secs\n";
		std::cout << "without FFT convolve= " << num_data << ", elapsed time = "
				<< end_time - start_time << "secs\n";
	}
}

/*
 * Test [With FFT] Performance of Convolve1D
 * RESULT: elapsed time will be output
 */
TEST_F(Convolve1DOperation , PerformanceTestWithFFT) {
	{ // [even], FFT, Gaussian Kernel Shape,input delta
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		SIMD_ALIGN
		float input_data[NUM_IN_MAX];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		input_data[num_data / 2] = 1.0; // center ch
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		SIMD_ALIGN
		float output_data[ELEMENTSOF(input_data)];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create;
		double start = sakura_GetCurrentTime();
		status_Create = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data,
				kernel_type, kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		double end = sakura_GetCurrentTime();
		LIBSAKURA_SYMBOL(Status) status_Convolve;
		double start_time = sakura_GetCurrentTime();
		for (size_t i = 0; i < LOOP_MAX; ++i) {
			status_Convolve =
			LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
					output_data);
		}
		double end_time = sakura_GetCurrentTime();
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);

		std::cout << "with FFT create = " << num_data << ", elapsed time = "
				<< end - start << "secs\n";
		std::cout << "with FFT convolve= " << num_data << ", elapsed time = "
				<< end_time - start_time << "secs\n";
	}
	{ // [odd], FFT, Gaussian Kernel Shape,input delta
		LIBSAKURA_SYMBOL(Convolve1DContext) *context = nullptr;
		SIMD_ALIGN
		float input_data[NUM_IN_MAX - 1];
		size_t const num_data(ELEMENTSOF(input_data));
		SIMD_ALIGN
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
		}
		input_data[num_data / 2] = 1.0; // center ch
		size_t const kernel_width(NUM_WIDTH);
		bool fftuse = true; // with FFT
		SIMD_ALIGN
		float output_data[ELEMENTSOF(input_data)];
		LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type =
		LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian);
		LIBSAKURA_SYMBOL(Status) status_Create;
		double start = sakura_GetCurrentTime();
		status_Create = LIBSAKURA_SYMBOL(CreateConvolve1DContext)(num_data,
				kernel_type, kernel_width, fftuse, &context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Create);
		double end = sakura_GetCurrentTime();
		LIBSAKURA_SYMBOL(Status) status_Convolve;
		double start_time = sakura_GetCurrentTime();
		for (size_t i = 0; i < LOOP_MAX; ++i) {
			status_Convolve =
			LIBSAKURA_SYMBOL(Convolve1D)(context, num_data, input_data,
					output_data);
		}
		double end_time = sakura_GetCurrentTime();
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Convolve);
		//verbose = true;
		if (verbose) {
			PrintArray("\n", num_data, output_data);
		}
		verbose = false;
		for (size_t i = 0; i < kernel_width - 1; ++i) {
		}
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Destroy);

		std::cout << "with FFT create = " << num_data << ", elapsed time = "
				<< end - start << "secs\n";
		std::cout << "with FFT convolve= " << num_data << ", elapsed time = "
				<< end_time - start_time << "secs\n";
	}
}

/*
 * Test
 * RESULT:
 */
//TEST_F(Convolve1DOperation , dummy) {
//}
