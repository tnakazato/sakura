/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2016
 * National Astronomical Observatory of Japan
 * 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
 *
 * This file is part of Sakura.
 *
 * Sakura is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * Sakura is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
 * @SAKURA_LICENSE_HEADER_END@
 */
#include <iostream>
#include <string>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <initializer_list>
#include <fftw3.h>
#include <climits>
#include <cfloat>
#include <memory>
#include <stdlib.h>
#include <vector>

#include <libsakura/localdef.h>
#include <libsakura/sakura.h>
#include "loginit.h"
#include "gtest/gtest.h"
#include "aligned_memory.h"

/* the number of elements in input/output array to test */
#define NUM_WIDTH 5
#define NUM_WIDTH_THIN 2
#define NUM_IN_EVEN 24
#define NUM_IN_ODD 25
#define NUM_IN_LARGE 64
#define NUM_IN_MAX 8192
#define LOOP_MAX 1000

#define TEST_GAUSS(Name) TEST(CreateGaussianKernelTest, Name)

#define BENCH "#x# benchmark "

extern "C" {
struct LIBSAKURA_SYMBOL(Convolve1DContextFloat) {
//	bool use_fft;
//	size_t num_data;
	size_t kernel_width;
	fftwf_plan plan_real_to_complex_float;
	fftwf_plan plan_complex_to_real_float;
	fftwf_complex *fft_applied_complex_kernel;
	float *real_array;
	void *real_array_work;
	float *real_kernel_array;
	void *real_kernel_array_work;
};
}
typedef enum {
	SpikeType_kleft,
	SpikeType_kcenter,
	SpikeType_kright,
	SpikeType_kleftright,
	SpikeType_kall,
	SpikeType_knegative
} SpikeType;

typedef enum {
	MaskType_ktrue, MaskType_kfalse, MaskType_knonzerofalse, MaskType_kleftright
} MaskType;

namespace {
void *DummyAlloc(size_t size) {
	return nullptr;
}

struct GaussianKernel {
	static LIBSAKURA_SYMBOL(Status) Generate(size_t num_kernel,
			float kernel[]) {
		return LIBSAKURA_SYMBOL(CreateGaussianKernelFloat)(FWHM, num_kernel,
				kernel);
	}

	static float FWHM;
};
float GaussianKernel::FWHM = 0.0;

void Deallocate(void *ptr) {
	free(ptr);
}

struct StandardInitializer {
	static void Initialize(size_t num_kernel, void *storage_ptr,
			float **kernel) {
		*kernel = reinterpret_cast<float *>(storage_ptr);
	}
};

struct NotAlignedInitializer {
	static void Initialize(size_t num_kernel, void *storage_ptr,
			float **kernel) {
		float *aligned_kernel = reinterpret_cast<float *>(storage_ptr);
		*kernel = aligned_kernel + 1;
	}
};

struct NullPointerInitializer {
	static void Initialize(size_t num_kernel, void *storage_ptr,
			float **kernel) {
		*kernel = nullptr;
	}
};

struct InvalidArgumentValidator {
	static void Validate(float kernel_width, size_t num_kernel,
			float const *kernel, sakura_Status const status) {
		// only validate returned status
		EXPECT_EQ(sakura_Status_kInvalidArgument, status);
	}
};

struct NotGoodValidator {
	static void Validate(float kernel_width, size_t num_kernel,
			float const *kernel, sakura_Status const status) {
		// only validate returned status
		EXPECT_EQ(sakura_Status_kNG, status);
	}
};

struct NumKernelOneValidator {
	static void Validate(float kernel_width, size_t num_kernel,
			float const *kernel, sakura_Status const status) {
		// num_kernel must be 1
		ASSERT_EQ(1, num_kernel);

		// execution must be successful
		EXPECT_EQ(sakura_Status_kOK, status);

		EXPECT_EQ(1.0f, kernel[0]);
	}
};

struct StandardValidator {
	static void Validate(float kernel_width, size_t num_kernel,
			float const *kernel, sakura_Status const status) {
		// execution must be successful
		EXPECT_EQ(sakura_Status_kOK, status);

		// validate each kernel elements
		std::unique_ptr<void, decltype(&Deallocate)> storage(
				malloc(num_kernel * sizeof(float)), Deallocate);
		float *expected_kernel = reinterpret_cast<float *>(storage.get());
		double sigma = kernel_width / (2.0 * sqrt(2.0 * log(2.0)));
		double peak_value = 1.0 / (sqrt(2.0 * M_PI) * sigma);
		size_t peak_location = num_kernel / 2;
		double expected_sum = 0.0;
		for (size_t i = 0; i < num_kernel; ++i) {
			auto separation = static_cast<double>(peak_location)
					- static_cast<double>(i);
			auto normalized_separation = separation / sigma;
			expected_kernel[i] = peak_value
					* exp(-normalized_separation * normalized_separation / 2.0);
			expected_sum += expected_kernel[i];
		}
		for (size_t i = 0; i < num_kernel; ++i) {
			float expected = expected_kernel[i] / expected_sum;
			EXPECT_TRUE(std::isfinite(kernel[i]));
			EXPECT_FLOAT_EQ(expected, kernel[i]);
		}

		// validate symmetry
		if (num_kernel > 1) {
			size_t num_tail = peak_location - ((num_kernel % 2 == 0) ? 1 : 0);
			for (size_t i = 0; i < num_tail; ++i) {
				// must match exactly so use EXPECT_EQ instead of EXPECT_FLOAT_EQ
				EXPECT_EQ(kernel[peak_location - i], kernel[peak_location + i]);
			}
		}

		// validate sum
		expected_sum = 1.0;
		double actual_sum = 0.0;
		for (size_t i = 0; i < num_kernel; ++i) {
			actual_sum += kernel[i];
		}
		EXPECT_FLOAT_EQ(expected_sum, actual_sum);
	}
};

struct HalfWidthValidator {
	static void Validate(float kernel_width, size_t num_kernel,
			float const *kernel, sakura_Status const status) {
		// kernel_width must be equal to num_kernel
		ASSERT_EQ(kernel_width, static_cast<float>(num_kernel));

		StandardValidator::Validate(kernel_width, num_kernel, kernel, status);

		// additional check: edge values should be half maximum
		size_t peak_location = num_kernel / 2;
		auto peak_value = kernel[peak_location];
		// correction factor for odd case:
		// in this case, FWHM is slightly shifted with respect to edge
		// so that correction must be needed for verification
		double half_width = kernel_width / 2.0f;
		double actual_width = static_cast<double>(peak_location);
		double sigma = kernel_width / (2.0 * sqrt(2.0 * log(2.0)));
		double correction_factor = 1.0;
		if (num_kernel % 2 > 0) {
			correction_factor = exp(
					-(actual_width * actual_width - half_width * half_width)
							/ (2.0 * sigma * sigma));
		}
		auto half_maximum = static_cast<double>(peak_value) / 2.0f
				* correction_factor;
		std::cout << "half_maximum=" << half_maximum << " (num_kernel "
				<< num_kernel << ", correction factor " << correction_factor
				<< ")" << std::endl;
		EXPECT_FLOAT_EQ(half_maximum, kernel[0]);
		// validate only when num_kernel is odd since even kernel is
		// not symmetric
		if (num_kernel % 2 != 0) {
			EXPECT_FLOAT_EQ(half_maximum, kernel[num_kernel - 1]);
		}
	}
};

struct NarrowKernelValidator {
	static void Validate(float kernel_width, size_t num_kernel,
			float const *kernel, sakura_Status const status) {
		// kernel_width must be equal to num_kernel
		ASSERT_EQ(FLT_MIN, kernel_width);

		StandardValidator::Validate(kernel_width, num_kernel, kernel, status);

		size_t peak_location = num_kernel / 2;
		std::cout << "peak: kernel[" << peak_location << "]="
				<< kernel[peak_location] << std::endl;
		std::cout << "off-peak: kernel[" << peak_location - 1 << "]="
				<< kernel[peak_location - 1] << std::endl;
		std::cout << "off-peak: kernel[" << peak_location + 1 << "]="
				<< kernel[peak_location + 1] << std::endl;

		// kernel is so narrow that value of the kernel except peak_location
		// is zero
		for (size_t i = 0; i < peak_location; ++i) {
			EXPECT_FLOAT_EQ(0.0f, kernel[i]);
		}
		for (size_t i = peak_location + 1; i < num_kernel; ++i) {
			EXPECT_FLOAT_EQ(0.0f, kernel[i]);
		}

	}
};

struct WideKernelValidator {
	static void Validate(float kernel_width, size_t num_kernel,
			float const *kernel, sakura_Status const status) {
		// kernel_width must be equal to num_kernel
		ASSERT_EQ(FLT_MAX, kernel_width);

		StandardValidator::Validate(kernel_width, num_kernel, kernel, status);

		size_t peak_location = num_kernel / 2;
		std::cout << "peak: kernel[" << peak_location << "]="
				<< kernel[peak_location] << std::endl;
		std::cout << "edge: kernel[" << 0 << "]=" << kernel[0] << std::endl;
		std::cout << "edge: kernel[" << num_kernel - 1 << "]="
				<< kernel[num_kernel - 1] << std::endl;

		// kernel is so wide that all the kernel value are the same
		for (size_t i = 0; i < num_kernel; ++i) {
			EXPECT_FLOAT_EQ(kernel[peak_location], kernel[i]);
		}
	}
};

struct NullLogger {
	static void PrintAllocatedSize(size_t storage_size) {
	}
	static void PrintElapsedTime(std::string const test_name,
			double elapsed_time) {
	}
};

struct PerformanceTestLogger {
	static void PrintAllocatedSize(size_t storage_size) {
		std::cout << "size to be allocated: "
				<< static_cast<double>(storage_size) * 1.0e-6 << " MB"
				<< std::endl;
	}
	static void PrintElapsedTime(std::string const test_name,
			double elapsed_time) {
		std::cout << BENCH << test_name << " " << elapsed_time << std::endl;
	}
};

template<typename Initializer, typename Validator, typename Logger = NullLogger>
inline void RunGaussianTest(float const kernel_width, size_t const num_kernel,
		std::string const test_name = "") {
	// prepare storage
	void *storage_ptr = nullptr;
	size_t const alignment = sakura_GetAlignment();
	size_t const storage_size = (num_kernel + 1) * sizeof(float); // margin for not-aligned test case
	Logger::PrintAllocatedSize(storage_size);
	int allocation_status = posix_memalign(&storage_ptr, alignment,
			storage_size);
	ASSERT_EQ(0, allocation_status);
	std::unique_ptr<void, decltype(&Deallocate)> storage(storage_ptr,
			Deallocate);

	// initialize kernel
	float *kernel = nullptr;
	Initializer::Initialize(num_kernel, storage_ptr, &kernel);

	// run
	double start_time = sakura_GetCurrentTime();
	sakura_Status status = sakura_CreateGaussianKernelFloat(kernel_width,
			num_kernel, kernel);
	double end_time = sakura_GetCurrentTime();
	Logger::PrintElapsedTime(test_name, end_time - start_time);

	// verification
	Validator::Validate(kernel_width, num_kernel, kernel, status);
}

} /* namespace */

using namespace std;

/*
 *  a struct to store reference data
 *  offset: the start index of data array to compare with reference,
 *          i.e. data[offset] is compared with reference[0].
 *  num_reference: the number of array elements to compare
 *  reference: an array of reference values
 */
template<typename DataType>
struct ReferenceData {
	size_t offset;
	vector<DataType> reference;

	ReferenceData(size_t offset_index = 0, initializer_list<DataType> ref = { }) :
			offset(offset_index), reference(ref) {
	}
};

template<typename Kernel>
struct ConvolveTestComponent {
	string name;
	Kernel kernel;bool use_fft;
	size_t num_data;
	SpikeType data_type;
	MaskType mask_type;
	vector<ReferenceData<float>> data_ref;
	vector<ReferenceData<float>> weight_ref;

	ConvolveTestComponent(string in_name, Kernel in_kernel, bool use_fft,
			size_t in_num_data, SpikeType in_dtype, MaskType in_mtype,
			initializer_list<ReferenceData<float>> data,
			initializer_list<ReferenceData<float>> weight) :
			name(in_name), use_fft(use_fft), num_data(in_num_data), data_type(
					in_dtype), mask_type(in_mtype), data_ref(data), weight_ref(
					weight) {
		kernel = in_kernel;
	}
};

struct ConvolveGaussKernel {
	float kernel_width;
	size_t num_kernel;
	//bool use_fft;
	void generate(float kernel[/*num_kernel*/]) {
		sakura_Status status = LIBSAKURA_SYMBOL(CreateGaussianKernelFloat)(
				kernel_width, num_kernel, kernel);
		ASSERT_EQ(sakura_Status_kOK, status);
	}
};

struct ConvolveRightAngledTriangleKernel {
	float kernel_width;
	size_t num_kernel;
	//bool use_fft;
	void generate(float kernel[/*num_kernel*/]) {
		float sum = 0.0;
		size_t const start = num_kernel / 2 - kernel_width / 2;
		for (size_t i = 0; i < num_kernel; ++i) {
			float value;
			if (i < start || i >= start + kernel_width) {
				value = 0.0;
			} else {
				value = static_cast<float>(i - start + 1);
			}
			sum += value;
			kernel[i] = value;
		}
		for (size_t i = 0; i < num_kernel; ++i) {
			kernel[i] = kernel[i] / sum;
		}
	}
};

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
		//PrintArray("in2", NUM_IN_EVEN, in2_);
	}
	void PrintArray(string const name, size_t num_in, float *in) {
		cout << name << " = (";
		for (size_t i = 0; i < num_in; ++i) {
			cout << setprecision(10) << in[i] << ", ";
		}
		cout << ")" << endl;
	}
	void PrintArrayCenter(string const name, size_t num_in, float *in) {
		size_t istart = num_in / 2 - 12;
		size_t iend = num_in / 2 + 12;
		cout << name << "[ " << istart << "~" << iend - 1 << " ] = (";
		for (size_t i = istart; i < iend; ++i) {
			cout << setprecision(10) << in[i] << ", ";
		}
		cout << ")" << endl;
	}
	template<typename Kernel>
	void RunBaseTest(size_t input_data_size, SpikeType input_spike_type,
			size_t num_data, bool use_dummy_num_data, size_t num_kernel,
			bool use_fft, float output_data[],
			LIBSAKURA_SYMBOL(Status) expected_status, bool align_check,
			bool verbose, size_t loop_max, double *elapsed_create = nullptr,
			double *elapsed_convolve = nullptr) {
		LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context = nullptr;
		SIMD_ALIGN
		float input_data[input_data_size];
		SIMD_ALIGN
		bool input_mask[ELEMENTSOF(input_data)];
		SIMD_ALIGN
		float output_weight[ELEMENTSOF(input_data)];
		InitializeDataAndMask(input_spike_type, MaskType_ktrue, input_data_size,
				input_data, input_mask);
		double start = sakura_GetCurrentTime();
		SIMD_ALIGN float kernel[num_kernel];
		LIBSAKURA_SYMBOL(Status) status = expected_status;
//		status = LIBSAKURA_SYMBOL(CreateGaussianKernelFloat)(kernel_width,
//				num_kernel, kernel);
		status = Kernel::Generate(num_kernel, kernel);
//		ASSERT_EQ(status, LIBSAKURA_SYMBOL(Status_kOK));
		ASSERT_EQ(status, expected_status);
		if (use_fft) {
			status =
			LIBSAKURA_SYMBOL(CreateConvolve1DContextFFTFloat)(num_kernel,
					kernel, &context);
			ASSERT_EQ(expected_status, status);
		}
		double end = sakura_GetCurrentTime();
		if (use_fft) {
			if (align_check) {
				ASSERT_TRUE(sakura_IsAligned(context->real_array))<< "real_array is not aligned";
				ASSERT_TRUE(sakura_IsAligned(context->real_kernel_array))<< "real_kernel_array is not aligned";
			}
		}
		size_t bad_num_data = num_data;
		if (use_dummy_num_data) {
			bad_num_data = num_data - 1;
			EXPECT_TRUE(num_data != bad_num_data)
					<< "In this test, num_data != bad_num_data is expected"; // assert
		}
		LIBSAKURA_SYMBOL(Status) status_Convolve;
		double start_time = sakura_GetCurrentTime();
		for (size_t i = 0; i < loop_max && expected_status == sakura_Status_kOK;
				++i) {
			if (use_fft) {
				status_Convolve = LIBSAKURA_SYMBOL(Convolve1DFFTFloat)(context,
						num_data, input_data, output_data);
			} else {
				status_Convolve = LIBSAKURA_SYMBOL(Convolve1DFloat)(num_kernel,
						kernel, num_data, input_data, input_mask, output_data,
						output_weight);

			}
			ASSERT_EQ(expected_status, status_Convolve);
		}
		double end_time = sakura_GetCurrentTime();
		if (verbose) {
			std::cout << "num_data : " << num_data << std::endl;
			std::cout << "use_fft : " << use_fft << std::endl;
			if (num_data <= NUM_IN_LARGE) {
				PrintArray("output_data", num_data, output_data);
			} else {
				PrintArrayCenter("output_data", num_data, output_data);
			}
		}
		verbose = false;
		if (use_fft) {
			LIBSAKURA_SYMBOL(Status) status_Destroy =
			LIBSAKURA_SYMBOL(DestroyConvolve1DContextFloat)(context);
			ASSERT_EQ(expected_status, status_Destroy);
		}
		if (loop_max > 1 && expected_status == sakura_Status_kOK) {
			std::cout << "use_fft : " << use_fft << std::endl;
			std::cout << "create elapsed time : " << end - start << "sec\n";
			std::cout << "convolve : " << num_data << "ch, loop = " << loop_max
					<< ", elapsed time : " << end_time - start_time << "sec\n";
		}
		if (elapsed_create) {
			*elapsed_create = end - start;
		}
		if (elapsed_convolve) {
			*elapsed_convolve = end_time - start_time;
		}
	}

	/*
	 * Invoke Convolution using a list of ConvolveTestComponent.
	 */
	template<typename ConvolveKernel>
	void RunConvolveTestComponentList(size_t num_test,
			ConvolveTestComponent<ConvolveKernel> const TestList[]) {
		for (size_t i = 0; i < num_test; ++i) {
			ConvolveTestComponent<ConvolveKernel> const test_component =
					TestList[i];
			cout << "[Test " << test_component.name << "]" << endl;
			size_t const num_data = test_component.num_data;
			SIMD_ALIGN float input_data[num_data];
			SIMD_ALIGN float output_data[ELEMENTSOF(input_data)];
			SIMD_ALIGN bool input_mask[ELEMENTSOF(input_data)];
			SIMD_ALIGN float output_weight[ELEMENTSOF(input_data)];
			RunConvolutionTest(test_component.kernel, test_component.use_fft,
					test_component.data_type, test_component.mask_type,
					test_component.num_data, input_data, input_mask,
					output_data, output_weight, LIBSAKURA_SYMBOL(Status_kOK), 1,
					test_component.data_ref, test_component.weight_ref);
		}
	}

	/*
	 * Run sakura_Convolve1DFloat and compare output data and mask with references.
	 */
	template<typename ConvolveKernel>
	void RunConvolutionTest(ConvolveKernel kernel_param, bool use_fft,
			SpikeType const spike_type, MaskType const mask_type,
			size_t const num_data, float *input_data, bool *input_mask,
			float *output_data, float *output_weight,
			LIBSAKURA_SYMBOL(Status) const convolution_status,
			size_t const num_repeat = 1,
			vector<ReferenceData<float>> const data_ref = { },
			vector<ReferenceData<float>> const weight_ref = { }) {
		// create kernel array
		size_t const num_kernel = kernel_param.num_kernel;
		SIMD_ALIGN float kernel[num_kernel];
		kernel_param.generate(kernel);
		// create context if FFT
		LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context = nullptr;
		if (use_fft) {
			LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL(CreateConvolve1DContextFFTFloat)(num_kernel,
					kernel, &context);
			ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		}
		// initialize data and mask
		InitializeDataAndMask(spike_type, mask_type, num_data, input_data,
				input_mask);
		// convolution
		if (num_repeat > 1) {
			cout << "Iterating convolution for " << num_repeat
					<< " loops. The length of arrays is " << num_data << endl;
		}
		double start = LIBSAKURA_SYMBOL(GetCurrentTime)();
		LIBSAKURA_SYMBOL(Status) exec_status;
		for (size_t i = 0; i < num_repeat; ++i) {
			if (use_fft) {
				exec_status = LIBSAKURA_SYMBOL(Convolve1DFFTFloat)(context,
						num_data, input_data, output_data);
			} else {
				exec_status = LIBSAKURA_SYMBOL(Convolve1DFloat)(num_kernel,
						kernel, num_data, input_data, input_mask, output_data,
						output_weight);
			}
		}
		double end = LIBSAKURA_SYMBOL(GetCurrentTime)();
		if (num_repeat > 1) {
			cout << "#x# benchmark Convolve1DFloat " << end - start << endl;
		}
		ASSERT_EQ(convolution_status, exec_status);
		// destroy context
		if (use_fft) {
			LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL(DestroyConvolve1DContextFloat)(context);
			ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		}
		// verify output data and mask
		if (exec_status == LIBSAKURA_SYMBOL(Status_kOK)) {
			for (size_t i = 0; i < data_ref.size(); ++i) {
				size_t const offset = data_ref[i].offset;
				size_t const num_reference = data_ref[i].reference.size();
				vector<float> const reference = data_ref[i].reference;
				ASSERT_TRUE(offset + num_reference <= num_data);
				for (size_t j = 0; j < num_reference; ++j) {
					EXPECT_FLOAT_EQ(reference[j], output_data[offset + j]);
				}
			}
			for (size_t i = 0; i < weight_ref.size(); ++i) {
				size_t const offset = weight_ref[i].offset;
				size_t const num_reference = weight_ref[i].reference.size();
				vector<float> const reference = weight_ref[i].reference;
				ASSERT_TRUE(offset + num_reference <= num_data);
				for (size_t j = 0; j < num_reference; ++j) {
					EXPECT_EQ(reference[j], output_weight[offset + j]);
				}
			}
		}
	}

	bool verbose;

private:
	void InitializeDataAndMask(SpikeType const spike_type,
			MaskType const mask_type, size_t const num_data, float data[],
			bool mask[]) {
		for (size_t i = 0; i < num_data; ++i) {
			data[i] = 0.0;
			mask[i] = (mask_type != MaskType_kfalse);
		}
		switch (spike_type) {
		case SpikeType_kleft:
			data[0] = 1.0;
			break;
		case SpikeType_kcenter:
			data[num_data / 2] = 1.0;
			break;
		case SpikeType_kright:
			data[num_data - 1] = 1.0;
			break;
		case SpikeType_kleftright:
			data[0] = 1.0;
			data[num_data - 1] = 1.0;
			break;
		case SpikeType_kall:
			data[0] = 1.0;
			data[num_data / 2] = 1.0;
			data[num_data - 1] = 1.0;
			break;
		case SpikeType_knegative:
			data[num_data / 2] = -1.0;
			break;
		}
		switch (mask_type) {
		case MaskType_ktrue:
		case MaskType_kfalse:
			break;
		case MaskType_kleftright:
			mask[0] = false;
			mask[num_data - 1] = false;
			break;
		case MaskType_knonzerofalse:
			for (size_t i = 0; i < num_data; ++i) {
				mask[i] = (data[i] == 0.0);
			}
			break;
		}
	}
}
;

/*
 * Test creating gaussian kernel by sakura_CreateConvolve1DContextFloat
 * RESULT:
 * if num_data = 0,kernel_width = 0, LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian) = UnknownType,
 * CreateConvolve1DContextFloat will return Status_kInvalidArgument)
 */
TEST_F(Convolve1DOperation ,InvalidArguments) {
	{ // num_kernel > INT_MAX
		std::cout << "num_kernel > INT_MAX" << std::endl;
		size_t input_data_size(NUM_IN_EVEN);
		size_t const num_kernel(size_t(INT_MAX) + 1);
		bool const use_dummy_num_data = false;
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		bool const align_check = false;
		bool const use_fft = true;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		float output_data[input_data_size];
		size_t num_data = NUM_IN_EVEN;
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_kernel, use_fft, output_data,
				sakura_Status_kInvalidArgument, align_check, verbose, loop_max);
	}
	{ // num_kernel == 0
		std::cout << "num_kernel == 0" << std::endl;
		size_t const num_kernel(0);
		bool const use_dummy_num_data = false;
		size_t input_data_size(NUM_IN_EVEN);
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		bool const align_check = false;
		bool const use_fft = true;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		size_t num_data = NUM_IN_EVEN;
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_kernel, use_fft, output_data,
				sakura_Status_kInvalidArgument, align_check, verbose, loop_max);
	}
}

/*
 * Test different number of data between context->num_data and num_data
 * with running Convolve1D
 * RESULT:
 * kInvalidArgument status will be returned and then
 * message "num_data does't equal to context->num_data" will be shown.
 * TODO KS memo: this test did not test what it is intended. Implement proper test here.
 */

/*
 * Test Failed malloc of context
 * RESULT:
 * Convolve1D will return Status_kInvalidArgument
 */
TEST(Convolve1DOperationFailed , FailedMallocContext) {
	{
		LIBSAKURA_SYMBOL(Status) status_Init = LIBSAKURA_SYMBOL(Initialize)(
				DummyAlloc, nullptr);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status_Init);
		LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context = nullptr;
		SIMD_ALIGN
		float input_data[NUM_IN_EVEN];
		SIMD_ALIGN
		bool input_mask[ELEMENTSOF(input_data)];
		size_t const num_data(ELEMENTSOF(input_data));
		for (size_t i = 0; i < ELEMENTSOF(input_data); ++i) {
			input_data[i] = 0.0;
			input_mask[i] = true;
		}
		bool use_fft = true; // with FFT
		constexpr size_t kNumKernel = NUM_IN_EVEN;
		SIMD_ALIGN float kernel[kNumKernel];
		LIBSAKURA_SYMBOL(Status) status_Create =
		LIBSAKURA_SYMBOL(CreateConvolve1DContextFFTFloat)(kNumKernel, kernel,
				&context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kNoMemory), status_Create);
		SIMD_ALIGN
		float output_data[num_data];
		SIMD_ALIGN
		bool output_mask[num_data];
		LIBSAKURA_SYMBOL(Status) status_Convolve1d =
		LIBSAKURA_SYMBOL(Convolve1DFFTFloat)(context, num_data, input_data,
				output_data);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Convolve1d);
		LIBSAKURA_SYMBOL(Status) status_Destroy =
		LIBSAKURA_SYMBOL(DestroyConvolve1DContextFloat)(context);
		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status_Destroy);
		LIBSAKURA_SYMBOL(CleanUp)();
	}
}

/*
 * Test creating gaussian kernel by sakura_CreateConvolve1DContextFloat
 * and checking kernel shape,value,even,odd.
 * RESULT:
 * intended gaussian kernel shape,value will be obtained.
 */
TEST_F(Convolve1DOperation , ValidateGaussianKernel) {
	{ // [even],FFT, Gaussian Kernel Shape,input only 1 spike at center
		float gaussian_kernel[11] = { 0.011742971, 0.031861119, 0.069249183,
				0.12056981, 0.16816399, 0.187887, 0.16816399, 0.12056981,
				0.069249183, 0.031861119, 0.011742971 }; // calculated data beforehand
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const kernel_width = static_cast<size_t>(GaussianKernel::FWHM);
		size_t const input_data_size(NUM_IN_EVEN);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = true;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
		for (size_t i = 0; i < kernel_width - 1; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
			EXPECT_FLOAT_EQ(gaussian_kernel[kernel_width + (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
		}
	}
	{ // [odd],FFT, Gaussian Kernel Shape,input only 1 spike at center
		float gaussian_kernel[11] = { 0.011742957, 0.031861119, 0.069249183,
				0.12056981, 0.16816399, 0.18788746, 0.16816399, 0.12056981,
				0.069249183, 0.031861119, 0.011742957 }; // calculated data beforehand
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const kernel_width = static_cast<size_t>(GaussianKernel::FWHM);
		size_t const input_data_size(NUM_IN_ODD);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = true;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
		for (size_t i = 0; i < kernel_width - 1; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
			EXPECT_FLOAT_EQ(gaussian_kernel[kernel_width + (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
		}
		EXPECT_FLOAT_EQ(gaussian_kernel[ELEMENTSOF(gaussian_kernel)/2],
				output_data[(num_data / 2)]);
	}
	{ // [even],without FFT, Gaussian Kernel Shape,input only 1 spike at center
		float gaussian_kernel[11] = { 0.011745105, 0.03186197, 0.069249394,
				0.12056985, 0.16816399, 0.18788747, 0.16816404, 0.1205702,
				0.069251029, 0.031866922, 0.011754746 }; // analytic value by emulating the code
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const kernel_width = static_cast<size_t>(GaussianKernel::FWHM);
		size_t const input_data_size(NUM_IN_EVEN);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = false;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
		for (size_t i = 0; i < kernel_width; ++i) {
			EXPECT_FLOAT_EQ(gaussian_kernel[kernel_width + (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
		}
	}
	{ // [odd],without FFT, Gaussian Kernel Shape,input only 1 spike at center
		float gaussian_kernel[11] = { 0.011745105, 0.031861969, 0.069249393,
				0.12056985, 0.16816399, 0.18788746, 0.16816399, 0.12056985,
				0.069249393, 0.031861969, 0.011745105 }; // analytic value by emulating the code
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const kernel_width = static_cast<size_t>(GaussianKernel::FWHM);
		size_t const input_data_size(NUM_IN_ODD);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = false;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
		for (size_t i = 0; i < kernel_width; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
			EXPECT_FLOAT_EQ(gaussian_kernel[kernel_width + (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
		}
		EXPECT_FLOAT_EQ(gaussian_kernel[ELEMENTSOF(gaussian_kernel)/2],
				output_data[(num_data / 2)]);
	}
}

/*
 * Test convolution against input spike
 * RESULT: same result will be yielded by with/without FFT
 */
TEST_F(Convolve1DOperation , OtherInputDataFFTonoff) {
	{ // [even](FFT) kernel_width > num_data
		size_t const input_data_size(NUM_IN_EVEN);
		GaussianKernel::FWHM = static_cast<float>(NUM_IN_EVEN + 1);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = true;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
		float gaussian_kernel[11] = { 0.0453697890043, 0.0472178347409,
				0.0487070940435, 0.0497995242476, 0.0504667051136,
				0.0506910830736, 0.0504667051136, 0.0497995242476,
				0.0487070940435, 0.0472178347409, 0.0453697890043 };
		for (size_t i = 0; i < 5; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);
			EXPECT_FLOAT_EQ(gaussian_kernel[5 + i],
					output_data[(num_data / 2) + i]);
		}
	}
	{ // [even](Without FFT) kernel_width > num_data
		size_t const input_data_size(NUM_IN_EVEN);
		GaussianKernel::FWHM = static_cast<float>(NUM_IN_EVEN + 1);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = false;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
		float gaussian_kernel[11] = { 0.052494097, 0.052322092, 0.051935662,
				0.051320437, 0.050466707, 0.052084927, 0.053482965, 0.054660225,
				0.05562174, 0.056377983, 0.056944642 }; // analytic values obtained by emulating code
		for (size_t i = 0; i < 5; ++i) {
			EXPECT_FLOAT_EQ(gaussian_kernel[5 + i],
					output_data[(num_data / 2) + i]);
		}
	}
	{ // [odd](FFT) kernel_width > num_data
		size_t const input_data_size(NUM_IN_ODD);
		GaussianKernel::FWHM = static_cast<float>(NUM_IN_ODD + 1);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = true;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
		float gaussian_kernel[11] = { 0.043915681541, 0.0455670394003,
				0.0468942411244, 0.047865845263, 0.0484584420919,
				0.0486575998366, 0.0484584420919, 0.047865845263,
				0.0468942411244, 0.0455670394003, 0.043915681541 };
		for (size_t i = 0; i < 5; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);EXPECT_FLOAT_EQ(
					gaussian_kernel[5 + i], output_data[(num_data / 2) + i]);
		}
	}
	{ // [odd](Without FFT) kernel_width > num_data
		size_t const input_data_size(NUM_IN_ODD);
		GaussianKernel::FWHM = static_cast<float>(NUM_IN_ODD + 1);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = false;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
		float gaussian_kernel[11] = { 0.052354915, 0.052003436, 0.051467945,
				0.050736419, 0.049800864, 0.048657602, 0.049800864, 0.050736419,
				0.051467945, 0.052003436, 0.052354915 }; // analytic values obtained by a script which emulates code

		for (size_t i = 0; i < 5; ++i) {
			EXPECT_FLOAT_EQ(output_data[(num_data / 2) - (i + 1)],
					output_data[(num_data / 2) + (i + 1)]);EXPECT_FLOAT_EQ(
					gaussian_kernel[5 + i], output_data[(num_data / 2) + i]);
		}
	}
	{ // [even],FFT, Gaussian Kernel Shape, 2 spikes
		size_t input_data_size(NUM_IN_EVEN);
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = true;
		size_t loop_max(1);
		bool const verbose = false;
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kleftright,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
	}
	{ // [even],without FFT, Gaussian Kernel Shape, 2 spikes
		size_t input_data_size(NUM_IN_EVEN);
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = false;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kleftright,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
	}
	{ // [odd],FFT, Gaussian Kernel Shape, 2 spikes
		size_t input_data_size(NUM_IN_ODD);
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = true;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kleftright,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
	}
	{ // [odd],without FFT, Gaussian Kernel Shape, 2 spike
		size_t input_data_size(NUM_IN_ODD);
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = false;
		bool const verbose = false;
		size_t loop_max(1);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kleftright,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
	}
}

/*
 * Test Compare result of with/without FFTW for gaussian kernel
 * by sakura_CreateConvolve1DContextFloat
 * RESULT:
 * each result will be equal between convolution with/without fft
 */
TEST_F(Convolve1DOperation , CompareResultWithFFTWithoutFFT) {
	auto get_num_kernel = [](size_t kernel_width) {
		double const six_sigma = 6.0 / sqrt(8.0 * log(2.0));
		size_t const threshold = static_cast<size_t>(kernel_width * six_sigma);
		return 2 * threshold + 1;
	};
	{ // simple compare shape at near the center
		size_t input_data_size(NUM_IN_EVEN);
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const kernel_width = static_cast<size_t>(GaussianKernel::FWHM);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool use_fft = false;
		bool verbose = false;
		size_t loop_max(1);
		size_t num_kernel = get_num_kernel(kernel_width);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_kernel, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
		SIMD_ALIGN
		float output_data_fft[input_data_size];
		use_fft = true;
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_data, use_fft,
				output_data_fft, sakura_Status_kOK, align_check, verbose,
				loop_max);
		for (size_t i = 0; i < input_data_size; ++i) {
			if ((i > num_data / 2 - kernel_width / 2)
					&& (i < num_data / 2 + kernel_width / 2)) {
				EXPECT_FLOAT_EQ(output_data_fft[i], output_data[i]);
			}
		}
	}
	{ // compare edge value with 2 spikes ( fft > without fft)
		size_t input_data_size(NUM_IN_EVEN);
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const kernel_width = static_cast<size_t>(GaussianKernel::FWHM);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool use_fft = false;
		bool verbose = false;
		size_t loop_max(1);
		size_t num_kernel = get_num_kernel(kernel_width);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kleftright,
				num_data, use_dummy_num_data, num_kernel, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
		SIMD_ALIGN
		float output_data_fft[input_data_size];
		use_fft = true;
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kleftright,
				num_data, use_dummy_num_data, num_data, use_fft,
				output_data_fft, sakura_Status_kOK, align_check, verbose,
				loop_max);
		bool gthan = false;
		if (output_data_fft[0] > output_data[0]
				&& output_data_fft[input_data_size - 1]
						> output_data[input_data_size - 1]) {
			gthan = true;
		}
		EXPECT_EQ(true, gthan);
		EXPECT_FLOAT_EQ(0.31633885, output_data[0]);
		EXPECT_FLOAT_EQ(0.31633885, output_data[input_data_size - 1]);
	}
	{ // compare edge value with 3 spikes ( fft > without fft)
		size_t input_data_size(NUM_IN_EVEN);
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const kernel_width = static_cast<size_t>(GaussianKernel::FWHM);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool use_fft = false;
		bool verbose = false;
		size_t loop_max(1);
		size_t num_kernel = get_num_kernel(kernel_width);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kall, num_data,
				use_dummy_num_data, num_kernel, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max);
		SIMD_ALIGN
		float output_data_fft[input_data_size];
		use_fft = true;
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kall, num_data,
				use_dummy_num_data, num_data, use_fft, output_data_fft,
				sakura_Status_kOK, align_check, verbose, loop_max);
		bool gthan = false;
		if (output_data_fft[0] > output_data[0]
				&& output_data_fft[input_data_size - 1]
						> output_data[input_data_size - 1]) {
			gthan = true;
		}
		EXPECT_EQ(true, gthan);
		EXPECT_FLOAT_EQ(0.316338854, output_data[0]);
		EXPECT_FLOAT_EQ(0.316339299, output_data[input_data_size - 1]);
		EXPECT_FLOAT_EQ(0.187887762, output_data[input_data_size / 2]);
		EXPECT_FLOAT_EQ(output_data_fft[input_data_size / 2],
				output_data[input_data_size / 2]);
	}
}

/*
 * Test [Without FFT] Performance of Convolve1D
 * RESULT: elapsed time will be output
 */
TEST_F(Convolve1DOperation , PerformanceTestWithoutFFT) {
	{ // [even],without FFT, Gaussian Kernel Shape,input delta
		auto get_num_kernel =
				[](size_t kernel_width) {
					double const six_sigma = 6.0 / sqrt(8.0 * log(2.0));
					size_t const threshold = static_cast<size_t>(kernel_width * six_sigma);
					return 2 * threshold + 1;
				};
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const kernel_width = static_cast<size_t>(GaussianKernel::FWHM);
		size_t const input_data_size(NUM_IN_MAX);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = false;
		bool const verbose = false;
		size_t loop_max(10000);
		size_t num_kernel = get_num_kernel(kernel_width);
		SIMD_ALIGN
		float output_data[input_data_size];
		double elapsed_create = 0.0;
		double elapsed_convolve = 0.0;
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_kernel, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max,
				&elapsed_create, &elapsed_convolve);
		char const *test_name = "Convolve1DOperation_PerformanceTestWithoutFFT";
		std::cout << BENCH << test_name << "_Create " << elapsed_create
				<< std::endl;
		std::cout << BENCH << test_name << "_Convolve " << elapsed_convolve
				<< std::endl;
	}
}

/*
 * Test [With FFT] Performance of Convolve1D
 * RESULT: elapsed time will be output
 */
TEST_F(Convolve1DOperation , PerformanceTestWithFFT) {
	double elapsed_create = 0.0;
	double elapsed_convolve = 0.0;
	char const *test_name = "Convolve1DOperation_PerformanceTestWithFFT";
	{ // [even], FFT, Gaussian Kernel Shape,input delta
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const input_data_size(NUM_IN_MAX);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = true;
		bool const verbose = false;
		size_t loop_max(10000);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max,
				&elapsed_create, &elapsed_convolve);

		std::cout << BENCH << test_name << "_Even_Create " << elapsed_create
				<< std::endl;
		std::cout << BENCH << test_name << "_Even_Convolve " << elapsed_convolve
				<< std::endl;
	}
	{ // [odd], FFT, Gaussian Kernel Shape,input delta
		GaussianKernel::FWHM = static_cast<float>(NUM_WIDTH);
		size_t const input_data_size(NUM_IN_MAX - 1);
		size_t const num_data(input_data_size);
		bool const use_dummy_num_data = false;
		bool const align_check = false;
		bool const use_fft = true;
		bool const verbose = false;
		size_t loop_max(10000);
		SIMD_ALIGN
		float output_data[input_data_size];
		RunBaseTest<GaussianKernel>(input_data_size, SpikeType_kcenter,
				num_data, use_dummy_num_data, num_data, use_fft, output_data,
				sakura_Status_kOK, align_check, verbose, loop_max,
				&elapsed_create, &elapsed_convolve);

		std::cout << BENCH << test_name << "_Odd_Create " << elapsed_create
				<< std::endl;
		std::cout << BENCH << test_name << "_Odd_Convolve " << elapsed_convolve
				<< std::endl;
	}
}

/*
 * Test user defined asymmetric kernel (right angled triangle)
 */

TEST_F(Convolve1DOperation, NumKernelWithoutFFTTriangle) {
	ConvolveTestComponent<ConvolveRightAngledTriangleKernel> TriangleNumKernelTest[] =
			{ { "num_kernel(odd) = num_data", { 4, NUM_IN_ODD }, false,
			NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { { 10, { 0.1, 0.2,
					0.3, 0.4 } } }, { } }, { "num_kernel(even) = num_data", { 4,
			NUM_IN_EVEN }, false,
			NUM_IN_EVEN, SpikeType_kcenter, MaskType_ktrue, { { 10, { 0.1, 0.2,
					0.3, 0.4 } } }, { } }, { "num_kernel(odd) > num_data", { 4,
					2 * NUM_IN_ODD + 1 }, false,
			NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { { 10, { 0.1, 0.2,
					0.3, 0.4 } } }, { } }, { "num_kernel(even) > num_data", { 4,
					2 * NUM_IN_ODD }, false,
			NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { { 10, { 0.1, 0.2,
					0.3, 0.4 } } }, { } }, { "num_kernel(odd) < num_data", { 4,
					13 }, false, NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue,
					{ { 10, { 0.1, 0.2, 0.3, 0.4 } } }, { } }, {
					"num_kernel(even) < num_data", { 4, 12 }, false,
					NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { { 10, {
							0.1, 0.2, 0.3, 0.4 } } }, { } } };
	RunConvolveTestComponentList<ConvolveRightAngledTriangleKernel>(
			ELEMENTSOF(TriangleNumKernelTest), TriangleNumKernelTest);
}

/**
 * Test for sakura_CreateGaussianKernel
 */
// T-001
TEST_GAUSS(NegativeKernelWidth) {
	RunGaussianTest<StandardInitializer, InvalidArgumentValidator>(-1.0f, 10);
}

// T-002
TEST_GAUSS(ZeroKernelWidth) {
	RunGaussianTest<StandardInitializer, InvalidArgumentValidator>(0.0f, 10);
}

// T-003
TEST_GAUSS(PositiveInfiniteKernelWidth) {
	constexpr float kOne = 1.0f;
	constexpr float kZero = 0.0f;
	float positive_inf = kOne / kZero;
	ASSERT_TRUE(std::isinf(positive_inf));
	ASSERT_GT(positive_inf, FLT_MAX);
	RunGaussianTest<StandardInitializer, InvalidArgumentValidator>(positive_inf,
			10);
}

// T-004
TEST_GAUSS(NegativeInfiniteKernelWidth) {
	constexpr float kOne = 1.0f;
	constexpr float kZero = 0.0f;
	float negative_inf = -kOne / kZero;
	ASSERT_TRUE(std::isinf(negative_inf));
	ASSERT_LT(negative_inf, FLT_MIN);
	RunGaussianTest<StandardInitializer, InvalidArgumentValidator>(negative_inf,
			10);
}

// T-005
TEST_GAUSS(NotANumberKernelWidth) {
	float nan_value = std::nan("");
	ASSERT_TRUE(std::isnan(nan_value));
	RunGaussianTest<StandardInitializer, InvalidArgumentValidator>(nan_value,
			10);
}

// T-006
TEST_GAUSS(KernelIsNull) {
	RunGaussianTest<NullPointerInitializer, InvalidArgumentValidator>(2.0f, 10);
}

// T-007
TEST_GAUSS(KernelIsNotAligned) {
	RunGaussianTest<NotAlignedInitializer, InvalidArgumentValidator>(2.0f, 10);
}

// T-008
TEST_GAUSS(NumKernelIsZero) {
	RunGaussianTest<StandardInitializer, InvalidArgumentValidator>(2.0f, 0);
}

// T-009
TEST_GAUSS(NumKernelIsOne) {
// result should be independent of kernel_width if num_kernel is 1
	RunGaussianTest<StandardInitializer, NumKernelOneValidator>(2.0f, 1);

	RunGaussianTest<StandardInitializer, NumKernelOneValidator>(256.0f, 1);
}

// T-010
TEST_GAUSS(NumKernelIsTwo) {
	RunGaussianTest<StandardInitializer, StandardValidator>(2.0f, 2);
}

// T-011
TEST_GAUSS(NumKernelIsThree) {
	RunGaussianTest<StandardInitializer, StandardValidator>(2.0f, 3);
}

// T-012
TEST_GAUSS(NumKernelIsFour) {
	RunGaussianTest<StandardInitializer, StandardValidator>(2.0f, 4);
}

// T-014
TEST_GAUSS(EvenNumkernel) {
	RunGaussianTest<StandardInitializer, StandardValidator>(24.0f, 64);
}

// T-013
TEST_GAUSS(OddNumkernel) {
	RunGaussianTest<StandardInitializer, StandardValidator>(24.0f, 65);
}

// T-016
TEST_GAUSS(NumkernelLessThanKernelWidth) {
// even
	RunGaussianTest<StandardInitializer, StandardValidator>(128.0f, 64);

	// odd
	RunGaussianTest<StandardInitializer, StandardValidator>(128.0f, 65);
}

// T-017
TEST_GAUSS(NumKernelGreaterThanKernelWidth) {
// even
	RunGaussianTest<StandardInitializer, StandardValidator>(64.0f, 128);

	// odd
	RunGaussianTest<StandardInitializer, StandardValidator>(64.0f, 129);
}

// T-018
TEST_GAUSS(HalfWidth) {
// even
	RunGaussianTest<StandardInitializer, HalfWidthValidator>(128.0f, 128);

	// odd
	RunGaussianTest<StandardInitializer, HalfWidthValidator>(129.0f, 129);
}

// T-019
TEST_GAUSS(NarrowKernelWidth) {
	RunGaussianTest<StandardInitializer, NarrowKernelValidator>(FLT_MIN, 128);
}

// T-020
TEST_GAUSS(WideKernelWidth) {
	RunGaussianTest<StandardInitializer, WideKernelValidator>(FLT_MAX, 128);
}

// T-021
TEST_GAUSS(PowerOfTwoNumKernelPerformance) {
	constexpr size_t kNumKernel = 268435456; // 2^28
	RunGaussianTest<StandardInitializer, StandardValidator,
			PerformanceTestLogger>(128.0f, kNumKernel,
			"CreateGaussian_PowerOfTwoNumKernelPerformanceTest");
}

// T-022
TEST_GAUSS(OddNumKernelPerformance) {
	constexpr size_t kNumKernel = 268435456 + 1; // 2^28 + 1
	RunGaussianTest<StandardInitializer, StandardValidator,
			PerformanceTestLogger>(128.0f, kNumKernel,
			"CreateGaussian_OddNumKernelPerformanceTest");
}

/**
 * Test for LIBSAKURA_SYMBOL(Convolve1DFloat)
 */
