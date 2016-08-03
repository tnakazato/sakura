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
#include <iomanip>
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

void Deallocate(void *ptr) {
	free(ptr);
}

struct GaussianKernel {
	float kernel_width;
	size_t num_kernel;
	//bool use_fft;
	void Generate(float kernel[/*num_kernel*/]) {
		sakura_Status status = LIBSAKURA_SYMBOL(CreateGaussianKernelFloat)(
				kernel_width, num_kernel, kernel);
		ASSERT_EQ(sakura_Status_kOK, status);
		if (false) {
			std::cout << "Gaussian kernel [" << num_kernel << "] = {";
			for (size_t i = 0; i < num_kernel; ++i) {
				std::cout << std::setprecision(8) << kernel[i] << ", ";
			}
			std::cout << "}" << std::endl;
		}
	}
};

struct RightAngledTriangleKernel {
	float kernel_width;
	size_t num_kernel;
	//bool use_fft;
	void Generate(float kernel[/*num_kernel*/]) {
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

struct SakuraOkInitializer {
	static LIBSAKURA_SYMBOL(Status) Initialize() {
		return LIBSAKURA_SYMBOL(Initialize)(nullptr, nullptr);
	}
};

struct SakuraNoMemoryInitializer {
	static LIBSAKURA_SYMBOL(Status) Initialize() {
		return LIBSAKURA_SYMBOL(Initialize)(DummyAlloc, nullptr);
	}
};

struct StandardArrayInitializer {
	static void Initialize(size_t /*num_kernel*/, void *storage_ptr,
			float **kernel) {
		*kernel = reinterpret_cast<float *>(storage_ptr);
	}
};

struct NotAlignedArrayInitializer {
	static void Initialize(size_t /*num_kernel*/, void *storage_ptr,
			float **kernel) {
		float *aligned_kernel = reinterpret_cast<float *>(storage_ptr);
		*kernel = aligned_kernel + 1;
	}
};

struct NullPointerArrayInitializer {
	static void Initialize(size_t /*num_kernel*/, void *storage_ptr,
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
	static void Validate(float /*kernel_width*/, size_t /*num_kernel*/,
			float const */*kernel*/, sakura_Status const status) {
		// only validate returned status
		EXPECT_EQ(sakura_Status_kNG, status);
	}
};

struct NoMemoryValidator {
	static void Validate(float /*kernel_width*/, size_t /*num_kernel*/,
			float const */*kernel*/, sakura_Status const status) {
		// only validate returned status
		EXPECT_EQ(sakura_Status_kNoMemory, status);
	}
};

struct StatusOKValidator {
	static void Validate(float /*kernel_width*/, size_t /*num_kernel*/,
			float const */*kernel*/, sakura_Status const status) {
		// execution must be successful
		EXPECT_EQ(sakura_Status_kOK, status);
	}
};

struct NumKernelOneValidator {
	static void Validate(float /*kernel_width*/, size_t num_kernel,
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

struct ContextValidator {
	static void Validate(size_t num_kernel, float const *kernel,
			LIBSAKURA_SYMBOL(Convolve1DContextFloat) const *context) {
		ASSERT_EQ(context->kernel_width, num_kernel);
//		float *flipped_kernel = nullptr;
//		std::unique_ptr<void, DefaultAlignedMemory> kernel_storage(
//				DefaultAlignedMemory::AlignedAllocateOrException(
//						sizeof(float) * num_kernel, &flipped_kernel));
//		size_t const nelem[1] = { num_kernel };
//		// KS TODO: workaround 1 element shift between FlipArrayFloat and FlipData1D
//		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(FlipArrayFloat)(
//				false, 1, nelem, kernel, flipped_kernel);
//		ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
//		// KS TODO: FFT flipped_kernel
//		for (size_t i = 0; i < num_kernel; ++i) {
//			EXPECT_FLOAT_EQ(flipped_kernel[i], context->real_kernel_array[i]);
//		}
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

/*
 * Test creation and destroy of context
 * Note: destroy_status is evaluated only when creation is successful
 */
template<typename SakuraInitializer, typename KernelInitializer,
		typename CreateValidator, typename DestroyValidator>
inline void RunContextTest(size_t const num_kernel,
		bool const valid_context_pointer = true) {
	//initialize sakura
	LIBSAKURA_SYMBOL(Status) init_status = SakuraInitializer::Initialize();
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), init_status);
	float *dummy_kernel = nullptr;
	// avoid allocating storage > INT_MAX
	size_t const num_dummy =
			(num_kernel <= INT_MAX) ? num_kernel + 1 : NUM_IN_ODD;
	std::unique_ptr<void, DefaultAlignedMemory> kernel_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(float) * num_dummy, &dummy_kernel));
	ASSERT_NE(dummy_kernel, nullptr);
	// set some values to kernel
	RightAngledTriangleKernel kernel_generator = {
			static_cast<float>(num_dummy), num_dummy };
	kernel_generator.Generate(dummy_kernel);
	// initialize kernel
	float *kernel = nullptr;
	KernelInitializer::Initialize(num_kernel, dummy_kernel, &kernel);
	// initialize context
	LIBSAKURA_SYMBOL(Convolve1DContextFloat) **context_ptr_ptr = nullptr;
	LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context_ptr = nullptr;
	if (valid_context_pointer) {
		context_ptr_ptr = &context_ptr;
	}
	// create context
	LIBSAKURA_SYMBOL(Status) create_status =
			LIBSAKURA_SYMBOL(CreateConvolve1DContextFFTFloat)(num_kernel,
					kernel, context_ptr_ptr);
	//std::cout << "create kernel status = " << create_status << std::endl;
	constexpr float kDummyKernelWidth = 0;
	CreateValidator::Validate(kDummyKernelWidth, num_kernel, kernel,
			create_status);
	// context should be nullptr when context creation failed
	if (create_status != LIBSAKURA_SYMBOL(Status_kOK)
			&& valid_context_pointer) {
		ASSERT_EQ(nullptr, *context_ptr_ptr);
	} else if (create_status == LIBSAKURA_SYMBOL(Status_kOK)) {
		ContextValidator::Validate(num_kernel, kernel, *context_ptr_ptr);
	}
	// destroy context
	if (context_ptr_ptr != nullptr) {
		LIBSAKURA_SYMBOL(Status) status =
				LIBSAKURA_SYMBOL(DestroyConvolve1DContextFloat)(
						*context_ptr_ptr);
		//std::cout << "destroy status = " << status << std::endl;
		DestroyValidator::Validate(kDummyKernelWidth, num_kernel, kernel,
				status);
	}
	// Clean-up sakura
	LIBSAKURA_SYMBOL(CleanUp)();
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
	Kernel kernel;
	bool use_fft;
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

	/*
	 * Invoke Convolution using a list of ConvolveTestComponent.
	 */
	template<typename ConvolveKernel>
	void RunConvolveTestComponentList(size_t num_test,
			ConvolveTestComponent<ConvolveKernel> const TestList[],
			size_t const num_repeat = 1) {
		for (size_t i = 0; i < num_test; ++i) {
			ConvolveTestComponent<ConvolveKernel> const test_component =
					TestList[i];
			if (num_test > 1) {
				cout << "[Test " << test_component.name << " w/"
						<< (test_component.use_fft ? "" : "o") << " FFT]" << endl;
			}
			size_t const num_data = test_component.num_data;
			SIMD_ALIGN float input_data[num_data];
			SIMD_ALIGN float output_data[ELEMENTSOF(input_data)];
			SIMD_ALIGN bool input_mask[ELEMENTSOF(input_data)];
			SIMD_ALIGN float output_weight[ELEMENTSOF(input_data)];
			RunConvolutionTest(test_component.kernel, test_component.use_fft,
					test_component.data_type, test_component.mask_type,
					test_component.num_data, input_data, input_mask,
					output_data, output_weight, LIBSAKURA_SYMBOL(Status_kOK),
					num_repeat, test_component.data_ref,
					test_component.weight_ref);
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
		kernel_param.Generate(kernel);
		// create context if FFT
		LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context = nullptr;
		if (use_fft) {
			LIBSAKURA_SYMBOL(Status) status =
					LIBSAKURA_SYMBOL(CreateConvolve1DContextFFTFloat)(
							num_kernel, kernel, &context);
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
		double start, end;
		LIBSAKURA_SYMBOL(Status) exec_status;
		if (use_fft) {
			start = LIBSAKURA_SYMBOL(GetCurrentTime)();
			for (size_t i = 0; i < num_repeat; ++i) {
				exec_status = LIBSAKURA_SYMBOL(Convolve1DFFTFloat)(context,
						num_data, input_data, output_data);
			}
			end = LIBSAKURA_SYMBOL(GetCurrentTime)();
		} else {
			start = LIBSAKURA_SYMBOL(GetCurrentTime)();
			for (size_t i = 0; i < num_repeat; ++i) {
				exec_status = LIBSAKURA_SYMBOL(Convolve1DFloat)(num_kernel,
						kernel, num_data, input_data, input_mask, output_data,
						output_weight);
			}
			end = LIBSAKURA_SYMBOL(GetCurrentTime)();
		}
		if (num_repeat > 1) {
			cout << "#x# benchmark Convolve1DFloat_With"
					<< (use_fft ? "" : "Out") << "FFT_Data"
					<< (num_data % 2 == 0 ? "Even" : "Odd") << " "
					<< end - start << endl;
		}
		ASSERT_EQ(convolution_status, exec_status);
		// destroy context
		if (use_fft) {
			LIBSAKURA_SYMBOL(Status) status =
					LIBSAKURA_SYMBOL(DestroyConvolve1DContextFloat)(context);
			ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		}
		// verify output data and weight
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
 * Test failure cases of context creation and destroy.
 * (context=nullptr after execution of create function)
 */
// create T-001, destroy T-001
TEST(ContextTest, NumKernelIsZero) {
	RunContextTest<SakuraOkInitializer, StandardArrayInitializer,
			InvalidArgumentValidator, InvalidArgumentValidator>(0);
}
//create T-002, destroy T-001
TEST(ContextTest, NumKernelIsGreaterThanIntMax) {
	RunContextTest<SakuraOkInitializer, StandardArrayInitializer,
			InvalidArgumentValidator, InvalidArgumentValidator>(
			(size_t) INT_MAX + 1);
}
//create T-003, destroy T-001
TEST(ContextTest, KernelIsNull) {
	RunContextTest<SakuraOkInitializer, NullPointerArrayInitializer,
			InvalidArgumentValidator, InvalidArgumentValidator>(25);
}
//create T-004, destroy T-001
TEST(ContextTest, KernelIsNotAligned) {
	RunContextTest<SakuraOkInitializer, NotAlignedArrayInitializer,
			InvalidArgumentValidator, InvalidArgumentValidator>(25);
}
//create T-005, destroy T-001
TEST(ContextTest, ContextIsNullPtr) {
	RunContextTest<SakuraOkInitializer, StandardArrayInitializer,
			InvalidArgumentValidator, InvalidArgumentValidator>(25, false);
}
//create T-006, destroy T-001
TEST(ContextTest, KernelNoMemory) {
	// context generation fails with Status_kNoMemory.
	RunContextTest<SakuraNoMemoryInitializer, StandardArrayInitializer,
			NoMemoryValidator, InvalidArgumentValidator>(25);
}

/*
 * Successful context creation and destroy
 * num_kernel = UB is not tested because of valgrind memory consumption issue.
 */
//create T-007, destroy T-002
TEST(ContextTest, NumKernelIsLB) {
	RunContextTest<SakuraOkInitializer, StandardArrayInitializer,
			StatusOKValidator, StatusOKValidator>(1);
}
//create T-009, destroy T-002
TEST(ContextTest, NumKernelIsIR) {
	RunContextTest<SakuraOkInitializer, StandardArrayInitializer,
			StatusOKValidator, StatusOKValidator>(25);
}

/**********************************************************************
 * Test Convolution operation with FFT
 *********************************************************************/
// convolve w/ FFT T-001
TEST_F(Convolve1DOperation , ContextIsNull) {
	LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context = nullptr;
	size_t const num_data = 1;
	SIMD_ALIGN float input_data[num_data];
	SIMD_ALIGN float output_data[ELEMENTSOF(input_data)];
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Convolve1DFFTFloat)(
			context, num_data, input_data, output_data);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
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
 * Convolution by narrow Gaussian kernel (kernel_width < num_data)
 */
// convolve w/ FFT T-012 (num_data=odd)
TEST_F(Convolve1DOperation , GaussWithFFTDataIsOdd) {
	initializer_list<float> gaussian_kernel = initializer_list<float>( {
			0.011742966, 0.031861115, 0.069249175, 0.12056981, 0.16816399,
			0.18788746, 0.16816399, 0.12056981, 0.069249175, 0.031861115,
			0.011742966 }); // Kernel values
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (odd)", { kernel_width, num_data }, true, num_data,
			SpikeType_kcenter, MaskType_ktrue, { { num_data / 2
					- gaussian_kernel.size() / 2, gaussian_kernel } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);
}
// convolve w/ FFT T-013 (num_data=even)
TEST_F(Convolve1DOperation , GaussWithFFTDataIsEven) {
	initializer_list<float> gaussian_kernel = initializer_list<float>( {
			0.011742966, 0.031861115, 0.069249175, 0.12056981, 0.16816399,
			0.18788746, 0.16816399, 0.12056981, 0.069249175, 0.031861115,
			0.011742966 }); // Kernel values
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_EVEN;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (even)", { kernel_width, num_data }, true,
			num_data, SpikeType_kcenter, MaskType_ktrue, { { num_data / 2
					- gaussian_kernel.size() / 2, gaussian_kernel } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);
}
/*
 * Convolution by wide Gaussian kernel (kernel_width > num_data)
 */
// convolve w/ FFT T-014 (num_data=odd)
TEST_F(Convolve1DOperation , WideGaussWithFFTDataIsOdd) {
	initializer_list<float> gaussian_kernel = initializer_list<float>( {
			0.043915681541, 0.0455670394003, 0.0468942411244, 0.047865845263,
			0.0484584420919, 0.0486575998366, 0.0484584420919, 0.047865845263,
			0.0468942411244, 0.0455670394003, 0.043915681541 }); // Kernel values
	float const kernel_width = NUM_IN_ODD + 1;
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"wide Gaussian kernel (odd)", { kernel_width, num_data }, true,
			num_data, SpikeType_kcenter, MaskType_ktrue, { { num_data / 2
					- gaussian_kernel.size() / 2, gaussian_kernel } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);
}
// convolve w/ FFT T-015 (num_data=even)
TEST_F(Convolve1DOperation , WideGaussWithFFTDataIsEven) {
	initializer_list<float> gaussian_kernel = initializer_list<float>( {
			0.0453697890043, 0.0472178347409, 0.0487070940435, 0.0497995242476,
			0.0504667051136, 0.0506910830736, 0.0504667051136, 0.0497995242476,
			0.0487070940435, 0.0472178347409, 0.0453697890043 }); // Kernel values
	float const kernel_width = NUM_IN_EVEN + 1;
	size_t const num_data = NUM_IN_EVEN;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"wide Gaussian kernel (even)", { kernel_width, num_data }, true,
			num_data, SpikeType_kcenter, MaskType_ktrue, { { num_data / 2
					- gaussian_kernel.size() / 2, gaussian_kernel } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);
}

/*
 * Convolution by Gaussian kernel. Data has 2 spikes at edges
 */
// convolve w/ FFT T-016 (num_data=odd)
TEST_F(Convolve1DOperation , TwoSpikesWithFFTDataIsOdd) {
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (odd) edge spikes data",
			{ kernel_width, num_data }, true, num_data, SpikeType_kleftright,
			MaskType_ktrue, { }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);
}
// convolve w/ FFT T-017 (num_data=even)
TEST_F(Convolve1DOperation , TwoSpikesWithFFTDataIsEven) {
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_EVEN;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (even) edge spikes data",
			{ kernel_width, num_data }, true, num_data, SpikeType_kleftright,
			MaskType_ktrue, { }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);
}
/*
 * Convolution by Gaussian kernel. Data has 3 spikes at center and edges
 */
// convolve w/ FFT T-019  (num_data = even)
TEST_F(Convolve1DOperation , ThreeSpikesWithFFTDataIsEven) {
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_EVEN;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (even) 3 spikes data", { kernel_width, num_data },
			true, num_data, SpikeType_kall, MaskType_ktrue, { { num_data / 2, {
					0.187887762 } } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);
}
/*
 * Test Performance of Convolve1D
 */
// convolve w/ FFT T-023 (num_data = even)
TEST_F(Convolve1DOperation , PerformanceTestWithFFTDataIsEven) {
	initializer_list<float> gaussian_kernel = initializer_list<float>( {
			0.011742966, 0.031861115, 0.069249175, 0.12056981, 0.16816399,
			0.18788746, 0.16816399, 0.12056981, 0.069249175, 0.031861115,
			0.011742966 }); // Kernel values
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_MAX;
	size_t num_repeat = 10000;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (even)", { kernel_width, num_data }, true,
			num_data, SpikeType_kcenter, MaskType_ktrue, { { num_data / 2
					- gaussian_kernel.size() / 2, gaussian_kernel } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest, num_repeat);
}

// convolve w/ FFT T-022 (num_data = odd)
TEST_F(Convolve1DOperation , PerformanceTestWithFFTDataIsOdd) {
	initializer_list<float> gaussian_kernel = initializer_list<float>( {
			0.011742966, 0.031861115, 0.069249175, 0.12056981, 0.16816399,
			0.18788746, 0.16816399, 0.12056981, 0.069249175, 0.031861115,
			0.011742966 }); // Kernel values
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_MAX - 1;
	size_t num_repeat = 10000;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (odd)", { kernel_width, num_data }, true, num_data,
			SpikeType_kcenter, MaskType_ktrue, { { num_data / 2
					- gaussian_kernel.size() / 2, gaussian_kernel } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest, num_repeat);

}

/**********************************************************************
 * Test Convolution operation withOUT FFT
 *********************************************************************/
/*
 * Convolution by narrow Gaussian kernel (kernel_width < num_data)
 */
// convolve w/o FFT T-018 (num_data = odd)
TEST_F(Convolve1DOperation , GaussWithOutFFTDataIsOdd) {
	initializer_list<float> gaussian_kernel = initializer_list<float>( {
			0.011745105, 0.031861969, 0.069249393, 0.12056985, 0.16816399,
			0.18788746, 0.16816399, 0.12056985, 0.069249393, 0.031861969,
			0.011745105 }); // analytic value by emulating the code
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (odd)", { kernel_width, num_data }, false,
			num_data, SpikeType_kcenter, MaskType_ktrue, { { num_data / 2
					- gaussian_kernel.size() / 2, gaussian_kernel } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);

}
// convolve w/o FFT T-019 (num_data = even)
TEST_F(Convolve1DOperation , GaussWithOutFFTDataIsEven) {
	initializer_list<float> gaussian_kernel = initializer_list<float>( {
			0.011745105, 0.03186197, 0.069249394, 0.12056985, 0.16816399,
			0.18788747, 0.16816404, 0.1205702, 0.069251029, 0.031866922,
			0.011754746 }); // analytic value by emulating the code
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_EVEN;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (even)", { kernel_width, num_data }, false,
			num_data, SpikeType_kcenter, MaskType_ktrue, { { num_data / 2
					- gaussian_kernel.size() / 2, gaussian_kernel } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);
}

/*
 * Convolution by wide Gaussian kernel (kernel_width > num_data)
 */
// convolve w/o FFT T-024 (num_data = odd)
TEST_F(Convolve1DOperation , WideGaussWithOutFFTDataIsOdd) {
	initializer_list<float> gaussian_kernel = initializer_list<float>( {
			0.052354915, 0.052003436, 0.051467945, 0.050736419, 0.049800864,
			0.048657602, 0.049800864, 0.050736419, 0.051467945, 0.052003436,
			0.052354915 }); // analytic values obtained by a script which emulates code
	float const kernel_width = NUM_IN_ODD + 1;
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"wide Gaussian kernel (odd)", { kernel_width, num_data }, false,
			num_data, SpikeType_kcenter, MaskType_ktrue, { { num_data / 2
					- gaussian_kernel.size() / 2, gaussian_kernel } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);
}
// convolve w/o FFT T-025 (num_data = even)
TEST_F(Convolve1DOperation , WideGaussWithOutFFTDataIsEven) {
	initializer_list<float> gaussian_kernel = initializer_list<float>( {
			0.052494097, 0.052322092, 0.051935662, 0.051320437, 0.050466707,
			0.052084927, 0.053482965, 0.054660225, 0.05562174, 0.056377983,
			0.056944642 }); // analytic values obtained by emulating code
	float const kernel_width = NUM_IN_EVEN + 1;
	size_t const num_data = NUM_IN_EVEN;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"wide Gaussian kernel (even)", { kernel_width, num_data }, false,
			num_data, SpikeType_kcenter, MaskType_ktrue, { { num_data / 2
					- gaussian_kernel.size() / 2, gaussian_kernel } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);
}
/*
 * Convolution by Gaussian kernel. Data has 2 spikes at edges
 */
// convolve w/o FFT T-026 (num_data = odd)
TEST_F(Convolve1DOperation , TwoSpikesWithOutFFTDataIsOdd) {
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (odd) edge spikes data",
			{ kernel_width, num_data }, false, num_data, SpikeType_kleftright,
			MaskType_ktrue, { }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);
}
// convolve w/o FFT T-027 (num_data = even)
TEST_F(Convolve1DOperation , TwoSpikesWithOutFFTDataIsEven) {
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_EVEN;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (even) edge spikes data",
			{ kernel_width, num_data }, false, num_data, SpikeType_kleftright,
			MaskType_ktrue, { { 0, { 0.31633885 } }, { num_data - 1, {
					0.31633885 } } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);
}

/*
 * Convolution by Gaussian kernel. Data has 3 spikes at center and edges
 */
// convolve w/o FFT T-029 (num_data = even)
TEST_F(Convolve1DOperation , ThreeSpikesWithoutFFTDataIsEven) {
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_EVEN;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (even) 3 spikes data", { kernel_width, num_data },
			false, num_data, SpikeType_kall, MaskType_ktrue, { { 0, {
					0.316338854 } }, { num_data / 2, { 0.187887762 } }, {
					num_data - 1, { 0.316339299 } } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest);

}

/*
 * Test Performance of Convolve1D
 */
// convolve w/o FFT T-036 (num_data = even)
TEST_F(Convolve1DOperation , PerformanceTestWithoutFFT) {
	initializer_list<float> gaussian_kernel = initializer_list<float>( {
			0.011742966, 0.031861115, 0.069249175, 0.12056981, 0.16816399,
			0.18788746, 0.16816399, 0.12056981, 0.069249175, 0.031861115,
			0.011742966 }); // Kernel values
	float const kernel_width = NUM_WIDTH;
	size_t const num_kernel = 25;
	size_t const num_data = NUM_IN_MAX;
	size_t num_repeat = 10000;
	ConvolveTestComponent<GaussianKernel> GaussianKernelTest[] = { {
			"Gaussian kernel (even)", { kernel_width, num_kernel }, false,
			num_data, SpikeType_kcenter, MaskType_ktrue, { { num_data / 2
					- gaussian_kernel.size() / 2, gaussian_kernel } }, { } }, };
	RunConvolveTestComponentList<GaussianKernel>(ELEMENTSOF(GaussianKernelTest),
			GaussianKernelTest, num_repeat);
}

/*
 * Test user defined asymmetric kernel (right angled triangle)
 */
// convolve w/o FFT T-037-042
TEST_F(Convolve1DOperation, NumKernelWithoutFFTTriangle) {
	ReferenceData<float> ref_data = { 10, { 0.1, 0.2, 0.3, 0.4 } };
	ConvolveTestComponent<RightAngledTriangleKernel> TriangleNumKernelTest[] = {
			// T-037
			{ "num_kernel(odd) = num_data", { 4, NUM_IN_ODD }, false,
			NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { ref_data }, { } },
			// T-038
			{ "num_kernel(even) = num_data", { 4, NUM_IN_EVEN }, false,
			NUM_IN_EVEN, SpikeType_kcenter, MaskType_ktrue, { ref_data }, { } },
			// T-039
			{ "num_kernel(odd) > num_data", { 4, 2 * NUM_IN_ODD + 1 },
			false, NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { ref_data }, { } },
			// T-040
			{ "num_kernel(even) > num_data", { 4, 2 * NUM_IN_ODD }, false,
			NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { ref_data }, { } },
			// T-041
			{ "num_kernel(odd) < num_data", { 4, 13 }, false, NUM_IN_ODD,
					SpikeType_kcenter, MaskType_ktrue, { ref_data }, { } },
			// T-042
			{ "num_kernel(even) < num_data", { 4, 12 }, false, NUM_IN_ODD,
					SpikeType_kcenter, MaskType_ktrue, { ref_data }, { } } };
	RunConvolveTestComponentList<RightAngledTriangleKernel>(
			ELEMENTSOF(TriangleNumKernelTest), TriangleNumKernelTest);
}

/**
 * Test for sakura_CreateGaussianKernel
 */
// T-001
TEST_GAUSS(NegativeKernelWidth) {
	RunGaussianTest<StandardArrayInitializer, InvalidArgumentValidator>(-1.0f,
			10);
}

// T-002
TEST_GAUSS(ZeroKernelWidth) {
	RunGaussianTest<StandardArrayInitializer, InvalidArgumentValidator>(0.0f,
			10);
}

// T-003
TEST_GAUSS(PositiveInfiniteKernelWidth) {
	constexpr float kOne = 1.0f;
	constexpr float kZero = 0.0f;
	float positive_inf = kOne / kZero;
	ASSERT_TRUE(std::isinf(positive_inf));
	ASSERT_GT(positive_inf, FLT_MAX);
	RunGaussianTest<StandardArrayInitializer, InvalidArgumentValidator>(
			positive_inf, 10);
}

// T-004
TEST_GAUSS(NegativeInfiniteKernelWidth) {
	constexpr float kOne = 1.0f;
	constexpr float kZero = 0.0f;
	float negative_inf = -kOne / kZero;
	ASSERT_TRUE(std::isinf(negative_inf));
	ASSERT_LT(negative_inf, FLT_MIN);
	RunGaussianTest<StandardArrayInitializer, InvalidArgumentValidator>(
			negative_inf, 10);
}

// T-005
TEST_GAUSS(NotANumberKernelWidth) {
	float nan_value = std::nan("");
	ASSERT_TRUE(std::isnan(nan_value));
	RunGaussianTest<StandardArrayInitializer, InvalidArgumentValidator>(
			nan_value, 10);
}

// T-006
TEST_GAUSS(KernelIsNull) {
	RunGaussianTest<NullPointerArrayInitializer, InvalidArgumentValidator>(2.0f,
			10);
}

// T-007
TEST_GAUSS(KernelIsNotAligned) {
	RunGaussianTest<NotAlignedArrayInitializer, InvalidArgumentValidator>(2.0f,
			10);
}

// T-008
TEST_GAUSS(NumKernelIsZero) {
	RunGaussianTest<StandardArrayInitializer, InvalidArgumentValidator>(2.0f,
			0);
}

// T-009
TEST_GAUSS(NumKernelIsOne) {
// result should be independent of kernel_width if num_kernel is 1
	RunGaussianTest<StandardArrayInitializer, NumKernelOneValidator>(2.0f, 1);

	RunGaussianTest<StandardArrayInitializer, NumKernelOneValidator>(256.0f, 1);
}

// T-010
TEST_GAUSS(NumKernelIsTwo) {
	RunGaussianTest<StandardArrayInitializer, StandardValidator>(2.0f, 2);
}

// T-011
TEST_GAUSS(NumKernelIsThree) {
	RunGaussianTest<StandardArrayInitializer, StandardValidator>(2.0f, 3);
}

// T-012
TEST_GAUSS(NumKernelIsFour) {
	RunGaussianTest<StandardArrayInitializer, StandardValidator>(2.0f, 4);
}

// T-014
TEST_GAUSS(EvenNumkernel) {
	RunGaussianTest<StandardArrayInitializer, StandardValidator>(24.0f, 64);
}

// T-013
TEST_GAUSS(OddNumkernel) {
	RunGaussianTest<StandardArrayInitializer, StandardValidator>(24.0f, 65);
}

// T-016
TEST_GAUSS(NumkernelLessThanKernelWidth) {
// even
	RunGaussianTest<StandardArrayInitializer, StandardValidator>(128.0f, 64);

	// odd
	RunGaussianTest<StandardArrayInitializer, StandardValidator>(128.0f, 65);
}

// T-017
TEST_GAUSS(NumKernelGreaterThanKernelWidth) {
// even
	RunGaussianTest<StandardArrayInitializer, StandardValidator>(64.0f, 128);

	// odd
	RunGaussianTest<StandardArrayInitializer, StandardValidator>(64.0f, 129);
}

// T-018
TEST_GAUSS(HalfWidth) {
// even
	RunGaussianTest<StandardArrayInitializer, HalfWidthValidator>(128.0f, 128);

	// odd
	RunGaussianTest<StandardArrayInitializer, HalfWidthValidator>(129.0f, 129);
}

// T-019
TEST_GAUSS(NarrowKernelWidth) {
	RunGaussianTest<StandardArrayInitializer, NarrowKernelValidator>(FLT_MIN,
			128);
}

// T-020
TEST_GAUSS(WideKernelWidth) {
	RunGaussianTest<StandardArrayInitializer, WideKernelValidator>(FLT_MAX,
			128);
}

// T-021
TEST_GAUSS(PowerOfTwoNumKernelPerformance) {
	constexpr size_t kNumKernel = 268435456; // 2^28
	RunGaussianTest<StandardArrayInitializer, StandardValidator,
			PerformanceTestLogger>(128.0f, kNumKernel,
			"CreateGaussian_PowerOfTwoNumKernelPerformanceTest");
}

// T-022
TEST_GAUSS(OddNumKernelPerformance) {
	constexpr size_t kNumKernel = 268435456 + 1; // 2^28 + 1
	RunGaussianTest<StandardArrayInitializer, StandardValidator,
			PerformanceTestLogger>(128.0f, kNumKernel,
			"CreateGaussian_OddNumKernelPerformanceTest");
}

/**
 * Test for LIBSAKURA_SYMBOL(Convolve1DFloat)
 */
