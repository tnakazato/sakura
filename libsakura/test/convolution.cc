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
#include "testutil.h"

/* the number of elements in input/output array to test */
#define NUM_WIDTH 5
#define NUM_IN_EVEN 24
#define NUM_IN_ODD 25
#define NUM_IN_LARGE 8192

#define TEST_GAUSS(Name) TEST(CreateGaussianKernelTest, Name)

#define BENCH "#x# benchmark "

extern "C" {
struct LIBSAKURA_SYMBOL(Convolve1DContextFloat) {
	size_t num_kernel;
	fftw_plan plan_r2c;
	fftw_plan plan_c2r;
	fftw_complex *ffted_kernel;
	fftw_complex *ffted_input_data;
	fftw_complex *ffted_convolve_data;
	double *real_array;
	void *real_array_work;
};
}
typedef enum {
	SpikeType_kleft,
	SpikeType_kcenter,
	SpikeType_kright,
	SpikeType_kleftright,
	SpikeType_kall,
	SpikeType_knegative,
	SpikeType_kedgenaninf
} SpikeType;

typedef enum {
	MaskType_ktrue,
	MaskType_kfalse,
	MaskType_knonzerofalse,
	MaskType_kleftright,
	MaskType_knone
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
		float const peak_location = static_cast<float>(num_kernel / 2);
//		sakura_Status status = LIBSAKURA_SYMBOL(CreateGaussianKernelFloat)(
//				peak_location, kernel_width, num_kernel, kernel);
//		ASSERT_EQ(sakura_Status_kOK, status);
		double kernel_sum = 0.0;
		auto const sigma = kernel_width / (2.0 * sqrt(2.0 * log(2.0)));
		for (size_t i = 0; i < num_kernel; ++i) {
			auto const offset = static_cast<float>(i) - peak_location;
			kernel[i] = exp(-offset * offset / (2.0 * sigma * sigma));
			kernel_sum += kernel[i];
		}
		for (size_t i = 0; i < num_kernel; ++i) {
			kernel[i] /= kernel_sum;
		}
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
		size_t const start = std::max(
				static_cast<int>(num_kernel / 2 - kernel_width / 2), 0);
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
	template<typename DataType>
	static void Initialize(void *storage_ptr, DataType **kernel) {
		*kernel = reinterpret_cast<DataType *>(storage_ptr);
	}
};

struct NotAlignedArrayInitializer {
	template<typename DataType>
	static void Initialize(void *storage_ptr, DataType **kernel) {
		DataType *aligned_kernel = reinterpret_cast<DataType *>(storage_ptr);
		*kernel = aligned_kernel + 1;
	}
};

struct NullPointerArrayInitializer {
	template<typename DataType>
	static void Initialize(void *storage_ptr, DataType **kernel) {
		*kernel = nullptr;
	}
};

struct InPlaceArrayInitializer { // dummy initializer for in-place operation
	template<typename DataType>
	static void Initialize(void *storage_ptr, DataType **kernel) {
		// should not be called
		FAIL();
	}
};

struct InvalidArgumentValidator {
	static void Validate(float /*peak_location*/, float kernel_width,
			size_t num_kernel, float const *kernel,
			sakura_Status const status) {
		// only validate returned status
		EXPECT_EQ(sakura_Status_kInvalidArgument, status);
	}
};

struct NotGoodValidator {
	static void Validate(float /*peak_location*/, float /*kernel_width*/,
			size_t /*num_kernel*/, float const */*kernel*/,
			sakura_Status const status) {
		// only validate returned status
		EXPECT_EQ(sakura_Status_kNG, status);
	}
};

struct NoMemoryValidator {
	static void Validate(float /*peak_location*/, float /*kernel_width*/,
			size_t /*num_kernel*/, float const */*kernel*/,
			sakura_Status const status) {
		// only validate returned status
		EXPECT_EQ(sakura_Status_kNoMemory, status);
	}
};

struct StatusOKValidator {
	static void Validate(float /*peak_location*/, float /*kernel_width*/,
			size_t /*num_kernel*/, float const */*kernel*/,
			sakura_Status const status) {
		// execution must be successful
		EXPECT_EQ(sakura_Status_kOK, status);
	}
};

struct NumKernelOneValidator {
	static void Validate(float /*peak_location*/, float /*kernel_width*/,
			size_t num_kernel, float const *kernel,
			sakura_Status const status) {
		// num_kernel must be 1
		ASSERT_EQ(1u, num_kernel);

		// execution must be successful
		EXPECT_EQ(sakura_Status_kOK, status);

		EXPECT_EQ(1.0f, kernel[0]);
	}
};

struct BasicGaussKernelValidator {
	template<typename T>
	static T GetGaussianIntegration(T mu, T sigma, T a, T b) {
		auto an = (a - mu) / (sqrt(2.0) * sigma);
		auto bn = (b - mu) / (sqrt(2.0) * sigma);
		return 0.5 * (std::erf(bn) - std::erf(an));
	}

	static void Validate(float peak_location, float kernel_width,
			size_t num_kernel, float const *kernel,
			sakura_Status const status) {
		// execution must be successful
		EXPECT_EQ(sakura_Status_kOK, status);

		// validate each kernel elements
		std::unique_ptr<void, decltype(&Deallocate)> storage(
				malloc(num_kernel * sizeof(float)), Deallocate);
		float *expected_kernel = reinterpret_cast<float *>(storage.get());
		double sigma = kernel_width / (2.0 * sqrt(2.0 * log(2.0)));
		double expected_sum = 0.0;
		for (size_t i = 0; i < num_kernel; ++i) {
			double a = static_cast<double>(i) - 0.5;
			double b = static_cast<double>(i) + 0.5;
			expected_kernel[i] = GetGaussianIntegration<double>(peak_location,
					sigma, a, b);
			expected_sum += expected_kernel[i];
		}
		for (size_t i = 0; i < num_kernel; ++i) {
			float expected =
					(expected_sum == 0.0) ?
							0.0 : expected_kernel[i] / expected_sum;
			ASSERT_TRUE(std::isfinite(kernel[i]));
			EXPECT_FLOAT_EQ(expected, kernel[i]);
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

struct SymmetricKernelValidator {
	static void Validate(float peak_location, float kernel_width,
			size_t num_kernel, float const *kernel,
			sakura_Status const status) {
		// the additional symmetry check makes sense only if peak is in kernel.
		// otherwise, use BasicGaussKernelValidator
		ASSERT_LE(0.0, peak_location);
		ASSERT_GT(num_kernel, peak_location);
		// peak_location should be x.0 or x.5 for symmetry
		size_t peak_index = static_cast<size_t>(peak_location);
		ASSERT_TRUE(
				fabs(peak_location - peak_index) < 1.0e-6
						|| fabs(peak_location - peak_index - 0.5f) < 1.0e-6);

		BasicGaussKernelValidator::Validate(peak_location, kernel_width,
				num_kernel, kernel, status);
		// validate symmetry
		if (num_kernel > 1) {
//			size_t num_tail = peak_index - ((num_kernel % 2 == 0) ? 1 : 0);
			size_t offset = fabs(peak_location - peak_index) < 1.0e-6 ? 0 : 1;
			size_t num_tail = std::min(peak_index,
					num_kernel - 1 - peak_index - offset);
			for (size_t i = 0; i < num_tail; ++i) {
				// must match exactly so use EXPECT_EQ instead of EXPECT_FLOAT_EQ
				EXPECT_EQ(kernel[peak_index - i],
						kernel[peak_index + offset + i]);
			}
		}
	}
};

struct HalfWidthValidator {
	static void Validate(float peak_location, float kernel_width,
			size_t num_kernel, float const *kernel,
			sakura_Status const status) {
		// kernel_width must be equal to num_kernel and peak must be equal to num_kernel/2
		ASSERT_EQ(kernel_width, static_cast<float>(num_kernel));
		ASSERT_FLOAT_EQ(peak_location, static_cast<float>(num_kernel / 2));

		SymmetricKernelValidator::Validate(peak_location, kernel_width,
				num_kernel, kernel, status);

		// additional check: edge values should be half maximum
		size_t peak_index = num_kernel / 2;
		auto peak_value = kernel[peak_index];
		// correction factor for odd case:
		// in this case, FWHM is slightly shifted with respect to edge
		// so that correction must be needed for verification
		double half_width = kernel_width / 2.0f;
		double actual_width = static_cast<double>(peak_index);
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
	static void Validate(float peak_location, float kernel_width,
			size_t num_kernel, float const *kernel,
			sakura_Status const status) {
		// kernel_width must be equal to num_kernel
		ASSERT_EQ(FLT_MIN, kernel_width);

		SymmetricKernelValidator::Validate(peak_location, kernel_width,
				num_kernel, kernel, status);

		size_t peak_index = static_cast<size_t>(peak_location);
		std::cout << "peak: kernel[" << peak_index << "]=" << kernel[peak_index]
				<< std::endl;
		std::cout << "off-peak: kernel[" << peak_index - 1 << "]="
				<< kernel[peak_index - 1] << std::endl;
		std::cout << "off-peak: kernel[" << peak_index + 1 << "]="
				<< kernel[peak_index + 1] << std::endl;

		// kernel is so narrow that value of the kernel except peak_location
		// is zero
		for (size_t i = 0; i < peak_index; ++i) {
			EXPECT_FLOAT_EQ(0.0f, kernel[i]);
		}
		for (size_t i = peak_index + 1; i < num_kernel; ++i) {
			EXPECT_FLOAT_EQ(0.0f, kernel[i]);
		}

	}
};

struct WideKernelValidator {
	static void Validate(float peak_location, float kernel_width,
			size_t num_kernel, float const *kernel,
			sakura_Status const status) {
		// kernel_width must be equal to FLT_MAX and peak must be equal to num_kernel/2
		ASSERT_EQ(FLT_MAX, kernel_width);
		ASSERT_FLOAT_EQ(peak_location, static_cast<float>(num_kernel / 2));

		SymmetricKernelValidator::Validate(peak_location, kernel_width,
				num_kernel, kernel, status);

		size_t peak_index = num_kernel / 2;
		std::cout << "peak: kernel[" << peak_index << "]=" << kernel[peak_index]
				<< std::endl;
		std::cout << "edge: kernel[" << 0 << "]=" << kernel[0] << std::endl;
		std::cout << "edge: kernel[" << num_kernel - 1 << "]="
				<< kernel[num_kernel - 1] << std::endl;

		// kernel is so wide that all the kernel value are the same
		for (size_t i = 0; i < num_kernel; ++i) {
			EXPECT_FLOAT_EQ(kernel[peak_index], kernel[i]);
		}
	}
};

struct ContextValidator {
	static void Validate(size_t num_kernel, float const *kernel,
	LIBSAKURA_SYMBOL(Convolve1DContextFloat) const *context) {
		ASSERT_EQ(context->num_kernel, num_kernel);
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
inline void RunGaussianTest(float const peak_location, float const kernel_width,
		size_t const num_kernel, std::string const test_name = "") {
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
	Initializer::Initialize(storage_ptr, &kernel);

	// run
	double start_time = GetCurrentTime();
	sakura_Status status = sakura_CreateGaussianKernelFloat(peak_location,
			kernel_width, num_kernel, kernel);
	double end_time = GetCurrentTime();
	Logger::PrintElapsedTime(test_name, end_time - start_time);

	// verification
	Validator::Validate(peak_location, kernel_width, num_kernel, kernel,
			status);
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
	// set some values to kernelsetprecision(8) <<
	RightAngledTriangleKernel kernel_generator = {
			static_cast<float>(num_dummy), num_dummy };
	kernel_generator.Generate(dummy_kernel);
	// initialize kernel
	float *kernel = nullptr;
	KernelInitializer::Initialize(dummy_kernel, &kernel);
	// initialize context
	LIBSAKURA_SYMBOL(Convolve1DContextFloat) **context_ptr_ptr = nullptr;
	LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context_ptr = nullptr;
	if (valid_context_pointer) {
		context_ptr_ptr = &context_ptr;
	}
	// create context
	LIBSAKURA_SYMBOL(Status) create_status =
	LIBSAKURA_SYMBOL(CreateConvolve1DContextFFTFloat)(num_kernel, kernel,
			context_ptr_ptr);
	//std::cout << "create kernel status = " << create_status << std::endl;
	constexpr float kDummyKernelValue = 0.0f;
	CreateValidator::Validate(kDummyKernelValue, kDummyKernelValue, num_kernel,
			kernel, create_status);
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
		LIBSAKURA_SYMBOL(DestroyConvolve1DContextFloat)(*context_ptr_ptr);
		//std::cout << "destroy status = " << status << std::endl;
		DestroyValidator::Validate(kDummyKernelValue, kDummyKernelValue,
				num_kernel, kernel, status);
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
 *  reference: a vector of reference values
 * see CompareResultWithReference for how to use it in tests.
 */
template<typename DataType>
struct ReferenceData {
	size_t offset;
	vector<DataType> reference;

	ReferenceData(size_t offset_index = 0, initializer_list<DataType> ref = { }) :
			offset(offset_index), reference(ref) {
	}
	ReferenceData(size_t offset_index = 0,
			vector<DataType> ref = vector<DataType>()) :
			offset(offset_index), reference(ref) {
	}
};
// Compare result array with a vector of reference data.
inline void CompareResultWithReference(size_t num_result, float const *result,
		vector<ReferenceData<float>> const &reference_list) {
	for (size_t i = 0; i < reference_list.size(); ++i) {
		size_t const offset = reference_list[i].offset;
		size_t const num_reference = reference_list[i].reference.size();
		vector<float> const reference = reference_list[i].reference;
		ASSERT_TRUE(offset + num_reference <= num_result);
		for (size_t j = 0; j < num_reference; ++j) {
			EXPECT_FLOAT_EQ(reference[j], result[offset + j]);
		}
	}
}

/*
 * A class to test convolution with or without FFT by invoking either
 * sakura_Convolve1DFFTFloat or sakura_Convolve1DFloat
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

	template<typename DataType>
	void PrintArray(char const *name, size_t num_data, DataType *data_array) {
		constexpr size_t kMaxLength(30);
		cout << name << " = ";
		if (data_array == nullptr) { // array is nullptr
			cout << "NULL" << endl;
		} else if (num_data > kMaxLength) { // long array (just show the length)
			cout << num_data << " elements" << endl;
		} else { // normal array
			cout << "[ ";
			if (num_data > 0) {
				for (size_t i = 0; i < num_data - 1; ++i)
					cout << setprecision(9) << data_array[i] << ", ";
				cout << setprecision(9) << data_array[num_data - 1];
			}
			cout << " ]" << endl;
		}
	}

	/*
	 *  a struct to store test parameters and expected results (ReferenceData) to
	 *  test convolution. 
	 *  name: a test name to print
	 *  kernel: kernel parameters. kernel width and size of kernel array is stored
	 *          in a struct Kernel
	 *  use_fft: if true, convolution is done with FFT
	 *  num_data: the size of data (and mask) arrays
	 *  data_type: describes how to initialize data array by an enum, SpikeType
	 *  mask_type: describes how to initialize mask array by an enum, SpikeType
	 *  data_ref: a vector of expected results (ReferenceData) of convolution.
	 *            The values are compared with output_data in the test.
	 *  weight_ref: a vector of expected results (ReferenceData) of convolution.
	 *              The values are compared with output_weight in the test.
	 *  see Convolve1DOperation::RunConvolveTestComponentList for how to use it in tests.
	 */
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
				initializer_list<ReferenceData<float>> data = { },
				initializer_list<ReferenceData<float>> weight = { }) :
				name(in_name), use_fft(use_fft), num_data(in_num_data), data_type(
						in_dtype), mask_type(in_mtype), data_ref(data), weight_ref(
						weight) {
			kernel = in_kernel;
		}
	};
	/*
	 * Invoke Convolution using a list of ConvolveTestComponent.
	 * this function loops over the list, allocate arrays for each ConvolveTestComponent
	 * and invoke RunConvolutionTest function in each loop.
	 */
	template<typename ConvolveKernel, typename KernelInitializer,
			typename InDataInitializer, typename InMaskInitializer,
			typename OutDataInitializer, typename WeightInitializer,
			typename ConvolveValidator>
	void RunConvolveTestComponentList(size_t num_test,
			ConvolveTestComponent<ConvolveKernel> const TestList[],
			size_t const num_repeat = 1) {
		for (size_t i = 0; i < num_test; ++i) {
			ConvolveTestComponent<ConvolveKernel> const test_component =
					TestList[i];
//			if (num_test > 1)
			{
				cout << "[Test " << test_component.name << " w/"
						<< (test_component.use_fft ? "" : "o") << " FFT]"
						<< endl;
			}
			// force setting mask type to MaskType_knone to avoid initialization of mask
			MaskType const mask_type =
					test_component.use_fft ?
							MaskType_knone : test_component.mask_type;
			// avoid allocating storage > INT_MAX
			size_t const num_dummy =
					(test_component.num_data <= INT_MAX) ?
							test_component.num_data + 1 : 1;
			// don't need allocate input_mask and output_weight if use_fft = true
			size_t const num_maskweight =
					(test_component.use_fft || test_component.num_data > INT_MAX) ?
							1 : test_component.num_data;
			size_t const num_maskweight_dummy = num_maskweight + 1;
			// allocate storage for input/output data arrays
			void *dummy_input_data = nullptr;
			std::unique_ptr<void, DefaultAlignedMemory> in_storage(
					DefaultAlignedMemory::AlignedAllocateOrException(
							sizeof(float) * num_dummy, &dummy_input_data));
			ASSERT_NE(dummy_input_data, nullptr);
			void *dummy_output_data = nullptr;
			std::unique_ptr<void, DefaultAlignedMemory> out_storage(
					DefaultAlignedMemory::AlignedAllocateOrException(
							sizeof(float) * num_dummy, &dummy_output_data));
			ASSERT_NE(dummy_output_data, nullptr);
			// allocate storage for weight and mask
			void *dummy_output_weight = nullptr;
			std::unique_ptr<void, DefaultAlignedMemory> weight_storage(
					DefaultAlignedMemory::AlignedAllocateOrException(
							sizeof(float) * num_maskweight_dummy,
							&dummy_output_weight));
			ASSERT_NE(dummy_output_weight, nullptr);
			void * dummy_input_mask = nullptr;
			std::unique_ptr<void, DefaultAlignedMemory> mask_storage(
					DefaultAlignedMemory::AlignedAllocateOrException(
							sizeof(bool) * num_maskweight_dummy,
							&dummy_input_mask));
			ASSERT_NE(dummy_input_mask, nullptr);
			// initialize arrays
			float *input_data = nullptr;
			InDataInitializer::Initialize(dummy_input_data, &input_data);
			bool *input_mask = nullptr;
			InMaskInitializer::Initialize(dummy_input_mask, &input_mask);
			float *output_data = nullptr;
			if (typeid(OutDataInitializer) == typeid(InPlaceArrayInitializer)) {
				output_data = input_data;
			} else {
				OutDataInitializer::Initialize(dummy_output_data, &output_data);
			}
			float *output_weight = nullptr;
			WeightInitializer::Initialize(dummy_output_weight, &output_weight);
			// Invoke a test
			RunConvolutionTest<ConvolveKernel, KernelInitializer,
					ConvolveValidator>(test_component.kernel,
					test_component.use_fft, test_component.data_type, mask_type,
					test_component.num_data, input_data, input_mask,
					output_data, output_weight, num_repeat,
					test_component.data_ref, test_component.weight_ref);
		}
	}

	/*
	 * Run sakura_Convolve1D[FFT]Float and compare output data and weight
	 * with references.
	 * In case of FFT, create and destroy of context should always succeed.
	 */
	template<typename ConvolveKernel, typename KernelInitializer,
			typename ConvolveValidator>
	void RunConvolutionTest(ConvolveKernel kernel_param, bool use_fft,
			SpikeType const spike_type, MaskType const mask_type,
			size_t const num_data, float *input_data, bool *input_mask,
			float *output_data, float *output_weight, size_t const num_repeat =
					1, vector<ReferenceData<float>> const data_ref = { },
			vector<ReferenceData<float>> const weight_ref = { }) {
		// create kernel array
		size_t const num_kernel = kernel_param.num_kernel; // nominal kernel array size
		size_t const num_dummy = (num_kernel <= INT_MAX) ? num_kernel + 1 : 1;
		float *dummy_kernel = nullptr;
		std::unique_ptr<void, DefaultAlignedMemory> kernel_storage(
				DefaultAlignedMemory::AlignedAllocateOrException(
						sizeof(float) * num_dummy, &dummy_kernel));
		ASSERT_NE(dummy_kernel, nullptr);
		// modify num_kernel in ConvolveKernel if necessary
		if (num_kernel < 1 || num_kernel > INT_MAX) {
			kernel_param.num_kernel = num_dummy;
		}
		kernel_param.Generate(dummy_kernel);
		float *kernel = nullptr;
		KernelInitializer::Initialize(dummy_kernel, &kernel);
		size_t const num_kernel_actural = kernel_param.num_kernel; // actual kernel array size
		// create context if FFT
		LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context = nullptr;
		bool in_place_operation = false;
		if (use_fft) {
			LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL(CreateConvolve1DContextFFTFloat)(
					num_kernel_actural, kernel, &context);
			ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		}
		if (input_data == output_data && input_data != nullptr) {
			in_place_operation = true;
			cout << "Invoking in-place convolution." << endl;
		}
		// initialize data and mask (only if the size of data and mask
		// arrays is actually num_data)
		if (input_data
				!= nullptr&& input_mask != nullptr && num_data <= INT_MAX) {
			if (use_fft) { // set only data
				InitializeData(spike_type, num_data, input_data);
			} else { // set data and mask
				InitializeDataAndMask(spike_type, mask_type, num_data,
						input_data, input_mask);
			}
		}
		auto reinitialize_function =
				in_place_operation ?
						InPlaceAction::reinitialize :
						OutOfPlaceAction::reinitialize;
		if (verbose) {
			PrintArray("kernel", num_kernel_actural, kernel);
			PrintArray("data (in)", num_data, input_data);
			if (!use_fft)
				PrintArray("mask (in)", num_data, input_mask);
		}
		// convolution
		if (num_repeat > 1) {
			cout << "Iterating convolution for " << num_repeat
					<< " loops. The length of arrays is " << num_data << endl;
		}
		LIBSAKURA_SYMBOL(Status) exec_status = LIBSAKURA_SYMBOL(Status_kOK);
		double start = GetCurrentTime();
		if (use_fft) {
			for (size_t i = 0; i < num_repeat; ++i) {
				reinitialize_function(spike_type, num_data, input_data);
				exec_status = LIBSAKURA_SYMBOL(Convolve1DFFTFloat)(context,
						num_data, input_data, output_data);
			}
		} else {
			for (size_t i = 0; i < num_repeat; ++i) {
				exec_status = LIBSAKURA_SYMBOL(Convolve1DFloat)(num_kernel,
						kernel, num_data, input_data, input_mask, output_data,
						output_weight);
			}
		}
		double end = GetCurrentTime();
		if (num_repeat > 1) {
			cout << "#x# benchmark Convolve1DFloat_With"
					<< (use_fft ? "" : "Out") << "FFT_Data"
					<< (num_data % 2 == 0 ? "Even" : "Odd") << " "
					<< end - start << endl;
		}
		ConvolveValidator::Validate(kernel_param.kernel_width / 2.,
				kernel_param.kernel_width, num_kernel_actural, kernel,
				exec_status);
		// destroy context
		if (use_fft) {
			LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL(DestroyConvolve1DContextFloat)(context);
			ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		}
		// verify output data and weight
		if (exec_status == LIBSAKURA_SYMBOL(Status_kOK)) {
			if (verbose) {
				PrintArray("data (out)", num_data, output_data);
				if (!use_fft)
					PrintArray("weight (out)", num_data, output_weight);
			}
			CompareResultWithReference(num_data, output_data, data_ref);
			CompareResultWithReference(num_data, output_weight, weight_ref);
		}
	}

	bool verbose;

	// representative references
	// Gaussian (kernel_width=5) 
	vector<float> const gaussian_ref_fft { 0.031861115, 0.069249175, 0.12056981,
			0.16816399, 0.18788746, 0.16816399, 0.12056981, 0.069249175,
			0.031861115 }; // Kernel values
	vector<float> const gauss_weight_ref_L { 0.59394373002975898,
			0.76210771502975894, 0.88267752502975894, 0.95192670052975892,
			0.98378781572975893 };
	vector<float> const gauss_weight_ref_R_odd { 0.98378781572975893,
			0.95192670052975892, 0.88267752502975894, 0.76210771502975894,
			0.59394373002975898 }; // analytic value
	vector<float> const gauss_weight_ref_R_even { 0.98378779394430194,
			0.95192667874430192, 0.88267750324430194, 0.76210769324430194,
			0.59394370824430198 };
	// right-angled triangle (kernel_width = 4)
	vector<float> const triangle_ref { 0.1, 0.2, 0.3, 0.4 }; // Kernel values
	vector<float> const triangle_weight_ref_L { 0.6, 1.0 };
	vector<float> const triangle_weight_ref_R { 1.0, 0.9, 0.7 };

private:
	// set data and mask arrays based on SpikeType and MaskType
	void InitializeDataAndMask(SpikeType const spike_type,
			MaskType const mask_type, size_t const num_data, float data[],
			bool mask[]) {
		// data
		InitializeData(spike_type, num_data, data);
		// mask
		if (mask_type != MaskType_knone) {
			for (size_t i = 0; i < num_data; ++i) {
				mask[i] = (mask_type != MaskType_kfalse);
			}
		}
		switch (mask_type) {
		case MaskType_ktrue:
		case MaskType_kfalse:
		case MaskType_knone:
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
	// set data array based on SpikeType
	static void InitializeData(SpikeType const spike_type,
			size_t const num_data, float data[]) {
		for (size_t i = 0; i < num_data; ++i) {
			data[i] = 0.0;
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
		case SpikeType_kedgenaninf:
			data[0] = (0.0f / 0.0f); // NaN (not a number) in float
			data[num_data / 2] = 1.0;
			data[num_data - 1] = (1.0f / 0.0f); // positive infinity in float
			break;
		}
	}
	// Re-initialization of input data and mask for the case of in-place operation.
	struct InPlaceAction {
		static void reinitialize(SpikeType const spike_type,
				size_t const num_data, float *data) {
			InitializeData(spike_type, num_data, data);
		}
	};
	// Dummy function (does nothing) for re-initialization for the case of out-of-place operation.
	struct OutOfPlaceAction {
		static void reinitialize(SpikeType const spike_type,
				size_t const num_data, float data[]) {
			// do nothing
		}
	};
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
/*
 * Failure cases of Convolution
 */
// convolve w/ FFT T-001
TEST_F(Convolve1DOperation , WithFFTContextIsNull) {
	LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context = nullptr;
	size_t const num_data = 1;
	SIMD_ALIGN float input_data[num_data];
	SIMD_ALIGN float output_data[ELEMENTSOF(input_data)];
	LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Convolve1DFFTFloat)(
			context, num_data, input_data, output_data);
	ASSERT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}
// convolve w/ FFT T-002 ~ 004
TEST_F(Convolve1DOperation , WithFFTNumDataIsInvalid) {
	ConvolveTestComponent<GaussianKernel> TestList[] = {
	// T-002 (num_data = 0)
			{ "num_data = 0", { NUM_WIDTH, 1 }, true, 0, SpikeType_kcenter,
					MaskType_knone },
			// T-003 (num_data = INT_MAX)
			{ "num_data > INT_MAX", { NUM_WIDTH, 1 }, true,
					(static_cast<size_t>(INT_MAX) + 1), SpikeType_kcenter,
					MaskType_knone },
			// T-004 (num_data != num_kernel)
			{ "num_data != num_kernel", { NUM_WIDTH, NUM_IN_EVEN }, true,
			NUM_IN_ODD, SpikeType_kcenter, MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/ FFT T-005
TEST_F(Convolve1DOperation , WithFFTInDataIsNull) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { { "in_data = nullptr",
			{ NUM_WIDTH, num_data }, true, num_data, SpikeType_kcenter,
			MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			NullPointerArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/ FFT T-006
TEST_F(Convolve1DOperation , WithFFTInDataIsNotAligned) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { {
			"not aligned in_data", { NUM_WIDTH, num_data }, true, num_data,
			SpikeType_kcenter, MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			NotAlignedArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/ FFT T-007
TEST_F(Convolve1DOperation , WithFFTOutDataIsNull) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { { "out_data = nullptr",
			{ NUM_WIDTH, num_data }, true, num_data, SpikeType_kcenter,
			MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			NullPointerArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/ FFT T-008
TEST_F(Convolve1DOperation , WithFFTOutDataIsNotAligned) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { {
			"Not aligned out_data", { NUM_WIDTH, num_data }, true, num_data,
			SpikeType_kcenter, MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			NotAlignedArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
/*
 * In-Place convolution
 */
// convolve w/ FFT T-009
TEST_F(Convolve1DOperation , WithFFTInPlaceConvolution) {
	ConvolveTestComponent<GaussianKernel> TestList[] =
			{ { "In-place convolution (&output_data == &input_data)", {
			NUM_WIDTH, NUM_IN_ODD }, true, NUM_IN_ODD, SpikeType_kcenter,
					MaskType_knone, { {
					NUM_IN_ODD / 2 - gaussian_ref_fft.size() / 2,
							gaussian_ref_fft } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			InPlaceArrayInitializer, StandardArrayInitializer, StatusOKValidator>(
			ELEMENTSOF(TestList), TestList);
}
/*
 * Convolution by narrow Gaussian kernel (kernel_width < num_data)
 */
// convolve w/ FFT T-010, 012 ~ 013
TEST_F(Convolve1DOperation , GaussWithFFTVariousNumData) {
	ConvolveTestComponent<GaussianKernel> TestList[] =
			{
			// T-010 (num_data = 1: LB)
					{ "Gaussian kernel (num_data = 1: LB)", { NUM_WIDTH, 1 },
					true, 1, SpikeType_kcenter, MaskType_knone,
							{ { 0, { 1.0 } } } },
					// T-012 (num_data=odd)
					{ "Gaussian kernel (num_data = odd)", { NUM_WIDTH,
					NUM_IN_ODD },
					true, NUM_IN_ODD, SpikeType_kcenter, MaskType_knone, { {
					NUM_IN_ODD / 2 - gaussian_ref_fft.size() / 2,
							gaussian_ref_fft } } },
					// T-013 (num_data=even)
					{ "Gaussian kernel (num_data = even)", { NUM_WIDTH,
					NUM_IN_EVEN },
					true, NUM_IN_EVEN, SpikeType_kcenter, MaskType_knone, { {
					NUM_IN_EVEN / 2 - gaussian_ref_fft.size() / 2,
							gaussian_ref_fft } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList);
}
/*
 * Convolution by wide Gaussian kernel (kernel_width > num_data)
 */
// convolve w/ FFT T-014 ~ 015
TEST_F(Convolve1DOperation , WideGaussWithFFT) {
	size_t const num_data_even = NUM_IN_EVEN;
	size_t const num_data_odd = NUM_IN_ODD;
	vector<float> gaussian_odd { 0.043915681541, 0.0455670394003,
			0.0468942411244, 0.047865845263, 0.0484584420919, 0.0486575998366,
			0.0484584420919, 0.047865845263, 0.0468942411244, 0.0455670394003,
			0.043915681541 }; // analytic value by emulating the code
	vector<float> gaussian_even { 0.0453697890043, 0.0472178347409,
			0.0487070940435, 0.0497995242476, 0.0504667051136, 0.0506910830736,
			0.0504667051136, 0.0497995242476, 0.0487070940435, 0.0472178347409,
			0.0453697890043 }; // analytic value by emulating the code
	ReferenceData<float> ref_odd = { num_data_odd / 2 - gaussian_odd.size() / 2,
			gaussian_odd };
	ReferenceData<float> ref_even = { num_data_even / 2
			- gaussian_even.size() / 2, gaussian_even };
	ConvolveTestComponent<GaussianKernel> TestList[] = {
	// T-014 (num_data=odd)
			{ "wide Gaussian kernel (num_data = odd)", { num_data_odd + 1,
					num_data_odd }, true, num_data_odd, SpikeType_kcenter,
					MaskType_knone, { { ref_odd } } },
			// T-015 (num_data=even)
			{ "wide Gaussian kernel (num_data = even)", { num_data_even + 1,
					num_data_even }, true, num_data_even, SpikeType_kcenter,
					MaskType_knone, { { ref_even } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList);
}
/*
 * Convolution by Gaussian kernel. Data has 2 spikes at edges
 */
// convolve w/ FFT T-016 ~ 017
TEST_F(Convolve1DOperation , GaussWithFFTTwoSpikes) {
	vector<float> gaussian_ref_ledge { 0.356051445, 0.288733795, 0.189818985,
			0.101110291, 0.0436040814 };
	vector<float> gaussian_ref_redge { 0.0436040814, 0.101110291, 0.189818985,
			0.288733795, 0.356051445 };
	ConvolveTestComponent<GaussianKernel> TestList[] =
			{
			// T-016 (num_data=odd)
					{ "Gaussian kernel (num_data = odd) edge spikes in data", {
					NUM_WIDTH, NUM_IN_ODD }, true, NUM_IN_ODD,
							SpikeType_kleftright, MaskType_knone, { { 0,
									gaussian_ref_ledge }, { NUM_IN_ODD
									- gaussian_ref_redge.size(),
									gaussian_ref_redge } } },
					// T-017 (num_data=even)
					{ "Gaussian kernel (num_data = even) edge spikes in data", {
					NUM_WIDTH, NUM_IN_EVEN }, true, NUM_IN_EVEN,
							SpikeType_kleftright, MaskType_knone, { { 0,
									gaussian_ref_ledge }, { NUM_IN_EVEN
									- gaussian_ref_redge.size(),
									gaussian_ref_redge } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList);
}
/*
 * Convolution by Gaussian kernel. Data has 3 spikes at center and edges
 */
// convolve w/ FFT T-018 ~ 019
TEST_F(Convolve1DOperation , GaussWithFFTThreeSpikes) {
	vector<float> gaussian_ref_L { 0.35605147, 0.28873407, 0.18982185,
			0.10113387 };
	vector<float> gaussian_ref_C_odd { 0.12057296, 0.16816429, 0.18788746,
			0.16816429, 0.12057296 };
	vector<float> gaussian_ref_R_odd { 0.10113387, 0.18982185, 0.28873407,
			0.35605147 };
	vector<float> gaussian_ref_C_even { 0.12057296, 0.16816429, 0.18788776,
			0.16816713, 0.12059626 };
	vector<float> gaussian_ref_R_even { 0.10126565, 0.18984257, 0.28873666,
			0.35605172 };
	ConvolveTestComponent<GaussianKernel> TestList[] = {
	// T-018  (num_data = odd)
			{ "Gaussian kernel (num_data = odd) 3 spikes in data", {
			NUM_WIDTH, NUM_IN_ODD }, true, NUM_IN_ODD, SpikeType_kall,
					MaskType_knone, { { 0, gaussian_ref_L },
							{ NUM_IN_ODD / 2 - gaussian_ref_C_odd.size() / 2,
									gaussian_ref_C_odd }, {
							NUM_IN_ODD - gaussian_ref_R_odd.size(),
									gaussian_ref_R_odd } } },
			// T-019  (num_data = even)
			{ "Gaussian kernel (num_data = even) 3 spikes in data", {
			NUM_WIDTH, NUM_IN_EVEN }, true, NUM_IN_EVEN, SpikeType_kall,
					MaskType_knone, { { 0, gaussian_ref_L }, { NUM_IN_EVEN / 2
							- gaussian_ref_C_even.size() / 2,
							gaussian_ref_C_even }, {
					NUM_IN_EVEN - gaussian_ref_R_even.size(),
							gaussian_ref_R_even } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList);
}
/*
 * Convolution by Gaussian kernel. Data has negative spikes at center
 */
// convolve w/ FFT T-020 ~ 021
TEST_F(Convolve1DOperation , GaussWithFFTNegativeSpike) {
	vector<float> data_ref { -0.031861115, -0.069249175, -0.12056981,
			-0.16816399, -0.18788746, -0.16816399, -0.12056981, -0.069249175,
			-0.031861115 };
	ConvolveTestComponent<GaussianKernel> TestList[] = {
	// T-020 (num_data=odd)
			{ "Gaussian kernel (num_data = odd)", { NUM_WIDTH, NUM_IN_ODD },
			true, NUM_IN_ODD, SpikeType_knegative, MaskType_knone, { {
			NUM_IN_ODD / 2 - data_ref.size() / 2, data_ref } } },
			// T-021 (num_data=even)
			{ "Gaussian kernel (num_data = even)", { NUM_WIDTH, NUM_IN_EVEN },
			true, NUM_IN_EVEN, SpikeType_knegative, MaskType_knone, { {
			NUM_IN_EVEN / 2 - data_ref.size() / 2, data_ref } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList);
}
/*
 * Test Performance of Convolve1D
 */
// convolve w/ FFT T-022 (num_data = odd)
TEST_F(Convolve1DOperation , PerformanceTestWithFFTDataIsOdd) {
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_LARGE - 1;
	size_t num_repeat = 10000;
	ConvolveTestComponent<GaussianKernel> TestList[] = { {
			"Gaussian kernel (num_data = odd)", { kernel_width, num_data },
			true, num_data, SpikeType_kcenter, MaskType_knone, { { num_data / 2
					- gaussian_ref_fft.size() / 2, gaussian_ref_fft } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList, num_repeat);

}
// convolve w/ FFT T-023 (num_data = even)
TEST_F(Convolve1DOperation , PerformanceTestWithFFTDataIsEven) {
	float const kernel_width = NUM_WIDTH;
	size_t const num_data = NUM_IN_LARGE;
	size_t num_repeat = 10000;
	ConvolveTestComponent<GaussianKernel> TestList[] = { {
			"Gaussian kernel (num_data = even)", { kernel_width, num_data },
			true, num_data, SpikeType_kcenter, MaskType_knone, { { num_data / 2
					- gaussian_ref_fft.size() / 2, gaussian_ref_fft } } }, };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList, num_repeat);
}
/*
 * Convolution by user defined asymmetric kernel (right angled triangle)
 */
// convolve w/ FFT T-024 ~ 023
TEST_F(Convolve1DOperation , TriangleWithFFT) {
	float const kernel_width = 4;
	size_t const num_data_odd = NUM_IN_ODD;
	size_t const num_data_even = NUM_IN_EVEN;
	// offset of reference data
	size_t const offset_odd = num_data_odd / 2 - triangle_ref.size() / 2;
	size_t const offset_even = num_data_even / 2 - triangle_ref.size() / 2;
	// reference data
	vector<float> triangle_ref_L { 0.7, 0.4 };
	vector<float> triangle_ref_R { 0.1, 0.3, 0.5 };
	vector<float> wide_ref_odd { 0.015384615, 0.01846154, 0.02153846,
			0.02461538, 0.02769231, 0.03076923, 0.03384615, 0.03692308, 0.04,
			0.04307692, 0.04615385, 0.04923077, 0.05230769, 0.05538461,
			0.05846154, 0.06153846, 0.06461538, 0.06769231, 0.07076923,
			0.07384615, 0.07692308 };
	vector<float> wide_ref_even { 0.01, 0.01333333, 0.016666667, 0.02,
			0.02333333, 0.026666667, 0.03, 0.03333334, 0.036666667, 0.04,
			0.04333333, 0.046666667, 0.05, 0.05333333, 0.056666667, 0.06,
			0.06333333, 0.066666667, 0.07, 0.07333333, 0.076666667, 0.08 };
	ConvolveTestComponent<RightAngledTriangleKernel> TestList[] = {
	// T-024 (num_data = odd)
			{ "Asynmetric kernel (center spike, num_data = odd)", {
					kernel_width, num_data_odd },
			true, num_data_odd, SpikeType_kcenter, MaskType_knone, { {
					offset_odd, triangle_ref } } },
			// T-025 (num_data = even)
			{ "Asynmetric kernel (center spike, num_data = even)", {
					kernel_width, num_data_even },
			true, num_data_even, SpikeType_kcenter, MaskType_knone, { {
					offset_even, triangle_ref } } },
			// T-026 (num_data = odd)
			{ "Asynmetric kernel (kernel_width > num_data, num_data = odd)", {
					num_data_odd + 1, num_data_odd }, true, num_data_odd,
					SpikeType_kcenter, MaskType_knone, { { num_data_odd
							- wide_ref_odd.size(), wide_ref_odd } } },
			// T-027 (num_data = even)
			{ "Asynmetric kernel (kernel_width > num_data, num_data = even)", {
					num_data_odd + 1, num_data_even },
			true, num_data_even, SpikeType_kcenter, MaskType_knone, { {
					num_data_even - wide_ref_even.size(), wide_ref_even } } },
			// T-028 (num_data = odd)
			{ "Asynmetric kernel (edge spikes, num_data = odd)", { kernel_width,
					num_data_odd },
			true, num_data_odd, SpikeType_kleftright, MaskType_knone, { { 0,
					triangle_ref_L }, { num_data_odd - triangle_ref_R.size(),
					triangle_ref_R } } },
			// T-029 (num_data = even)
			{ "Asynmetric kernel (edge spikes, num_data = even)", {
					kernel_width, num_data_even },
			true, num_data_even, SpikeType_kleftright, MaskType_knone, { { 0,
					triangle_ref_L }, { num_data_even - triangle_ref_R.size(),
					triangle_ref_R } } },
			// T-030 (num_data = odd)
			{ "Asynmetric kernel (3 spikes, num_data = odd)", { kernel_width,
					num_data_odd },
			true, num_data_odd, SpikeType_kall, MaskType_knone, { { 0,
					triangle_ref_L }, { offset_odd, triangle_ref }, {
					num_data_odd - triangle_ref_R.size(), triangle_ref_R } } },
			// T-031 (num_data = even)
			{ "Asynmetric kernel (3 spikes, num_data = even)", { kernel_width,
					num_data_even },
			true, num_data_even, SpikeType_kall, MaskType_knone, { { 0,
					triangle_ref_L }, { offset_even, triangle_ref }, {
					num_data_even - triangle_ref_R.size(), triangle_ref_R } } },
			// T-032 (num_data = odd)
			{ "Asynmetric kernel (negative spike, num_data = odd)", {
					kernel_width, num_data_odd },
			true, num_data_odd, SpikeType_knegative, MaskType_knone, { {
					offset_odd, { -0.1, -0.2, -0.3, -0.4 } } } },
			// T-033 (num_data = even)
			{ "Asynmetric kernel (negative spike, num_data = even)", {
					kernel_width, num_data_even },
			true, num_data_even, SpikeType_knegative, MaskType_knone, { {
					offset_even, { -0.1, -0.2, -0.3, -0.4 } } } } };
	RunConvolveTestComponentList<RightAngledTriangleKernel,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StatusOKValidator>(ELEMENTSOF(TestList),
			TestList);

}

/**********************************************************************
 * Test Convolution operation withOUT FFT
 *********************************************************************/
/*
 * Failure cases of Convolution
 */
// convolve w/o FFT T-001 ~ 004
TEST_F(Convolve1DOperation , WithOutFFTNumDataIsInvalid) {
	ConvolveTestComponent<GaussianKernel> TestList[] = {
	// T-001 (num_kernel = 0)
			{ "num_kernel = 0", { NUM_WIDTH, 0 }, false, NUM_IN_ODD,
					SpikeType_kcenter, MaskType_knone },
			// T-002 (num_kernel > INT_MAX)
			{ "num_kernel > INT_MAX", { NUM_WIDTH, static_cast<size_t>(INT_MAX)
					+ 1 }, false,
			NUM_IN_ODD, SpikeType_kcenter, MaskType_knone },
			// T-003 (num_data = 0)
			{ "num_data = 0", { NUM_WIDTH, 1 }, false, 0, SpikeType_kcenter,
					MaskType_knone },
			// T-004 (num_data > INT_MAX)
			{ "num_data > INT_MAX", { NUM_WIDTH, 1 }, false,
					(static_cast<size_t>(INT_MAX) + 1), SpikeType_kcenter,
					MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/o FFT T-005
TEST_F(Convolve1DOperation , WithOutFFTKernelIsNull) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] =
			{ { "kernel = nullptr", {
			NUM_WIDTH, num_data }, false, num_data, SpikeType_kcenter,
					MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, NullPointerArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/o FFT T-006
TEST_F(Convolve1DOperation , WithOutFFTKernelIsNotAligned) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { { "not aligned kernel",
			{ NUM_WIDTH, num_data }, false, num_data, SpikeType_kcenter,
			MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, NotAlignedArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/o FFT T-007
TEST_F(Convolve1DOperation , WithOutFFTInDataIsNull) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { {
			"input_data = nullptr", { NUM_WIDTH, num_data }, false, num_data,
			SpikeType_kcenter, MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			NullPointerArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/o FFT T-008
TEST_F(Convolve1DOperation , WithOutFFTInDataIsNotAligned) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { {
			"not aligned input_data", { NUM_WIDTH, num_data }, false, num_data,
			SpikeType_kcenter, MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			NotAlignedArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/o FFT T-009
TEST_F(Convolve1DOperation , WithOutFFTInMaskIsNull) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { {
			"input_mask = nullptr", { NUM_WIDTH, num_data }, false, num_data,
			SpikeType_kcenter, MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, NullPointerArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/o FFT T-010
TEST_F(Convolve1DOperation , WithOutFFTInMaskIsNotAligned) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { {
			"Not aligned input_mask", { NUM_WIDTH, num_data }, false, num_data,
			SpikeType_kcenter, MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, NotAlignedArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/o FFT T-011
TEST_F(Convolve1DOperation , WithOutFFTOutDataIsNull) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { { "out_data = nullptr",
			{ NUM_WIDTH, num_data }, false, num_data, SpikeType_kcenter,
			MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			NullPointerArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/o FFT T-012
TEST_F(Convolve1DOperation , WithOutFFTOutDataIsNotAligned) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { {
			"Not aligned out_data", { NUM_WIDTH, num_data }, false, num_data,
			SpikeType_kcenter, MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			NotAlignedArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/o FFT T-013
TEST_F(Convolve1DOperation , WithOutFFTOutWeightIsNull) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { {
			"out_weight = nullptr", { NUM_WIDTH, num_data }, false, num_data,
			SpikeType_kcenter, MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, NullPointerArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/o FFT T-014
TEST_F(Convolve1DOperation , WithOutFFTOutWeightIsNotAligned) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { {
			"Not aligned out_weight", { NUM_WIDTH, num_data }, false, num_data,
			SpikeType_kcenter, MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, NotAlignedArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
// convolve w/o FFT T-015
TEST_F(Convolve1DOperation , WithOutFFTInPlaceConvolution) {
	size_t const num_data = NUM_IN_ODD;
	ConvolveTestComponent<GaussianKernel> TestList[] = { {
			"In-place convolution (&output_data == &input_data)", {
			NUM_WIDTH, NUM_IN_ODD }, false, num_data, SpikeType_kcenter,
			MaskType_knone } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			InPlaceArrayInitializer, StandardArrayInitializer,
			InvalidArgumentValidator>(ELEMENTSOF(TestList), TestList);
}
/*
 * Convolution by narrow Gaussian kernel (kernel_width < num_data = num_kernel)
 */
// convolve w/o FFT T-016, 018 ~ 019
TEST_F(Convolve1DOperation , GaussWithOutFFTVariousNumData) {
	vector<float> gaussian_ref_odd { 0.031861969, 0.069249393, 0.12056985,
			0.16816399, 0.18788746, 0.16816399, 0.12056985, 0.069249393,
			0.031861969 }; // analytic value by emulating the code
	vector<float> gaussian_ref_even { 0.011745105, 0.03186197, 0.069249394,
			0.12056985, 0.16816399, 0.18788747, 0.16816404, 0.1205702,
			0.069251029, 0.031866922, 0.011754746 }; // analytic value by emulating the code
	ConvolveTestComponent<GaussianKernel> TestList[] = {
	// T-016 (num_data = 1: LB)
			{ "Gaussian kernel (num_data = 1: LB)", { NUM_WIDTH, 1 },
			false, 1, SpikeType_kcenter, MaskType_ktrue, { { 0, { 1.0 } } }, { {
					0, { 1.0 } } } },
			// T-018 (num_data=odd)
			{ "Gaussian kernel (num_data = odd)", { NUM_WIDTH, NUM_IN_ODD },
			false, NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { {
			NUM_IN_ODD / 2 - gaussian_ref_odd.size() / 2, gaussian_ref_odd } },
					{ { 0, gauss_weight_ref_L }, { NUM_IN_ODD
							- gauss_weight_ref_R_odd.size(),
							gauss_weight_ref_R_odd } } },
			// T-019 (num_data=even)
			{ "Gaussian kernel (num_data = even)", { NUM_WIDTH, NUM_IN_EVEN },
			false, NUM_IN_EVEN, SpikeType_kcenter, MaskType_ktrue,
					{ {
					NUM_IN_EVEN / 2 - gaussian_ref_even.size() / 2,
							gaussian_ref_even } }, { { 0, gauss_weight_ref_L },
							{ NUM_IN_EVEN - gauss_weight_ref_R_even.size(),
									gauss_weight_ref_R_even } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList);

}
/*
 * num_kernel != num_data
 */
// convolve w/o FFT T-020-023
TEST_F(Convolve1DOperation, GaussWithoutFFTNumKernelVSNumData) {
	vector<float> long_ref_odd { 0.0318619674, 0.0692493947, 0.120569846,
			0.168163988, 0.187887460, 0.168163988, 0.120569846, 0.06924939477,
			0.0318619674 }; // analytic value by emulating the code
	vector<float> long_ref_even { 0.031861967, 0.069249395, 0.12056985,
			0.16816399, 0.18788746, 0.16816399, 0.12056985, 0.069249395,
			0.031861967 }; // analytic value by emulating the code
	vector<float> short_ref_odd { 0.031925101, 0.069388248, 0.12081194,
			0.16850170, 0.18826478, 0.16850170, 0.120811941, 0.069388248,
			0.031925101 };
	vector<float> short_ref_even { 0.032036397, 0.069630146, 0.12123312,
			0.16908912, 0.18892111, 0.16908912, 0.12123312, 0.069630146,
			0.032036397 };
	vector<float> short_weight_L_odd { 0.59413238, 0.76263409, 0.88344604,
			0.95283428, 0.98475938 };
	vector<float> short_weight_R_odd { 0.98475938, 0.95283428, 0.88344604,
			0.76263409, 0.59413238 };
	vector<float> short_weight_L_even { 0.59620363, 0.76529275, 0.88652587,
			0.95615602, 0.98819241 };
	vector<float> short_weight_R_even { 0.98470625, 0.95266985, 0.88303971,
			0.76180659, 0.59271746 };
	ConvolveTestComponent<GaussianKernel> TestList[] =
			{
			// T-020
					{ "num_kernel(odd) > num_data", { NUM_WIDTH, 2 * NUM_IN_ODD
							+ 1 },
					false, NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { {
					NUM_IN_ODD / 2 - long_ref_odd.size() / 2, long_ref_odd } },
							{ { 0, gauss_weight_ref_L }, { NUM_IN_ODD
									- gauss_weight_ref_R_odd.size(),
									gauss_weight_ref_R_odd } } },
					// T-021
					{ "num_kernel(even) > num_data",
							{ NUM_WIDTH, 2 * NUM_IN_ODD }, false,
							NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { {
							NUM_IN_ODD / 2 - long_ref_even.size() / 2,
									long_ref_even } }, {
									{ 0, gauss_weight_ref_L }, { NUM_IN_ODD
											- gauss_weight_ref_R_even.size(),
											gauss_weight_ref_R_even } } },
					// T-022
					{ "num_kernel(odd) < num_data", { NUM_WIDTH, 13 }, false,
					NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue,
							{ {
							NUM_IN_ODD / 2 - short_ref_odd.size() / 2,
									short_ref_odd } }, {
									{ 0, short_weight_L_odd }, { NUM_IN_ODD
											- short_weight_R_odd.size(),
											short_weight_R_odd } } },
					// T-023
					{ "num_kernel(even) < num_data", { NUM_WIDTH, 12 }, false,
					NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue,
							{ {
							NUM_IN_ODD / 2 - short_ref_even.size() / 2,
									short_ref_even } }, { { 0,
									short_weight_L_even }, { NUM_IN_ODD
									- short_weight_R_even.size(),
									short_weight_R_even } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList);
}
/*
 * Convolution by wide Gaussian kernel (kernel_width > num_data)
 */
// convolve w/o FFT T-024 ~ 025
TEST_F(Convolve1DOperation , WideGaussWithOutFFT) {
	vector<float> gaussian_ref_odd { 0.052354915, 0.052003436, 0.051467945,
			0.050736419, 0.049800864, 0.048657602, 0.049800864, 0.050736419,
			0.051467945, 0.052003436, 0.052354915 }; // analytic values obtained by a script which emulates code
	vector<float> gaussian_ref_even { 0.052494097, 0.052322092, 0.051935662,
			0.051320437, 0.050466707, 0.052084927, 0.053482965, 0.054660225,
			0.05562174, 0.056377983, 0.056944642 }; // analytic values obtained by emulating code
	vector<float> weight_ref_odd { 0.52432880, 0.5727872420, 0.62065308730,
			0.66754732840, 0.71311436780, 0.7570300493, 0.7990084570,
			0.838807242, 0.8762313258, 0.9111349117, 0.9434218202, 0.9730442278,
			1.0, 0.9730442278, 0.9434218202, 0.9111349117, 0.8762313258,
			0.8388072420, 0.7990084570, 0.7570300493, 0.7131143678,
			0.6675473284, 0.62065308730, 0.572787242, 0.52432879990 };
	vector<float> weight_ref_even { 0.5387260371, 0.5891927422, 0.63899226640,
			0.68769936040, 0.73491719510, 0.78028698410, 0.82349598920,
			0.8642836264, 0.90244549650, 0.9378352551, 0.9703643782,
			0.999999992, 0.97323899290, 0.94360337910, 0.9110742560,
			0.8756844974, 0.8375226273, 0.7967349901, 0.753525985, 0.708156196,
			0.6609383613, 0.6122312673, 0.5624317431, 0.5119650380 };
	size_t const num_data_odd = NUM_IN_ODD;
	size_t const num_data_even = NUM_IN_EVEN;
	ConvolveTestComponent<GaussianKernel> TestList[] = {
	// T-024 (num_data = odd)
			{ "wide Gaussian kernel (odd)", { NUM_IN_ODD + 1, num_data_odd },
			false, num_data_odd, SpikeType_kcenter, MaskType_ktrue, { {
					num_data_odd / 2 - gaussian_ref_odd.size() / 2,
					gaussian_ref_odd } }, { { 0, weight_ref_odd } } },
			// T-025 (num_data = even)
			{ "wide Gaussian kernel (even)", { NUM_IN_EVEN + 1, num_data_even },
			false, num_data_even, SpikeType_kcenter, MaskType_ktrue, { {
					num_data_even / 2 - gaussian_ref_even.size() / 2,
					gaussian_ref_even } }, { { 0, weight_ref_even } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList);
}
/*
 * Convolution by Gaussian kernel. Data has 2 spikes at edges
 */
// convolve w/o FFT T-026 ~ 027
TEST_F(Convolve1DOperation, GaussWithoutFFTTwoSpikes) {
	vector<float> data_ref_L { 0.316338822, 0.220656452, 0.136595536,
			0.0727463316 };
	vector<float> data_ref_R_odd { 0.072746331688628943, 0.13659553640038027,
			0.22065645273441889, 0.31633882218200382 };
	vector<float> data_ref_R_even { 0.072746333353475745, 0.13659553977170918,
			0.2206564590420598, 0.31633883378509969 };
	ConvolveTestComponent<GaussianKernel> TestList[] = {
	// T-026 (num_data = odd)
			{ "Gaussian kernel (odd) edge spikes data", { NUM_WIDTH,
			NUM_IN_ODD }, false, NUM_IN_ODD, SpikeType_kleftright,
					MaskType_ktrue, { { 0, data_ref_L }, { NUM_IN_ODD
							- data_ref_R_odd.size(), data_ref_R_odd } }, { { 0,
							gauss_weight_ref_L }, {
					NUM_IN_ODD - gauss_weight_ref_R_odd.size(),
							gauss_weight_ref_R_odd } } },
			// T-027 (num_data = even)
			{ "Gaussian kernel (even) edge spikes data", { NUM_WIDTH,
			NUM_IN_EVEN }, false, NUM_IN_EVEN, SpikeType_kleftright,
					MaskType_ktrue, { { 0, data_ref_L }, { NUM_IN_EVEN
							- data_ref_R_even.size(), data_ref_R_even } }, { {
							0, gauss_weight_ref_L }, {
					NUM_IN_EVEN - gauss_weight_ref_R_even.size(),
							gauss_weight_ref_R_even } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList);
}
/*
 * Convolution by Gaussian kernel. Data has 3 spikes at center and edges
 */
// convolve w/o FFT T-028 ~ 029
TEST_F(Convolve1DOperation, GaussWithoutFFTThreeSpikes) {
	vector<float> gaussian_ref_L { 0.316338859, 0.220656819, 0.136598784,
			0.0727711028 };
	vector<float> gaussian_ref_C_odd { 0.120572713, 0.168164267, 0.187887503,
			0.168164267, 0.120572713 };
	vector<float> gaussian_ref_R_odd { 0.0727711028, 0.136598784, 0.220656819,
			0.316338858 };
	vector<float> gaussian_ref_C_even { 0.120572713, 0.168164289, 0.187887747,
			0.16816690, 0.120593775 };
	vector<float> gaussian_ref_R_even { 0.072909543, 0.136622254, 0.220660221,
			0.316339304 };
	ConvolveTestComponent<GaussianKernel> TestList[] = {
	// T-028 (num_data = odd)
			{ "Gaussian kernel (odd) edge spikes data", { NUM_WIDTH,
			NUM_IN_ODD }, false, NUM_IN_ODD, SpikeType_kall, MaskType_ktrue, { {
					0, gaussian_ref_L }, { NUM_IN_ODD / 2
					- gaussian_ref_C_odd.size() / 2, gaussian_ref_C_odd }, {
			NUM_IN_ODD - gaussian_ref_R_odd.size(), gaussian_ref_R_odd } },
					{ { 0, gauss_weight_ref_L }, {
					NUM_IN_ODD - gauss_weight_ref_R_odd.size(),
							gauss_weight_ref_R_odd } } },
			// T-029 (num_data = even)
			{ "Gaussian kernel (even) 3 spikes data",
					{ NUM_WIDTH, NUM_IN_EVEN },
					false, NUM_IN_EVEN, SpikeType_kall, MaskType_ktrue, { { 0,
							gaussian_ref_L }, { NUM_IN_EVEN / 2
							- gaussian_ref_C_even.size() / 2,
							gaussian_ref_C_even }, {
					NUM_IN_EVEN - gaussian_ref_R_even.size(),
							gaussian_ref_R_even } }, {
							{ 0, gauss_weight_ref_L }, {
							NUM_IN_EVEN - gauss_weight_ref_R_even.size(),
									gauss_weight_ref_R_even } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList);

}
/*
 * Convolution by Gaussian kernel. Data has negative spikes at center
 */
// convolve w/ FFT T-030 ~ 031
TEST_F(Convolve1DOperation , GaussWithOutFFTNegativeSpike) {
	vector<float> data_odd { -0.031861969, -0.069249393, -0.12056985,
			-0.16816399, -0.18788746, -0.16816399, -0.12056985, -0.069249393,
			-0.031861969 }; // analytic value by emulating the code
	vector<float> data_even { -0.011745105, -0.03186197, -0.069249394,
			-0.12056985, -0.16816399, -0.18788747, -0.16816404, -0.1205702,
			-0.069251029, -0.031866922, -0.011754746 }; // analytic value by emulating the code
	ReferenceData<float> ref_odd = { NUM_IN_ODD / 2 - data_odd.size() / 2,
			data_odd };
	ReferenceData<float> ref_even = { NUM_IN_EVEN / 2 - data_even.size() / 2,
			data_even };
	ConvolveTestComponent<GaussianKernel> TestList[] = {
	// T-030 (num_data=odd)
			{ "Gaussian kernel (num_data = odd)", { NUM_WIDTH,
			NUM_IN_ODD },
			false, NUM_IN_ODD, SpikeType_knegative, MaskType_ktrue, {
					{ ref_odd } },
					{ { 0, gauss_weight_ref_L }, {
					NUM_IN_ODD - gauss_weight_ref_R_odd.size(),
							gauss_weight_ref_R_odd } } },
			// T-031 (num_data=even)
			{ "Gaussian kernel (num_data = even)", { NUM_WIDTH,
			NUM_IN_EVEN },
			false, NUM_IN_EVEN, SpikeType_knegative, MaskType_ktrue, { {
					ref_even } }, { { 0, gauss_weight_ref_L }, {
			NUM_IN_EVEN - gauss_weight_ref_R_even.size(),
					gauss_weight_ref_R_even } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList);
}
/*
 * Convolution by Gaussian kernel with mask
 */
// convolve w/o FFT T-032 ~ 035
TEST_F(Convolve1DOperation, GaussWithoutFFTVariousMask) {
	vector<float> gaussian_ref { 0.031866919, 0.069251028, 0.12057019,
			0.16816404, 0.18788747, 0.16816404, 0.120570191, 0.069251028,
			0.031866919 }; // analytic value by emulating the code
	vector<float> maskall_weight_ref { 0.40605625, 0.59394345, 0.76210485,
			0.88265394, 0.95177134, 0.98296780 };
	vector<float> maskedge_weight_ref { 0.40605627, 0.59394373, 0.76210772,
			0.88267754, 0.95192670, 0.98378782 };
	ReferenceData<float> const maskall_weight_L = { 0, maskall_weight_ref };
	reverse(maskall_weight_ref.begin(), maskall_weight_ref.end());
	ReferenceData<float> const maskall_weight_R = { NUM_IN_ODD
			- maskall_weight_ref.size(), maskall_weight_ref };
	ReferenceData<float> const maskedge_weight_L = { 0, maskedge_weight_ref };
	reverse(maskedge_weight_ref.begin(), maskedge_weight_ref.end());
	ReferenceData<float> const maskedge_weight_R = { NUM_IN_ODD
			- maskedge_weight_ref.size(), maskedge_weight_ref };
	ReferenceData<float> const zeros = { 0, vector<float>(NUM_IN_ODD, 0.0) };
	ConvolveTestComponent<GaussianKernel> TestList[] = {
	// T-032 (mask non-zero elements)
			{ "mask out non-zero elements", { NUM_WIDTH,
			NUM_IN_ODD }, false, NUM_IN_ODD, SpikeType_kall,
					MaskType_knonzerofalse, { { zeros } }, { maskall_weight_L,
							maskall_weight_R } },
			// T-033 (mask all)
			{ "mask out all elements", { NUM_WIDTH,
			NUM_IN_ODD }, false, NUM_IN_ODD, SpikeType_kall, MaskType_kfalse, {
					{ zeros } }, { { zeros } } },
			// T-034 (mask edge spikes)
			{ "mask out spikes at edges", { NUM_WIDTH,
			NUM_IN_ODD }, false, NUM_IN_ODD, SpikeType_kall,
					MaskType_kleftright, { { NUM_IN_ODD / 2
							- gaussian_ref.size() / 2, gaussian_ref } }, {
							maskedge_weight_L, maskedge_weight_R } },
			// T-035 (mask nan and inf elements at edges)
			{ "mask out nan and inf elements", { NUM_WIDTH,
			NUM_IN_ODD }, false, NUM_IN_ODD, SpikeType_kedgenaninf,
					MaskType_kleftright, { { NUM_IN_ODD / 2
							- gaussian_ref.size() / 2, gaussian_ref } }, {
							maskedge_weight_L, maskedge_weight_R } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList);
}
/*
 * Test Performance of Convolve1D
 */
// convolve w/o FFT T-036 (num_data = even)
TEST_F(Convolve1DOperation, PerformanceTestWithoutFFT) {
	float const kernel_width = NUM_WIDTH;
	size_t const num_kernel = 25;
	size_t const num_data = NUM_IN_LARGE;
	size_t num_repeat = 10000;
	ConvolveTestComponent<GaussianKernel> TestList[] = {
			{ "Gaussian kernel (even)", { kernel_width, num_kernel }, false,
					num_data, SpikeType_kcenter, MaskType_ktrue, {
							{ num_data / 2 - gaussian_ref_fft.size() / 2,
									gaussian_ref_fft } }, { { 0,
							gauss_weight_ref_L }, { num_data
							- gauss_weight_ref_R_odd.size(),
							gauss_weight_ref_R_odd } } } };
	RunConvolveTestComponentList<GaussianKernel, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StatusOKValidator>(ELEMENTSOF(TestList), TestList, num_repeat);
}
/*
 * Test user defined asymmetric kernel (right angled triangle)
 */
// convolve w/o FFT T-037 ~ 042
TEST_F(Convolve1DOperation, TriangleWithoutFFTNumKernelVSNumData) {
	constexpr size_t kKernelWidth = 4;
	ReferenceData<float> ref_data = { 10, triangle_ref };
	ConvolveTestComponent<RightAngledTriangleKernel> TriangleNumKernelTest[] =
			{
			// T-037
					{ "num_kernel(odd) = num_data",
							{ kKernelWidth, NUM_IN_ODD }, false,
							NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, {
									ref_data }, { { 0, triangle_weight_ref_L },
									{ NUM_IN_ODD - triangle_weight_ref_R.size(),
											triangle_weight_ref_R } } },
					// T-038
					{ "num_kernel(even) = num_data",
							{ kKernelWidth, NUM_IN_EVEN },
							false,
							NUM_IN_EVEN, SpikeType_kcenter, MaskType_ktrue, {
									ref_data }, { { 0, triangle_weight_ref_L },
									{ NUM_IN_EVEN
											- triangle_weight_ref_R.size(),
											triangle_weight_ref_R } } },
					// T-039
					{ "num_kernel(odd) > num_data", { kKernelWidth, 2
							* NUM_IN_ODD + 1 },
					false, NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, {
							ref_data }, { { 0, triangle_weight_ref_L }, {
					NUM_IN_ODD - triangle_weight_ref_R.size(),
							triangle_weight_ref_R } } },
					// T-040
					{ "num_kernel(even) > num_data", { kKernelWidth, 2
							* NUM_IN_ODD },
					false,
					NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { ref_data },
							{ { 0, triangle_weight_ref_L }, { NUM_IN_ODD
									- triangle_weight_ref_R.size(),
									triangle_weight_ref_R } } },
					// T-041
					{ "num_kernel(odd) < num_data", { kKernelWidth, 13 }, false,
					NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { ref_data },
							{ { 0, triangle_weight_ref_L }, { NUM_IN_ODD
									- triangle_weight_ref_R.size(),
									triangle_weight_ref_R } } },
					// T-042
					{ "num_kernel(even) < num_data", { kKernelWidth, 12 },
					false,
					NUM_IN_ODD, SpikeType_kcenter, MaskType_ktrue, { ref_data },
							{ { 0, triangle_weight_ref_L }, { NUM_IN_ODD
									- triangle_weight_ref_R.size(),
									triangle_weight_ref_R } } } };
	RunConvolveTestComponentList<RightAngledTriangleKernel,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StatusOKValidator>(
			ELEMENTSOF(TriangleNumKernelTest), TriangleNumKernelTest);
}
// convolve w/o FFT T-043 ~ 044
TEST_F(Convolve1DOperation, WideTriangleWithoutFFT) {
	size_t const num_data_odd = NUM_IN_ODD;
	size_t const num_data_even = NUM_IN_EVEN;
	// reference data
	vector<float> data_odd { 0.010989011, 0.019047619, 0.025, 0.0294117649,
			0.032679739, 0.035087720, 0.036842104, 0.038095238, 0.038961038,
			0.039525692, 0.039855074, 0.04, 0.039999999, 0.0432098749,
			0.046583852, 0.050156740, 0.053968253, 0.058064514, 0.0625,
			0.067340068, 0.072664359, 0.078571431, 0.085185182, 0.092664093,
			0.10121458 };
	vector<float> data_even { 0.01098901, 0.019047620, 0.024999999, 0.029411765,
			0.032679740, 0.035087718, 0.036842105, 0.038095239, 0.038961038296,
			0.039525694, 0.039855071, 0.0399999990, 0.043478260, 0.047138047,
			0.051020409, 0.055172415, 0.059649125, 0.064516128, 0.069852940,
			0.075757580, 0.082352941, 0.089795915, 0.098290600, 0.108108112 };
	vector<float> weight_L_odd { 0.280, 0.32307692, 0.3692308, 0.41846154 };
	vector<float> weight_R_odd { 0.86153846, 0.83076923, 0.79692307, 0.760 };
	vector<float> weight_L_even { 0.30333333, 0.35, 0.4, 0.453333345 };
	vector<float> weight_R_even { 0.85, 0.81666667, 0.78, 0.74 };
	ConvolveTestComponent<RightAngledTriangleKernel> TestList[] = {
	// T-043 (num_data = odd)
			{ "Asynmetric kernel (kernel_width > num_data, num_data = odd)", {
					num_data_odd + 1, num_data_odd }, false, num_data_odd,
					SpikeType_kcenter, MaskType_ktrue, { { 0, data_odd } }, { {
							0, weight_L_odd }, { num_data_odd
							- weight_R_odd.size(), weight_R_odd } } },
			// T-044 (num_data = even)
			{ "Asynmetric kernel (kernel_width > num_data, num_data = even)", {
					num_data_odd + 1, num_data_even },
			false, num_data_even, SpikeType_kcenter, MaskType_ktrue, { { 0,
					data_even } }, { { 0, weight_L_even }, { num_data_even
					- weight_R_even.size(), weight_R_even } } } };
	RunConvolveTestComponentList<RightAngledTriangleKernel,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StatusOKValidator>(ELEMENTSOF(TestList),
			TestList);
}

// convolve w/o FFT T-045 ~ 050
TEST_F(Convolve1DOperation, TriangleWithoutFFTVariousInData) {
	constexpr size_t kKernelWidth = 4;
	size_t const num_data_odd = NUM_IN_ODD;
	size_t const num_data_even = NUM_IN_EVEN;
	ReferenceData<float> ref_data_C = { 10, triangle_ref };
	ReferenceData<float> ref_deta_negative = { 10, { -0.1, -0.2, -0.3, -0.4 } };
	ReferenceData<float> ref_data_L = { 0, { 0.5, 0.4, 0.0 } };
	vector<float> data_R { 0.0, 0.10, 0.22222222, 0.42857143 };
	ReferenceData<float> ref_data_R_odd =
			{ num_data_odd - data_R.size(), data_R };
	ReferenceData<float> ref_data_R_even = { num_data_even - data_R.size(),
			data_R };
	ReferenceData<float> ref_weight_L = { 0, triangle_weight_ref_L };
	ReferenceData<float> ref_weight_R_odd = { num_data_odd
			- triangle_weight_ref_R.size(), triangle_weight_ref_R };
	ReferenceData<float> ref_weight_R_even = { num_data_even
			- triangle_weight_ref_R.size(), triangle_weight_ref_R };
	ConvolveTestComponent<RightAngledTriangleKernel> TestList[] =
			{
// T-045 (num_data = odd)
					{ "Asynmetric kernel (2 spikes at edges, num_data = odd)", {
							kKernelWidth, num_data_odd }, false, num_data_odd,
							SpikeType_kleftright, MaskType_ktrue, { ref_data_L,
									ref_data_R_odd }, { ref_weight_L,
									ref_weight_R_odd } },
					// T-046 (num_data = even)
					{ "Asynmetric kernel (2 spikes at edges, num_data = even)",
							{ kKernelWidth, num_data_even },
							false, num_data_even, SpikeType_kleftright,
							MaskType_ktrue, { ref_data_L, ref_data_R_even }, {
									ref_weight_L, ref_weight_R_even } },
					// T-047 (num_data = odd)
					{
							"Asynmetric kernel (3 spikes at center and edges, num_data = odd)",
							{ kKernelWidth, num_data_odd }, false, num_data_odd,
							SpikeType_kall, MaskType_ktrue, { ref_data_L,
									ref_data_C, ref_data_R_odd }, {
									ref_weight_L, ref_weight_R_odd } },
					// T-048 (num_data = even)
					{
							"Asynmetric kernel (3 spikes at center and edges, num_data = even)",
							{ kKernelWidth, num_data_even },
							false, num_data_even, SpikeType_kall,
							MaskType_ktrue, { ref_data_L, ref_data_C,
									ref_data_R_even }, { ref_weight_L,
									ref_weight_R_even } },
					// T-049 (num_data = odd)
					{
							"Asynmetric kernel (negative spike at center, num_data = odd)",
							{ kKernelWidth, num_data_odd }, false, num_data_odd,
							SpikeType_knegative, MaskType_ktrue, {
									ref_deta_negative }, { ref_weight_L,
									ref_weight_R_odd } },
					// T-050 (num_data = even)
					{
							"Asynmetric kernel (negative spike at center, num_data = even)",
							{ kKernelWidth, num_data_even },
							false, num_data_even, SpikeType_knegative,
							MaskType_ktrue, { ref_deta_negative }, {
									ref_weight_L, ref_weight_R_even } }, };
	RunConvolveTestComponentList<RightAngledTriangleKernel,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StatusOKValidator>(ELEMENTSOF(TestList),
			TestList);
}
// convolve w/o FFT T-051 ~ 054
TEST_F(Convolve1DOperation, TriangleWithoutFFTVariousMask) {
	constexpr size_t kKernelWidth = 4;
	ReferenceData<float> const zeros = { 0, vector<float>(NUM_IN_ODD, 0.0) };
	ReferenceData<float> const maskedge_weight_L = { 0, { 0.3, 0.6, 1.0 } };
	ReferenceData<float> const maskcenter_weight = { 9, { 1.0, 0.9, 0.8, 0.7,
			0.6, 1.0 } };
	ReferenceData<float> const maskedge_weight_R =
			{ 21, { 1.0, 0.9, 0.7, 0.4 } };
	ConvolveTestComponent<RightAngledTriangleKernel> TestList[] = {
// T-051 (mask non-zero elements)
			{ "Asynmetric kernel (mask out non-zero elements)", { kKernelWidth,
			NUM_IN_ODD }, false, NUM_IN_ODD, SpikeType_kall,
					MaskType_knonzerofalse, { zeros }, { maskedge_weight_L,
							maskcenter_weight, maskedge_weight_R } },
			// T-052 (mask all)
			{ "Asynmetric kernel (mask out non-zero elements)", { kKernelWidth,
			NUM_IN_ODD }, false, NUM_IN_ODD, SpikeType_kall, MaskType_kfalse, {
					zeros }, { zeros } },
			// T-053 (mask edge spikes)
			{ "Asynmetric kernel (mask out spikes at edges)", { kKernelWidth,
			NUM_IN_ODD }, false, NUM_IN_ODD, SpikeType_kall,
					MaskType_kleftright, { { 10, triangle_ref } }, {
							maskedge_weight_L, maskedge_weight_R } },
			// T-054 (mask nan and inf elements at edges)
			{ "Asynmetric kernel (mask out nan and inf elements)", {
					kKernelWidth, NUM_IN_ODD }, false, NUM_IN_ODD,
					SpikeType_kedgenaninf, MaskType_kleftright, { { 10,
							triangle_ref } }, { maskedge_weight_L,
							maskedge_weight_R } } };
	RunConvolveTestComponentList<RightAngledTriangleKernel,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StandardArrayInitializer,
			StandardArrayInitializer, StatusOKValidator>(ELEMENTSOF(TestList),
			TestList);
}

/**
 * Test for sakura_CreateGaussianKernel
 */
// T-001
TEST_GAUSS(NegativeKernelWidth) {
	RunGaussianTest<StandardArrayInitializer, InvalidArgumentValidator>(5.0f,
			-1.0f, 10);
}

// T-002
TEST_GAUSS(ZeroKernelWidth) {
	RunGaussianTest<StandardArrayInitializer, InvalidArgumentValidator>(5.0f,
			0.0f, 10);
}

// T-003
TEST_GAUSS(PositiveInfiniteKernelWidth) {
	constexpr float kOne = 1.0f;
	constexpr float kZero = 0.0f;
	float positive_inf = kOne / kZero;
	ASSERT_TRUE(std::isinf(positive_inf));
	ASSERT_GT(positive_inf, FLT_MAX);
	RunGaussianTest<StandardArrayInitializer, InvalidArgumentValidator>(5.0f,
			positive_inf, 10);
}

// T-004
TEST_GAUSS(NegativeInfiniteKernelWidth) {
	constexpr float kOne = 1.0f;
	constexpr float kZero = 0.0f;
	float negative_inf = -kOne / kZero;
	ASSERT_TRUE(std::isinf(negative_inf));
	ASSERT_LT(negative_inf, FLT_MIN);
	RunGaussianTest<StandardArrayInitializer, InvalidArgumentValidator>(5.0f,
			negative_inf, 10);
}

// T-005
TEST_GAUSS(NotANumberKernelWidth) {
	float nan_value = std::nan("");
	ASSERT_TRUE(std::isnan(nan_value));
	RunGaussianTest<StandardArrayInitializer, InvalidArgumentValidator>(5.0f,
			nan_value, 10);
}

// T-006
TEST_GAUSS(KernelIsNull) {
	RunGaussianTest<NullPointerArrayInitializer, InvalidArgumentValidator>(5.0f,
			2.0f, 10);
}

// T-007
TEST_GAUSS(KernelIsNotAligned) {
	RunGaussianTest<NotAlignedArrayInitializer, InvalidArgumentValidator>(5.0f,
			2.0f, 10);
}

// T-008
TEST_GAUSS(NumKernelIsZero) {
	RunGaussianTest<StandardArrayInitializer, InvalidArgumentValidator>(0.0f,
			2.0f, 0);
}

// T-009
TEST_GAUSS(NumKernelIsOne) {
// result should be independent of kernel_width if num_kernel is 1
	RunGaussianTest<StandardArrayInitializer, NumKernelOneValidator>(0.0f, 2.0f,
			1);

	RunGaussianTest<StandardArrayInitializer, NumKernelOneValidator>(0.0f,
			256.0f, 1);
}

// T-010
TEST_GAUSS(NumKernelIsTwo) {
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator>(1.0f,
			2.0f, 2);
}

// T-011
TEST_GAUSS(NumKernelIsThree) {
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator>(1.0f,
			2.0f, 3);
}

// T-012
TEST_GAUSS(NumKernelIsFour) {
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator>(2.0f,
			2.0f, 4);
}

// T-014
TEST_GAUSS(EvenNumkernel) {
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator>(32.0f,
			24.0f, 64);
}

// T-013
TEST_GAUSS(OddNumkernel) {
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator>(32.0f,
			24.0f, 65);
}

// T-016
TEST_GAUSS(NumkernelLessThanKernelWidth) {
// even
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator>(32.0f,
			128.0f, 64);

	// odd
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator>(32.0f,
			128.0f, 65);
}

// T-017
TEST_GAUSS(NumKernelGreaterThanKernelWidth) {
// even
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator>(64.0f,
			64.0f, 128);

	// odd
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator>(64.0f,
			64.0f, 129);
}

// T-018
// 2016/09/13 TN
// Disable this test since now calculation of Gaussian at k is changed to definite
// integral between k - 1/2 and k + 1/2 so idea of validation in HalfWidthValidator
// doesn't make sense anymore.
//TEST_GAUSS(HalfWidth) {
//	// even
//	RunGaussianTest<StandardArrayInitializer, HalfWidthValidator>(64.0f, 128.0f,
//			128);
//
//	// odd
//	RunGaussianTest<StandardArrayInitializer, HalfWidthValidator>(64.0f, 129.0f,
//			129);
//}

// T-019
TEST_GAUSS(NarrowKernelWidth) {
	RunGaussianTest<StandardArrayInitializer, NarrowKernelValidator>(64.0f,
	FLT_MIN, 128);
}

// T-020
TEST_GAUSS(WideKernelWidth) {
	RunGaussianTest<StandardArrayInitializer, WideKernelValidator>(64.0f,
	FLT_MAX, 128);
}

// T-021
TEST_GAUSS(PowerOfTwoNumKernelPerformance) {
	constexpr size_t kNumKernel = 268435456; // 2^28
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator,
			PerformanceTestLogger>(static_cast<float>(kNumKernel / 2), 128.0f,
			kNumKernel, "CreateGaussian_PowerOfTwoNumKernelPerformanceTest");
}

// T-022
TEST_GAUSS(OddNumKernelPerformance) {
	constexpr size_t kNumKernel = 268435456 + 1; // 2^28 + 1
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator,
			PerformanceTestLogger>(static_cast<float>(kNumKernel / 2), 128.0f,
			kNumKernel, "CreateGaussian_OddNumKernelPerformanceTest");
}

// T-023
TEST_GAUSS(EvenNumkernelSymmetric) {
// even
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator>(31.5f,
			24.0f, 64);
}

// T-024
TEST_GAUSS(PeaklocationIsNegative) {
	RunGaussianTest<StandardArrayInitializer, BasicGaussKernelValidator>(-1.0f,
			24.0f, 64);
}

// T-025
TEST_GAUSS(PeaklocationFarFromCenterOffPeakSymmetric) {
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator>(21.5f,
			24.0f, 64);
}

// T-026
TEST_GAUSS(PeaklocationFarFromCenterOnPeakSymmetric) {
	RunGaussianTest<StandardArrayInitializer, SymmetricKernelValidator>(42.0f,
			24.0f, 64);
}

// T-027
TEST_GAUSS(PeaklocationIsOffPixel) {
	RunGaussianTest<StandardArrayInitializer, BasicGaussKernelValidator>(32.3f,
			24.0f, 64);
}

// T-028
TEST_GAUSS(PeaklocationIsGraterThanNumkernel) {
	RunGaussianTest<StandardArrayInitializer, BasicGaussKernelValidator>(70.0f,
			24.0f, 64);
}

/**
 * Test for LIBSAKURA_SYMBOL(Convolve1DFloat)
 */
