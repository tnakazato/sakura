#include <iostream>
#include <string>
#include <sys/time.h>
#include <stdarg.h>
#include <cmath>
#include <cfloat>
#include <memory>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "aligned_memory.h"
#include "gtest/gtest.h"

void InitializeFloatArray(size_t num_array, float array[], ...) {
	va_list arguments_list;
	va_start(arguments_list, array);
	for (size_t i = 0; i < num_array; ++i) {
		array[i] = static_cast<float>(va_arg(arguments_list, double));
	}
}

class ApplyCalibrationTest: public ::testing::Test {
public:
	virtual void SetUp() {
		// Initialize sakura
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Initialize)(nullptr,
				nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}
	virtual void TearDown() {
		LIBSAKURA_SYMBOL(CleanUp)();
	}
	virtual void PerformTest(LIBSAKURA_SYMBOL(Status) expected_status,
			size_t num_scaling_factor, float const scaling_factor[],
			size_t num_data, float const target[], float const reference[],
			float result[], float const expected[], bool check_result) {
		std::string message =
				(expected_status == sakura_Status_kOK) ?
						"ApplyCalibration had any problems during execution." :
						"ApplyCalibration should fail!";
		LIBSAKURA_SYMBOL(Status) result_status =
		LIBSAKURA_SYMBOL(ApplyPositionSwitchCalibration)(num_scaling_factor,
				scaling_factor, num_data, target, reference, result);

		EXPECT_EQ(expected_status, result_status) << message;

		if (check_result && (expected_status == result_status)) {
			// Value check
			for (size_t index = 0; index < num_data; ++index) {
				std::cout << "Expected value at index " << index << ": "
						<< expected[index] << std::endl;
				EXPECT_FLOAT_EQ(expected[index], result[index])
						<< "interpolated value differs from expected value at "
						<< index << ": " << expected[index] << ", "
						<< result[index];
			}
		}

	}
};
TEST_F(ApplyCalibrationTest, ZeroLengthData) {
	size_t const num_scaling_factor = 1;
	size_t const num_data = 0;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	InitializeFloatArray(num_scaling_factor, scaling_factor, 1.0);
	SIMD_ALIGN
	float target[num_data];

	PerformTest(LIBSAKURA_SYMBOL(Status_kOK), num_scaling_factor,
			scaling_factor, num_data, target, target, target, nullptr, true);
}

TEST_F(ApplyCalibrationTest, ZeroLengthScalingFactor) {
	size_t const num_scaling_factor = 0;
	size_t const num_data = 1;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	SIMD_ALIGN
	float target[num_data];
	InitializeFloatArray(num_data, target, 1.0);

	PerformTest(LIBSAKURA_SYMBOL(Status_kInvalidArgument), num_scaling_factor,
			scaling_factor, num_data, target, target, target, nullptr, false);
}

TEST_F(ApplyCalibrationTest, NullPointer) {
	size_t const num_scaling_factor = 1;
	size_t const num_data = 1;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	InitializeFloatArray(num_scaling_factor, scaling_factor, 1.0);
	SIMD_ALIGN
	float target[num_data];
	InitializeFloatArray(num_data, target, 1.0);

	PerformTest(LIBSAKURA_SYMBOL(Status_kInvalidArgument), num_scaling_factor,
			scaling_factor, num_data, target, nullptr, target, nullptr, false);
}

TEST_F(ApplyCalibrationTest, InputArrayNotAligned) {
	size_t const num_scaling_factor = 1;
	size_t const num_data = 1;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	InitializeFloatArray(num_scaling_factor, scaling_factor, 1.0);
	SIMD_ALIGN
	float target[num_data + 1];
	InitializeFloatArray(num_data + 1, target, 1.0, 1.0);

	PerformTest(LIBSAKURA_SYMBOL(Status_kInvalidArgument), num_scaling_factor,
			scaling_factor, num_data, &target[1], target, target, nullptr,
			false);
}

TEST_F(ApplyCalibrationTest, InvaidNumberOfScalingFactor) {
	size_t const num_scaling_factor = 2;
	size_t const num_data = 3;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	InitializeFloatArray(num_scaling_factor, scaling_factor, 1.0, 1.0);
	SIMD_ALIGN
	float target[num_data];
	InitializeFloatArray(num_data, target, 1.0, 1.0, 1.0);

	PerformTest(LIBSAKURA_SYMBOL(Status_kInvalidArgument), num_scaling_factor,
			scaling_factor, num_data, target, target, target, nullptr, false);
}

TEST_F(ApplyCalibrationTest, ZeroDivision) {
	size_t const num_scaling_factor = 1;
	size_t const num_data = 1;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	InitializeFloatArray(num_scaling_factor, scaling_factor, 1.0);
	SIMD_ALIGN
	float target[num_data];
	InitializeFloatArray(num_data, target, 1.0);
	SIMD_ALIGN
	float reference[num_data];
	InitializeFloatArray(num_data, reference, 0.0);
	SIMD_ALIGN
	float result[num_data];

	PerformTest(LIBSAKURA_SYMBOL(Status_kOK), num_scaling_factor,
			scaling_factor, num_data, target, reference, result, nullptr,
			false);

	// Check result is inf or not.
	// Since isinf sometimes didn't work, additional constraint is added.
	EXPECT_TRUE(isinf(result[0]) || result[0] > FLT_MAX)
			<< "result must be inf! (" << result[0] << ")";
}

TEST_F(ApplyCalibrationTest, BasicTest) {
	size_t const num_scaling_factor = 2;
	size_t const num_data = 2;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	InitializeFloatArray(num_scaling_factor, scaling_factor, 0.5, 1.0);
	SIMD_ALIGN
	float target[num_data];
	InitializeFloatArray(num_data, target, 1.0, 1.0);
	SIMD_ALIGN
	float reference[num_data];
	InitializeFloatArray(num_data, reference, 1.0, 0.5);
	SIMD_ALIGN
	float result[num_data];
	SIMD_ALIGN
	float expected[num_data];
	InitializeFloatArray(num_data, expected, 0.0, 1.0);

	PerformTest(LIBSAKURA_SYMBOL(Status_kOK), num_scaling_factor,
			scaling_factor, num_data, target, reference, result, expected,
			true);
}

TEST_F(ApplyCalibrationTest, ShareInputOutputStorage) {
	size_t const num_scaling_factor = 2;
	size_t const num_data = 2;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	InitializeFloatArray(num_scaling_factor, scaling_factor, 0.5, 2.0);
	SIMD_ALIGN
	float target[num_data];
	InitializeFloatArray(num_data, target, 1.0, 1.0);
	SIMD_ALIGN
	float reference[num_data];
	InitializeFloatArray(num_data, reference, 1.0, 0.5);
	SIMD_ALIGN
	float expected[num_data];
	InitializeFloatArray(num_data, expected, 0.0, 2.0);
	SIMD_ALIGN
	float original_target[num_data];
	for (size_t i = 0; i < num_data; ++i)
		original_target[i] = target[i];

	PerformTest(LIBSAKURA_SYMBOL(Status_kOK), num_scaling_factor,
			scaling_factor, num_data, target, reference, target, expected,
			true);

	// additional test:
	// make sure that input storage is overwritten
	for (size_t i = 0; i < num_data; ++i) {
		std::cout << "Original value at " << i << ": " << original_target[i]
				<< std::endl;
		EXPECT_NE(original_target[i], target[i])
				<< "storage must be overwritten (" << i << ": "
				<< original_target[i] << ", " << target[i] << ")";
	}
}

TEST_F(ApplyCalibrationTest, SingleScalingFactor) {
	size_t const num_scaling_factor = 1;
	size_t const num_data = 2;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	InitializeFloatArray(num_scaling_factor, scaling_factor, 0.5);
	SIMD_ALIGN
	float target[num_data];
	InitializeFloatArray(num_data, target, 1.0, 1.0);
	SIMD_ALIGN
	float reference[num_data];
	InitializeFloatArray(num_data, reference, 1.0, 0.5);
	SIMD_ALIGN
	float result[num_data];
	SIMD_ALIGN
	float expected[num_data];
	InitializeFloatArray(num_data, expected, 0.0, 0.5);

	PerformTest(LIBSAKURA_SYMBOL(Status_kOK), num_scaling_factor,
			scaling_factor, num_data, target, reference, result, expected,
			true);
}

TEST_F(ApplyCalibrationTest, TooManyScalingFactor) {
	size_t const num_scaling_factor = 3;
	size_t const num_data = 2;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	InitializeFloatArray(num_scaling_factor, scaling_factor, 0.5, 1.0, -5.0);
	SIMD_ALIGN
	float target[num_data];
	InitializeFloatArray(num_data, target, 1.0, 1.0);
	SIMD_ALIGN
	float reference[num_data];
	InitializeFloatArray(num_data, reference, 1.0, 0.5);
	SIMD_ALIGN
	float result[num_data];
	SIMD_ALIGN
	float expected[num_data];
	InitializeFloatArray(num_data, expected, 0.0, 1.0);

	PerformTest(LIBSAKURA_SYMBOL(Status_kOK), num_scaling_factor,
			scaling_factor, num_data, target, reference, result, expected,
			true);
}

TEST_F(ApplyCalibrationTest, PerformanceTestAllAtOnce) {
	size_t const iteration = 1000;
	size_t const length = 400000;
	size_t const num_scaling_factor = length * iteration;
	size_t const num_data = length * iteration;
	size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
	size_t num_arena = num_scaling_factor + sakura_alignment - 1;
	std::unique_ptr<float[]> storage_for_scaling_factor(new float[num_arena]);
	float *scaling_factor = sakura_AlignFloat(num_arena,
			storage_for_scaling_factor.get(), num_scaling_factor);
	for (size_t i = 0; i < num_scaling_factor; ++i)
		scaling_factor[i] = 1.0;
	std::unique_ptr<float[]> storage_for_target(new float[num_arena]);
	float *target = sakura_AlignFloat(num_arena, storage_for_target.get(),
			num_data);
	std::unique_ptr<float[]> storage_for_reference(new float[num_arena]);
	float *reference = sakura_AlignFloat(num_arena, storage_for_reference.get(),
			num_data);
	std::unique_ptr<float[]> storage_for_result(new float[num_arena]);
	float *result = sakura_AlignFloat(num_arena, storage_for_result.get(),
			num_data);
	for (size_t i = 0; i < num_data; ++i) {
		target[i] = 2.0;;
		reference[i] = 1.0;
	}

	PerformTest(LIBSAKURA_SYMBOL(Status_kOK), num_scaling_factor,
			scaling_factor, num_data, target, reference, result, nullptr,
			false);
}

TEST_F(ApplyCalibrationTest, PerformanceTestIndividual) {
	size_t const iteration = 1000;
	size_t const length = 400000;
	size_t const num_scaling_factor = length * iteration;
	size_t const num_data = length * iteration;
	size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
	size_t num_arena = num_scaling_factor + sakura_alignment - 1;
	std::unique_ptr<float[]> storage_for_scaling_factor(new float[num_arena]);
	float *scaling_factor = sakura_AlignFloat(num_arena,
			storage_for_scaling_factor.get(), num_scaling_factor);
	for (size_t i = 0; i < num_scaling_factor; ++i)
		scaling_factor[i] = 1.0;
	std::unique_ptr<float[]> storage_for_target(new float[num_arena]);
	float *target = sakura_AlignFloat(num_arena, storage_for_target.get(),
			num_data);
	std::unique_ptr<float[]> storage_for_reference(new float[num_arena]);
	float *reference = sakura_AlignFloat(num_arena, storage_for_reference.get(),
			num_data);
	std::unique_ptr<float[]> storage_for_result(new float[num_arena]);
	float *result = sakura_AlignFloat(num_arena, storage_for_result.get(),
			num_data);
	for (size_t i = 0; i < num_data; ++i) {
		target[i] = 2.0;
		reference[i] = 1.0;
	}

	for (size_t i = 0; i < iteration; ++i) {
		float const *scaling_factor_local = &scaling_factor[i * length];
		float const *target_local = &target[i * length];
		float const *reference_local = &reference[i * length];
		float *result_local = &result[i * length];
		PerformTest(LIBSAKURA_SYMBOL(Status_kOK), length, scaling_factor_local,
				length, target_local, reference_local, result_local, nullptr,
				false);
	}
}

TEST_F(ApplyCalibrationTest, PerformanceTestSingleScalingFactor) {
	size_t const iteration = 1000;
	size_t const length = 400000;
	size_t const num_scaling_factor = length * iteration;
	size_t const num_data = length * iteration;
	size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
	size_t num_arena = num_scaling_factor + sakura_alignment - 1;
	std::unique_ptr<float[]> storage_for_scaling_factor(new float[num_arena]);
	float *scaling_factor = sakura_AlignFloat(num_arena,
			storage_for_scaling_factor.get(), num_scaling_factor);
	for (size_t i = 0; i < num_scaling_factor; ++i)
		scaling_factor[i] = 1.0;
	std::unique_ptr<float[]> storage_for_target(new float[num_arena]);
	float *target = sakura_AlignFloat(num_arena, storage_for_target.get(),
			num_data);
	std::unique_ptr<float[]> storage_for_reference(new float[num_arena]);
	float *reference = sakura_AlignFloat(num_arena, storage_for_reference.get(),
			num_data);
	std::unique_ptr<float[]> storage_for_result(new float[num_arena]);
	float *result = sakura_AlignFloat(num_arena, storage_for_result.get(),
			num_data);
	for (size_t i = 0; i < num_data; ++i) {
		target[i] = 2.0;
		reference[i] = 1.0;
	}

	PerformTest(LIBSAKURA_SYMBOL(Status_kOK), 1, scaling_factor, num_data,
			target, reference, result, nullptr,
			false);
}

TEST_F(ApplyCalibrationTest, PerformanceTestShareInputOutput) {
	size_t const iteration = 1000;
	size_t const length = 400000;
	size_t const num_scaling_factor = length * iteration;
	size_t const num_data = length * iteration;
	size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
	size_t num_arena = num_scaling_factor + sakura_alignment - 1;
	std::unique_ptr<float[]> storage_for_scaling_factor(new float[num_arena]);
	float *scaling_factor = sakura_AlignFloat(num_arena,
			storage_for_scaling_factor.get(), num_scaling_factor);
	for (size_t i = 0; i < num_scaling_factor; ++i)
		scaling_factor[i] = 1.0;
	std::unique_ptr<float[]> storage_for_target(new float[num_arena]);
	float *target = sakura_AlignFloat(num_arena, storage_for_target.get(),
			num_data);
	std::unique_ptr<float[]> storage_for_reference(new float[num_arena]);
	float *reference = sakura_AlignFloat(num_arena, storage_for_reference.get(),
			num_data);
	for (size_t i = 0; i < num_data; ++i) {
		target[i] = 2.0;
		reference[i] = 1.0;
	}

	PerformTest(LIBSAKURA_SYMBOL(Status_kOK), num_scaling_factor,
			scaling_factor, num_data, target, reference, target, nullptr,
			false);
}
