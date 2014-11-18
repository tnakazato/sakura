/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2014
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
#include <sys/time.h>
#include <stdarg.h>
#include <math.h>
#include <cfloat>
#include <memory>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "loginit.h"
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
			size_t num_data, size_t num_segments, float const target[],
			float const reference[], float result[], float const expected[],
			bool check_result, size_t iteration = 1) {
		size_t memory_size = (num_scaling_factor + num_data * num_segments * ((target == result) ? 2 : 3)) * 4; // bytes
		std::cout << "memory usage: " << (double)memory_size / 1.0e9 << "GB" << std::endl;
		std::string message =
				(expected_status == sakura_Status_kOK) ?
						"ApplyCalibration had any problems during execution." :
						"ApplyCalibration should fail!";
		float *target_saved = nullptr;
		std::unique_ptr<float[]> storage_for_saved;
		if (target == result) {
			size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
			size_t num_arena = num_data + sakura_alignment - 1;
			storage_for_saved.reset(new float[num_arena]);
			target_saved = sakura_AlignFloat(num_arena, storage_for_saved.get(),
					num_data);
			for (size_t i = 0; i < num_data; ++i) {
				target_saved[i] = target[i];
			}
		}
		double elapsed_time = 0.0;
		for (size_t iter = 0; iter < iteration; ++iter) {
			if (target == result) {
				float *tmp = const_cast<float *>(target);
				for (size_t i = 0; i < num_data; ++i) {
					tmp[i] = target_saved[i];
				}
			}
			for (size_t i = 0; i < num_segments; ++i) {
				float const *ws =
						(num_scaling_factor == 1) ?
								scaling_factor : &scaling_factor[i * num_data];
				float const *wt = &target[i * num_data];
				float const *wf = &reference[i * num_data];
				float *wr = &result[i * num_data];
				double start = sakura_GetCurrentTime();
				LIBSAKURA_SYMBOL(Status) result_status =
				LIBSAKURA_SYMBOL(ApplyPositionSwitchCalibration)(
						num_scaling_factor, ws, num_data, wt, wf, wr);
				double end = sakura_GetCurrentTime();
				elapsed_time += end - start;
				EXPECT_EQ(expected_status, result_status) << message;

				if (check_result && (expected_status == result_status)) {
					// Value check
					for (size_t index = i * num_data; index < (i + 1) * num_data; ++index) {
						std::cout << "Expected value at index " << index << ": "
								<< expected[index] << std::endl;
						EXPECT_FLOAT_EQ(expected[index], result[index])
								<< "calibrated value differs from expected value at "
								<< index << ": " << expected[index] << ", "
								<< result[index];
					}
				}
			}
		}
		std::cout << "Elapsed time (" << num_segments << " segments with "
				<< iteration << " iterations) " << elapsed_time << " sec"
				<< std::endl;
	}
	template<size_t ITER, size_t NUM_ARRAY, size_t NUM_CHAN, size_t NUM_SCALING, bool IN_PLACE>
	void RunPerformanceTest() {
		size_t const num_data = NUM_CHAN * NUM_ARRAY;
		size_t const num_scaling_factor = get_num_scaling_factor<NUM_SCALING>(num_data);
		size_t sakura_alignment = LIBSAKURA_SYMBOL(GetAlignment)();
		size_t num_arena = num_scaling_factor + sakura_alignment - 1;
		std::unique_ptr<float[]> storage_for_scaling_factor(new float[num_arena]);
		float *scaling_factor = sakura_AlignFloat(num_arena,
				storage_for_scaling_factor.get(), num_scaling_factor);
		for (size_t i = 0; i < num_scaling_factor; ++i)
			scaling_factor[i] = 1.0;
                num_arena = num_data + sakura_alignment - 1;
		std::unique_ptr<float[]> storage_for_target(new float[num_arena]);
		float *target = sakura_AlignFloat(num_arena, storage_for_target.get(),
				num_data);
		std::unique_ptr<float[]> storage_for_reference(new float[num_arena]);
		float *reference = sakura_AlignFloat(num_arena, storage_for_reference.get(),
				num_data);
		float *storage = nullptr;
		float *result = get_result_pointer<IN_PLACE>(num_data, target, num_arena, &storage);
		std::unique_ptr<float[]> storage_for_result(storage);
		for (size_t i = 0; i < num_data; ++i) {
			target[i] = 2.0;
			;
			reference[i] = 1.0;
		}

		PerformTest(LIBSAKURA_SYMBOL(Status_kOK), NUM_SCALING,
				scaling_factor, NUM_CHAN, NUM_ARRAY, target, reference, result, nullptr,
				false, ITER);
	}
private:
	template<size_t NUM_SCALING>
	size_t get_num_scaling_factor(size_t num_data) {return num_data;}
	template<bool IN_PLACE>
	float *get_result_pointer(size_t num_data, float const *target, size_t num_storage, float **storage) {
		return const_cast<float *>(target);
	}
};

template<>
size_t ApplyCalibrationTest::get_num_scaling_factor<1>(size_t num_data) {return 1;}
template<>
float *ApplyCalibrationTest::get_result_pointer<false>(size_t num_data, float const *target, size_t num_storage, float **storage) {
	*storage = new float[num_storage];
	return sakura_AlignFloat(num_storage, *storage, num_data);
}

TEST_F(ApplyCalibrationTest, ZeroLengthData) {
	size_t const num_scaling_factor = 1;
	size_t const num_data = 0;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	InitializeFloatArray(num_scaling_factor, scaling_factor, 1.0);
	SIMD_ALIGN
	float target[num_data];

	PerformTest(LIBSAKURA_SYMBOL(Status_kOK), num_scaling_factor,
			scaling_factor, num_data, 1, target, target, target, nullptr, true);
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
			scaling_factor, num_data, 1, target, target, target, nullptr,
			false);
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
			scaling_factor, num_data, 1, target, nullptr, target, nullptr,
			false);
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
			scaling_factor, num_data, 1, &target[1], target, target, nullptr,
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
			scaling_factor, num_data, 1, target, target, target, nullptr,
			false);
}

TEST_F(ApplyCalibrationTest, ZeroDivision) {
	size_t const num_scaling_factor = 1;
	size_t const num_data = 2;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	InitializeFloatArray(num_scaling_factor, scaling_factor, 1.0);
	SIMD_ALIGN
	float target[num_data];
	InitializeFloatArray(num_data, target, 1.0, -1.0);
	SIMD_ALIGN
	float reference[num_data];
	InitializeFloatArray(num_data, reference, 0.0, 0.0);
	SIMD_ALIGN
	float result[num_data];

	PerformTest(LIBSAKURA_SYMBOL(Status_kOK), num_scaling_factor,
			scaling_factor, num_data, 1, target, reference, result, nullptr,
			false);

	// Check result is inf or not.
	// Since isinf sometimes didn't work, additional constraint is added.
	EXPECT_TRUE(isinf(result[0])) << "result must be inf! (" << result[0]
			<< ")";
	EXPECT_TRUE(isinf(result[1])) << "result must be -inf! (" << result[1]
			<< ")";
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
			scaling_factor, num_data, 1, target, reference, result, expected,
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
			scaling_factor, num_data, 1, target, reference, target, expected,
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
			scaling_factor, num_data, 1, target, reference, result, expected,
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
			scaling_factor, num_data, 1, target, reference, result, expected,
			true);
}
TEST_F(ApplyCalibrationTest, IntrinsicsTest) {
	size_t const num_scaling_factor = 10;
	size_t const num_data = 10;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	SIMD_ALIGN
	float target[num_data];
	SIMD_ALIGN
	float reference[num_data];
	SIMD_ALIGN
	float result[num_data];
	SIMD_ALIGN
	float expected[num_data];
	for (size_t i = 0; i < num_data; ++i) {
		scaling_factor[i] = 10.0;
		target[i] = 2.0;
		reference[i] = 1.0;
		expected[i] = 10.0;
	}

	PerformTest(LIBSAKURA_SYMBOL(Status_kOK), num_scaling_factor,
			scaling_factor, num_data, 1, target, reference, result, expected,
			true);
}

TEST_F(ApplyCalibrationTest, PerformanceTestAllAtOnce) {
	RunPerformanceTest<10, 1, 40000000, 40000000, false>();
}

TEST_F(ApplyCalibrationTest, PerformanceTestIndividual) {
	RunPerformanceTest<10, 1000, 40000, 40000, false>();
}

TEST_F(ApplyCalibrationTest, PerformanceTestSingleScalingFactor) {
	RunPerformanceTest<10, 1, 40000000, 1, false>();
}

TEST_F(ApplyCalibrationTest, PerformanceTestShareInputOutput) {
	RunPerformanceTest<50, 1000, 40000, 40000, true>();
}

TEST_F(ApplyCalibrationTest, PerformanceTestShareInputOutputSingleScalingFactor) {
	RunPerformanceTest<50, 1000, 40000, 1, true>();
}
