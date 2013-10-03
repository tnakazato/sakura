#include <iostream>
#include <string>
#include <sys/time.h>
#include <math.h>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "aligned_memory.h"
#include "gtest/gtest.h"

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
};
TEST_F(ApplyCalibrationTest, ZeroLengthData) {
	size_t const num_scaling_factor = 1;
	size_t const num_data = 0;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	scaling_factor[0] = 1.0;
	SIMD_ALIGN
	float target[num_data];

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(ApplyPositionSwitchCalibration)(num_scaling_factor,
			scaling_factor, num_data, target, target, target);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
}

TEST_F(ApplyCalibrationTest, ZeroLengthScalingFactor) {
	size_t const num_scaling_factor = 0;
	size_t const num_data = 1;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	SIMD_ALIGN
	float target[num_data];
	target[0] = 1.0;

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(ApplyPositionSwitchCalibration)(num_scaling_factor,
			scaling_factor, num_data, target, target, target);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(ApplyCalibrationTest, NullPointer) {
	size_t const num_scaling_factor = 1;
	size_t const num_data = 1;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	scaling_factor[0] = 1.0;
	SIMD_ALIGN
	float target[num_data];
	target[0] = 1.0;

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(ApplyPositionSwitchCalibration)(num_scaling_factor,
			scaling_factor, num_data, target, nullptr, target);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(ApplyCalibrationTest, InputArrayNotAligned) {
	size_t const num_scaling_factor = 1;
	size_t const num_data = 1;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	scaling_factor[0] = 1.0;
	SIMD_ALIGN
	float target[num_data + 1];
	target[0] = 1.0;
	target[1] = 1.0;

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(ApplyPositionSwitchCalibration)(num_scaling_factor,
			scaling_factor, num_data, &target[1], target, target);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(ApplyCalibrationTest, InvaidNumberOfScalingFactor) {
	size_t const num_scaling_factor = 2;
	size_t const num_data = 3;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	scaling_factor[0] = 1.0;
	scaling_factor[1] = 1.0;
	SIMD_ALIGN
	float target[num_data];
	target[0] = 1.0;
	target[1] = 1.0;
	target[2] = 1.0;

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(ApplyPositionSwitchCalibration)(num_scaling_factor,
			scaling_factor, num_data, target, target, target);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), status);
}

TEST_F(ApplyCalibrationTest, ZeroDivision) {
	size_t const num_scaling_factor = 1;
	size_t const num_data = 1;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	scaling_factor[0] = 1.0;
	SIMD_ALIGN
	float target[num_data];
	target[0] = 1.0;
	SIMD_ALIGN
	float reference[num_data];
	reference[0] = 0.0;
	SIMD_ALIGN
	float result[num_data];

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(ApplyPositionSwitchCalibration)(num_scaling_factor,
			scaling_factor, num_data, target, reference, result);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	EXPECT_TRUE(isinf(result[0]));
}

TEST_F(ApplyCalibrationTest, BasicTest) {
	size_t const num_scaling_factor = 1;
	size_t const num_data = 1;
	SIMD_ALIGN
	float scaling_factor[num_scaling_factor];
	scaling_factor[0] = 1.0;
	SIMD_ALIGN
	float target[num_data];
	target[0] = 1.0;
	float expected[num_data];
	expected[0] = 0.0;

	LIBSAKURA_SYMBOL(Status) status =
	LIBSAKURA_SYMBOL(ApplyPositionSwitchCalibration)(num_scaling_factor,
			scaling_factor, num_data, target, target, target);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0; i < num_data; ++i) {
		EXPECT_EQ(expected[0], target[0]);
	}
}
