#include <iostream>

#include <libsakura/sakura.h>

#include "gtest/gtest.h"

#include "aligned_memory.h"

class Interpolate1dFloatTest: public ::testing::Test {
protected:
	virtual void SetUp() {
		initialize_result_ = sakura_Initialize();
	}
	virtual void TearDown() {
		sakura_CleanUp();
	}
	sakura_Status initialize_result_;
};

TEST_F(Interpolate1dFloatTest, InvalidType) {
	EXPECT_EQ(initialize_result_, sakura_Status_kOK);

	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	SIMD_ALIGN
	double const x_base[num_base] = { 0.0, 1.0 };
	SIMD_ALIGN
	float const y_base[num_base] = { 1.0, -1.0 };
	SIMD_ALIGN
	double const x_interpolated[num_interpolated] = { -1.0, 0.0, 0.1, 0.7, 1.5 };
	SIMD_ALIGN
	float y_interpolated[num_interpolated];

	// check alignment
	ASSERT_TRUE(sakura_IsAligned(x_base))<< "x_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_base))<< "y_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(x_interpolated))<< "x_interpolated is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_interpolated))<< "y_interpolated is not aligned";

	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNumMethod, 0, num_base, x_base, y_base,
			num_interpolated, x_interpolated, y_interpolated);

	// Should return InvalidArgument status
	EXPECT_EQ(result, sakura_Status_kInvalidArgument)
			<< "Interpolate1dFloat should fail!";
}

TEST_F(Interpolate1dFloatTest, Nearest) {
	EXPECT_EQ(initialize_result_, sakura_Status_kOK);

	size_t const num_base = 2;
	size_t const num_interpolated = 6;
	SIMD_ALIGN
	double const x_base[num_base] = { 0.0, 1.0 };
	SIMD_ALIGN
	float const y_base[num_base] = { 1.0, -1.0 };
	SIMD_ALIGN
	double const x_interpolated[num_interpolated] = { -1.0, 0.0, 0.1, 0.5, 0.7,
			1.5 };
	SIMD_ALIGN
	float y_interpolated[num_interpolated];

	// check alignment
	ASSERT_TRUE(sakura_IsAligned(x_base))<< "x_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_base))<< "y_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(x_interpolated))<< "x_interpolated is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_interpolated))<< "y_interpolated is not aligned";

	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNearest, 0, num_base, x_base, y_base,
			num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(result, sakura_Status_kOK)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		float reference;
		if (x_interpolated[index] <= 0.5 * (x_base[0] + x_base[1])) {
			reference = y_base[0];
		} else {
			reference = y_base[1];
		}
		std::cout << "Expected value at index " << index << ": " << reference
				<< std::endl;
		EXPECT_EQ(y_interpolated[index], reference)
				<< "interpolated value differs from expected value at " << index
				<< ": " << reference << ", " << y_interpolated[index];
	}
}

TEST_F(Interpolate1dFloatTest, NearestSingleBase) {
	EXPECT_EQ(initialize_result_, sakura_Status_kOK);

	size_t const num_base = 1;
	size_t const num_interpolated = 3;
	SIMD_ALIGN
	double const x_base[num_base] = { 0.0 };
	SIMD_ALIGN
	float const y_base[num_base] = { 1.0 };
	SIMD_ALIGN
	double const x_interpolated[num_interpolated] = { -1.0, 0.0, 0.1 };
	SIMD_ALIGN
	float y_interpolated[num_interpolated];

	// check alignment
	ASSERT_TRUE(sakura_IsAligned(x_base))<< "x_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_base))<< "y_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(x_interpolated))<< "x_interpolated is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_interpolated))<< "y_interpolated is not aligned";

	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNearest, 0, num_base, x_base, y_base,
			num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(result, sakura_Status_kOK)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		float reference = y_base[0];
		std::cout << "Expected value at index " << index << ": " << reference
				<< std::endl;
		EXPECT_EQ(y_interpolated[index], reference)
				<< "interpolated value differs from expected value at " << index
				<< ": " << reference << ", " << y_interpolated[index];
	}
}

TEST_F(Interpolate1dFloatTest, Linear) {
	EXPECT_EQ(initialize_result_, sakura_Status_kOK);

	size_t const num_base = 2;
	size_t const num_interpolated = 6;
	SIMD_ALIGN
	double const x_base[num_base] = { 0.0, 1.0 };
	SIMD_ALIGN
	float const y_base[num_base] = { 1.0, -1.0 };
	SIMD_ALIGN
	double const x_interpolated[num_interpolated] = { -1.0, 0.0, 0.1, 0.5, 0.7,
			1.5 };
	SIMD_ALIGN
	float y_interpolated[num_interpolated];
	SIMD_ALIGN
	float const y_expected[num_interpolated] =
			{ 1.0, 1.0, 0.8, 0.0, -0.4, -1.0 };

	// check alignment
	ASSERT_TRUE(sakura_IsAligned(x_base))<< "x_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_base))<< "y_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(x_interpolated))<< "x_interpolated is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_interpolated))<< "y_interpolated is not aligned";

	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kLinear, 0, num_base, x_base, y_base,
			num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(result, sakura_Status_kOK)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected[index] << std::endl;
		EXPECT_EQ(y_interpolated[index], y_expected[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected[index] << ", " << y_interpolated[index];
	}
}

TEST_F(Interpolate1dFloatTest, LinearSingleBase) {
	EXPECT_EQ(initialize_result_, sakura_Status_kOK);

	size_t const num_base = 1;
	size_t const num_interpolated = 3;
	SIMD_ALIGN
	double x_base[num_base] = { 0.0 };
	SIMD_ALIGN
	float y_base[num_base] = { 1.0 };
	SIMD_ALIGN
	double x_interpolated[num_interpolated] = { -1.0, 0.0, 0.1 };
	SIMD_ALIGN
	float y_interpolated[num_interpolated];

	// check alignment
	ASSERT_TRUE(sakura_IsAligned(x_base))<< "x_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_base))<< "y_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(x_interpolated))<< "x_interpolated is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_interpolated))<< "y_interpolated is not aligned";

	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kLinear, 0, num_base, x_base, y_base,
			num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(result, sakura_Status_kOK)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		float reference = y_base[0];
		std::cout << "Expected value at index " << index << ": " << reference
				<< std::endl;
		EXPECT_EQ(y_interpolated[index], reference)
				<< "interpolated value differs from expected value at " << index
				<< ": " << reference << ", " << y_interpolated[index];
	}
}
