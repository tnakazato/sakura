#include <iostream>

#include <libsakura/sakura.h>

#include "gtest/gtest.h"

#include "aligned_memory.h"

class Interpolate1dFloatTest: public ::testing::Test {
protected:
	virtual void SetUp() {
		initialize_result_ = sakura_Initialize();
		polynomial_order_ = 0;
	}
	virtual void TearDown() {
		sakura_CleanUp();
	}
	sakura_Status initialize_result_;
	int polynomial_order_;
};

TEST_F(Interpolate1dFloatTest, InvalidType) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

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
			sakura_InterpolationMethod_kNumMethod, polynomial_order_, num_base,
			x_base, y_base, num_interpolated, x_interpolated, y_interpolated);

	// Should return InvalidArgument status
	EXPECT_EQ(sakura_Status_kInvalidArgument, result)
			<< "Interpolate1dFloat should fail!";
}

TEST_F(Interpolate1dFloatTest, Nearest) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

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
			sakura_InterpolationMethod_kNearest, polynomial_order_, num_base,
			x_base, y_base, num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
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
		EXPECT_EQ(reference, y_interpolated[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << reference << ", " << y_interpolated[index];
	}
}

TEST_F(Interpolate1dFloatTest, NearestSingleBase) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

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
			sakura_InterpolationMethod_kNearest, polynomial_order_, num_base,
			x_base, y_base, num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		float reference = y_base[0];
		std::cout << "Expected value at index " << index << ": " << reference
				<< std::endl;
		EXPECT_EQ(reference, y_interpolated[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << reference << ", " << y_interpolated[index];
	}
}

TEST_F(Interpolate1dFloatTest, Linear) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

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
	ASSERT_TRUE(sakura_IsAligned(y_expected))<< "y_expected is not aligned";

	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kLinear, polynomial_order_, num_base,
			x_base, y_base, num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected[index], y_interpolated[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected[index] << ", " << y_interpolated[index];
	}
}

TEST_F(Interpolate1dFloatTest, LinearSingleBase) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

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
			sakura_InterpolationMethod_kLinear, polynomial_order_, num_base,
			x_base, y_base, num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		float reference = y_base[0];
		std::cout << "Expected value at index " << index << ": " << reference
				<< std::endl;
		EXPECT_EQ(reference, y_interpolated[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << reference << ", " << y_interpolated[index];
	}
}

TEST_F(Interpolate1dFloatTest, PolynomialOrder0) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 0;

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
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base, y_base, num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 0-th order polynomial interpolation acts like NearestInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		float reference;
		if (x_interpolated[index] <= 0.5 * (x_base[0] + x_base[1])) {
			reference = y_base[0];
		} else {
			reference = y_base[1];
		}
		std::cout << "Expected value at index " << index << ": " << reference
				<< std::endl;
		EXPECT_EQ(reference, y_interpolated[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << reference << ", " << y_interpolated[index];
	}
}

TEST_F(Interpolate1dFloatTest, PolynomialOrder1) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 1;

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
	ASSERT_TRUE(sakura_IsAligned(y_expected))<< "y_expected is not aligned";

	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base, y_base, num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 1-st order polynomial interpolation acts like LinearInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected[index], y_interpolated[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected[index] << ", " << y_interpolated[index];
	}
}

TEST_F(Interpolate1dFloatTest, PolynomialOrder2Full) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 2;

	size_t const num_base = 3;
	size_t const num_interpolated = 6;
	SIMD_ALIGN
	double const x_base[num_base] = { 0.0, 1.0, 2.0 };
	SIMD_ALIGN
	float const y_base[num_base] = { 1.0, -1.0, 0.0 };
	SIMD_ALIGN
	double const x_interpolated[num_interpolated] = { -1.0, 0.0, 0.1, 0.5, 0.7,
			1.5 };
	SIMD_ALIGN
	float y_interpolated[num_interpolated];
	SIMD_ALIGN
	float y_expected[num_interpolated];

	// expected value can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
	y_expected[0] = 1.0; // out of range
	for (size_t i = 1; i < num_interpolated; ++i) {
		y_expected[i] = (1.5 * x_interpolated[i] - 3.5) * x_interpolated[i]
				+ 1.0;
	}

	// check alignment
	ASSERT_TRUE(sakura_IsAligned(x_base))<< "x_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_base))<< "y_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(x_interpolated))<< "x_interpolated is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_interpolated))<< "y_interpolated is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_expected))<< "y_expected is not aligned";

	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base, y_base, num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 1-st order polynomial interpolation acts like LinearInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected[index], y_interpolated[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected[index] << ", " << y_interpolated[index];
	}
}

TEST_F(Interpolate1dFloatTest, PolynomialOrder1Sub) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 1;

	size_t const num_base = 3;
	size_t const num_interpolated = 6;
	SIMD_ALIGN
	double const x_base[num_base] = { 0.0, 1.0, 2.0 };
	SIMD_ALIGN
	float const y_base[num_base] = { 1.0, -1.0, 0.0 };
	SIMD_ALIGN
	double const x_interpolated[num_interpolated] = { -1.0, 0.0, 0.1, 0.5, 0.7,
			1.5 };
	SIMD_ALIGN
	float y_interpolated[num_interpolated];
	SIMD_ALIGN
	float const y_expected[num_interpolated] =
			{ 1.0, 1.0, 0.8, 0.0, -0.4, -0.5 };

	// check alignment
	ASSERT_TRUE(sakura_IsAligned(x_base))<< "x_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_base))<< "y_base is not aligned";
	ASSERT_TRUE(sakura_IsAligned(x_interpolated))<< "x_interpolated is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_interpolated))<< "y_interpolated is not aligned";
	ASSERT_TRUE(sakura_IsAligned(y_expected))<< "y_expected is not aligned";

	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base, y_base, num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 1-st order polynomial interpolation acts like LinearInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected[index], y_interpolated[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected[index] << ", " << y_interpolated[index];
	}
}

TEST_F(Interpolate1dFloatTest, PolynomialSingleBase) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 1;

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
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base, y_base, num_interpolated, x_interpolated, y_interpolated);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		float reference = y_base[0];
		std::cout << "Expected value at index " << index << ": " << reference
				<< std::endl;
		EXPECT_EQ(reference, y_interpolated[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << reference << ", " << y_interpolated[index];
	}
}
