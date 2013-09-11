#include <iostream>
#include <memory>
#include <cmath>

#include <libsakura/sakura.h>

#include "gtest/gtest.h"

#include "aligned_memory.h"

class Interpolate1dFloatTest: public ::testing::Test {
protected:
	virtual void SetUp() {
		initialize_result_ = sakura_Initialize();
		polynomial_order_ = 0;
		sakura_alignment_ = sakura_GetAlignment();
	}
	virtual void TearDown() {
		sakura_CleanUp();
	}
	virtual void AllocateMemory(size_t num_base, size_t num_interpolated) {
		size_t num_arena_base = num_base + sakura_alignment_ - 1;
		size_t num_arena_interpolated = num_interpolated + sakura_alignment_
				- 1;
		storage_for_x_base_.reset(new double[num_arena_base]);
		x_base_ = sakura_AlignDouble(num_arena_base, storage_for_x_base_.get(),
				num_base);
		storage_for_x_interpolated_.reset(new double[num_arena_interpolated]);
		x_interpolated_ = sakura_AlignDouble(num_arena_interpolated,
				storage_for_x_interpolated_.get(), num_interpolated);
		storage_for_y_base_.reset(new float[num_arena_base]);
		y_base_ = sakura_AlignFloat(num_arena_base, storage_for_y_base_.get(),
				num_base);
		storage_for_y_interpolated_.reset(new float[num_arena_interpolated]);
		y_interpolated_ = sakura_AlignFloat(num_arena_interpolated,
				storage_for_y_interpolated_.get(), num_interpolated);
		storage_for_y_expected_.reset(new float[num_arena_interpolated]);
		y_expected_ = sakura_AlignFloat(num_arena_interpolated,
				storage_for_y_expected_.get(), num_interpolated);

		// check alignment
		ASSERT_TRUE(x_base_ != nullptr)<< "x_base_ is null";
		ASSERT_TRUE(sakura_IsAligned(x_base_))<< "x_base_ is not aligned";
		ASSERT_TRUE(sakura_IsAligned(y_base_))<< "y_base_ is not aligned";
		ASSERT_TRUE(sakura_IsAligned(x_interpolated_))<< "x_interpolated_ is not aligned";
		ASSERT_TRUE(sakura_IsAligned(y_interpolated_))<< "y_interpolated_ is not aligned";
	}

	sakura_Status initialize_result_;
	size_t sakura_alignment_;
	int polynomial_order_;

	std::unique_ptr<double[]> storage_for_x_base_;
	std::unique_ptr<float[]> storage_for_y_base_;
	std::unique_ptr<double[]> storage_for_x_interpolated_;
	std::unique_ptr<float[]> storage_for_y_interpolated_;
	std::unique_ptr<float[]> storage_for_y_expected_;
	double *x_base_;
	double *x_interpolated_;
	float *y_base_;
	float *y_interpolated_;
	float *y_expected_;
};

TEST_F(Interpolate1dFloatTest, InvalidType) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 1.0;
	y_base_[0] = 1.0;
	y_base_[1] = -1.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.7;
	x_interpolated_[4] = 1.5;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNumMethod, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Should return InvalidArgument status
	EXPECT_EQ(sakura_Status_kInvalidArgument, result)
			<< "Interpolate1dFloat should fail!";
}

TEST_F(Interpolate1dFloatTest, ZeroLengthBaseArray) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	// initial setup
	size_t const num_base = 0;
	size_t const num_interpolated = 5;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNearest, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Should return InvalidArgument status
	EXPECT_EQ(sakura_Status_kInvalidArgument, result)
			<< "Interpolate1dFloat should fail!";
}

TEST_F(Interpolate1dFloatTest, NegativePolynomialOrder) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	polynomial_order_ = -1;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Should return InvalidArgument status
	EXPECT_EQ(sakura_Status_kInvalidArgument, result)
			<< "Interpolate1dFloat should fail!";
}
TEST_F(Interpolate1dFloatTest, NegativePolynomialOrderButNearest) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 2;
	size_t const num_interpolated = 6;
	polynomial_order_ = -1;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 1.0;
	y_base_[0] = 1.0;
	y_base_[1] = -1.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.5;
	x_interpolated_[4] = 0.7;
	x_interpolated_[5] = 1.5;
	y_expected_[0] = 1.0;
	y_expected_[1] = 1.0;
	y_expected_[2] = 1.0;
	y_expected_[3] = 1.0;
	y_expected_[4] = -1.0;
	y_expected_[5] = -1.0;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNearest, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, Nearest) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 2;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 1.0;
	y_base_[0] = 1.0;
	y_base_[1] = -1.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.5;
	x_interpolated_[4] = 0.7;
	x_interpolated_[5] = 1.5;
	y_expected_[0] = 1.0;
	y_expected_[1] = 1.0;
	y_expected_[2] = 1.0;
	y_expected_[3] = 1.0;
	y_expected_[4] = -1.0;
	y_expected_[5] = -1.0;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNearest, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, NearestDescending) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 2;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[1] = 0.0;
	x_base_[0] = 1.0;
	y_base_[1] = 1.0;
	y_base_[0] = -1.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.5;
	x_interpolated_[4] = 0.7;
	x_interpolated_[5] = 1.5;
	y_expected_[0] = 1.0;
	y_expected_[1] = 1.0;
	y_expected_[2] = 1.0;
	y_expected_[3] = 1.0;
	y_expected_[4] = -1.0;
	y_expected_[5] = -1.0;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNearest, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, NearestOpposite) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 2;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[1] = 0.0;
	x_base_[0] = 1.0;
	y_base_[1] = 1.0;
	y_base_[0] = -1.0;
	x_interpolated_[5] = -1.0;
	x_interpolated_[4] = 0.0;
	x_interpolated_[3] = 0.1;
	x_interpolated_[2] = 0.5;
	x_interpolated_[1] = 0.7;
	x_interpolated_[0] = 1.5;
	y_expected_[5] = 1.0;
	y_expected_[4] = 1.0;
	y_expected_[3] = 1.0;
	y_expected_[2] = 1.0;
	y_expected_[1] = -1.0;
	y_expected_[0] = -1.0;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNearest, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

//TEST_F(Interpolate1dFloatTest, NearestSingleBase) {
TEST_F(Interpolate1dFloatTest, SingleBase) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 1;
	size_t const num_interpolated = 3;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	y_base_[0] = 1.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNearest, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		float reference = y_base_[0];
		std::cout << "Expected value at index " << index << ": " << reference
				<< std::endl;
		EXPECT_EQ(reference, y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << reference << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, Linear) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 2;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 1.0;
	y_base_[0] = 1.0;
	y_base_[1] = -1.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.5;
	x_interpolated_[4] = 0.7;
	x_interpolated_[5] = 1.5;
	y_expected_[0] = 1.0;
	y_expected_[1] = 1.0;
	y_expected_[2] = 0.8;
	y_expected_[3] = 0.0;
	y_expected_[4] = -0.4;
	y_expected_[5] = -1.0;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kLinear, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, LinearDescending) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 2;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[1] = 0.0;
	x_base_[0] = 1.0;
	y_base_[1] = 1.0;
	y_base_[0] = -1.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.5;
	x_interpolated_[4] = 0.7;
	x_interpolated_[5] = 1.5;
	y_expected_[0] = 1.0;
	y_expected_[1] = 1.0;
	y_expected_[2] = 0.8;
	y_expected_[3] = 0.0;
	y_expected_[4] = -0.4;
	y_expected_[5] = -1.0;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kLinear, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, LinearOpposite) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 2;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[1] = 0.0;
	x_base_[0] = 1.0;
	y_base_[1] = 1.0;
	y_base_[0] = -1.0;
	x_interpolated_[5] = -1.0;
	x_interpolated_[4] = 0.0;
	x_interpolated_[3] = 0.1;
	x_interpolated_[2] = 0.5;
	x_interpolated_[1] = 0.7;
	x_interpolated_[0] = 1.5;
	y_expected_[5] = 1.0;
	y_expected_[4] = 1.0;
	y_expected_[3] = 0.8;
	y_expected_[2] = 0.0;
	y_expected_[1] = -0.4;
	y_expected_[0] = -1.0;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kLinear, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

//TEST_F(Interpolate1dFloatTest, LinearSingleBase) {
//	EXPECT_EQ(sakura_Status_kOK, initialize_result_);
//
//	size_t const num_base = 1;
//	size_t const num_interpolated = 3;
//
//	// initial setup
//	AllocateMemory(num_base, num_interpolated);
//	x_base_[0] = 0.0;
//	y_base_[0] = 1.0;
//	x_interpolated_[0] = -1.0;
//	x_interpolated_[1] = 0.0;
//	x_interpolated_[2] = 0.1;
//
//	// execute interpolation
//	sakura_Status result = sakura_Interpolate1dFloat(
//			sakura_InterpolationMethod_kLinear, polynomial_order_, num_base,
//			x_base_, y_base_, num_interpolated, x_interpolated_,
//			y_interpolated_);
//
//	// Basic check whether function is completed or not
//	EXPECT_EQ(sakura_Status_kOK, result)
//			<< "Interpolate1dFloat had any problems during execution.";
//
//	// Value check
//	for (size_t index = 0; index < num_interpolated; ++index) {
//		float reference = y_base_[0];
//		std::cout << "Expected value at index " << index << ": " << reference
//				<< std::endl;
//		EXPECT_EQ(reference, y_interpolated_[index])
//				<< "interpolated value differs from expected value at " << index
//				<< ": " << reference << ", " << y_interpolated_[index];
//	}
//}

TEST_F(Interpolate1dFloatTest, PolynomialOrder0) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 0;

	size_t const num_base = 2;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 1.0;
	y_base_[0] = 1.0;
	y_base_[1] = -1.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.5;
	x_interpolated_[4] = 0.7;
	x_interpolated_[5] = 1.5;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 0-th order polynomial interpolation acts like NearestInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		float reference;
		if (x_interpolated_[index] <= 0.5 * (x_base_[0] + x_base_[1])) {
			reference = y_base_[0];
		} else {
			reference = y_base_[1];
		}
		std::cout << "Expected value at index " << index << ": " << reference
				<< std::endl;
		EXPECT_EQ(reference, y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << reference << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, PolynomialOrder1) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 1;

	size_t const num_base = 2;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 1.0;
	y_base_[0] = 1.0;
	y_base_[1] = -1.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.5;
	x_interpolated_[4] = 0.7;
	x_interpolated_[5] = 1.5;
	y_expected_[0] = 1.0;
	y_expected_[1] = 1.0;
	y_expected_[2] = 0.8;
	y_expected_[3] = 0.0;
	y_expected_[4] = -0.4;
	y_expected_[5] = -1.0;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 1-st order polynomial interpolation acts like LinearInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, PolynomialOrder2Full) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 2;

	size_t const num_base = 3;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 1.0;
	x_base_[2] = 2.0;
	y_base_[0] = 1.0;
	y_base_[1] = -1.0;
	y_base_[2] = 0.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.5;
	x_interpolated_[4] = 0.7;
	x_interpolated_[5] = 1.5;

	// expected value can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
	y_expected_[0] = 1.0; // out of range
	for (size_t i = 1; i < num_interpolated; ++i) {
		y_expected_[i] = (1.5 * x_interpolated_[i] - 3.5) * x_interpolated_[i]
				+ 1.0;
	}

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 1-st order polynomial interpolation acts like LinearInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, PolynomialOrder2FullDescending) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 2;

	size_t const num_base = 3;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[2] = 0.0;
	x_base_[1] = 1.0;
	x_base_[0] = 2.0;
	y_base_[2] = 1.0;
	y_base_[1] = -1.0;
	y_base_[0] = 0.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.5;
	x_interpolated_[4] = 0.7;
	x_interpolated_[5] = 1.5;

	// expected value can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
	y_expected_[0] = 1.0; // out of range
	for (size_t i = 1; i < num_interpolated; ++i) {
		y_expected_[i] = (1.5 * x_interpolated_[i] - 3.5) * x_interpolated_[i]
				+ 1.0;
	}

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 1-st order polynomial interpolation acts like LinearInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, PolynomialOrder2FullOpposite) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 2;

	size_t const num_base = 3;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[2] = 0.0;
	x_base_[1] = 1.0;
	x_base_[0] = 2.0;
	y_base_[2] = 1.0;
	y_base_[1] = -1.0;
	y_base_[0] = 0.0;
	x_interpolated_[5] = -1.0;
	x_interpolated_[4] = 0.0;
	x_interpolated_[3] = 0.1;
	x_interpolated_[2] = 0.5;
	x_interpolated_[1] = 0.7;
	x_interpolated_[0] = 1.5;

	// expected value can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
	y_expected_[5] = 1.0; // out of range
	for (size_t i = 0; i < num_interpolated - 1; ++i) {
		y_expected_[i] = (1.5 * x_interpolated_[i] - 3.5) * x_interpolated_[i]
				+ 1.0;
	}

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 1-st order polynomial interpolation acts like LinearInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, PolynomialOrder1Sub) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 1;

	size_t const num_base = 3;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 1.0;
	x_base_[2] = 2.0;
	y_base_[0] = 1.0;
	y_base_[1] = -1.0;
	y_base_[2] = 0.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.5;
	x_interpolated_[4] = 0.7;
	x_interpolated_[5] = 1.5;
	y_expected_[0] = 1.0;
	y_expected_[1] = 1.0;
	y_expected_[2] = 0.8;
	y_expected_[3] = 0.0;
	y_expected_[4] = -0.4;
	y_expected_[5] = -0.5;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 1-st order polynomial interpolation acts like LinearInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

//TEST_F(Interpolate1dFloatTest, PolynomialSingleBase) {
//	EXPECT_EQ(sakura_Status_kOK, initialize_result_);
//
//	polynomial_order_ = 1;
//
//	size_t const num_base = 1;
//	size_t const num_interpolated = 3;
//
//	// initial setup
//	AllocateMemory(num_base, num_interpolated);
//	x_base_[0] = 0.0;
//	y_base_[0] = 1.0;
//	x_interpolated_[0] = -1.0;
//	x_interpolated_[1] = 0.0;
//	x_interpolated_[2] = 0.1;
//
//	// execute interpolation
//	sakura_Status result = sakura_Interpolate1dFloat(
//			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
//			x_base_, y_base_, num_interpolated, x_interpolated_,
//			y_interpolated_);
//
//	// Basic check whether function is completed or not
//	EXPECT_EQ(sakura_Status_kOK, result)
//			<< "Interpolate1dFloat had any problems during execution.";
//
//	// Value check
//	for (size_t index = 0; index < num_interpolated; ++index) {
//		float reference = y_base_[0];
//		std::cout << "Expected value at index " << index << ": " << reference
//				<< std::endl;
//		EXPECT_EQ(reference, y_interpolated_[index])
//				<< "interpolated value differs from expected value at " << index
//				<< ": " << reference << ", " << y_interpolated_[index];
//	}
//}

TEST_F(Interpolate1dFloatTest, Spline) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 3;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 1.0;
	x_base_[2] = 2.0;
	y_base_[0] = 1.0;
	y_base_[1] = -1.0;
	y_base_[2] = 0.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.5;
	x_interpolated_[4] = 0.7;
	x_interpolated_[5] = 1.5;
	y_expected_[0] = 1.0;
	y_expected_[1] = 1.0;
	y_expected_[2] = 0.72575;
	y_expected_[3] = -0.28125;
	y_expected_[4] = -0.66775;
	y_expected_[5] = -0.78125;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kSpline, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 1-st order polynomial interpolation acts like LinearInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, SplineDescending) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 3;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[2] = 0.0;
	x_base_[1] = 1.0;
	x_base_[0] = 2.0;
	y_base_[2] = 1.0;
	y_base_[1] = -1.0;
	y_base_[0] = 0.0;
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.0;
	x_interpolated_[2] = 0.1;
	x_interpolated_[3] = 0.5;
	x_interpolated_[4] = 0.7;
	x_interpolated_[5] = 1.5;
	y_expected_[0] = 1.0;
	y_expected_[1] = 1.0;
	y_expected_[2] = 0.72575;
	y_expected_[3] = -0.28125;
	y_expected_[4] = -0.66775;
	y_expected_[5] = -0.78125;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kSpline, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 1-st order polynomial interpolation acts like LinearInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

TEST_F(Interpolate1dFloatTest, SplineOpposite) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 3;
	size_t const num_interpolated = 6;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[2] = 0.0;
	x_base_[1] = 1.0;
	x_base_[0] = 2.0;
	y_base_[2] = 1.0;
	y_base_[1] = -1.0;
	y_base_[0] = 0.0;
	x_interpolated_[5] = -1.0;
	x_interpolated_[4] = 0.0;
	x_interpolated_[3] = 0.1;
	x_interpolated_[2] = 0.5;
	x_interpolated_[1] = 0.7;
	x_interpolated_[0] = 1.5;
	y_expected_[5] = 1.0;
	y_expected_[4] = 1.0;
	y_expected_[3] = 0.72575;
	y_expected_[2] = -0.28125;
	y_expected_[1] = -0.66775;
	y_expected_[0] = -0.78125;

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kSpline, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";

	// Value check
	// 1-st order polynomial interpolation acts like LinearInterpolation
	for (size_t index = 0; index < num_interpolated; ++index) {
		std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
		EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
	}
}

//TEST_F(Interpolate1dFloatTest, SplineSingleBase) {
//	EXPECT_EQ(sakura_Status_kOK, initialize_result_);
//
//	size_t const num_base = 1;
//	size_t const num_interpolated = 3;
//
//	// initial setup
//	AllocateMemory(num_base, num_interpolated);
//	x_base_[0] = 0.0;
//	y_base_[0] = 1.0;
//	x_interpolated_[0] = -1.0;
//	x_interpolated_[1] = 0.0;
//	x_interpolated_[2] = 0.1;
//
//	// execute interpolation
//	sakura_Status result = sakura_Interpolate1dFloat(
//			sakura_InterpolationMethod_kSpline, polynomial_order_, num_base,
//			x_base_, y_base_, num_interpolated, x_interpolated_,
//			y_interpolated_);
//
//	// Basic check whether function is completed or not
//	EXPECT_EQ(sakura_Status_kOK, result)
//			<< "Interpolate1dFloat had any problems during execution.";
//
//	// Value check
//	for (size_t index = 0; index < num_interpolated; ++index) {
//		float reference = y_base_[0];
//		std::cout << "Expected value at index " << index << ": " << reference
//				<< std::endl;
//		EXPECT_EQ(reference, y_interpolated_[index])
//				<< "interpolated value differs from expected value at " << index
//				<< ": " << reference << ", " << y_interpolated_[index];
//	}
//}

TEST_F(Interpolate1dFloatTest, SingleBasePerformance) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 1;
	size_t const num_interpolated = 200000000;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	y_base_[0] = 1.0;
	double dx = 1.0 / static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = -0.5 + dx * static_cast<double>(i);
	}

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNearest, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";
}

TEST_F(Interpolate1dFloatTest, NearestPerformance) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 2;
	size_t const num_interpolated = 200000000;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 100.0;
	y_base_[0] = 1.0;
	y_base_[1] = 0.0;
	double dx = fabs(x_base_[1] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
	}

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNearest, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";
}

TEST_F(Interpolate1dFloatTest, NearestOppositePerformance) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 2;
	size_t const num_interpolated = 200000000;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[1] = 0.0;
	x_base_[0] = 100.0;
	y_base_[1] = 1.0;
	y_base_[0] = 0.0;
	double dx = fabs(x_base_[1] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] - dx * static_cast<double>(i);
	}

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kNearest, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";
}

TEST_F(Interpolate1dFloatTest, LinearPerformance) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 2;
	size_t const num_interpolated = 200000000;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 100.0;
	y_base_[0] = 1.0;
	y_base_[1] = 0.0;
	double dx = fabs(x_base_[1] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
	}

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kLinear, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";
}

TEST_F(Interpolate1dFloatTest, LinearOppositePerformance) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 2;
	size_t const num_interpolated = 200000000;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[1] = 0.0;
	x_base_[0] = 100.0;
	y_base_[1] = 1.0;
	y_base_[0] = 0.0;
	double dx = fabs(x_base_[1] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] - dx * static_cast<double>(i);
	}

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kLinear, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";
}

TEST_F(Interpolate1dFloatTest, PolynomialOrder2FullPerformance) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 2;

	size_t const num_base = 3;
	size_t const num_interpolated = 20000000; // 1/10 of Nearest and Linear

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 100.0;
	x_base_[2] = 200.0;
	y_base_[0] = 1.0;
	y_base_[1] = -1.0;
	y_base_[2] = 0.0;
	double dx = fabs(x_base_[2] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
	}

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";
}

TEST_F(Interpolate1dFloatTest, PolynomialOrder2FullOppositePerformance) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	polynomial_order_ = 2;

	size_t const num_base = 3;
	size_t const num_interpolated = 20000000; // 1/10 of Nearest and Linear

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[2] = 0.0;
	x_base_[1] = 100.0;
	x_base_[0] = 200.0;
	y_base_[2] = 1.0;
	y_base_[1] = -1.0;
	y_base_[0] = 0.0;
	double dx = fabs(x_base_[2] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] - dx * static_cast<double>(i);
	}

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";
}

TEST_F(Interpolate1dFloatTest, SplinePerformance) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 3;
	size_t const num_interpolated = 200000000;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[0] = 0.0;
	x_base_[1] = 100.0;
	x_base_[2] = 200.0;
	y_base_[0] = 1.0;
	y_base_[1] = -1.0;
	y_base_[2] = 0.0;
	double dx = fabs(x_base_[2] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
	}

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kSpline, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";
}

TEST_F(Interpolate1dFloatTest, SplineOppositePerformance) {
	EXPECT_EQ(sakura_Status_kOK, initialize_result_);

	size_t const num_base = 3;
	size_t const num_interpolated = 200000000;

	// initial setup
	AllocateMemory(num_base, num_interpolated);
	x_base_[2] = 0.0;
	x_base_[1] = 100.0;
	x_base_[0] = 200.0;
	y_base_[2] = 1.0;
	y_base_[1] = -1.0;
	y_base_[0] = 0.0;
	double dx = fabs(x_base_[2] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] - dx * static_cast<double>(i);
	}

	// execute interpolation
	sakura_Status result = sakura_Interpolate1dFloat(
			sakura_InterpolationMethod_kSpline, polynomial_order_, num_base,
			x_base_, y_base_, num_interpolated, x_interpolated_,
			y_interpolated_);

	// Basic check whether function is completed or not
	EXPECT_EQ(sakura_Status_kOK, result)
			<< "Interpolate1dFloat had any problems during execution.";
}
