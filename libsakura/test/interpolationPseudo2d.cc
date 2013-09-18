#include <iostream>
#include <memory>
#include <cmath>
#include <string>
#include <stdarg.h>

#include <libsakura/sakura.h>

#include "gtest/gtest.h"
#include "interpolation.h"

//#include "CubicSplineInterpolator1D.h"

class InterpolatePseudo2dFloatTest: public ::testing::Test {
protected:
	virtual void SetUp() {
		initialize_result_ = sakura_Initialize(nullptr, nullptr);
		polynomial_order_ = 0;
		sakura_alignment_ = sakura_GetAlignment();
	}
	virtual void TearDown() {
		sakura_CleanUp();
	}
	virtual void AllocateMemory(size_t num_base, size_t num_interpolated) {
		size_t num_arena_xbase = num_base + sakura_alignment_ - 1;
		size_t num_arena_ybase = num_base * num_interpolated
				+ sakura_alignment_;
		size_t num_arena_interpolated = num_interpolated + sakura_alignment_
				- 1;
		storage_for_x_base_.reset(new double[num_arena_xbase]);
		x_base_ = sakura_AlignDouble(num_arena_xbase, storage_for_x_base_.get(),
				num_base);
		storage_for_y_base_.reset(new float[num_arena_ybase]);
		y_base_ = sakura_AlignFloat(num_arena_ybase, storage_for_y_base_.get(),
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
		ASSERT_TRUE(sakura_IsAligned(y_interpolated_))<< "y_interpolated_ is not aligned";
	}

	virtual void RunInterpolatePseudo2d(sakura_InterpolationMethod interpolation_method,
			double x_interpolated, size_t num_base, size_t num_interpolated,
			sakura_Status expected_status, bool check_result) {
		// sakura must be properly initialized
		ASSERT_EQ(sakura_Status_kOK, initialize_result_)
		<< "sakura must be properly initialized!";

		// execute interpolation
		sakura_Status result = sakura_InterpolatePseudo2dFloat(
				interpolation_method, polynomial_order_, x_interpolated, num_base,
				x_base_, num_interpolated, y_base_, y_interpolated_);

		// Should return InvalidArgument status
		std::string message = (expected_status == sakura_Status_kOK) ?
		"InterpolatePseudo2dFloat had any problems during execution." :
		"InterpolatePseudo2dFloat should fail!";
		EXPECT_EQ(expected_status, result) << message;

		if (check_result && (expected_status == result)) {
			// Value check
			for (size_t index = 0; index < num_interpolated; ++index) {
				std::cout << "Expected value at index " << index << ": "
				<< y_expected_[index] << std::endl;
				EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
				<< "interpolated value differs from expected value at " << index
				<< ": " << y_expected_[index] << ", " << y_interpolated_[index];
			}
		}
	}

	sakura_Status initialize_result_;
	size_t sakura_alignment_;
	int polynomial_order_;

	std::unique_ptr<double[]> storage_for_x_base_;
	std::unique_ptr<float[]> storage_for_y_base_;
	std::unique_ptr<float[]> storage_for_y_interpolated_;
	std::unique_ptr<float[]> storage_for_y_expected_;
	double *x_base_;
	float *y_base_;
	float *y_interpolated_;
	float *y_expected_;
};

TEST_F(InterpolatePseudo2dFloatTest, InvalidType) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	AllocateMemory(num_base, num_interpolated);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base, y_base_, 1.0, -1.0);

	// execute interpolation
	// Should return InvalidArgument status
	RunInterpolatePseudo2d(sakura_InterpolationMethod_kNumMethod, 0.0, num_base,
			num_interpolated, sakura_Status_kInvalidArgument, false);
}

TEST_F(InterpolatePseudo2dFloatTest, ZeroLengthBaseArray) {
	// initial setup
	size_t const num_base = 0;
	size_t const num_interpolated = 5;

	// execute interpolation
	// Should return InvalidArgument status
	RunInterpolatePseudo2d(sakura_InterpolationMethod_kNearest, 0.0, num_base,
			num_interpolated, sakura_Status_kInvalidArgument, false);
}

TEST_F(InterpolatePseudo2dFloatTest, NegativePolynomialOrder) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	polynomial_order_ = -1;

	// execute interpolation
	// Should return InvalidArgument status
	RunInterpolatePseudo2d(sakura_InterpolationMethod_kPolynomial, 0.0, num_base,
			num_interpolated, sakura_Status_kInvalidArgument, false);
}

//TEST_F(InterpolatePseudo2dFloatTest, NegativePolynomialOrderButNearest) {
//	// initial setup
//	size_t const num_base = 2;
//	size_t const num_interpolated = 6;
//	polynomial_order_ = -1;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
//	InitializeFloatArray(num_base, y_base_, 1.0, -1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
//			0.5, 0.7, 1.5);
//	InitializeFloatArray(num_interpolated, y_expected_, 1.0, 1.0, 1.0, 1.0,
//			-1.0, -1.0);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kNearest, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, InputArrayNotAligned) {
//	// initial setup
//	size_t const num_base = 2;
//	size_t const num_interpolated = 5;
//	AllocateMemory(num_base, num_interpolated);
//
//	// execute interpolation
//	x_base_ = &(x_base_[1]);
//	RunInterpolate1d(sakura_InterpolationMethod_kPolynomial, 1,
//			num_interpolated, sakura_Status_kInvalidArgument, false);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, Nearest) {
//	// initial setup
//	size_t const num_base = 2;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
//	InitializeFloatArray(num_base, y_base_, 1.0, -1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
//			0.5, 0.7, 1.5);
//	InitializeFloatArray(num_interpolated, y_expected_, 1.0, 1.0, 1.0, 1.0,
//			-1.0, -1.0);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kNearest, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, NearestDescending) {
//	// initial setup
//	size_t const num_base = 2;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 1.0, 0.0);
//	InitializeFloatArray(num_base, y_base_, -1.0, 1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
//			0.5, 0.7, 1.5);
//	InitializeFloatArray(num_interpolated, y_expected_, 1.0, 1.0, 1.0, 1.0,
//			-1.0, -1.0);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kNearest, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, NearestOpposite) {
//	// initial setup
//	size_t const num_base = 2;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 1.0, 0.0);
//	InitializeFloatArray(num_base, y_base_, -1.0, 1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, 1.5, 0.7, 0.5, 0.1,
//			0.0, -1.0);
//	InitializeFloatArray(num_interpolated, y_expected_, -1.0, -1.0, 1.0, 1.0,
//			1.0, 1.0);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kNearest, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, SingleBase) {
//	// initial setup
//	size_t const num_base = 1;
//	size_t const num_interpolated = 3;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0);
//	InitializeFloatArray(num_base, y_base_, 1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1);
//	InitializeFloatArray(num_interpolated, y_expected_, 1.0, 1.0, 1.0);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kNearest, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, Linear) {
//	// initial setup
//	size_t const num_base = 2;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
//	InitializeFloatArray(num_base, y_base_, 1.0, -1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
//			0.5, 0.7, 1.5);
//	InitializeFloatArray(num_interpolated, y_expected_, 1.0, 1.0, 0.8, 0.0,
//			-0.4, -1.0);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kLinear, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, LinearDescending) {
//	// initial setup
//	size_t const num_base = 2;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 1.0, 0.0);
//	InitializeFloatArray(num_base, y_base_, -1.0, 1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
//			0.5, 0.7, 1.5);
//	InitializeFloatArray(num_interpolated, y_expected_, 1.0, 1.0, 0.8, 0.0,
//			-0.4, -1.0);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kLinear, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, LinearOpposite) {
//	// initial setup
//	size_t const num_base = 2;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 1.0, 0.0);
//	InitializeFloatArray(num_base, y_base_, -1.0, 1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, 1.5, 0.7, 0.5, 0.1,
//			0.0, -1.0);
//	InitializeFloatArray(num_interpolated, y_expected_, -1.0, -0.4, 0.0, 0.8,
//			1.0, 1.0);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kLinear, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, PolynomialOrder0) {
//	// initial setup
//	polynomial_order_ = 0;
//	size_t const num_base = 2;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
//	InitializeFloatArray(num_base, y_base_, 1.0, -1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
//			0.5, 0.7, 1.5);
//	InitializeFloatArray(num_interpolated, y_expected_, 1.0, 1.0, 1.0, 1.0,
//			-1.0, -1.0);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kPolynomial, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, PolynomialOrder1) {
//	// initial setup
//	polynomial_order_ = 1;
//	size_t const num_base = 2;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
//	InitializeFloatArray(num_base, y_base_, 1.0, -1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
//			0.5, 0.7, 1.5);
//	InitializeFloatArray(num_interpolated, y_expected_, 1.0, 1.0, 0.8, 0.0,
//			-0.4, -1.0);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kPolynomial, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, PolynomialOrder2Full) {
//	// initial setup
//	polynomial_order_ = 2;
//	size_t const num_base = 3;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
//	InitializeFloatArray(num_base, y_base_, 1.0, -1.0, 0.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
//			0.5, 0.7, 1.5);
//
//	// expected value can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
//	y_expected_[0] = 1.0; // out of range
//	for (size_t i = 1; i < num_interpolated; ++i) {
//		y_expected_[i] = (1.5 * x_interpolated_[i] - 3.5) * x_interpolated_[i]
//				+ 1.0;
//	}
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kPolynomial, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, PolynomialOrder2FullDescending) {
//	// initial setup
//	polynomial_order_ = 2;
//	size_t const num_base = 3;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
//	InitializeFloatArray(num_base, y_base_, 0.0, -1.0, 1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
//			0.5, 0.7, 1.5);
//
//	// expected value can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
//	y_expected_[0] = 1.0; // out of range
//	for (size_t i = 1; i < num_interpolated; ++i) {
//		y_expected_[i] = (1.5 * x_interpolated_[i] - 3.5) * x_interpolated_[i]
//				+ 1.0;
//	}
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kPolynomial, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, PolynomialOrder2FullOpposite) {
//	// initial setup
//	polynomial_order_ = 2;
//	size_t const num_base = 3;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
//	InitializeFloatArray(num_base, y_base_, 0.0, -1.0, 1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, 1.5, 0.7, 0.5, 0.1,
//			0.0, -1.0);
//
//	// expected value can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
//	y_expected_[5] = 1.0; // out of range
//	for (size_t i = 0; i < num_interpolated - 1; ++i) {
//		y_expected_[i] = (1.5 * x_interpolated_[i] - 3.5) * x_interpolated_[i]
//				+ 1.0;
//	}
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kPolynomial, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, PolynomialOrder1Sub) {
//	// initial setup
//	polynomial_order_ = 1;
//	size_t const num_base = 3;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
//	InitializeFloatArray(num_base, y_base_, 1.0, -1.0, 0.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
//			0.5, 0.7, 1.5);
//	InitializeFloatArray(num_interpolated, y_expected_, 1.0, 1.0, 0.8, 0.0,
//			-0.4, -0.5);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kPolynomial, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, Spline) {
//	// initial setup
//	size_t const num_base = 3;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
//	InitializeFloatArray(num_base, y_base_, 1.0, -1.0, 0.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
//			0.5, 0.7, 1.5);
//	InitializeFloatArray(num_interpolated, y_expected_, 1.0, 1.0, 0.72575,
//			-0.28125, -0.66775, -0.78125);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kSpline, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, SplineDescending) {
//	// initial setup
//	size_t const num_base = 3;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
//	InitializeFloatArray(num_base, y_base_, 0.0, -1.0, 1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
//			0.5, 0.7, 1.5);
//	InitializeFloatArray(num_interpolated, y_expected_, 1.0, 1.0, 0.72575,
//			-0.28125, -0.66775, -0.78125);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kSpline, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, SplineOpposite) {
//	// initial setup
//	size_t const num_base = 3;
//	size_t const num_interpolated = 6;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
//	InitializeFloatArray(num_base, y_base_, 0.0, -1.0, 1.0);
//	InitializeDoubleArray(num_interpolated, x_interpolated_, 1.5, 0.7, 0.5, 0.1,
//			0.0, -1.0);
//	InitializeFloatArray(num_interpolated, y_expected_, -0.78125, -0.66775,
//			-0.28125, 0.72575, 1.0, 1.0);
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kSpline, num_base,
//			num_interpolated, sakura_Status_kOK, true);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, SingleBasePerformance) {
//	// initial setup
//	size_t const num_base = 1;
//	size_t const num_interpolated = 200000000;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0);
//	InitializeFloatArray(num_base, y_base_, 1.0);
//	double dx = 1.0 / static_cast<double>(num_interpolated - 1);
//	for (size_t i = 0; i < num_interpolated; ++i) {
//		x_interpolated_[i] = -0.5 + dx * static_cast<double>(i);
//	}
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kNearest, num_base,
//			num_interpolated, sakura_Status_kOK, false);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, NearestPerformance) {
//	// initial setup
//	size_t const num_base = 2;
//	size_t const num_interpolated = 200000000;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 100.0);
//	InitializeFloatArray(num_base, y_base_, 1.0, 0.0);
//	double dx = fabs(x_base_[1] - x_base_[0])
//			/ static_cast<double>(num_interpolated - 1);
//	for (size_t i = 0; i < num_interpolated; ++i) {
//		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
//	}
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kNearest, num_base,
//			num_interpolated, sakura_Status_kOK, false);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, NearestOppositePerformance) {
//	// initial setup
//	size_t const num_base = 2;
//	size_t const num_interpolated = 200000000;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 100.0, 0.0);
//	InitializeFloatArray(num_base, y_base_, 0.0, 1.0);
//	double dx = fabs(x_base_[1] - x_base_[0])
//			/ static_cast<double>(num_interpolated - 1);
//	for (size_t i = 0; i < num_interpolated; ++i) {
//		x_interpolated_[i] = x_base_[0] - dx * static_cast<double>(i);
//	}
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kNearest, num_base,
//			num_interpolated, sakura_Status_kOK, false);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, LinearPerformance) {
//	// initial setup
//	size_t const num_base = 2;
//	size_t const num_interpolated = 200000000;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 100.0);
//	InitializeFloatArray(num_base, y_base_, 1.0, 0.0);
//	double dx = fabs(x_base_[1] - x_base_[0])
//			/ static_cast<double>(num_interpolated - 1);
//	for (size_t i = 0; i < num_interpolated; ++i) {
//		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
//	}
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kLinear, num_base,
//			num_interpolated, sakura_Status_kOK, false);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, LinearOppositePerformance) {
//	// initial setup
//	size_t const num_base = 2;
//	size_t const num_interpolated = 200000000;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 100.0, 0.0);
//	InitializeFloatArray(num_base, y_base_, 0.0, 1.0);
//	double dx = fabs(x_base_[1] - x_base_[0])
//			/ static_cast<double>(num_interpolated - 1);
//	for (size_t i = 0; i < num_interpolated; ++i) {
//		x_interpolated_[i] = x_base_[0] - dx * static_cast<double>(i);
//	}
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kLinear, num_base,
//			num_interpolated, sakura_Status_kOK, false);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, PolynomialOrder2FullPerformance) {
//	// initial setup
//	polynomial_order_ = 2;
//	size_t const num_base = 3;
//	size_t const num_interpolated = 100000000; // 1/2 of Nearest and Linear
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 100.0, 200.0);
//	InitializeFloatArray(num_base, y_base_, 1.0, -1.0, 0.0);
//	double dx = fabs(x_base_[2] - x_base_[0])
//			/ static_cast<double>(num_interpolated - 1);
//	for (size_t i = 0; i < num_interpolated; ++i) {
//		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
//	}
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kPolynomial, num_base,
//			num_interpolated, sakura_Status_kOK, false);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, PolynomialOrder2FullOppositePerformance) {
//	// initial setup
//	polynomial_order_ = 2;
//	size_t const num_base = 3;
//	size_t const num_interpolated = 100000000; // 1/2 of Nearest and Linear
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 200.0, 100.0, 0.0);
//	InitializeFloatArray(num_base, y_base_, 0.0, -1.0, 1.0);
//	double dx = fabs(x_base_[2] - x_base_[0])
//			/ static_cast<double>(num_interpolated - 1);
//	for (size_t i = 0; i < num_interpolated; ++i) {
//		x_interpolated_[i] = x_base_[0] - dx * static_cast<double>(i);
//	}
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kPolynomial, num_base,
//			num_interpolated, sakura_Status_kOK, false);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, SplinePerformance) {
//	// initial setup
//	size_t const num_base = 3;
//	size_t const num_interpolated = 200000000;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 100.0, 200.0);
//	InitializeFloatArray(num_base, y_base_, 1.0, -1.0, 0.0);
//	double dx = fabs(x_base_[2] - x_base_[0])
//			/ static_cast<double>(num_interpolated - 1);
//	for (size_t i = 0; i < num_interpolated; ++i) {
//		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
//	}
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kSpline, num_base,
//			num_interpolated, sakura_Status_kOK, false);
//}
//
//TEST_F(InterpolatePseudo2dFloatTest, SplineOppositePerformance) {
//	// initial setup
//	size_t const num_base = 3;
//	size_t const num_interpolated = 200000000;
//	AllocateMemory(num_base, num_interpolated);
//	InitializeDoubleArray(num_base, x_base_, 200.0, 100.0, 0.0);
//	InitializeFloatArray(num_base, y_base_, 0.0, -1.0, 1.0);
//	double dx = fabs(x_base_[2] - x_base_[0])
//			/ static_cast<double>(num_interpolated - 1);
//	for (size_t i = 0; i < num_interpolated; ++i) {
//		x_interpolated_[i] = x_base_[0] - dx * static_cast<double>(i);
//	}
//
//	// execute interpolation
//	RunInterpolate1d(sakura_InterpolationMethod_kSpline, num_base,
//			num_interpolated, sakura_Status_kOK, false);
//}

