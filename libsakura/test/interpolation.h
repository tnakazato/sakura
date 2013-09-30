/*
 * interpolation.h
 *
 *  Created on: Sep 18, 2013
 *      Author: nakazato
 */

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include <stdarg.h>

void InitializeDoubleArray(size_t num_array, double array[], ...) {
	va_list arguments_list;
	va_start(arguments_list, array);
	for (size_t i = 0; i < num_array; ++i) {
		array[i] = va_arg(arguments_list, double);
	}
}

void InitializeFloatArray(size_t num_array, float array[], ...) {
	va_list arguments_list;
	va_start(arguments_list, array);
	for (size_t i = 0; i < num_array; ++i) {
		array[i] = static_cast<float>(va_arg(arguments_list, double));
	}
}

class InterpolateFloatTestBase: public ::testing::Test {
protected:
	virtual void SetUp() {
		initialize_result_ = sakura_Initialize(nullptr, nullptr);
		polynomial_order_ = 0;
		sakura_alignment_ = sakura_GetAlignment();
	}
	virtual void TearDown() {
		sakura_CleanUp();
	}
	virtual void AllocateMemory(size_t num_base, size_t num_interpolated,
			size_t num_array) {
		size_t num_arena_xbase = num_base + sakura_alignment_ - 1;
		size_t num_arena_ybase = num_base * num_array + sakura_alignment_ - 1;
		size_t num_arena_xinterpolated = num_interpolated + sakura_alignment_
				- 1;
		size_t num_arena_yinterpolated = num_interpolated * num_array
				+ sakura_alignment_ - 1;
		storage_for_x_base_.reset(new double[num_arena_xbase]);
		x_base_ = sakura_AlignDouble(num_arena_xbase, storage_for_x_base_.get(),
				num_base);
		storage_for_y_base_.reset(new float[num_arena_ybase]);
		y_base_ = sakura_AlignFloat(num_arena_ybase, storage_for_y_base_.get(),
				num_base * num_array);
		storage_for_x_interpolated_.reset(new double[num_arena_xinterpolated]);
		x_interpolated_ = sakura_AlignDouble(num_arena_xinterpolated,
				storage_for_x_interpolated_.get(), num_interpolated);
		storage_for_y_interpolated_.reset(new float[num_arena_yinterpolated]);
		y_interpolated_ = sakura_AlignFloat(num_arena_yinterpolated,
				storage_for_y_interpolated_.get(),
				num_interpolated * num_array);
		storage_for_y_expected_.reset(new float[num_arena_yinterpolated]);
		y_expected_ = sakura_AlignFloat(num_arena_yinterpolated,
				storage_for_y_expected_.get(), num_interpolated * num_array);

		// check alignment
		ASSERT_TRUE(x_base_ != nullptr)<< "x_base_ is null";
		ASSERT_TRUE(sakura_IsAligned(x_base_))<< "x_base_ is not aligned";
		ASSERT_TRUE(sakura_IsAligned(y_base_))<< "y_base_ is not aligned";
		ASSERT_TRUE(sakura_IsAligned(y_interpolated_))<< "y_interpolated_ is not aligned";
	}

	virtual void InspectResult(sakura_Status expected_status, sakura_Status result_status,
			size_t num_interpolated, size_t num_array, bool check_result) {
		// Should return InvalidArgument status
		std::string message = (expected_status == sakura_Status_kOK) ?
		"InterpolateArray1DFloat had any problems during execution." :
		"InterpolateArray1DFloat should fail!";
		EXPECT_EQ(expected_status, result_status) << message;

		if (check_result && (expected_status == result_status)) {
			// Value check
			for (size_t index = 0; index < num_interpolated * num_array; ++index) {
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
	uint8_t polynomial_order_;

	std::unique_ptr<double[]> storage_for_x_base_;
	std::unique_ptr<float[]> storage_for_y_base_;
	std::unique_ptr<double[]> storage_for_x_interpolated_;
	std::unique_ptr<float[]> storage_for_y_interpolated_;
	std::unique_ptr<float[]> storage_for_y_expected_;
	double *x_base_;
	float *y_base_;
	double *x_interpolated_;
	float *y_interpolated_;
	float *y_expected_;
};

#endif /* INTERPOLATION_H_ */
