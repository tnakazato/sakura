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
/*
 * interpolation.h
 *
 *  Created on: Sep 18, 2013
 *      Author: nakazato
 */

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include <stdarg.h>
#include <vector>

#include "testutil.h"

inline void InitializeDoubleArray(size_t num_array, double array[], ...) {
	va_list arguments_list;
	va_start(arguments_list, array);
	for (size_t i = 0; i < num_array; ++i) {
		array[i] = va_arg(arguments_list, double);
	}
}

inline void InitializeFloatArray(size_t num_array, float array[], ...) {
	va_list arguments_list;
	va_start(arguments_list, array);
	for (size_t i = 0; i < num_array; ++i) {
		array[i] = static_cast<float>(va_arg(arguments_list, double));
	}
}

template<typename T>
inline void EquallySpacedGrid(size_t const num_grid, T const start, T const end,
		T grid[]) {
	T const incr = (end - start) / static_cast<T>(num_grid - 1);
	for (size_t i = 0; i < num_grid; ++i) {
		grid[i] = start + i * incr;
	}
}

void LogElapsed(const char *test_name, double elapsed) {
	std::cout << "#x# benchmark " << test_name << " " << elapsed << std::endl;
}

template<typename T>
class InterpolateFloatTestBase: public ::testing::Test {
protected:
	virtual void SetUp() {
		initialize_result_ = sakura_Initialize(nullptr, nullptr);
		polynomial_order_ = 0;
		sakura_alignment_ = sakura_GetAlignment();
		InitializePointers();
	}
	virtual void TearDown() {
		sakura_CleanUp();
	}
	double RunInterpolateArray1D(sakura_InterpolationMethod interpolation_method,
			size_t num_base, size_t num_interpolated, size_t num_array,
			sakura_Status expected_status, bool check_result, bool check_mask =
					false, size_t iteration = 1) {
		// sakura must be properly initialized
		EXPECT_EQ(sakura_Status_kOK, initialize_result_)
				<< "sakura must be properly initialized!";

		// execute interpolation
		double elapsed = 0.0;
		for (size_t iter = 0; iter < iteration; ++iter) {
			double start = GetCurrentTime();
			sakura_Status result = T::Run(interpolation_method,
					polynomial_order_, num_base, x_base_, num_array, y_base_,
					mask_base_, num_interpolated, x_interpolated_,
					y_interpolated_, mask_interpolated_);
			double end = GetCurrentTime();
			elapsed += end - start;
			InspectResult(expected_status, result, num_interpolated, num_array,
					check_result, check_mask);
		}
		//std::cout << "Elapsed time " << elapsed << " sec" << std::endl;
		return elapsed;
	}
	void InitializePointers() {
		x_base_ = nullptr;
		y_base_ = nullptr;
		x_interpolated_ = nullptr;
		y_interpolated_ = nullptr;
		y_expected_ = nullptr;
		mask_base_ = nullptr;
		mask_interpolated_ = nullptr;
		mask_expected_ = nullptr;
	}
	virtual void AllocateMemory(size_t num_base, size_t num_interpolated,
			size_t num_array) {
		size_t margin_double = ((sakura_alignment_ - 1) / sizeof(double) + 1);
		size_t margin_float = ((sakura_alignment_ - 1) / sizeof(float) + 1);
		size_t margin = sakura_alignment_ - 1;
		size_t num_arena_xbase = num_base + margin_double;
		size_t num_arena_ybase = num_base * num_array + margin_float;
		size_t num_arena_xinterpolated = num_interpolated + margin_double;
		size_t num_arena_yinterpolated = num_interpolated * num_array
				+ margin_float;
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
		storage_mask_.resize(3);
		size_t margin_bool = ((sakura_alignment_ - 1) / sizeof(bool) + 1);
		size_t required_mask = num_base * num_array * sizeof(bool);
		storage_mask_[0].reset(new bool[required_mask + margin_bool]);
		mask_base_ =
				reinterpret_cast<bool *>(sakura_AlignAny(required_mask + margin,
						reinterpret_cast<void *>(storage_mask_[0].get()),
						required_mask));
		required_mask = num_interpolated * num_array * sizeof(bool);
		storage_mask_[1].reset(new bool[required_mask + margin_bool]);
		storage_mask_[2].reset(new bool[required_mask + margin_bool]);
		mask_interpolated_ =
				reinterpret_cast<bool *>(sakura_AlignAny(required_mask + margin,
						reinterpret_cast<void *>(storage_mask_[1].get()),
						required_mask));
		mask_expected_ =
				reinterpret_cast<bool *>(sakura_AlignAny(required_mask + margin,
						reinterpret_cast<void *>(storage_mask_[2].get()),
						required_mask));

		float mem_bytes = 4.0
				* (num_arena_ybase + num_arena_yinterpolated
						+ num_arena_yinterpolated)
				+ 8.0 * (num_arena_xbase + num_arena_xinterpolated);
		std::cout << "memory size: " << mem_bytes / 1.0e9 << "GB" << std::endl;

		// check alignment
		ASSERT_TRUE(x_base_ != nullptr) << "x_base_ is null";
		ASSERT_TRUE(y_base_ != nullptr) << "y_base_ is null";
		ASSERT_TRUE(x_interpolated_ != nullptr) << "x_interpolated_ is null";
		ASSERT_TRUE(y_interpolated_ != nullptr) << "y_interpolated_ is null";
		ASSERT_TRUE(mask_base_ != nullptr) << "mask_base_ is null";
		ASSERT_TRUE(mask_interpolated_ != nullptr)
				<< "mask_interpoalted_ is null";

		ASSERT_TRUE(sakura_IsAligned(x_base_)) << "x_base_ is not aligned";
		ASSERT_TRUE(sakura_IsAligned(y_base_)) << "y_base_ is not aligned";
		ASSERT_TRUE(sakura_IsAligned(y_interpolated_))
				<< "y_interpolated_ is not aligned";

		// mask is initialized so that all the elements are valid
		for (size_t i = 0; i < num_base * num_array; ++i) {
			mask_base_[i] = true;
		}
		for (size_t i = 0; i < num_interpolated * num_array; ++i) {
			mask_interpolated_[i] = true;
			mask_expected_[i] = true;
		}
	}

	virtual void InspectResult(sakura_Status expected_status,
			sakura_Status result_status, size_t num_interpolated,
			size_t num_array, bool check_result, bool check_mask) {
		// Should return InvalidArgument status
		std::string message =
				(expected_status == sakura_Status_kOK) ?
						"There was any problems during execution." :
						"The execution should fail!";
		EXPECT_EQ(expected_status, result_status) << message;

		if (check_result && (expected_status == result_status)) {
			// Value check
			for (size_t index = 0; index < num_interpolated * num_array;
					++index) {
				//std::cout << "Expected value at index " << index << ": "
				//		<< y_expected_[index] << std::endl;
				EXPECT_FLOAT_EQ(y_expected_[index], y_interpolated_[index])
						<< "interpolated value differs from expected value at "
						<< index << ": " << y_expected_[index] << ", "
						<< y_interpolated_[index];
			}
		}
		// Mask check
		if (check_mask && (expected_status == result_status)) {
			for (size_t index = 0; index < num_interpolated * num_array;
					++index) {
				//std::cout << "Expected mask at index " << index << ": "
				//		<< mask_expected_[index] << std::endl;
				EXPECT_EQ(mask_expected_[index], mask_interpolated_[index])
						<< "interpolated mask differs from expected value at "
						<< index << ": " << mask_expected_[index] << ", "
						<< mask_interpolated_[index];
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
	std::vector<std::unique_ptr<bool[]> > storage_mask_;
	double *x_base_;
	float *y_base_;
	double *x_interpolated_;
	float *y_interpolated_;
	float *y_expected_;
	bool *mask_base_;
	bool *mask_interpolated_;
	bool *mask_expected_;
};


#endif /* INTERPOLATION_H_ */
