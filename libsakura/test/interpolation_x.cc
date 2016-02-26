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
#include <libsakura/sakura.h>

#include <stdarg.h>

#include <iostream>
#include <memory>
#include <cmath>
#include <string>

#include "gtest/gtest.h"
#include "asap/CubicSplineInterpolator1D.h"

#include "aligned_memory.h"
#include "interpolation.h"
#include "loginit.h"

#define TEST_INTERP_X(name) TEST_F(InterpolateArray1DFloatTest,name)

struct Runner {
	static sakura_Status Run(sakura_InterpolationMethod method,
			uint8_t polynomial_order, size_t num_base, double const x_base[],
			size_t num_array, float const y_base[], bool const mask_base[],
			size_t num_interpolated, double const x_interpolated[],
			float y_interpolated[], bool mask_interpolated[]) {
		return sakura_InterpolateXAxisFloat(method, polynomial_order, num_base,
				x_base, num_array, y_base, mask_base, num_interpolated,
				x_interpolated, y_interpolated, mask_interpolated);
	}
};
typedef InterpolateFloatTestBase<Runner> InterpolateArray1DFloatTest;

TEST_INTERP_X(InvalidType) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 1;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, 1);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 0.0);

	// execute interpolation
	// Should return InvalidArgument status
	RunInterpolateArray1D(sakura_InterpolationMethod_kNumElements, num_base,
			num_interpolated, num_array, sakura_Status_kInvalidArgument, false);
}

TEST_INTERP_X(ZeroLengthBaseArray) {
	// initial setup
	size_t const num_base = 1;
	size_t const num_interpolated = 1;
	size_t const num_array = 5;
	AllocateMemory(num_base, num_interpolated, num_array);

	// execute interpolation
	// Should return InvalidArgument status
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, 0,
			num_interpolated, num_array, sakura_Status_kInvalidArgument, false);
}

TEST_INTERP_X(ZeroLengthInterpolatedArray) {
	// initial setup
	size_t const num_base = 1;
	size_t const num_interpolated = 1;
	size_t const num_array = 5;
	AllocateMemory(num_base, num_interpolated, num_array);

	// execute interpolation
	// Should return InvalidArgument status
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base, 0,
			num_array, sakura_Status_kOK, false);
}

TEST_INTERP_X(InputArrayNotAligned) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 1;
	size_t const num_array = 5;
	AllocateMemory(num_base, num_interpolated, num_array);

	// execute interpolation
	x_base_ = &(x_base_[1]);
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, 1,
			num_interpolated, num_array, sakura_Status_kInvalidArgument,
			false);
}

namespace {
void *AlwaysFailNew(size_t size) {
	return nullptr;
}
#define ELEMENTSOF(x) (sizeof(x) / sizeof((x)[0]))
}

TEST(InterpolateArray1DFloatStandaloneTest, BadAlloc) {
	sakura_Status status = sakura_Initialize(AlwaysFailNew, nullptr);
	EXPECT_EQ(status, sakura_Status_kOK);
	uint8_t polynomial_order = 2;
	constexpr size_t num_array = 1;
	SIMD_ALIGN
	double base_position[] = { 0.0, 1.0, 2.0 };
	constexpr size_t num_base = ELEMENTSOF(base_position);
	SIMD_ALIGN
	double interpolated_position[] = { -1.0, 0.5, 1.5, 2.5 };
	constexpr size_t num_interpolated = ELEMENTSOF(interpolated_position);
	SIMD_ALIGN
	float base_data[] = { 5.0, 4.0, 3.0 };
	SIMD_ALIGN
	float interpolated_data[num_interpolated * num_array];
	SIMD_ALIGN
	bool base_mask[] = { true, true, true };
	SIMD_ALIGN
	bool interpolated_mask[num_interpolated * num_array];
	static_assert(ELEMENTSOF(base_data) == num_base * num_array, "");
	static_assert(ELEMENTSOF(base_mask) == num_base * num_array, "");
	static_assert(ELEMENTSOF(interpolated_data) == num_interpolated * num_array, "");
	static_assert(ELEMENTSOF(interpolated_mask) == num_interpolated * num_array, "");
	status = sakura_InterpolateXAxisFloat(
			sakura_InterpolationMethod_kPolynomial, polynomial_order, num_base,
			base_position, num_array, base_data, base_mask, num_interpolated,
			interpolated_position, interpolated_data, interpolated_mask);
	EXPECT_EQ(status, sakura_Status_kNoMemory);
}

#if 0
// Comment out since these tests will fail Release build
TEST_INTERP_X(DuplicateBasePosition) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 1;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 0.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 1.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kInvalidArgument,
			false);
}

TEST_INTERP_X(BasePositionNotSorted) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 1;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 0.2, 0.1);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 1.0, 1.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kInvalidArgument,
			false);
}

TEST_INTERP_X(DuplicateInterpolatedPosition) {
	// initial setup
	size_t const num_base = 1;
	size_t const num_interpolated = 2;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 1.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kInvalidArgument,
			false);
}

TEST_INTERP_X(InterpolatedPositionNotSorted) {
	// initial setup
	size_t const num_base = 1;
	size_t const num_interpolated = 3;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 1.0, 0.0, 0.5);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kInvalidArgument,
			false);
}
#endif

TEST_INTERP_X(SingleBase) {
	// initial setup
	size_t const num_base = 1;
	size_t const num_interpolated = 1;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 5.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 0.0);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 0.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(OutOfRangeLeft) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 2;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -2.0, -1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 0.0, 3.0, 0.0);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			3.0, 3.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(OutOfRangeRight) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 2;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 2.0, 3.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, 1.0, 0.0, 3.0);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			3.0, 3.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(AllDataMasked) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 2;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 2.0, 3.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, 1.0, 0.0, 3.0);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			3.0, 3.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		mask_base_[i] = false;
	}
	for (size_t i = 0; i < num_interpolated * num_array; ++i) {
		mask_expected_[i] = false;
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false, true);
}

TEST_INTERP_X(Nearest) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, -0.5, 0.2);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.1, 0.5,
			0.9, 1.2);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			1.0, -1.0, -1.0, -0.5, -0.5, -0.5, 0.2, 0.2);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(NearestOffset) {
	// initial setup
	size_t const num_base = 4;
	size_t const num_interpolated = 5;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, -3.0, -2.1, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 3.0, 2.0, 1.0, -1.0,
			-2.0, -1.0, -0.5, 0.2);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.1, 0.5,
			0.9, 1.2);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			1.0, -1.0, -1.0, -0.5, -0.5, -0.5, 0.2, 0.2);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(Nearest1D) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 6;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
			0.5, 0.7, 1.5);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			1.0, 1.0, -1.0, -1.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(NearestDescending) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, -1.0, 1.0, 0.2, -0.5);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.1, 0.5,
			0.9, 1.2);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			1.0, -1.0, -1.0, -0.5, -0.5, -0.5, 0.2, 0.2);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(NearestOpposite) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, -1.0, 1.0, 0.2, -0.5);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 1.2, 0.9, 0.5, 0.1,
			-1.0);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, -1.0, -1.0,
			1.0, 1.0, 1.0, 0.2, 0.2, -0.5, -0.5, -0.5);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(NearestMask) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 5;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 0.5, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 1.0e10, -1.0, -0.5,
			1.0, 0.2);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.1, 0.5,
			0.9, 1.2);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			1.0, -1.0, -1.0, -0.5, -0.5, 1.0, 0.2, 0.2);
	mask_base_[1] = false;

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false, true);
}

TEST_INTERP_X(Linear) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, 0.0, 0.5);
	x_interpolated_[0] = -1.0;
	y_expected_[0] = 1.0;
	y_expected_[num_interpolated] = 0.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated * num_array - 1] = 0.5;
	double dx = x_base_[num_base - 1] - x_base_[0];
	EquallySpacedGrid(num_interpolated - 2, x_base_[0], x_base_[num_base - 1],
			&x_interpolated_[1]);

	for (size_t i = 0; i < num_array; ++i) {
		for (size_t j = 1; j < num_interpolated - 1; ++j) {
			double fraction = (x_interpolated_[j] - x_base_[0]) / dx;
			size_t start_position = i * num_base;
			size_t end_position = start_position + num_base - 1;
			size_t index = i * num_interpolated + j;
			y_expected_[index] = y_base_[start_position]
					+ fraction
							* (y_base_[end_position] - y_base_[start_position]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(LinearOffset) {
	// initial setup
	size_t const num_base = 4;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, -3.0, -2.0, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 3.0, 2.0, 1.0, -1.0,
			-2.0, -1.0, 0.0, 0.5);
	x_interpolated_[0] = -1.0;
	y_expected_[0] = 1.5;
	y_expected_[num_interpolated] = -0.5;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated * num_array - 1] = 0.5;
	double dx = x_base_[num_base - 1] - x_base_[2];
	EquallySpacedGrid(num_interpolated - 2, x_base_[2], x_base_[num_base - 1],
			&x_interpolated_[1]);
	for (size_t i = 0; i < num_array; ++i) {
		for (size_t j = 1; j < num_interpolated - 1; ++j) {
			double fraction = (x_interpolated_[j] - x_base_[2]) / dx;
			size_t start_position = i * num_base + 2;
			size_t end_position = start_position + num_base - 3;
			size_t index = i * num_interpolated + j;
			y_expected_[index] = y_base_[start_position]
					+ fraction
							* (y_base_[end_position] - y_base_[start_position]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(Linear1D) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 6;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
			0.5, 0.7, 1.5);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			0.8, 0.0, -0.4, -1.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(Linear1DSingleBaseMask) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 6;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 0.5, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 0.0, -1.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
			0.5, 0.7, 1.5);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0);
	mask_base_[0] = false;
	mask_base_[2] = false;

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true, true);
}

TEST_INTERP_X(SimpleLinear) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 1;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, 0.0, 0.5);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 0.4);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 0.2, 0.2);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(LinearDescending) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, -1.0, 1.0, 0.5, 0.0);

	x_interpolated_[0] = -1.0;
	y_expected_[0] = 1.0;
	y_expected_[num_interpolated] = 0.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated * num_array - 1] = 0.5;
	double dx = x_base_[0] - x_base_[num_base - 1];
	EquallySpacedGrid(num_interpolated - 2, x_base_[num_base - 1], x_base_[0],
			&x_interpolated_[1]);
	for (size_t i = 0; i < num_array; ++i) {
		for (size_t j = 1; j < num_interpolated - 1; ++j) {
			double fraction = -(x_interpolated_[j] - x_base_[0]) / dx;
			size_t start_position = i * num_base;
			size_t end_position = start_position + num_base - 1;
			size_t index = i * num_interpolated + j;
			y_expected_[index] = y_base_[start_position]
					+ fraction
							* (y_base_[end_position] - y_base_[start_position]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(LinearOpposite) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, -1.0, 1.0, 0.5, 0.0);

	x_interpolated_[0] = 10.0;
	y_expected_[0] = -1.0;
	y_expected_[num_interpolated] = 0.5;
	x_interpolated_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated - 1] = 1.0;
	y_expected_[num_interpolated * num_array - 1] = 0.0;
	double dx = x_base_[num_base - 1] - x_base_[0];
	EquallySpacedGrid(num_interpolated - 2, x_base_[0], x_base_[num_base - 1],
			&x_interpolated_[1]);
	for (size_t i = 0; i < num_array; ++i) {
		for (size_t j = 1; j < num_interpolated - 1; ++j) {
			double fraction = (x_interpolated_[j] - x_base_[0]) / dx;
			size_t start_position = i * num_base;
			size_t end_position = start_position + num_base - 1;
			size_t index = i * num_interpolated + j;
			y_expected_[index] = y_base_[start_position]
					+ fraction
							* (y_base_[end_position] - y_base_[start_position]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(LinearMask) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 0.5, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 0.0, -1.0, 0.0,
			-1.0e10, 0.5);
	mask_base_[4] = false;
	x_interpolated_[0] = -1.0;
	y_expected_[0] = 1.0;
	y_expected_[num_interpolated] = 0.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated * num_array - 1] = 0.5;
	double dx = x_base_[num_base - 1] - x_base_[0];
	EquallySpacedGrid(num_interpolated - 2, x_base_[0], x_base_[num_base - 1],
			&x_interpolated_[1]);
	for (size_t i = 0; i < num_array; ++i) {
		for (size_t j = 1; j < num_interpolated - 1; ++j) {
			double fraction = (x_interpolated_[j] - x_base_[0]) / dx;
			size_t start_position = i * num_base;
			size_t end_position = start_position + num_base - 1;
			size_t index = i * num_interpolated + j;
			y_expected_[index] = y_base_[start_position]
					+ fraction
							* (y_base_[end_position] - y_base_[start_position]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true, true);
}

TEST_INTERP_X(LinearMaskEdge) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 0.5, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0e10, 0.0, -1.0, 0.0,
			0.5, -1.0e10);
	mask_base_[0] = false;
	mask_base_[5] = false;
	x_interpolated_[0] = -1.0;
	y_expected_[0] = 0.0;
	y_expected_[num_interpolated] = 0.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated * num_array - 1] = 0.5;
	EquallySpacedGrid(num_interpolated - 2, x_base_[0], x_base_[num_base - 1],
			&x_interpolated_[1]);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, -0.2, -0.4, -0.6, -0.8, -1.0, -1.0, 0.0,
			0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true, true);
}

TEST_INTERP_X(LinearMaskArray) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 0.5, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 0.0, -1.0, 0.0,
			-1.0e10, 0.5);
	mask_base_[3] = false;
	mask_base_[4] = false;
	mask_base_[5] = false;
	for (size_t i = num_interpolated; i < num_interpolated * num_array; ++i) {
		mask_expected_[i] = false;
	}
	x_interpolated_[0] = -1.0;
	y_expected_[0] = 1.0;
	y_expected_[num_interpolated] = 0.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated * num_array - 1] = 0.5;
	double dx = x_base_[num_base - 1] - x_base_[0];
	EquallySpacedGrid(num_interpolated - 2, x_base_[0], x_base_[num_base - 1],
			&x_interpolated_[1]);
	for (size_t i = 0; i < num_array; ++i) {
		for (size_t j = 1; j < num_interpolated - 1; ++j) {
			double fraction = (x_interpolated_[j] - x_base_[0]) / dx;
			size_t start_position = i * num_base;
			size_t end_position = start_position + num_base - 1;
			size_t index = i * num_interpolated + j;
			y_expected_[index] = y_base_[start_position]
					+ fraction
							* (y_base_[end_position] - y_base_[start_position]);
		}
	}

	// workaround for avoiding unexpected test failure
	for (size_t i = 0; i < num_interpolated * num_array; ++i) {
		y_interpolated_[i] = y_expected_[i];
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true, true);
}

TEST_INTERP_X(PolynomialOrder0) {
	// initial setup
	polynomial_order_ = 0;
	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, -0.5, 0.2);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.1, 0.5,
			0.9, 1.2);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			1.0, -1.0, -1.0, -0.5, -0.5, -0.5, 0.2, 0.2);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(PolynomialOrder1) {
	// initial setup
	polynomial_order_ = 1;
	size_t const num_base = 2;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, 0.0, 0.5);

	x_interpolated_[0] = -1.0;
	y_expected_[0] = 1.0;
	y_expected_[num_interpolated] = 0.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated * num_array - 1] = 0.5;
	double dx = x_base_[num_base - 1] - x_base_[0];
	EquallySpacedGrid(num_interpolated - 2, x_base_[0], x_base_[num_base - 1],
			&x_interpolated_[1]);
	for (size_t i = 0; i < num_array; ++i) {
		for (size_t j = 1; j < num_interpolated - 1; ++j) {
			double fraction = (x_interpolated_[j] - x_base_[0]) / dx;
			size_t start_position = i * num_base;
			size_t end_position = start_position + num_base - 1;
			size_t index = i * num_interpolated + j;
			y_expected_[index] = y_base_[start_position]
					+ fraction
							* (y_base_[end_position] - y_base_[start_position]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(PolynomialOrder1Offset) {
	// initial setup
	polynomial_order_ = 1;
	size_t const num_base = 4;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, -3.0, -2.0, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 3.0, 2.0, 1.0, -1.0,
			-2.0, -1.0, 0.0, 0.5);
	x_interpolated_[0] = -1.0;
	y_expected_[0] = 1.5;
	y_expected_[num_interpolated] = -0.5;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated * num_array - 1] = 0.5;
	double dx = x_base_[num_base - 1] - x_base_[2];
	EquallySpacedGrid(num_interpolated - 2, x_base_[2], x_base_[num_base - 1],
			&x_interpolated_[1]);
	for (size_t i = 0; i < num_array; ++i) {
		for (size_t j = 1; j < num_interpolated - 1; ++j) {
			double fraction = (x_interpolated_[j] - x_base_[2]) / dx;
			size_t start_position = i * num_base + 2;
			size_t end_position = start_position + num_base - 3;
			size_t index = i * num_interpolated + j;
			y_expected_[index] = y_base_[start_position]
					+ fraction
							* (y_base_[end_position] - y_base_[start_position]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(Polynomial1DOrder0) {
	// initial setup
	polynomial_order_ = 0;
	size_t const num_base = 2;
	size_t const num_interpolated = 6;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
			0.5, 0.7, 1.5);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			1.0, 1.0, -1.0, -1.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(Polynomial1DOrder1) {
	// initial setup
	polynomial_order_ = 1;
	size_t const num_base = 2;
	size_t const num_interpolated = 6;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
			0.5, 0.7, 1.5);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			0.8, 0.0, -0.4, -1.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(PolynomialOrder2Full) {
	// initial setup
	polynomial_order_ = 2;
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, 0.0, 5.0,
			6.0, 9.5);

	x_interpolated_[0] = -1.0;
	y_expected_[0] = 1.0;
	y_expected_[num_interpolated] = 5.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = 0.0;
	y_expected_[num_interpolated * num_array - 1] = 9.5;
	EquallySpacedGrid(num_interpolated - 2, x_base_[0], x_base_[num_base - 1],
			&x_interpolated_[1]);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		// expected value for first array can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
		y_expected_[i] = (1.5 * x_interpolated_[i] - 3.5) * x_interpolated_[i]
				+ 1.0;
		// expected value for second array can be calculated by y = 1.25 x^2 - 0.25 x + 5.0
		y_expected_[num_interpolated + i] = (1.25 * x_interpolated_[i] - 0.25)
				* x_interpolated_[i] + 5.0;
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(PolynomialOrder2FullDescending) {
	// initial setup
	polynomial_order_ = 2;
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, -1.0, 1.0, 9.5,
			6.0, 5.0);

	x_interpolated_[0] = -1.0;
	y_expected_[0] = 1.0;
	y_expected_[num_interpolated] = 5.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = 0.0;
	y_expected_[num_interpolated * num_array - 1] = 9.5;
	EquallySpacedGrid(num_interpolated - 2, x_base_[num_base - 1], x_base_[0],
			&x_interpolated_[1]);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		// expected value for first array can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
		y_expected_[i] = (1.5 * x_interpolated_[i] - 3.5) * x_interpolated_[i]
				+ 1.0;
		// expected value for second array can be calculated by y = 1.25 x^2 - 0.25 x + 5.0
		y_expected_[num_interpolated + i] = (1.25 * x_interpolated_[i] - 0.25)
				* x_interpolated_[i] + 5.0;
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(PolynomialOrder2FullOpposite) {
	// initial setup
	polynomial_order_ = 2;
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, -1.0, 1.0, 9.5,
			6.0, 5.0);

	x_interpolated_[0] = 10.0;
	y_expected_[0] = 0.0;
	y_expected_[num_interpolated] = 9.5;
	x_interpolated_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated - 1] = 1.0;
	y_expected_[num_interpolated * num_array - 1] = 5.0;
	EquallySpacedGrid(num_interpolated - 2, x_base_[0], x_base_[num_base - 1],
			&x_interpolated_[1]);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		// expected value for first array can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
		y_expected_[i] = (1.5 * x_interpolated_[i] - 3.5) * x_interpolated_[i]
				+ 1.0;
		// expected value for second array can be calculated by y = 1.25 x^2 - 0.25 x + 5.0
		y_expected_[num_interpolated + i] = (1.25 * x_interpolated_[i] - 0.25)
				* x_interpolated_[i] + 5.0;
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(PolynomialOrder1Sub) {
	// initial setup
	polynomial_order_ = 1;
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, 0.0, 5.0,
			6.0, 9.5);

	x_interpolated_[0] = -1.0;
	y_expected_[0] = 1.0;
	y_expected_[num_interpolated] = 5.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = 0.0;
	y_expected_[num_interpolated * num_array - 1] = 9.5;
	EquallySpacedGrid(num_interpolated - 2, x_base_[0], x_base_[num_base - 1],
			&x_interpolated_[1]);

	for (size_t i = 0; i < num_array; ++i) {
		for (size_t j = 1; j < num_interpolated - 1; ++j) {
			// polynomial interpolation is a linear interpolation in the sub-region
			// x_base_[offset] ~ x_base_[offset+1]
			size_t offset = (x_interpolated_[j] < x_base_[1]) ? 0 : 1;
			double fraction = (x_interpolated_[j] - x_base_[offset])
					/ (x_base_[offset + 1] - x_base_[offset]);
			size_t start_position = i * num_base + offset;
			size_t end_position = start_position + 1;
			y_expected_[num_interpolated * i + j] = y_base_[start_position]
					+ fraction
							* (y_base_[end_position] - y_base_[start_position]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(ThreePointSpline) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, 0.0, 5.0,
			6.0, 9.5);
	x_interpolated_[0] = -1.0;
	EquallySpacedGrid(num_interpolated - 2, x_base_[0], x_base_[num_base - 1],
			&x_interpolated_[1]);
	x_interpolated_[num_interpolated - 1] = 10.0;
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			0.456, -0.052, -0.488, -0.816, -1.0, -1.016, -0.888, -0.652, -0.344,
			0.0, 0.0, 5.0, 5.0, 5.08, 5.19, 5.36, 5.62, 6.0, 6.52, 7.16, 7.89,
			8.68, 9.5, 9.5);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(FourPointSpline) {
	// initial setup
	constexpr size_t num_base = 4;
	size_t const num_interpolated = 27;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 5.0, 10.0, 15.0, 20.0);
	InitializeFloatArray(num_base, y_base_, 7.0, 10.0, -5.0, 8.0);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = i + 0.5;
	}
	// second derivatives calculated by hand
	SIMD_ALIGN
	float d2ydx2[] = { 0.0, -1.6, 2.08, 0.0 };
	size_t const left_edge = static_cast<size_t>(x_base_[0]);
	size_t const right_edge = static_cast<size_t>(x_base_[num_base - 1]);
	for (size_t i = 0; i < left_edge; ++i) {
		y_expected_[i] = y_base_[0];
	}
	for (size_t i = left_edge; i < right_edge; ++i) {
		size_t base_index = i / 5 - 1;
		auto const dx = x_base_[base_index + 1] - x_base_[base_index];
		auto const a = (x_base_[base_index + 1] - x_interpolated_[i]) / dx;
		auto const b = 1.0 - a;
		y_expected_[i] = a * y_base_[base_index] + b * y_base_[base_index + 1]
				+ ((a * a * a - a) * d2ydx2[base_index]
						+ (b * b * b - b) * d2ydx2[base_index + 1]) * dx * dx
						/ 6.0;
	}
	for (size_t i = right_edge; i < num_interpolated; ++i) {
		y_expected_[i] = y_base_[num_base - 1];
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(SparseInterpolation) {
	// initial setup
	size_t const num_base = 4;
	size_t const num_interpolated = 3;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0, 3.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, 0.0, 3.0);
	x_interpolated_[0] = -1.0;
	x_interpolated_[1] = 0.5;
	x_interpolated_[2] = 2.5;
	y_expected_[0] = 1.0;
	y_expected_[1] = 0.0;
	y_expected_[2] = 1.5;

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(OnePointInterpolation) {
	// initial setup
	size_t const num_base = 4;
	size_t const num_interpolated = 1;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0, 3.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, 0.0, 3.0);
	x_interpolated_[0] = 2.2;
	y_expected_[0] = 0.6;

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(Spline1D) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 6;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, 0.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.0, 0.1,
			0.5, 0.7, 1.5);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			0.72575, -0.28125, -0.66775, -0.78125);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(SplineForLinear) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 0.0, -1.0, 0.0,
			0.25, 0.5);
	x_interpolated_[0] = -1.0;
	y_expected_[0] = 1.0;
	y_expected_[num_interpolated] = 0.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated * num_array - 1] = 0.5;
	double dx = x_base_[num_base - 1] - x_base_[0];
	EquallySpacedGrid(num_interpolated - 2, x_base_[0], x_base_[num_base - 1],
			&x_interpolated_[1]);
	for (size_t i = 0; i < num_array; ++i) {
		for (size_t j = 1; j < num_interpolated - 1; ++j) {
			double fraction = (x_interpolated_[j] - x_base_[0]) / dx;
			size_t start_position = i * num_base;
			size_t end_position = start_position + num_base - 1;
			size_t index = i * num_interpolated + j;
			y_expected_[index] = y_base_[start_position]
					+ fraction
							* (y_base_[end_position] - y_base_[start_position]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(SplineDescending) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, -1.0, 1.0, 9.5,
			6.0, 5.0);
	x_interpolated_[0] = -1.0;
	EquallySpacedGrid(num_interpolated - 2, x_base_[num_base - 1], x_base_[0],
			&x_interpolated_[1]);
	x_interpolated_[num_interpolated - 1] = 10.0;
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			0.456, -0.052, -0.488, -0.816, -1.0, -1.016, -0.888, -0.652, -0.344,
			0.0, 0.0, 5.0, 5.0, 5.08, 5.19, 5.36, 5.62, 6.0, 6.52, 7.16, 7.89,
			8.68, 9.5, 9.5);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(SplineOpposite) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, -1.0, 1.0, 9.5,
			6.0, 5.0);
	x_interpolated_[0] = 10.0;
	EquallySpacedGrid(num_interpolated - 2, x_base_[0], x_base_[num_base - 1],
			&x_interpolated_[1]);
	x_interpolated_[num_interpolated - 1] = -1.0;
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 0.0, 0.0,
			-0.344, -0.652, -0.888, -1.016, -1.0, -0.816, -0.488, -0.052, 0.456,
			1.0, 1.0, 9.5, 9.5, 8.68, 7.89, 7.16, 6.52, 6.0, 5.62, 5.36, 5.19,
			5.08, 5.0, 5.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_INTERP_X(SingleBasePerformance) {
	// initial setup
	size_t const num_base = 1;
	size_t const num_interpolated = 4096;
	size_t const num_array = 30000;
	size_t const iter = 20;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	EquallySpacedGrid(num_interpolated, -1.0, 1.0, x_interpolated_);

	// execute interpolation
	double elapsed = RunInterpolateArray1D(sakura_InterpolationMethod_kNearest,
			num_base, num_interpolated, num_array, sakura_Status_kOK, false,
			false, iter);
	LogElapsed("InterpolationX_SingleBasePerformance", elapsed);
}

TEST_INTERP_X(SingleBase1DPerformance) {
	// initial setup
	size_t const num_base = 1;
	size_t const num_interpolated = 20000000;
	size_t const num_array = 1;
	size_t const iter = 120;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0);
	EquallySpacedGrid(num_interpolated, -0.5, 0.5, x_interpolated_);

	// execute interpolation
	double elapsed = RunInterpolateArray1D(sakura_InterpolationMethod_kNearest,
			num_base, num_interpolated, num_array, sakura_Status_kOK, false,
			false, iter);
	LogElapsed("InterpolationX_SingleBase1DPerformance", elapsed);
}

#define PERFORMANCE_TEST(NAME,METHOD,NUM_BASE,NUM_INTERP,NUM_ARRAY,LEFT,RIGHT,START_POS,INCREMENT,ITER) \
	TEST_INTERP_X(NAME) { \
	polynomial_order_ = 2; \
	AllocateMemory((NUM_BASE), (NUM_INTERP), (NUM_ARRAY)); \
	EquallySpacedGrid((NUM_BASE), (LEFT), (RIGHT), x_base_); \
	for (size_t i = 0; i < (NUM_BASE) * (NUM_ARRAY); ++i) { \
		y_base_[i] = static_cast<float>(i); \
	} \
	EquallySpacedGrid((NUM_INTERP), (START_POS), (INCREMENT) + (START_POS), x_interpolated_);\
	double elapsed = RunInterpolateArray1D((METHOD), (NUM_BASE), \
			(NUM_INTERP), (NUM_ARRAY), sakura_Status_kOK, false, false, (ITER)); \
    LogElapsed("InterpolationX_" #NAME, elapsed); \
}

// Nearest
PERFORMANCE_TEST(NearestPerformance, sakura_InterpolationMethod_kNearest, 2,
		4096, 30000, 0.0, 100.0, 0.0, 100.0, 10)
PERFORMANCE_TEST(Nearest1DPerformance, sakura_InterpolationMethod_kNearest, 2,
		20000000, 1, 0.0, 100.0, 0.0, 100.0, 40)
PERFORMANCE_TEST(NearestDescendingPerformance,
		sakura_InterpolationMethod_kNearest, 2, 4096, 30000, 100.0, 0.0, 0.0,
		100.0, 10)
PERFORMANCE_TEST(NearestOppositePerformance,
		sakura_InterpolationMethod_kNearest, 2, 4096, 30000, 100.0, 0.0, 100.0,
		-100.0, 10)

// Linear
PERFORMANCE_TEST(LinearPerformance, sakura_InterpolationMethod_kLinear, 2, 4096,
		30000, 0.0, 100.0, 0.0, 100.0, 20)
PERFORMANCE_TEST(Linear1DPerformance, sakura_InterpolationMethod_kLinear, 2,
		20000000, 1, 0.0, 100.0, 0.0, 100.0, 40)
PERFORMANCE_TEST(LinearDescendingPerformance,
		sakura_InterpolationMethod_kNearest, 2, 4096, 30000, 100.0, 0.0, 0.0,
		100.0, 10)
PERFORMANCE_TEST(LinearOppositePerformance, sakura_InterpolationMethod_kNearest,
		2, 4096, 30000, 100.0, 0.0, 100.0, -100.0, 10)

// Polynomial
PERFORMANCE_TEST(PolynomialOrder2FullPerformance,
		sakura_InterpolationMethod_kPolynomial, 3, 4096, 10000, 0.0, 200.0, 0.0,
		200.0, 2)
PERFORMANCE_TEST(Polynomial1DOrder2FullPerformance,
		sakura_InterpolationMethod_kPolynomial, 3, 10000000, 1, 0.0, 200.0, 0.0,
		200.0, 8)
PERFORMANCE_TEST(PolynomialOrder2FullDescendingPerformance,
		sakura_InterpolationMethod_kPolynomial, 3, 4096, 10000, 200.0, 0.0, 0.0,
		200.0, 2)
PERFORMANCE_TEST(PolynomialOrder2FullOppositePerformance,
		sakura_InterpolationMethod_kPolynomial, 3, 4096, 10000, 200.0, 0.0,
		200.0, -200.0, 2)

// Spline
PERFORMANCE_TEST(SplinePerformance, sakura_InterpolationMethod_kSpline, 2, 4096,
		30000, 0.0, 100.0, 0.0, 100.0, 5)
PERFORMANCE_TEST(Spline1DPerformance, sakura_InterpolationMethod_kSpline, 2,
		20000000, 1, 0.0, 100.0, 0.0, 100.0, 20)
PERFORMANCE_TEST(SplineDescendingPerformance,
		sakura_InterpolationMethod_kSpline, 2, 4096, 30000, 100.0, 0.0, 0.0,
		100.0, 5)
PERFORMANCE_TEST(SplineOppositePerformance, sakura_InterpolationMethod_kSpline,
		2, 4096, 30000, 100.0, 0.0, 100.0, -100.0, 5)

TEST_INTERP_X(SplineComparisonWithAsap) {
	// initial setup
	size_t const num_base = 128;
	size_t const num_interpolated = 500;
	size_t const num_array = 1;
	AllocateMemory(num_base, num_interpolated, num_array);
	double refpix_base = 64.0;
	double refval_base = 115.0e9;
	double increment_base = 14.0625e6;
	for (size_t i = 0; i < num_base; ++i) {
		x_base_[i] = (static_cast<double>(i) - refpix_base) * increment_base
				+ refval_base;
	}
	float y_base[num_base] = { -0.30260521, -0.27442384, -0.24624248,
			-0.21806112, -0.18987976, -0.1616984, -0.13351703, -0.10533567,
			-0.07715431, -0.04897295, -0.02079158, 0.00738978, 0.03557114,
			0.0637525, 0.09193387, 0.12011523, 0.14829659, 0.17647795,
			0.20465931, 0.23284069, 0.26102203, 0.28920341, 0.31738478,
			0.34556612, 0.3737475, 0.40192887, 0.43011022, 0.45829159,
			0.48647293, 0.51465434, 0.54283565, 0.57101703, 0.5991984,
			0.62737978, 0.65556115, 0.68374246, 0.71192384, 0.74010521,
			0.76828659, 0.79646796, 0.82464927, 0.85283065, 0.88101202,
			0.9091934, 0.93737477, 0.96555609, 0.99373746, 1.02191889,
			1.05010021, 1.07828152, 1.10646296, 1.13464427, 1.1628257,
			1.19100702, 1.21918833, 1.24736977, 1.27555108, 1.30373251,
			1.33191383, 1.36009514, 1.38827658, 1.41645789, 1.44463933,
			1.47282064, 1.50100195, 1.52918339, 1.5573647, 1.58554614,
			1.61372745, 1.64190876, 1.6700902, 1.69827151, 1.72645295,
			1.75463426, 1.78281558, 1.81099701, 1.83917832, 1.86735976,
			1.89554107, 1.92372239, 1.95190382, 1.98008513, 2.00826645,
			2.036448, 2.06462932, 2.09281063, 2.12099195, 2.14917326,
			2.17735481, 2.20553613, 2.23371744, 2.26189876, 2.29008007,
			2.31826162, 2.34644294, 2.37462425, 2.40280557, 2.43098688,
			2.45916843, 2.48734975, 2.51553106, 2.54371238, 2.57189369,
			2.60007524, 2.62825656, 2.65643787, 2.68461919, 2.7128005,
			2.74098206, 2.76916337, 2.79734468, 2.825526, 2.85370731,
			2.88188887, 2.91007018, 2.9382515, 2.96643281, 2.99461412,
			3.02279568, 3.05097699, 3.07915831, 3.10733962, 3.13552094,
			3.16370249, 3.1918838, 3.22006512, 3.24824643, 3.27642775 };
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = y_base[i];
	}
	double refpix = 250.0;
	double refval = 115.0e9;
	double increment = 1.0e6;
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = (static_cast<double>(i) - refpix) * increment
				+ refval;
	}

	// execute asap interpolation
	std::unique_ptr<asap::CubicSplineInterpolator1D<double, float> > interpolator(
			new asap::CubicSplineInterpolator1D<double, float>());
	interpolator->setX(x_base_, num_base);
	for (size_t j = 0; j < num_array; ++j) {
		float *y_base_work = &(y_base_[j * num_base]);
		float *y_interpolated_work = &(y_interpolated_[j * num_interpolated]);
		interpolator->setY(y_base_work, num_base);
		for (size_t i = 0; i < num_interpolated; ++i) {
			y_interpolated_work[i] = interpolator->interpolate(
					x_interpolated_[i]);
		}
	}

	for (size_t i = 0; i < num_interpolated * num_array; ++i) {
		y_expected_[i] = y_interpolated_[i];
	}

	// execute sakura interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

//TEST_INTERP_X(AsapSplinePerformance) {
//	// initial setup
//	size_t const num_base = 3;
//	size_t const num_interpolated = 4096;
//	size_t const num_array = 300000;
//	AllocateMemory(num_base, num_interpolated, num_array);
//	InitializeDoubleArray(num_base, x_base_, 0.0, 100.0, 200.0);
//	for (size_t i = 0; i < num_base * num_array; ++i) {
//		y_base_[i] = static_cast<float>(i);
//	}
//	double dx = fabs(x_base_[2] - x_base_[0])
//			/ static_cast<double>(num_interpolated - 1);
//	for (size_t i = 0; i < num_interpolated; ++i) {
//		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
//	}
//
//	// execute interpolation
//	std::unique_ptr<asap::CubicSplineInterpolator1D<double, float> > interpolator(
//			new asap::CubicSplineInterpolator1D<double, float>());
//	interpolator->setX(x_base_, num_base);
//	for (size_t j = 0; j < num_array; ++j) {
//		float *y_base_work = &(y_base_[j*num_base]);
//		float *y_interpolated_work = &(y_interpolated_[j*num_interpolated]);
//		interpolator->setY(y_base_work, num_base);
//		for (size_t i = 0; i < num_interpolated; ++i) {
//			y_interpolated_work[i] = interpolator->interpolate(x_interpolated_[i]);
//		}
//	}
//}

