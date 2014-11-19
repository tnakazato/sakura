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
#include <memory>
#include <cmath>
#include <string>
#include <stdarg.h>

#include <libsakura/sakura.h>

#include "loginit.h"
#include "gtest/gtest.h"
#include "interpolation.h"

//#include "asap/CubicSplineInterpolator1D.h"

struct Runner {
	static sakura_Status Run(sakura_InterpolationMethod method,
			uint8_t polynomial_order, size_t num_base, double const x_base[],
			size_t num_array, float const y_base[], size_t num_interpolated,
			double const x_interpolated[], float y_interpolated[]) {
		return sakura_InterpolateYAxisFloat(method, polynomial_order, num_base,
				x_base, num_array, y_base, num_interpolated, x_interpolated,
				y_interpolated);
	}
};
typedef InterpolateFloatTestBase<Runner> InterpolateArray1DFloatTest;

TEST_F(InterpolateArray1DFloatTest, InvalidType) {
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
	RunInterpolateArray1D(sakura_InterpolationMethod_kNumMethod, num_base,
			num_interpolated, num_array, sakura_Status_kInvalidArgument, false);
}

TEST_F(InterpolateArray1DFloatTest, ZeroLengthBaseArray) {
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

TEST_F(InterpolateArray1DFloatTest, ZeroLengthInterpolatedArray) {
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

TEST_F(InterpolateArray1DFloatTest, InputArrayNotAligned) {
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

TEST_F(InterpolateArray1DFloatTest, SingleBase) {
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

TEST_F(InterpolateArray1DFloatTest, OutOfRangeLeft) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 2;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -2.0, -1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 3.0, 0.0, 0.0);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 3.0,
			1.0, 3.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, OutOfRangeRight) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 2;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 2.0, 3.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, 0.0, 1.0, 3.0);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 3.0,
			1.0, 3.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, Nearest) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -0.5, -1.0, 0.2);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.1, 0.5,
			0.9, 1.2);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, -0.5,
			1.0, -0.5, 1.0, -0.5, -1.0, 0.2, -1.0, 0.2);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, NearestOffset) {
	// initial setup
	size_t const num_base = 4;
	size_t const num_interpolated = 5;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, -3.0, -2.1, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 3.0, -2.0, 2.0, -1.0,
			1.0, -0.5, -1.0, 0.2);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.1, 0.5,
			0.9, 1.2);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, -0.5,
			1.0, -0.5, 1.0, -0.5, -1.0, 0.2, -1.0, 0.2);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, NearestDescending) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, -1.0, 0.2, 1.0, -0.5);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.1, 0.5,
			0.9, 1.2);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, -0.5,
			1.0, -0.5, 1.0, -0.5, -1.0, 0.2, -1.0, 0.2);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, NearestOpposite) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, -1.0, 0.2, 1.0, -0.5);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 1.2, 0.9, 0.5, 0.1,
			-1.0);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, -1.0, 0.2,
			-1.0, 0.2, 1.0, -0.5, 1.0, -0.5, 1.0, -0.5);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, Linear) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 0.0, -1.0, 0.5,
			0.0, -0.2);
	x_interpolated_[0] = -1.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	for (size_t i = 0; i < num_array; ++i) {
		y_expected_[i] = y_base_[i];
		y_expected_[(num_interpolated - 1) * num_array + i] = y_base_[(num_base
				- 1) * num_array + i];
	}
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		size_t left_index = (x_interpolated_[i] < x_base_[1]) ? 0 : 1;
		size_t right_index = left_index + 1;
		double fraction = (x_interpolated_[i] - x_base_[left_index])
				/ (x_base_[right_index] - x_base_[left_index]);
		for (size_t j = 0; j < num_array; ++j) {
			y_expected_[i * num_array + j] = y_base_[left_index * num_array + j]
					+ fraction
							* (y_base_[right_index * num_array + j]
									- y_base_[left_index * num_array + j]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, LinearOffset) {
	// initial setup
	size_t const num_base = 5;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, -3.0, -2.0, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 3.0, -2.0, 2.0, -1.0,
			1.0, 0.0, -1.0, 0.5, 0.0, -0.2);
	x_interpolated_[0] = -1.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[0] = 1.5;
	y_expected_[1] = -0.5;
	for (size_t i = 0; i < num_array; ++i) {
		//y_expected_[i] = y_base_[i];
		y_expected_[(num_interpolated - 1) * num_array + i] = y_base_[(num_base
				- 1) * num_array + i];
	}
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[2];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[2] + (i - 1) * increment;
	}
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		size_t left_index = (x_interpolated_[i] < x_base_[3]) ? 2 : 3;
		size_t right_index = left_index + 1;
		double fraction = (x_interpolated_[i] - x_base_[left_index])
				/ (x_base_[right_index] - x_base_[left_index]);
		for (size_t j = 0; j < num_array; ++j) {
			y_expected_[i * num_array + j] = y_base_[left_index * num_array + j]
					+ fraction
							* (y_base_[right_index * num_array + j]
									- y_base_[left_index * num_array + j]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, SimpleLinear) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 1;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 0.0, -1.0, 0.5);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 0.4);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 0.2, 0.2);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, LinearDescending) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, -0.2, -1.0, 0.5,
			1.0, 0.0);
	x_interpolated_[0] = -1.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	for (size_t i = 0; i < num_array; ++i) {
		y_expected_[i] = y_base_[(num_base - 1) * num_array + i];
		y_expected_[(num_interpolated - 1) * num_array + i] = y_base_[i];
	}
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[0] - x_base_[num_base - 1];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[num_base - 1] + (i - 1) * increment;
	}
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		size_t left_index = (x_interpolated_[i] > x_base_[1]) ? 0 : 1;
		size_t right_index = left_index + 1;
		double fraction = (x_interpolated_[i] - x_base_[left_index])
				/ (x_base_[right_index] - x_base_[left_index]);
		for (size_t j = 0; j < num_array; ++j) {
			y_expected_[i * num_array + j] = y_base_[left_index * num_array + j]
					+ fraction
							* (y_base_[right_index * num_array + j]
									- y_base_[left_index * num_array + j]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, LinearOpposite) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, -0.2, -1.0, 0.5,
			1.0, 0.0);
	x_interpolated_[0] = 10.0;
	x_interpolated_[num_interpolated - 1] = -1.0;
	for (size_t i = 0; i < num_array; ++i) {
		y_expected_[i] = y_base_[i];
		y_expected_[(num_interpolated - 1) * num_array + i] = y_base_[(num_base
				- 1) * num_array + i];
	}
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		size_t left_index = (x_interpolated_[i] > x_base_[1]) ? 0 : 1;
		size_t right_index = left_index + 1;
		double fraction = (x_interpolated_[i] - x_base_[left_index])
				/ (x_base_[right_index] - x_base_[left_index]);
		for (size_t j = 0; j < num_array; ++j) {
			y_expected_[i * num_array + j] = y_base_[left_index * num_array + j]
					+ fraction
							* (y_base_[right_index * num_array + j]
									- y_base_[left_index * num_array + j]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, PolynomialOrder0) {
	// initial setup
	polynomial_order_ = 0;
	size_t const num_base = 2;
	size_t const num_interpolated = 5;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -0.5, -1.0, 0.2);
	InitializeDoubleArray(num_interpolated, x_interpolated_, -1.0, 0.1, 0.5,
			0.9, 1.2);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, -0.5,
			1.0, -0.5, 1.0, -0.5, -1.0, 0.2, -1.0, 0.2);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, PolynomialOrder1) {
	// initial setup
	polynomial_order_ = 1;
	size_t const num_base = 2;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 0.0, -1.0, 0.5);

	x_interpolated_[0] = -1.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	for (size_t i = 0; i < num_array; ++i) {
		y_expected_[i] = y_base_[i];
		y_expected_[(num_interpolated - 1) * num_array + i] = y_base_[(num_base
				- 1) * num_array + i];
	}
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
	for (size_t j = 1; j < num_interpolated - 1; ++j) {
		double fraction = (x_interpolated_[j] - x_base_[0]) / dx;
		size_t start = j * num_array;
		for (size_t i = 0; i < num_array; ++i) {
			y_expected_[start + i] = y_base_[i]
					+ fraction * (y_base_[i + num_array] - y_base_[i]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, PolynomialOrder1Offset) {
	// initial setup
	polynomial_order_ = 1;
	size_t const num_base = 5;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, -3.0, -2.0, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 3.0, -2.0, 2.0, -1.0,
			1.0, 0.0, -1.0, 0.5, 0.0, -0.2);
	x_interpolated_[0] = -1.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[0] = 1.5;
	y_expected_[1] = -0.5;
	for (size_t i = 0; i < num_array; ++i) {
		//y_expected_[i] = y_base_[i];
		y_expected_[(num_interpolated - 1) * num_array + i] = y_base_[(num_base
				- 1) * num_array + i];
	}
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[2];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[2] + (i - 1) * increment;
	}
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		size_t left_index = (x_interpolated_[i] < x_base_[3]) ? 2 : 3;
		size_t right_index = left_index + 1;
		double fraction = (x_interpolated_[i] - x_base_[left_index])
				/ (x_base_[right_index] - x_base_[left_index]);
		for (size_t j = 0; j < num_array; ++j) {
			y_expected_[i * num_array + j] = y_base_[left_index * num_array + j]
					+ fraction
							* (y_base_[right_index * num_array + j]
									- y_base_[left_index * num_array + j]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, PolynomialOrder2Full) {
	// initial setup
	polynomial_order_ = 2;
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 5.0, -1.0, 6.0,
			0.0, 9.5);

	x_interpolated_[0] = -1.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	for (size_t i = 0; i < num_array; ++i) {
		y_expected_[i] = y_base_[i];
		y_expected_[(num_interpolated - 1) * num_array + i] = y_base_[(num_base
				- 1) * num_array + i];
	}
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		size_t offset = i * num_array;
		// expected value for first array can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
		y_expected_[offset] = (1.5 * x_interpolated_[i] - 3.5)
				* x_interpolated_[i] + 1.0;
		// expected value for second array can be calculated by y = 1.25 x^2 - 0.25 x + 5.0
		y_expected_[offset + 1] = (1.25 * x_interpolated_[i] - 0.25)
				* x_interpolated_[i] + 5.0;
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, PolynomialOrder2FullDescending) {
	// initial setup
	polynomial_order_ = 2;
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, 9.5, -1.0, 6.0,
			1.0, 5.0);

	x_interpolated_[0] = -1.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	for (size_t i = 0; i < num_array; ++i) {
		y_expected_[i] = y_base_[(num_base - 1) * num_array + i];
		y_expected_[(num_interpolated - 1) * num_array + i] = y_base_[i];
	}
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[0] - x_base_[num_base - 1];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[num_base - 1] + (i - 1) * increment;
	}
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		size_t offset = i * num_array;
		// expected value for first array can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
		y_expected_[offset] = (1.5 * x_interpolated_[i] - 3.5)
				* x_interpolated_[i] + 1.0;
		// expected value for second array can be calculated by y = 1.25 x^2 - 0.25 x + 5.0
		y_expected_[offset + 1] = (1.25 * x_interpolated_[i] - 0.25)
				* x_interpolated_[i] + 5.0;
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, PolynomialOrder2FullOpposite) {
	// initial setup
	polynomial_order_ = 2;
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, 9.5, -1.0, 6.0,
			1.0, 5.0);

	x_interpolated_[0] = 10.0;
	x_interpolated_[num_interpolated - 1] = -1.0;
	for (size_t i = 0; i < num_array; ++i) {
		y_expected_[i] = y_base_[i];
		y_expected_[(num_interpolated - 1) * num_array + i] = y_base_[(num_base
				- 1) * num_array + i];
	}
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		size_t offset = i * num_array;
		// expected value for first array can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
		y_expected_[offset] = (1.5 * x_interpolated_[i] - 3.5)
				* x_interpolated_[i] + 1.0;
		// expected value for second array can be calculated by y = 1.25 x^2 - 0.25 x + 5.0
		y_expected_[offset + 1] = (1.25 * x_interpolated_[i] - 0.25)
				* x_interpolated_[i] + 5.0;
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, PolynomialOrder1Sub) {
	// initial setup
	polynomial_order_ = 1;
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 5.0, -1.0, 6.0,
			0.0, 9.5);

	x_interpolated_[0] = -1.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	for (size_t i = 0; i < num_array; ++i) {
		y_expected_[i] = y_base_[i];
		y_expected_[(num_interpolated - 1) * num_array + i] = y_base_[(num_base
				- 1) * num_array + i];
	}
	size_t const num_segments = 10;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}

	for (size_t j = 1; j < num_interpolated - 1; ++j) {
		// polynomial interpolation is a linear interpolation in the sub-region
		// x_base_[offset] ~ x_base_[offset+1]
		size_t offset = (x_interpolated_[j] < x_base_[1]) ? 0 : 1;
		double fraction = (x_interpolated_[j] - x_base_[offset])
				/ (x_base_[offset + 1] - x_base_[offset]);
		for (size_t i = 0; i < num_array; ++i) {
			y_expected_[j * num_array + i] = y_base_[offset * num_array + i]
					+ fraction
							* (y_base_[(offset + 1) * num_array + i]
									- y_base_[offset * num_array + i]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, Spline) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, 5.0, -1.0, 6.0,
			0.0, 9.5);
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	x_interpolated_[0] = -1.0;
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
	x_interpolated_[num_interpolated - 1] = 10.0;
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 5.0,
			1.0, 5.0, 0.456, 5.08, -0.052, 5.19, -0.488, 5.36, -0.816, 5.62,
			-1.0, 6.0, -1.016, 6.52, -0.888, 7.16, -0.652, 7.89, -0.344, 8.68,
			0.0, 9.5, 0.0, 9.5);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, SplineOffset) {
	// initial setup
	size_t const num_base = 5;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, -3.0, -2.0, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 4.0, -2.0, 3.0, -1.0,
			1.0, 1.0, 0.0, 2.0, -1.0, 3.0);
	x_interpolated_[0] = -1.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[0] = 2.0;
	y_expected_[1] = 0.0;
	for (size_t i = 0; i < num_array; ++i) {
		//y_expected_[i] = y_base_[i];
		y_expected_[(num_interpolated - 1) * num_array + i] = y_base_[(num_base
				- 1) * num_array + i];
	}
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[2];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[2] + (i - 1) * increment;
	}
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		size_t left_index = (x_interpolated_[i] < x_base_[3]) ? 2 : 3;
		size_t right_index = left_index + 1;
		double fraction = (x_interpolated_[i] - x_base_[left_index])
				/ (x_base_[right_index] - x_base_[left_index]);
		for (size_t j = 0; j < num_array; ++j) {
			y_expected_[i * num_array + j] = y_base_[left_index * num_array + j]
					+ fraction
							* (y_base_[right_index * num_array + j]
									- y_base_[left_index * num_array + j]);
		}
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, SplineDescending) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, 9.5, -1.0, 6.0,
			1.0, 5.0);
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[0] - x_base_[num_base - 1];
	double increment = dx / static_cast<double>(num_segments);
	x_interpolated_[0] = -1.0;
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[num_base - 1] + (i - 1) * increment;
	}
	x_interpolated_[num_interpolated - 1] = 10.0;
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 5.0,
			1.0, 5.0, 0.456, 5.08, -0.052, 5.19, -0.488, 5.36, -0.816, 5.62,
			-1.0, 6.0, -1.016, 6.52, -0.888, 7.16, -0.652, 7.89, -0.344, 8.68,
			0.0, 9.5, 0.0, 9.5);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, SplineOpposite) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, 9.5, -1.0, 6.0,
			1.0, 5.0);
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	x_interpolated_[0] = 10.0;
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
	x_interpolated_[num_interpolated - 1] = -1.0;
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 0.0, 9.5,
			0.0, 9.5, -0.344, 8.68, -0.652, 7.89, -0.888, 7.16, -1.016, 6.52,
			-1.0, 6.0, -0.816, 5.62, -0.488, 5.36, -0.052, 5.19, 0.456, 5.08,
			1.0, 5.0, 1.0, 5.0);

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1DFloatTest, SingleBasePerformance) {
	// initial setup
	size_t const num_base = 1;
	size_t const num_interpolated = 4096;
	size_t const num_array = 30000;
	size_t const iter = 5;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = -1.0
				+ static_cast<double>(i) * 2.0
						/ static_cast<double>(num_interpolated);
	}

	// execute interpolation
	RunInterpolateArray1D(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false, iter);
}

#define PERFORMANCE_TEST(NAME,METHOD,NUM_BASE,NUM_INTERP,NUM_ARRAY,LEFT,RIGHT,START_POS,INCREMENT,ITER) \
	TEST_F(InterpolateArray1DFloatTest, NAME) { \
	polynomial_order_ = 2; \
	AllocateMemory(NUM_BASE, NUM_INTERP, NUM_ARRAY); \
	double dx_base = (RIGHT - LEFT) / static_cast<double>(NUM_BASE - 1); \
	for (size_t i = 0; i < NUM_BASE; ++i) { \
		x_base_[i] = LEFT + dx_base * static_cast<double>(i); \
	} \
	for (size_t i = 0; i < NUM_BASE * NUM_ARRAY; ++i) { \
		y_base_[i] = static_cast<float>(i); \
	} \
	double dx = INCREMENT \
			/ static_cast<double>(NUM_INTERP - 1); \
	for (size_t i = 0; i < NUM_INTERP; ++i) { \
		x_interpolated_[i] = START_POS + dx * static_cast<double>(i); \
	} \
	RunInterpolateArray1D(METHOD, NUM_BASE, \
			NUM_INTERP, NUM_ARRAY, sakura_Status_kOK, false, ITER); \
}

// Nearest
PERFORMANCE_TEST(NearestPerformance, sakura_InterpolationMethod_kNearest, 2,
		4096, 30000, 0.0, 100.0, 0.0, 100.0, 10)
PERFORMANCE_TEST(Nearest1DPerformance, sakura_InterpolationMethod_kNearest, 2,
		20000000, 1, 0.0, 100.0, 0.0, 100.0, 10)
PERFORMANCE_TEST(NearestDescendingPerformance,
		sakura_InterpolationMethod_kNearest, 2, 4096, 30000, 100.0, 0.0, 0.0,
		100.0, 10)
PERFORMANCE_TEST(NearestOpppsitePerformance,
		sakura_InterpolationMethod_kNearest, 2, 4096, 30000, 100.0, 0.0, 100.0,
		-100.0, 3)

// Linear
PERFORMANCE_TEST(LinearPerformance, sakura_InterpolationMethod_kLinear, 2, 4096,
		30000, 0.0, 100.0, 0.0, 100.0, 5)
PERFORMANCE_TEST(Linear1DPerformance, sakura_InterpolationMethod_kLinear, 2,
		20000000, 1, 0.0, 100.0, 0.0, 100.0, 10)
PERFORMANCE_TEST(LinearDescendingPerformance,
		sakura_InterpolationMethod_kNearest, 2, 4096, 30000, 100.0, 0.0, 0.0,
		100.0, 10)
PERFORMANCE_TEST(LinearOpppsitePerformance, sakura_InterpolationMethod_kNearest,
		2, 4096, 30000, 100.0, 0.0, 100.0, -100.0, 3)

// Polynomial
PERFORMANCE_TEST(PolynomialOrder2FullPerformance,
		sakura_InterpolationMethod_kPolynomial, 3, 4096, 10000, 0.0, 200.0, 0.0,
		200.0, 2)
PERFORMANCE_TEST(Polynomial1DOrder2FullPerformance,
		sakura_InterpolationMethod_kPolynomial, 3, 10000000, 1, 0.0, 200.0, 0.0,
		200.0, 3)
PERFORMANCE_TEST(PolynomialOrder2FullDescendingPerformance,
		sakura_InterpolationMethod_kPolynomial, 3, 4096, 10000, 200.0, 0.0, 0.0,
		200.0, 2)
PERFORMANCE_TEST(PolynomialOrder2FullOppositePerformance,
		sakura_InterpolationMethod_kPolynomial, 3, 4096, 10000, 200.0, 0.0,
		200.0, -200.0, 2)

// Spline
PERFORMANCE_TEST(SplinePerformance, sakura_InterpolationMethod_kSpline, 2, 4096,
		30000, 0.0, 100.0, 0.0, 100.0, 3)
PERFORMANCE_TEST(Spline1DPerformance, sakura_InterpolationMethod_kSpline, 2,
		20000000, 1, 0.0, 100.0, 0.0, 100.0, 4)
PERFORMANCE_TEST(SplineDescendingPerformance,
		sakura_InterpolationMethod_kSpline, 2, 4096, 30000, 100.0, 0.0, 0.0,
		100.0, 3)
PERFORMANCE_TEST(SplineOpppsitePerformance, sakura_InterpolationMethod_kSpline,
		2, 4096, 30000, 100.0, 0.0, 100.0, -100.0, 2)

// AasapSplinePerformance
// Test performance of asap spline interpolator.
// To enable this test, you have to follow the procedure below:
//
// 1. Uncomment the line at the top of this file that includes
//    "CubicSplineInterpolator1D.h".
//
// 2. Copy the following files to test/asap directory from asap/src:
//        - Interpolator1D.h, Interpolator1D.tcc
//        - CubicSplineInterpolator1D.h, CubicSplineInterpolator1D.tcc
//        - Locator.h, Locator.tcc
//        - BisectionLocator.h, BisectionLocator.tcc
//    These files should be located at test directory.
//
// 3. Comment out all casa related lines from these files.
//
// 4. Change Interpolator1D::createLocator() so that it always uses
//    BisectionLocator.
//
//TEST_F(InterpolateArray1DFloatTest, AsapSplinePerformance) {
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

