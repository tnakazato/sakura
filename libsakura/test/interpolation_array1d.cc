#include <iostream>
#include <memory>
#include <cmath>
#include <string>
#include <stdarg.h>

#include <libsakura/sakura.h>

#include "gtest/gtest.h"
#include "interpolation.h"

//#include "asap/CubicSplineInterpolator1D.h"

class InterpolateArray1dFloatTest: public InterpolateFloatTestBase {
protected:
	virtual void RunInterpolateArray1d(
			sakura_InterpolationMethod interpolation_method, size_t num_base,
			size_t num_interpolated, size_t num_array,
			sakura_Status expected_status, bool check_result) {
		// sakura must be properly initialized
		ASSERT_EQ(sakura_Status_kOK, initialize_result_)<< "sakura must be properly initialized!";

		// execute interpolation
		sakura_Status result = sakura_InterpolateArray1dFloat(
				interpolation_method, polynomial_order_, num_base,
				x_base_, num_array, y_base_, num_interpolated, x_interpolated_, y_interpolated_);

		InspectResult(expected_status, result, num_interpolated, num_array, check_result);
	}
};

TEST_F(InterpolateArray1dFloatTest, InvalidType) {
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
	RunInterpolateArray1d(sakura_InterpolationMethod_kNumMethod, num_base,
			num_interpolated, num_array, sakura_Status_kInvalidArgument, false);
}

TEST_F(InterpolateArray1dFloatTest, ZeroLengthBaseArray) {
	// initial setup
	size_t const num_base = 0;
	size_t const num_interpolated = 1;
	size_t const num_array = 5;

	// execute interpolation
	// Should return InvalidArgument status
	RunInterpolateArray1d(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kInvalidArgument, false);
}

TEST_F(InterpolateArray1dFloatTest, ZeroLengthInterpolatedArray) {
	// initial setup
	size_t const num_base = 1;
	size_t const num_interpolated = 0;
	size_t const num_array = 5;

	// execute interpolation
	// Should return InvalidArgument status
	RunInterpolateArray1d(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, NegativePolynomialOrder) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 1;
	size_t const num_array = 5;
	polynomial_order_ = -1;

	// execute interpolation
	// Should return InvalidArgument status
	RunInterpolateArray1d(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kInvalidArgument, false);
}

TEST_F(InterpolateArray1dFloatTest, NegativePolynomialOrderButNearest) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 1;
	size_t const num_array = 2;
	polynomial_order_ = -1;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0);
	InitializeDoubleArray(num_interpolated, x_interpolated_, 0.1);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, -0.5, 0.2);
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, -0.5);

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, InputArrayNotAligned) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 1;
	size_t const num_array = 5;
	AllocateMemory(num_base, num_interpolated, num_array);

	// execute interpolation
	x_base_ = &(x_base_[1]);
	RunInterpolateArray1d(sakura_InterpolationMethod_kPolynomial, 1,
			num_interpolated, num_array, sakura_Status_kInvalidArgument,
			false);
}

TEST_F(InterpolateArray1dFloatTest, SingleBase) {
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
	RunInterpolateArray1d(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, Nearest) {
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
	RunInterpolateArray1d(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, NearestDescending) {
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
	RunInterpolateArray1d(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, NearestOpposite) {
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
	RunInterpolateArray1d(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, Linear) {
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
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
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
	RunInterpolateArray1d(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, LinearDescending) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 1.0, 0.0);
	InitializeFloatArray(num_base * num_interpolated, y_base_, -1.0, 1.0, 0.5,
			0.0);

	x_interpolated_[0] = -1.0;
	y_expected_[0] = 1.0;
	y_expected_[num_interpolated] = 0.0;
	x_interpolated_[num_interpolated - 1] = 10.0;
	y_expected_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated * num_array - 1] = 0.5;
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[0] - x_base_[num_base - 1];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[num_base - 1] + (i - 1) * increment;
	}
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
	RunInterpolateArray1d(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, LinearOpposite) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 1.0, 0.0);
	InitializeFloatArray(num_base * num_interpolated, y_base_, -1.0, 1.0, 0.5,
			0.0);

	x_interpolated_[0] = 10.0;
	y_expected_[0] = -1.0;
	y_expected_[num_interpolated] = 0.5;
	x_interpolated_[num_interpolated - 1] = -1.0;
	y_expected_[num_interpolated - 1] = 1.0;
	y_expected_[num_interpolated * num_array - 1] = 0.0;
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
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
	RunInterpolateArray1d(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, PolynomialOrder0) {
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
	RunInterpolateArray1d(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, PolynomialOrder1) {
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
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
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
	RunInterpolateArray1d(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, PolynomialOrder2Full) {
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
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		// expected value for first array can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
		y_expected_[i] = (1.5 * x_interpolated_[i] - 3.5) * x_interpolated_[i]
				+ 1.0;
		// expected value for second array can be calculated by y = 1.25 x^2 - 0.25 x + 5.0
		y_expected_[num_interpolated + i] = (1.25 * x_interpolated_[i] - 0.25)
				* x_interpolated_[i] + 5.0;
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, PolynomialOrder2FullDescending) {
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
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[0] - x_base_[num_base - 1];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[num_base - 1] + (i - 1) * increment;
	}
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		// expected value for first array can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
		y_expected_[i] = (1.5 * x_interpolated_[i] - 3.5) * x_interpolated_[i]
				+ 1.0;
		// expected value for second array can be calculated by y = 1.25 x^2 - 0.25 x + 5.0
		y_expected_[num_interpolated + i] = (1.25 * x_interpolated_[i] - 0.25)
				* x_interpolated_[i] + 5.0;
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, PolynomialOrder2FullOpposite) {
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
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
		std::cout << x_interpolated_[i] << std::endl;
	}
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		// expected value for first array can be calculated by y = 1.5 x^2 - 3.5 x + 1.0
		y_expected_[i] = (1.5 * x_interpolated_[i] - 3.5) * x_interpolated_[i]
				+ 1.0;
		// expected value for second array can be calculated by y = 1.25 x^2 - 0.25 x + 5.0
		y_expected_[num_interpolated + i] = (1.25 * x_interpolated_[i] - 0.25)
				* x_interpolated_[i] + 5.0;
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, PolynomialOrder1Sub) {
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
	size_t const num_segments = 10;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}

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
	RunInterpolateArray1d(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, Spline) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 1.0, 2.0);
	InitializeFloatArray(num_base * num_array, y_base_, 1.0, -1.0, 0.0, 5.0,
			6.0, 9.5);
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	x_interpolated_[0] = -1.0;
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
	x_interpolated_[num_interpolated - 1] = 10.0;
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			0.456, -0.052, -0.488, -0.816, -1.0, -1.016, -0.888, -0.652, -0.344,
			0.0, 0.0, 5.0, 5.0, 5.08, 5.19, 5.36, 5.62, 6.0, 6.52, 7.16, 7.89,
			8.68, 9.5, 9.5);

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, SplineDescending) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, -1.0, 1.0, 9.5,
			6.0, 5.0);
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[0] - x_base_[num_base - 1];
	double increment = dx / static_cast<double>(num_segments);
	x_interpolated_[0] = -1.0;
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[num_base - 1] + (i - 1) * increment;
	}
	x_interpolated_[num_interpolated - 1] = 10.0;
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 1.0, 1.0,
			0.456, -0.052, -0.488, -0.816, -1.0, -1.016, -0.888, -0.652, -0.344,
			0.0, 0.0, 5.0, 5.0, 5.08, 5.19, 5.36, 5.62, 6.0, 6.52, 7.16, 7.89,
			8.68, 9.5, 9.5);

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, SplineOpposite) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 13;
	size_t const num_array = 2;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 2.0, 1.0, 0.0);
	InitializeFloatArray(num_base * num_array, y_base_, 0.0, -1.0, 1.0, 9.5,
			6.0, 5.0);
	size_t const num_segments = num_interpolated - 3;
	double dx = x_base_[num_base - 1] - x_base_[0];
	double increment = dx / static_cast<double>(num_segments);
	x_interpolated_[0] = 10.0;
	for (size_t i = 1; i < num_interpolated - 1; ++i) {
		x_interpolated_[i] = x_base_[0] + (i - 1) * increment;
	}
	x_interpolated_[num_interpolated - 1] = -1.0;
	InitializeFloatArray(num_interpolated * num_array, y_expected_, 0.0, 0.0,
			-0.344, -0.652, -0.888, -1.016, -1.0, -0.816, -0.488, -0.052, 0.456,
			1.0, 1.0, 9.5, 9.5, 8.68, 7.89, 7.16, 6.52, 6.0, 5.62, 5.36, 5.19,
			5.08, 5.0, 5.0);

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, true);
}

TEST_F(InterpolateArray1dFloatTest, SingleBasePerformance) {
	// initial setup
	size_t const num_base = 1;
	size_t const num_interpolated = 4096;
	size_t const num_array = 300000;
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
	RunInterpolateArray1d(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, NearestPerformance) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 4096;
	size_t const num_array = 300000;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 100.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	double dx = fabs(x_base_[1] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, NearestDescendingPerformance) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 4096;
	size_t const num_array = 300000;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 100.0, 0.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	double dx = fabs(x_base_[0] - x_base_[1])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[1] + dx * static_cast<double>(i);
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, NearestOppositePerformance) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 4096;
	size_t const num_array = 300000;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 100.0, 0.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	double dx = fabs(x_base_[1] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] - dx * static_cast<double>(i);
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kNearest, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, LinearPerformance) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 4096;
	size_t const num_array = 300000;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 100.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	double dx = fabs(x_base_[1] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, LinearDescendingPerformance) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 4096;
	size_t const num_array = 300000;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 100.0, 0.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	double dx = fabs(x_base_[1] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[1] + dx * static_cast<double>(i);
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, LinearOppositePerformance) {
	// initial setup
	size_t const num_base = 2;
	size_t const num_interpolated = 4096;
	size_t const num_array = 300000;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 100.0, 0.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	double dx = fabs(x_base_[1] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] - dx * static_cast<double>(i);
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kLinear, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, PolynomialOrder2FullPerformance) {
	// initial setup
	polynomial_order_ = 2;
	size_t const num_base = 3;
	size_t const num_interpolated = 4096;
	size_t const num_array = 10000; // 1/3
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 100.0, 200.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	double dx = fabs(x_base_[2] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, PolynomialOrder2FullDescendingPerformance) {
	// initial setup
	polynomial_order_ = 2;
	size_t const num_base = 3;
	size_t const num_interpolated = 4096;
	size_t const num_array = 10000; // 1/3
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 200.0, 100.0, 0.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	double dx = fabs(x_base_[2] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[2] + dx * static_cast<double>(i);
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, PolynomialOrder2FullOppositePerformance) {
	// initial setup
	polynomial_order_ = 2;
	size_t const num_base = 3;
	size_t const num_interpolated = 4096;
	size_t const num_array = 10000; // 1/3
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 200.0, 100.0, 0.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	double dx = fabs(x_base_[2] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] - dx * static_cast<double>(i);
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kPolynomial, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, SplinePerformance) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 4096;
	size_t const num_array = 300000;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 0.0, 100.0, 200.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	double dx = fabs(x_base_[2] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] + dx * static_cast<double>(i);
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, SplineDescendingPerformance) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 4096;
	size_t const num_array = 300000;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 200.0, 100.0, 0.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	double dx = fabs(x_base_[2] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[2] + dx * static_cast<double>(i);
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

TEST_F(InterpolateArray1dFloatTest, SplineOppositePerformance) {
	// initial setup
	size_t const num_base = 3;
	size_t const num_interpolated = 4096;
	size_t const num_array = 300000;
	AllocateMemory(num_base, num_interpolated, num_array);
	InitializeDoubleArray(num_base, x_base_, 200.0, 100.0, 0.0);
	for (size_t i = 0; i < num_base * num_array; ++i) {
		y_base_[i] = static_cast<float>(i);
	}
	double dx = fabs(x_base_[2] - x_base_[0])
			/ static_cast<double>(num_interpolated - 1);
	for (size_t i = 0; i < num_interpolated; ++i) {
		x_interpolated_[i] = x_base_[0] - dx * static_cast<double>(i);
	}

	// execute interpolation
	RunInterpolateArray1d(sakura_InterpolationMethod_kSpline, num_base,
			num_interpolated, num_array, sakura_Status_kOK, false);
}

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
//TEST_F(InterpolateArray1dFloatTest, AsapSplinePerformance) {
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

