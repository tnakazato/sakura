/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2015
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

#include <aligned_memory.h>
#include <cmath>
#include <cstdint>
#include <gtest/gtest.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

template<typename Initializer, typename Printer, typename Checker>
void RunTest(size_t num_data, float fraction, double pixel_scale,
		double const *blc_x = NULL, double const *blc_y = NULL,
		double const *trc_x = NULL, double const *trc_y = NULL,
		std::string const test_name = "") {
	// Initialize Sakura
	sakura_Status init_status = sakura_Initialize(nullptr, nullptr);
	ASSERT_EQ(sakura_Status_kOK, init_status);

	SIMD_ALIGN
	double x_storage[num_data + 1];
	SIMD_ALIGN
	double y_storage[num_data + 1];
	SIMD_ALIGN
	bool mask_storage[num_data + 1];
	SIMD_ALIGN
	bool mask_expected[num_data];
	double *x = x_storage;
	double *y = y_storage;
	bool *mask = mask_storage;
	sakura_Status status_expected = sakura_Status_kOK;
	Initializer::Initialize(num_data, fraction, x_storage, y_storage,
			mask_storage, &x, &y, &mask, mask_expected, &status_expected);

	double start_time = sakura_GetCurrentTime();
	sakura_Status status = sakura_CreateMaskNearEdgeDouble(fraction,
			pixel_scale, num_data, x, y, blc_x, blc_y, trc_x, trc_y, mask);
	double end_time = sakura_GetCurrentTime();
	Printer::Print(test_name, end_time - start_time);

	ASSERT_EQ(status_expected, status);
	Checker::Check(num_data, mask, mask_expected);

	sakura_CleanUp();
}

#define TEST_MASK(Name) TEST(CreateMaskNearEdgeTest, Name)

struct InvalidArgumentInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[],
			bool const mask_in[], double **x_out, double **y_out,
			bool **mask_out,
			bool mask_expected[], sakura_Status *status_expected) {
		*status_expected = sakura_Status_kInvalidArgument;
	}
};

struct NotAlignedXInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[],
			bool const mask_in[], double **x_out, double **y_out,
			bool **mask_out,
			bool mask_expected[], sakura_Status *status_expected) {
		InvalidArgumentInitializer::Initialize(num_data, fraction, x_in, y_in,
				mask_in, x_out, y_out, mask_out, mask_expected,
				status_expected);
		*x_out = const_cast<double *>(x_in + 1);
	}
};

struct NotAlignedYInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[],
			bool const mask_in[], double **x_out, double **y_out,
			bool **mask_out,
			bool mask_expected[], sakura_Status *status_expected) {
		InvalidArgumentInitializer::Initialize(num_data, fraction, x_in, y_in,
				mask_in, x_out, y_out, mask_out, mask_expected,
				status_expected);
		*y_out = const_cast<double *>(y_in + 1);
	}
};

struct NotAlignedMaskInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[],
			bool const mask_in[], double **x_out, double **y_out,
			bool **mask_out,
			bool mask_expected[], sakura_Status *status_expected) {
		InvalidArgumentInitializer::Initialize(num_data, fraction, x_in, y_in,
				mask_in, x_out, y_out, mask_out, mask_expected,
				status_expected);
		*mask_out = const_cast<bool *>(mask_in + 1);
	}
};

struct NullXInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[],
			bool const mask_in[], double **x_out, double **y_out,
			bool **mask_out,
			bool mask_expected[], sakura_Status *status_expected) {
		InvalidArgumentInitializer::Initialize(num_data, fraction, x_in, y_in,
				mask_in, x_out, y_out, mask_out, mask_expected,
				status_expected);
		*x_out = nullptr;
	}
};

struct NullYInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[],
			bool const mask_in[], double **x_out, double **y_out,
			bool **mask_out,
			bool mask_expected[], sakura_Status *status_expected) {
		InvalidArgumentInitializer::Initialize(num_data, fraction, x_in, y_in,
				mask_in, x_out, y_out, mask_out, mask_expected,
				status_expected);
		*y_out = nullptr;
	}
};

struct NullMaskInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[],
			bool const mask_in[], double **x_out, double **y_out,
			bool **mask_out,
			bool mask_expected[], sakura_Status *status_expected) {
		InvalidArgumentInitializer::Initialize(num_data, fraction, x_in, y_in,
				mask_in, x_out, y_out, mask_out, mask_expected,
				status_expected);
		*mask_out = nullptr;
	}
};

struct BaseInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[],
			bool const mask_in[], double **x_out, double **y_out,
			bool **mask_out,
			bool mask_expected[], sakura_Status *status_expected) {
		*mask_out = const_cast<bool *>(mask_in);
		for (size_t i = 0; i < num_data; ++i) {
			(*mask_out)[i] = false;
			mask_expected[i] = false;
		}
		*status_expected = sakura_Status_kOK;
	}
};

struct StandardSquareParamSet {
	static constexpr double kOuterLength = 6.0;
	static constexpr double kMiddleLength = 4.0;
	static constexpr double kInnerLength = 2.0;
};

struct WiderSquareParamSet {
	static constexpr double kOuterLength = 250.0;
	static constexpr double kMiddleLength = 4.0;
	static constexpr double kInnerLength = 2.0;
};

template<typename ParamSet>
struct SquareShapeInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[],
			bool const mask_in[], double **x_out, double **y_out,
			bool **mask_out,
			bool mask_expected[], sakura_Status *status_expected) {
		BaseInitializer::Initialize(num_data, fraction, x_in, y_in, mask_in,
				x_out, y_out, mask_out, mask_expected, status_expected);

		*x_out = const_cast<double *>(x_in);
		*y_out = const_cast<double *>(y_in);

		constexpr size_t kNumSquare = 10;
		ASSERT_GT(num_data, kNumSquare);
		ASSERT_EQ(0, num_data % kNumSquare);
		size_t const num_data_per_square = num_data / kNumSquare;

		// outer square
		constexpr double kOuterLength = ParamSet::kOuterLength;
		DrawSquare(0, num_data_per_square, kOuterLength, *x_out, *y_out);

		// middle square
		constexpr double kMiddleLength = ParamSet::kMiddleLength;
		DrawSquare(num_data_per_square, num_data_per_square, kMiddleLength,
				*x_out, *y_out);

		// inner square
		constexpr double kInnerLength = ParamSet::kInnerLength;
		for (size_t i = 0; i < kNumSquare - 2; ++i) {
			size_t start_index = (i + 2) * num_data_per_square;
			DrawSquare(start_index, num_data_per_square, kInnerLength, *x_out,
					*y_out);
		}

		// expected mask
		size_t mask_length = fraction * static_cast<float>(num_data);
		//std::cout << "SquareShapeInitializer: mask_length = " << mask_length
		//		<< std::endl;
		for (size_t i = 0; i < mask_length; ++i) {
			mask_expected[i] = true;
		}
	}

private:
	static void DrawSide(size_t start_index, size_t length, double x0,
			double dx, double y0, double dy, double x[], double y[]) {
		size_t k = 0;
		size_t end_index = start_index + length;
		for (size_t i = start_index; i < end_index; ++i) {
			x[i] = x0 + k * dx;
			y[i] = y0 + k * dy;
			++k;
		}
	}

	static void DrawSquare(size_t start_index, size_t num_points,
			double half_length, double x[], double y[]) {
		size_t num_points_side = num_points / 4l;
		double separation = 2.0 * half_length
				/ static_cast<double>(num_points_side);

		// from top-left to top-right
		size_t start_index_local = start_index;
		DrawSide(start_index_local, num_points_side, -half_length, separation,
				half_length, 0.0, x, y);

		// from top-right to bottom-right
		start_index_local += num_points_side;
		DrawSide(start_index_local, num_points_side, half_length, 0.0,
				half_length, -separation, x, y);

		// from bottom-right to bottom-left
		start_index_local += num_points_side;
		DrawSide(start_index_local, num_points_side, half_length, -separation,
				-half_length, 0.0, x, y);

		// from bottom-left to top-right
		start_index_local += num_points_side;
		DrawSide(start_index_local, num_points_side, -half_length, 0.0,
				-half_length, separation, x, y);
	}
};

template<typename ParamSet>
struct FailedSquareShapeInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[],
			bool const mask_in[], double **x_out, double **y_out,
			bool **mask_out,
			bool mask_expected[], sakura_Status *status_expected) {
		SquareShapeInitializer<ParamSet>::Initialize(num_data, fraction, x_in,
				y_in, mask_in, x_out, y_out, mask_out, mask_expected,
				status_expected);
		*status_expected = sakura_Status_kInvalidArgument;
	}
};

template<typename ParamSet>
struct UserDefinedRangeSquareShapeInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[],
			bool const mask_in[], double **x_out, double **y_out,
			bool **mask_out,
			bool mask_expected[], sakura_Status *status_expected) {
		SquareShapeInitializer<ParamSet>::Initialize(num_data, fraction, x_in,
				y_in, mask_in, x_out, y_out, mask_out, mask_expected,
				status_expected);
		for (size_t i = 0; i < num_data; ++i) {
			mask_expected[i] = false;
		}
		// outermost square is excluded by the user-defined range
		size_t mask_length = fraction * static_cast<float>(num_data);
		for (size_t i = 0; i < mask_length; ++i) {
			mask_expected[i + 100] = true;
		}
	}
};

struct StandardCircularParamSet {
	static constexpr double kOutermostRadius = 6.4;
	static constexpr double kSecondOutermostRadius = 6.0;
	static constexpr double kInnermostRadius = 3.2;
	static constexpr double kIncrement = 0.2;
};

template<typename ParamSet>
struct CircularShapeInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[],
			bool const mask_in[], double **x_out, double **y_out,
			bool **mask_out,
			bool mask_expected[], sakura_Status *status_expected) {
		BaseInitializer::Initialize(num_data, fraction, x_in, y_in, mask_in,
				x_out, y_out, mask_out, mask_expected, status_expected);

		*x_out = const_cast<double *>(x_in);
		*y_out = const_cast<double *>(y_in);

		// PI
		double const const_pi = std::atan(1.0) * 4.0;

		constexpr size_t kNumAngle = 100;
		ASSERT_GT(num_data, kNumAngle);
		ASSERT_EQ(0, num_data % kNumAngle);
		size_t const num_radius = num_data / kNumAngle;
		SIMD_ALIGN
		double angle[kNumAngle];
		double const angle_increment = 2.0 * const_pi
				/ static_cast<double>(kNumAngle);
		for (size_t i = 0; i < kNumAngle; ++i) {
			angle[i] = static_cast<double>(i) * angle_increment;
		}
		SIMD_ALIGN
		double radius[num_radius];
		for (size_t i = 0; i < num_radius - 2; ++i) {
			radius[i] = ParamSet::kInnermostRadius
					+ ParamSet::kIncrement * static_cast<double>(i);
		}
		radius[num_radius - 2] = ParamSet::kSecondOutermostRadius;
		radius[num_radius - 1] = ParamSet::kOutermostRadius;
		size_t k = 0;
		for (size_t i = 0; i < kNumAngle; ++i) {
			double s = std::sin(angle[i]);
			double c = std::cos(angle[i]);
			for (size_t j = 0; j < num_radius; ++j) {
				double r = radius[j];
				(*x_out)[k] = r * s;
				(*y_out)[k] = r * c;
				++k;
			}
		}

		// mask
		size_t threshold = static_cast<size_t>(std::round(
				(1.0 - fraction) * static_cast<float>(num_radius)));
		//std::cout << "index mod threshold for masking: " << threshold
		//		<< std::endl;
		for (size_t i = 0; i < num_data; ++i) {
			size_t mod = i % num_radius;
			if (mod >= threshold) {
				mask_expected[i] = true;
			}
		}
	}
};

typedef SquareShapeInitializer<StandardSquareParamSet> StandardSquare;
typedef SquareShapeInitializer<WiderSquareParamSet> WiderSquare;
typedef FailedSquareShapeInitializer<StandardSquareParamSet> FailedSquare;
typedef UserDefinedRangeSquareShapeInitializer<WiderSquareParamSet> UserDefinedRangeSquare;
typedef CircularShapeInitializer<StandardCircularParamSet> StandardCircle;

struct NullChecker {
	static void Check(size_t num_data, bool const mask[],
	bool const mask_expected[]) {
	}
};

struct StandardChecker {
	static void Check(size_t num_data, bool const mask[],
	bool const mask_expected[]) {
		for (size_t i = 0; i < num_data; ++i) {
			EXPECT_EQ(mask_expected[i], mask[i]) << i;
		}
	}
};

struct PerformanceTestOutput {
	static void Print(std::string const test_name, double elapsed_time) {
		std::cout << "#x# benchmark " << test_name << " " << elapsed_time
				<< std::endl;
	}
};

struct NullOutput {
	static void Print(std::string const test_name, double elapsed_time) {
	}
};

// FAILURE CASES
TEST_MASK(NagativeFraction) {
	RunTest<InvalidArgumentInitializer, NullOutput, NullChecker>(10, -0.1f,
			0.5);
}

TEST_MASK(FractionLargerThanOne) {
	RunTest<InvalidArgumentInitializer, NullOutput, NullChecker>(10, 1.5f, 0.5);
}

TEST_MASK(NegativePixelScale) {
	RunTest<InvalidArgumentInitializer, NullOutput, NullChecker>(10, 0.1f,
			-0.5);
}

TEST_MASK(ZeroPixelScale) {
	RunTest<InvalidArgumentInitializer, NullOutput, NullChecker>(10, 0.1f, 0.0);
}

TEST_MASK(ArrayNotAligned) {
	// x is not aligned
	RunTest<NotAlignedXInitializer, NullOutput, NullChecker>(10, 0.1f, 0.5);

	// y is not aligned
	RunTest<NotAlignedYInitializer, NullOutput, NullChecker>(10, 0.1f, 0.5);

	// mask is not aligned
	RunTest<NotAlignedMaskInitializer, NullOutput, NullChecker>(10, 0.1f, 0.5);
}

TEST_MASK(ArrayIsNull) {
	// x is nullptr
	RunTest<NullXInitializer, NullOutput, NullChecker>(10, 0.1f, 0.5);

	// y is nullptr
	RunTest<NullYInitializer, NullOutput, NullChecker>(10, 0.1f, 0.5);

	// mask is nullptr
	RunTest<NullMaskInitializer, NullOutput, NullChecker>(10, 0.1f, 0.5);
}

TEST_MASK(InvalidUserDefinedRange) {
	// all data are located at left side of user defined bottom left corner
	double blc_x = 11.0;
	RunTest<FailedSquare, NullOutput, NullChecker>(1000, 0.1f, 0.5, &blc_x,
	NULL, NULL, NULL);

	// all data are located at right side of user defined top right corner
	double trc_x = -11.0;
	RunTest<FailedSquare, NullOutput, NullChecker>(1000, 0.1f, 0.5, NULL,
	NULL, &trc_x, NULL);

	// all data are located at lower side of user defined bottom left corner
	double blc_y = 11.0;
	RunTest<FailedSquare, NullOutput, NullChecker>(1000, 0.1f, 0.5, NULL,
			&blc_y, NULL, NULL);

	// all data are located at lower side of user defined top right corner
	double trc_y = -11.0;
	RunTest<FailedSquare, NullOutput, NullChecker>(1000, 0.1f, 0.5, NULL,
	NULL, NULL, &trc_y);
}

// SUCCESSFUL CASES
TEST_MASK(ZeroFraction) {
	RunTest<BaseInitializer, NullOutput, StandardChecker>(10, 0.0f, 0.5);
}

TEST_MASK(ZeroDataSize) {
	RunTest<BaseInitializer, NullOutput, StandardChecker>(0, 0.1f, 0.5);
}

TEST_MASK(FractionEqualToOne) {
	RunTest<StandardSquare, NullOutput, StandardChecker>(1000, 1.0f, 0.5);
}

TEST_MASK(FractionTenPercent) {
	// square spiral pattern
	RunTest<StandardSquare, NullOutput, StandardChecker>(1000, 0.1f, 0.5);

	// circular zigzag pattern
	RunTest<StandardCircle, NullOutput, StandardChecker>(1000, 0.1f, 0.5);
}

TEST_MASK(FractionTwentyPercent) {
	// square spiral pattern
	RunTest<StandardSquare, NullOutput, StandardChecker>(1000, 0.2f, 0.5);

	// circular zigzag pattern
	RunTest<StandardCircle, NullOutput, StandardChecker>(1000, 0.2f, 0.5);
}

TEST_MASK(ValidUserDefinedRange) {
	double blc_x = -5.0;
	double blc_y = -5.0;
	double trc_x = 5.0;
	double trc_y = 5.0;
	RunTest<UserDefinedRangeSquare, NullOutput, StandardChecker>(1000, 0.1f,
			0.5, &blc_x, &blc_y, &trc_x, &trc_y);
}

// PERFORMANCE TEST
TEST_MASK(LargeNumDataPerformance) {
	RunTest<StandardSquare, PerformanceTestOutput, StandardChecker>(40000, 0.1f,
			0.5, NULL, NULL, NULL, NULL,
			"MaskDataNearEdge_LargeNumDataPerformance");
}

TEST_MASK(LargeAreaPerformance) {
	RunTest<WiderSquare, PerformanceTestOutput, StandardChecker>(1000, 0.1f,
			0.5, NULL, NULL, NULL, NULL,
			"MaskDataNearEdge_LargeAreaPerformance");
}

