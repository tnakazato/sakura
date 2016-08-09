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

#include <aligned_memory.h>
#include <cmath>
#include <cstdint>
#include <gtest/gtest.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "testutil.h"

#define TEST_MASK(Name) TEST(CreateMaskNearEdgeTest, Name)

template<typename XInitializer, typename YInitializer, typename MaskInitializer>
struct InvalidArgumentInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[], bool const mask_in[], double **x_out,
			double **y_out, bool **mask_out, bool mask_expected[],
			sakura_Status *status_expected) {
		*status_expected = sakura_Status_kInvalidArgument;
		XInitializer::Initialize(num_data, x_in, x_out);
		YInitializer::Initialize(num_data, y_in, y_out);
		MaskInitializer::Initialize(num_data, mask_in, mask_out);
	}
};

struct NotAlignedInitializer {
	template<typename T>
	static void Initialize(size_t num_data, T const in[], T **out) {
		*out = const_cast<T *>(in + 1);
	}
};

struct NullPointerInitializer {
	template<typename T>
	static void Initialize(size_t num_data, T const in[], T **out) {
		*out = nullptr;
	}
};

struct DefaultDataInitializer {
	template<typename T>
	static void Initialize(size_t num_data, T const in[], T **out) {
		*out = const_cast<T *>(in);
		for (size_t i = 0; i < num_data; ++i) {
			(*out)[i] = static_cast<T>(i);
		}
	}
};

struct DefaultMaskInitializer {
	static void Initialize(size_t num_data, bool const in[], bool **out) {
		*out = const_cast<bool *>(in);
		for (size_t i = 0; i < num_data; ++i) {
			(*out)[i] = false;
		}
	}
};

typedef InvalidArgumentInitializer<DefaultDataInitializer,
		DefaultDataInitializer, DefaultMaskInitializer> BasicInvalidArgumentInitializer;
typedef InvalidArgumentInitializer<NotAlignedInitializer,
		DefaultDataInitializer, DefaultDataInitializer> NotAlignedXInitializer;
typedef InvalidArgumentInitializer<DefaultDataInitializer,
		NotAlignedInitializer, DefaultMaskInitializer> NotAlignedYInitializer;
typedef InvalidArgumentInitializer<DefaultDataInitializer,
		DefaultDataInitializer, NotAlignedInitializer> NotAlignedMaskInitializer;
typedef InvalidArgumentInitializer<NullPointerInitializer,
		DefaultDataInitializer, DefaultMaskInitializer> NullXInitializer;
typedef InvalidArgumentInitializer<DefaultDataInitializer,
		NullPointerInitializer, DefaultMaskInitializer> NullYInitializer;
typedef InvalidArgumentInitializer<DefaultDataInitializer,
		DefaultDataInitializer, NullPointerInitializer> NullMaskInitializer;

struct BaseInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[], bool const mask_in[], double **x_out,
			double **y_out, bool **mask_out, bool mask_expected[],
			sakura_Status *status_expected) {
		*mask_out = const_cast<bool *>(mask_in);
		for (size_t i = 0; i < num_data; ++i) {
			(*mask_out)[i] = false;
			mask_expected[i] = false;
		}
		*x_out = const_cast<double *>(x_in);
		for (size_t i = 0; i < num_data; ++i) {
			(*x_out)[i] = 0.0;
		}
		*y_out = const_cast<double *>(y_in);
		for (size_t i = 0; i < num_data; ++i) {
			(*y_out)[i] = 0.0;
		}
		*status_expected = sakura_Status_kOK;
	}
};

struct MaskNotProperlyConfiguredInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[], bool const mask_in[], double **x_out,
			double **y_out, bool **mask_out, bool mask_expected[],
			sakura_Status *status_expected) {
		BaseInitializer::Initialize(num_data, fraction, x_in, y_in, mask_in,
				x_out, y_out, mask_out, mask_expected, status_expected);

		// mimic uninitialized mask
		for (size_t i = 0; i < num_data; ++i) {
			(*mask_out)[i] = true;
		}
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

template<typename ParamSet, typename Initializer = BaseInitializer>
struct SquareShapeInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[], bool const mask_in[], double **x_out,
			double **y_out, bool **mask_out, bool mask_expected[],
			sakura_Status *status_expected) {
		Initializer::Initialize(num_data, fraction, x_in, y_in, mask_in, x_out,
				y_out, mask_out, mask_expected, status_expected);

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
			double const y_in[], bool const mask_in[], double **x_out,
			double **y_out, bool **mask_out, bool mask_expected[],
			sakura_Status *status_expected) {
		SquareShapeInitializer<ParamSet>::Initialize(num_data, fraction, x_in,
				y_in, mask_in, x_out, y_out, mask_out, mask_expected,
				status_expected);
		*status_expected = sakura_Status_kInvalidArgument;
	}
};

template<typename T>
struct WrongValueNaN {
	static T Get() {
		T nan_value = static_cast<T>(0);
		DoGet(&nan_value);
		return nan_value;
	}
private:
	static void DoGet(T *value) {
		*value = std::nan("");
		ASSERT_TRUE(std::isnan(*value));
	}
};

template<typename T>
struct WrongValuePositiveInf {
	static T Get() {
		T inf_value = static_cast<T>(0);
		DoGet(&inf_value);
		return inf_value;
	}
private:
	static void DoGet(T *value) {
		constexpr T kOne = static_cast<T>(1);
		constexpr T kZero = static_cast<T>(0);
		*value = kOne / kZero;
		ASSERT_TRUE(std::isinf(*value));
		ASSERT_TRUE(*value > kZero);
	}
};

template<typename T>
struct WrongValueNegativeInf {
	static T Get() {
		T inf_value = static_cast<T>(0);
		DoGet(&inf_value);
		return inf_value;
	}
private:
	static void DoGet(T *value) {
		constexpr T kOne = static_cast<T>(1);
		constexpr T kZero = static_cast<T>(0);
		*value = -kOne / kZero;
		ASSERT_TRUE(std::isinf(*value));
		ASSERT_TRUE(*value < kZero);
	}
};

template<typename ParamSet, typename WrongValueProvider, bool X = true>
struct FailedSquareShapeInitializerWithWrongValue {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[], bool const mask_in[], double **x_out,
			double **y_out, bool **mask_out, bool mask_expected[],
			sakura_Status *status_expected) {
		SquareShapeInitializer<ParamSet>::Initialize(num_data, fraction, x_in,
				y_in, mask_in, x_out, y_out, mask_out, mask_expected,
				status_expected);
		const double wrong_value = WrongValueProvider::Get();
		if (X) {
			(*x_out)[num_data - 1] = wrong_value;
		} else {
			(*y_out)[num_data - 1] = wrong_value;
		}
		*status_expected = sakura_Status_kInvalidArgument;
	}
};

template<typename ParamSet>
struct UserDefinedRangeSquareShapeInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[], bool const mask_in[], double **x_out,
			double **y_out, bool **mask_out, bool mask_expected[],
			sakura_Status *status_expected) {
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
	static constexpr double kInnermostRadius = 0.0;
	static constexpr double kIncrement = 0.2;
	static constexpr double kStartAngle = 0.0;
	static constexpr double kEndAngle = 360.0;
};

struct DonutShapeParamSet {
	static constexpr double kOutermostRadius = 6.4;
	static constexpr double kSecondOutermostRadius = 6.0;
	static constexpr double kInnermostRadius = 3.2;
	static constexpr double kIncrement = 0.2;
	static constexpr double kStartAngle = 0.0;
	static constexpr double kEndAngle = 360.0;
};

struct UpperCShapeParamSet {
	static constexpr double kOutermostRadius = 6.4;
	static constexpr double kSecondOutermostRadius = 6.0;
	static constexpr double kInnermostRadius = 3.2;
	static constexpr double kIncrement = 0.2;
	static constexpr double kStartAngle = 110.0;
	static constexpr double kEndAngle = 430.0;
};

struct HorizontalCShapeParamSet {
	static constexpr double kOutermostRadius = 6.4;
	static constexpr double kSecondOutermostRadius = 6.0;
	static constexpr double kInnermostRadius = 3.2;
	static constexpr double kIncrement = 0.2;
	static constexpr double kStartAngle = 20.0;
	static constexpr double kEndAngle = 340.0;
};

struct InclinedCShapeParamSet {
	static constexpr double kOutermostRadius = 6.4;
	static constexpr double kSecondOutermostRadius = 6.0;
	static constexpr double kInnermostRadius = 3.2;
	static constexpr double kIncrement = 0.2;
	static constexpr double kStartAngle = 65.0;
	static constexpr double kEndAngle = 385.0;
};

template<typename ParamSet>
struct CircularShapeInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[], bool const mask_in[], double **x_out,
			double **y_out, bool **mask_out, bool mask_expected[],
			sakura_Status *status_expected) {
		BaseInitializer::Initialize(num_data, fraction, x_in, y_in, mask_in,
				x_out, y_out, mask_out, mask_expected, status_expected);

		*x_out = const_cast<double *>(x_in);
		*y_out = const_cast<double *>(y_in);

		// PI
		double const const_pi = std::atan(1.0) * 4.0;

		constexpr size_t kNumAngle = 100;
		double const deg2rad = 2.0 * const_pi / 360.0;
		ASSERT_GT(num_data, kNumAngle);
		ASSERT_EQ(0, num_data % kNumAngle);
		size_t const num_radius = num_data / kNumAngle;
		SIMD_ALIGN
		double angle[kNumAngle];
		static_assert(0.0 <= ParamSet::kStartAngle && ParamSet::kStartAngle <= 720.0, "Internal Error");
		static_assert(0.0 <= ParamSet::kEndAngle && ParamSet::kEndAngle <= 720.0, "Internal Error");
		static_assert(ParamSet::kStartAngle < ParamSet::kEndAngle, "Internal Error");
		double const angle_increment = (ParamSet::kEndAngle
				- ParamSet::kStartAngle) / 360.0 * 2.0 * const_pi
				/ static_cast<double>(kNumAngle);
		for (size_t i = 0; i < kNumAngle; ++i) {
			angle[i] = ParamSet::kStartAngle * deg2rad
					+ static_cast<double>(i) * angle_increment;
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
				(*x_out)[k] = r * c;
				(*y_out)[k] = r * s;
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

struct IdenticalValueInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[], bool const mask_in[], double **x_out,
			double **y_out, bool **mask_out, bool mask_expected[],
			sakura_Status *status_expected) {
		// all elements in x and y will be initialized to 0
		BaseInitializer::Initialize(num_data, fraction, x_in, y_in, mask_in,
				x_out, y_out, mask_out, mask_expected, status_expected);
		*status_expected = sakura_Status_kNG;
	}
};

struct NumDataTwoInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[], bool const mask_in[], double **x_out,
			double **y_out, bool **mask_out, bool mask_expected[],
			sakura_Status *status_expected) {
		ASSERT_EQ(num_data, 2);
		BaseInitializer::Initialize(num_data, fraction, x_in, y_in, mask_in,
				x_out, y_out, mask_out, mask_expected, status_expected);
		(*x_out)[0] = -5.0;
		(*x_out)[1] = 3.0;
		(*y_out)[0] = -4.7;
		(*y_out)[1] = 4.3;
		mask_expected[0] = true;
		mask_expected[1] = true;
		*status_expected = sakura_Status_kOK;
	}
};

struct NumDataThreeInitializer {
	static void Initialize(size_t num_data, float fraction, double const x_in[],
			double const y_in[], bool const mask_in[], double **x_out,
			double **y_out, bool **mask_out, bool mask_expected[],
			sakura_Status *status_expected) {
		ASSERT_EQ(num_data, 3);
		BaseInitializer::Initialize(num_data, fraction, x_in, y_in, mask_in,
				x_out, y_out, mask_out, mask_expected, status_expected);
		(*x_out)[0] = -5.0;
		(*x_out)[1] = 3.0;
		(*x_out)[2] = 109.9;
		(*y_out)[0] = -4.7;
		(*y_out)[1] = 4.3;
		(*y_out)[2] = -328.1;
		mask_expected[0] = false;
		mask_expected[1] = false;
		mask_expected[2] = true;
		*status_expected = sakura_Status_kOK;
	}

};

typedef SquareShapeInitializer<StandardSquareParamSet> StandardSquare;
typedef SquareShapeInitializer<WiderSquareParamSet> WiderSquare;
typedef FailedSquareShapeInitializer<StandardSquareParamSet> FailedSquare;
typedef SquareShapeInitializer<StandardSquareParamSet,
		MaskNotProperlyConfiguredInitializer> WrongMaskSquare;
typedef UserDefinedRangeSquareShapeInitializer<WiderSquareParamSet> UserDefinedRangeSquare;
typedef CircularShapeInitializer<StandardCircularParamSet> StandardCircle;
typedef CircularShapeInitializer<DonutShapeParamSet> StandardDonutShape;
typedef CircularShapeInitializer<UpperCShapeParamSet> StandardUpperCShape;
typedef CircularShapeInitializer<HorizontalCShapeParamSet> StandardHorizontalCShape;
typedef CircularShapeInitializer<InclinedCShapeParamSet> StandardInclinedCShape;

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

struct AllMaskedChecker {
	static void Check(size_t num_data, bool const mask[],
	bool const mask_expected[]) {
		for (size_t i = 0; i < num_data; ++i) {
			EXPECT_EQ(true, mask[i]) << i;
		}
	}
};

struct AllUnmaskedChecker {
	static void Check(size_t num_data, bool const mask[],
	bool const mask_expected[]) {
		for (size_t i = 0; i < num_data; ++i) {
			EXPECT_FALSE(mask[i]) << i;
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

// RunTest
template<typename Initializer, typename Checker = NullChecker,
		typename Printer = NullOutput>
void RunTest(size_t num_data, float fraction, double pixel_size,
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

	double start_time = GetCurrentTime();
	sakura_Status status = sakura_CreateMaskNearEdgeDouble(fraction, pixel_size,
			num_data, x, y, blc_x, blc_y, trc_x, trc_y, mask);
	double end_time = GetCurrentTime();
	Printer::Print(test_name, end_time - start_time);

	ASSERT_EQ(status_expected, status);
	Checker::Check(num_data, mask, mask_expected);

	sakura_CleanUp();
}

// FAILURE CASES
TEST_MASK(NagativeFraction) {
	RunTest<BasicInvalidArgumentInitializer>(10, -0.1f, 0.0);
}

TEST_MASK(FractionLargerThanOne) {
	RunTest<BasicInvalidArgumentInitializer>(10, 1.5f, 0.0);
}

TEST_MASK(FractionIsNaN) {
	float const nan_float = WrongValueNaN<float>::Get();
	RunTest<BasicInvalidArgumentInitializer>(10, nan_float, 0.0);
}

TEST_MASK(FractionIsInf) {
	float const inf_float = WrongValuePositiveInf<float>::Get();
	RunTest<BasicInvalidArgumentInitializer>(10, inf_float, 0.0);

	float const ninf_float = WrongValueNegativeInf<float>::Get();
	RunTest<BasicInvalidArgumentInitializer>(10, ninf_float, 0.0);
}

TEST_MASK(NegativePixelScale) {
	RunTest<BasicInvalidArgumentInitializer>(10, 0.1f, -0.5);
}

//TEST_MASK(ZeroPixelScale) {
//	RunTest<BasicInvalidArgumentInitializer>(10, 0.1f, 0.0);
//}

TEST_MASK(PixelScaleIsNaN) {
	double const nan_double = WrongValueNaN<double>::Get();
	RunTest<BasicInvalidArgumentInitializer>(10, 0.1f, nan_double);
}

TEST_MASK(PixelScaleIsInf) {
	double const inf_double = WrongValuePositiveInf<double>::Get();
	RunTest<BasicInvalidArgumentInitializer>(10, 0.1f, inf_double);

	double const ninf_double = WrongValueNegativeInf<double>::Get();
	RunTest<BasicInvalidArgumentInitializer>(10, 0.1f, ninf_double);
}

TEST_MASK(ArrayNotAligned) {
	// x is not aligned
	RunTest<NotAlignedXInitializer>(10, 0.1f, 0.0);

	// y is not aligned
	RunTest<NotAlignedYInitializer>(10, 0.1f, 0.0);

	// mask is not aligned
	RunTest<NotAlignedMaskInitializer>(10, 0.1f, 0.0);
}

TEST_MASK(ArrayIsNull) {
	// x is nullptr
	RunTest<NullXInitializer>(10, 0.1f, 0.0);

	// y is nullptr
	RunTest<NullYInitializer>(10, 0.1f, 0.0);

	// mask is nullptr
	RunTest<NullMaskInitializer>(10, 0.1f, 0.0);
}

TEST_MASK(ArrayHasNaN) {
	// x has NaN
	RunTest<
			FailedSquareShapeInitializerWithWrongValue<StandardSquareParamSet,
					WrongValueNaN<double>, true> >(100, 0.1f, 0.0);

	// y has NaN
	RunTest<
			FailedSquareShapeInitializerWithWrongValue<StandardSquareParamSet,
					WrongValueNaN<double>, false> >(100, 0.1f, 0.0);
}

TEST_MASK(ArrayHasPositiveInf) {
	// x has +Inf
	RunTest<
			FailedSquareShapeInitializerWithWrongValue<StandardSquareParamSet,
					WrongValuePositiveInf<double>, true> >(100, 0.1f, 0.0);

	// y has +Inf
	RunTest<
			FailedSquareShapeInitializerWithWrongValue<StandardSquareParamSet,
					WrongValuePositiveInf<double>, false> >(100, 0.1f, 0.0);

	// x has -Inf
	RunTest<
			FailedSquareShapeInitializerWithWrongValue<StandardSquareParamSet,
					WrongValueNegativeInf<double>, true> >(100, 0.1f, 0.0);

	// y has -Inf
	RunTest<
			FailedSquareShapeInitializerWithWrongValue<StandardSquareParamSet,
					WrongValueNegativeInf<double>, false> >(100, 0.1f, 0.0);
}

TEST_MASK(TooFewEffectiveNumberOfMaskedData) {
	RunTest<BaseInitializer, AllUnmaskedChecker>(5, 0.1f, 0.0);
}

TEST_MASK(NumDataIsOne) {
	// If pixel size is not manually set, it should fail
	RunTest<BasicInvalidArgumentInitializer>(1, 0.1f, 0.0);

	// If blc and trc are not manually set, it should fail in general
	// However, it passes since expected number of masked data
	// (num_data * fraction) is less than 1.
	RunTest<BaseInitializer, AllUnmaskedChecker>(1, 0.1f, 1.0);

	// The following case should fail since num_data * fraction is 1
	RunTest<BasicInvalidArgumentInitializer>(1, 1.0f, 1.0);

	// If both pixel size and (blc, trc) are given, it should work
	// resulting mask will be all false since effective
	// mask ratio 1 * 0.1 = 0.1 is less than 1
	double blc_x = -1.0;
	double blc_y = -1.0;
	double trc_x = 1.0;
	double trc_y = 1.0;
	RunTest<BaseInitializer, AllUnmaskedChecker>(1, 0.1f, 1.0, &blc_x, &blc_y,
			&trc_x, &trc_y);

	// resulting mask will be all true since effective
	// mask ratio 1 * 0.1 = 0.1 is less than 1
	RunTest<BaseInitializer, AllMaskedChecker>(1, 1.0f, 1.0, &blc_x, &blc_y,
			&trc_x, &trc_y);

}

TEST_MASK(InvalidUserDefinedRange) {
	// blc_x == trc_x
	double blc_x = 2.0;
	double trc_x = 2.0;
	RunTest<FailedSquare>(1000, 0.1f, 0.0, &blc_x, NULL, &trc_x, NULL);

	// blc_x > trc_x
	blc_x = 2.0;
	trc_x = -2.0;
	RunTest<FailedSquare>(1000, 0.1f, 0.0, &blc_x, NULL, &trc_x, NULL);

	// blc_y == trc_y
	double blc_y = 2.0;
	double trc_y = 2.0;
	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL, &blc_y, NULL, &trc_y);

	// blc_y > trc_y
	blc_y = 2.0;
	trc_y = -2.0;
	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL, &blc_y, NULL, &trc_y);

	// all data are located at left side of user defined bottom left corner
	blc_x = 11.0;
	RunTest<FailedSquare>(1000, 0.1f, 0.0, &blc_x,
	NULL, NULL, NULL);

	// all data are located at right side of user defined top right corner
	trc_x = -11.0;
	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL,
	NULL, &trc_x, NULL);

	// all data are located at lower side of user defined bottom left corner
	blc_y = 11.0;
	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL, &blc_y, NULL, NULL);

	// all data are located at lower side of user defined top right corner
	trc_y = -11.0;
	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL,
	NULL, NULL, &trc_y);

	// any user defined range is +Inf
	double positive_inf_value = WrongValuePositiveInf<double>::Get();
	RunTest<FailedSquare>(1000, 0.1f, 0.0, &positive_inf_value,
	NULL, NULL, NULL);

	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL,
	NULL, &positive_inf_value, NULL);

	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL, &positive_inf_value, NULL,
	NULL);

	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL,
	NULL, NULL, &positive_inf_value);

	// any user defined range is -Inf
	double negative_inf_value = WrongValueNegativeInf<double>::Get();
	RunTest<FailedSquare>(1000, 0.1f, 0.0, &negative_inf_value,
	NULL, NULL, NULL);

	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL,
	NULL, &negative_inf_value, NULL);

	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL, &negative_inf_value, NULL,
	NULL);

	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL,
	NULL, NULL, &negative_inf_value);

	// any user defined range is NaN
	double nan_value = WrongValueNaN<double>::Get();
	RunTest<FailedSquare>(1000, 0.1f, 0.0, &nan_value,
	NULL, NULL, NULL);

	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL,
	NULL, &nan_value, NULL);

	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL, &nan_value, NULL, NULL);

	RunTest<FailedSquare>(1000, 0.1f, 0.0, NULL,
	NULL, NULL, &nan_value);
}

TEST_MASK(IdenticalValue) {
	// x and y are filled with identical value (0.0)
	RunTest<IdenticalValueInitializer>(1000, 0.1f, 0.0);
}

// SUCCESSFUL CASES
TEST_MASK(ZeroFraction) {
	RunTest<BaseInitializer, StandardChecker>(10, 0.0f, 0.0);
}

TEST_MASK(ZeroDataSize) {
	RunTest<BaseInitializer, StandardChecker>(0, 0.1f, 0.0);
}

TEST_MASK(FractionEqualToOne) {
//	RunTest<StandardSquare, StandardChecker>(1000, 1.0f, 0.0);
	RunTest<StandardSquare, AllMaskedChecker>(1000, 1.0f, 0.0);
}

TEST_MASK(NumDataEqualToTwo) {
	RunTest<NumDataTwoInitializer>(2, 0.1f, 0.0);
}

TEST_MASK(NumDataEqualToThree) {
	RunTest<NumDataThreeInitializer>(3, 0.1f, 0.0);
}

TEST_MASK(FractionTenPercent) {
	float const fraction = 0.1f;
	double const pixel_size = 0.0;

	// square spiral pattern
	RunTest<StandardSquare, StandardChecker>(1000, fraction, pixel_size);

	// circular zigzag pattern
	RunTest<StandardCircle, StandardChecker>(1000, fraction, pixel_size);

	// donut shape
	RunTest<StandardDonutShape, StandardChecker>(1000, fraction, pixel_size);

	// upper C-shape
	RunTest<StandardUpperCShape, StandardChecker>(1000, fraction, pixel_size);

	// horizontal C-shape
	RunTest<StandardHorizontalCShape, StandardChecker>(1000, fraction,
			pixel_size);

	// inclined C-shape
	RunTest<StandardInclinedCShape, StandardChecker>(1000, fraction,
			pixel_size);

//	// test
//	size_t num_data = 1000;
//	float fraction = 0.1f;
//	SIMD_ALIGN
//	double x_in[num_data+1];
//	SIMD_ALIGN
//	double y_in[num_data+1];
//	SIMD_ALIGN
//	bool mask_in[num_data+1];
//	//bool mask_expected[num_data];
//	double *x = x;
//	double *y = y;
//	bool *mask = mask_in;
//	sakura_Status status;
//	std::ofstream ofs;
//	StandardCircle::Initialize(num_data, fraction, x_in, y_in, mask_in, &x, &y, &mask, mask, &status);
//	ofs.open("circule.dat", std::ofstream::out);
//	for (size_t i = 0; i < num_data; ++i) {
//		ofs << x[i] << " " << y[i] << std::endl;
//	}
//	ofs.close();
//
//	StandardDonutShape::Initialize(num_data, fraction, x_in, y_in, mask_in, &x, &y,
//			&mask, mask, &status);
//	ofs.open("donut.dat", std::ofstream::out);
//	for (size_t i = 0; i < num_data; ++i) {
//		ofs << x[i] << " " << y[i] << std::endl;
//	}
//	ofs.close();
//
//	StandardUpperCShape::Initialize(num_data, fraction, x_in, y_in, mask_in, &x, &y,
//			&mask, mask, &status);
//	ofs.open("cupper.dat", std::ofstream::out);
//	for (size_t i = 0; i < num_data; ++i) {
//		ofs << x[i] << " " << y[i] << std::endl;
//	}
//	ofs.close();
//
//	StandardHorizontalCShape::Initialize(num_data, fraction, x_in, y_in, mask_in, &x, &y,
//			&mask, mask, &status);
//	ofs.open("chorizontal.dat", std::ofstream::out);
//	for (size_t i = 0; i < num_data; ++i) {
//		ofs << x[i] << " " << y[i] << std::endl;
//	}
//	ofs.close();
//
//	StandardInclinedCShape::Initialize(num_data, fraction, x_in, y_in, mask_in, &x, &y,
//			&mask, mask, &status);
//	ofs.open("cinclined.dat", std::ofstream::out);
//	for (size_t i = 0; i < num_data; ++i) {
//		ofs << x[i] << " " << y[i] << std::endl;
//	}
//	ofs.close();
}

TEST_MASK(FractionTwentyPercent) {
	float const fraction = 0.2f;
	double const pixel_size = 0.0;

	// square spiral pattern
	RunTest<StandardSquare, StandardChecker>(1000, fraction, pixel_size);

	// circular zigzag pattern
	RunTest<StandardCircle, StandardChecker>(1000, fraction, pixel_size);

	// donut shape
	RunTest<StandardDonutShape, StandardChecker>(1000, fraction, pixel_size);

	// upper C-shape
	RunTest<StandardUpperCShape, StandardChecker>(1000, fraction, pixel_size);

	// horizontal C-shape
	RunTest<StandardHorizontalCShape, StandardChecker>(1000, fraction,
			pixel_size);

	// inclined C-shape
	RunTest<StandardInclinedCShape, StandardChecker>(1000, fraction,
			pixel_size);
}

TEST_MASK(OutputMaskNotInitialized) {
	float const fraction = 0.1f;
	double const pixel_size = 0.0;

	// square spiral pattern
	RunTest<WrongMaskSquare, StandardChecker>(1000, fraction, pixel_size);
}

TEST_MASK(ManualPixelSize) {
	float const fraction = 0.1f;
	double pixel_size = 0.08;

	// square spiral pattern
	RunTest<StandardSquare, StandardChecker>(1000, fraction, pixel_size);

	pixel_size = 0.1;

	// circular zigzag pattern
	RunTest<StandardCircle, StandardChecker>(1000, fraction, pixel_size);

	// donut shape
	RunTest<StandardDonutShape, StandardChecker>(1000, fraction, pixel_size);

	// upper C-shape
	RunTest<StandardUpperCShape, StandardChecker>(1000, fraction, pixel_size);

	// horizontal C-shape
	RunTest<StandardHorizontalCShape, StandardChecker>(1000, fraction,
			pixel_size);

	// inclined C-shape
	RunTest<StandardInclinedCShape, StandardChecker>(1000, fraction,
			pixel_size);

	// large pixel size causes to mask all data regardless of fraction
	pixel_size = 1.0e10;
	RunTest<StandardCircle, AllMaskedChecker>(1000, fraction, pixel_size);
}

TEST_MASK(ValidUserDefinedRange) {
	double blc_x = -5.0;
	double blc_y = -5.0;
	double trc_x = 5.0;
	double trc_y = 5.0;
	RunTest<UserDefinedRangeSquare, StandardChecker, NullOutput>(1000, 0.1f,
			0.0, &blc_x, &blc_y, &trc_x, &trc_y);
}

// PERFORMANCE TEST
TEST_MASK(LargeNumDataPerformance) {
	RunTest<StandardSquare, StandardChecker, PerformanceTestOutput>(40000, 0.1f,
			0.0, NULL, NULL, NULL, NULL,
			"MaskDataNearEdge_LargeNumDataPerformance");
}

TEST_MASK(LargeAreaPerformance) {
	RunTest<WiderSquare, StandardChecker, PerformanceTestOutput>(1000, 0.1f,
			0.0, NULL, NULL, NULL, NULL,
			"MaskDataNearEdge_LargeAreaPerformance");
}

