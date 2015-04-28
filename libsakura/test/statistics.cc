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
/*
 * statistics.cc
 *
 *  Created on: 2013/08/20
 *      Author: kohji
 */

#include <sys/time.h>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <array>
#include <algorithm>
#include <numeric>

#include <libsakura/sakura.h>

#include <libsakura/localdef.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"

namespace {
double Rms2(size_t elements, float const data[], bool const is_valid[]) {
	double result = 0;
	for (size_t i = 0; i < elements; ++i) {
		if (is_valid[i]) {
			double v = data[i];
			result += v * v;
		}
	}
	return result;
}

template<typename T>
void ExpectEQ(T const &ref, T const &result) {
	if (std::isnan(ref)) {
		EXPECT_TRUE(isnan(result));
	} else {
		EXPECT_EQ(ref, result);
	}
}

template<typename T>
void ExpectRealEQ(T const &ref, T const &result) {
	assert(false);
}

template<>
void ExpectRealEQ<double>(double const &ref, double const &result) {
	if (std::isnan(ref)) {
		EXPECT_TRUE(isnan(result));
	} else {
		EXPECT_DOUBLE_EQ(ref, result);
	}
}

void TestResult(LIBSAKURA_SYMBOL(StatisticsResultFloat) const &ref,
LIBSAKURA_SYMBOL(StatisticsResultFloat) const &result) {
	EXPECT_EQ(ref.count, result.count);
	EXPECT_EQ(ref.index_of_max, result.index_of_max);
	EXPECT_EQ(ref.index_of_min, result.index_of_min);
	ExpectEQ(ref.max, result.max);
	ExpectEQ(ref.min, result.min);
	//std::cout << "sum\n";
	ExpectRealEQ(ref.sum, result.sum);
	//std::cout << "mean\n";
	ExpectRealEQ(ref.mean, result.mean);
	//std::cout << "rms\n";
	ExpectRealEQ(ref.rms, result.rms);
	//std::cout << "stddev\n";
	ExpectRealEQ(ref.stddev, result.stddev);
}

template<typename T>
void DisplayDiff(char const str[], T const &ref, T const &fast,
		T const &accurate) {
	if (ref != 0) {
		std::cout << str << ": fast " << std::setprecision(7)
				<< std::abs(ref - fast) * 100 / ref << "%, accurate "
				<< std::abs(ref - accurate) * 100 / ref << "%\n";
	}
}

template<typename CompareFast, typename CompareAccurate>
bool CallAndTestResult(LIBSAKURA_SYMBOL(StatisticsResultFloat) const &ref,
		size_t num_data, float const data[], bool const is_valid[],
		CompareFast compare_fast, CompareAccurate compare_accurate) {
	//std::cout << "Fast\n";
	LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
	{
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(ComputeStatisticsFloat)(num_data, data, is_valid,
				&result);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		compare_fast(ref, result);

	}
	//std::cout << "Accurate\n";
	LIBSAKURA_SYMBOL(StatisticsResultFloat) result_accurate;
	{
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(ComputeAccurateStatisticsFloat)(num_data, data,
				is_valid, &result_accurate);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		compare_accurate(ref, result_accurate);
	}

	auto min_ptr = &decltype(result)::min;
	auto min_index_ptr = &decltype(result)::index_of_min;
	auto min_max_test =
			[&data, &result, &result_accurate](decltype(min_ptr) attr, decltype(min_index_ptr) index) {
				if (isnan(result.*attr)) {
					EXPECT_TRUE(isnan(result_accurate.*attr));
					EXPECT_EQ(-1, result.*index);
					EXPECT_EQ(-1, result_accurate.*index);
				} else {
					EXPECT_FALSE(isnan(result_accurate.*attr));
					EXPECT_EQ(result.*attr, result_accurate.*attr);
					EXPECT_LE(0, result.*index);
					EXPECT_LE(0, result_accurate.*index);
					EXPECT_EQ(result.*attr, data[result.*index]);
					EXPECT_EQ(data[result.*index], data[result_accurate.*index]);
				}
				if (result.*index < 0) {
					EXPECT_EQ(result.*index, result_accurate.*index);
				} else {
					EXPECT_LE(0, result.*index);
					EXPECT_EQ(result.*attr, data[result.*index]);
					EXPECT_EQ(data[result.*index], data[result_accurate.*index]);
				}
			};
	min_max_test(min_ptr, min_index_ptr);
	auto max_ptr = &decltype(result)::max;
	auto max_index_ptr = &decltype(result)::index_of_max;
	min_max_test(max_ptr, max_index_ptr);
	EXPECT_EQ(result.count, result_accurate.count);
	auto value_test =
			[&result, &result_accurate](decltype(&decltype(result)::sum) attr) {
				if (isnan(result.*attr)) {
					EXPECT_TRUE(isnan(result_accurate.*attr));
				} else {
					EXPECT_TRUE(!isnan(result_accurate.*attr));
				}
			};
	value_test(&decltype(result)::sum);
	value_test(&decltype(result)::rms);
	value_test(&decltype(result)::stddev);

	if ((!isnan(result.sum) && !isnan(result_accurate.sum)
			&& result.sum != result_accurate.sum)
			|| (!isnan(result.rms) && !isnan(result_accurate.rms)
					&& result.rms != result_accurate.rms)
			|| (!isnan(result.stddev) && !isnan(result_accurate.stddev)
					&& result.stddev != result_accurate.stddev)) {
		if (false) { // for debug
			std::cout << "error ratio: \n";
			DisplayDiff("sum", ref.sum, result.sum, result_accurate.sum);
			DisplayDiff("rms", ref.rms, result.rms, result_accurate.rms);
			DisplayDiff("stddev", ref.stddev, result.stddev,
					result_accurate.stddev);
		}
		return true;
	}
	return false;
}

void TestResultNop(LIBSAKURA_SYMBOL(StatisticsResultFloat) const &ref,
LIBSAKURA_SYMBOL(StatisticsResultFloat) const &result) {
}

bool CallAndIsDifferent(size_t num_data, float const data[],
bool const is_valid[]) {
	static LIBSAKURA_SYMBOL(StatisticsResultFloat) const ref = { 0 };
	return CallAndTestResult(ref, num_data, data, is_valid, TestResultNop,
			TestResultNop);
}

}

TEST(Statistics, SortValidValuesDensely) {
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
			nullptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

	{
		SIMD_ALIGN
		static float data[] = { 2.f, -2.f, 3.f, 0.f, -2.f, -1.f };
		SIMD_ALIGN
		static bool const is_valid[] = { true, true, true, true, true, true };
		STATIC_ASSERT(ELEMENTSOF(data) == ELEMENTSOF(is_valid));
		size_t new_elements = size_t(-1);
		result = LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(ELEMENTSOF(data),
				is_valid, data, &new_elements);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
		EXPECT_EQ(ELEMENTSOF(data), new_elements);
		EXPECT_EQ(-2.f, data[0]);
		EXPECT_EQ(-2.f, data[1]);
		EXPECT_EQ(-1.f, data[2]);
		EXPECT_EQ(0.f, data[3]);
		EXPECT_EQ(2.f, data[4]);
		EXPECT_EQ(3.f, data[5]);

		new_elements = size_t(-1);
		result = LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(0, is_valid,
				data, &new_elements);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
		EXPECT_EQ(0, new_elements);
	}
	{
		SIMD_ALIGN
		static float data[] = { 2.f, -2.f, 3.f, 0.f, -2.f, -1.f };
		SIMD_ALIGN
		static bool const is_valid[] = { true, false, true, false, true, true };
		STATIC_ASSERT(ELEMENTSOF(data) == ELEMENTSOF(is_valid));
		size_t new_elements = size_t(-1);
		result = LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(ELEMENTSOF(data),
				is_valid, data, &new_elements);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
		EXPECT_EQ(ELEMENTSOF(data) - 2, new_elements);
		EXPECT_EQ(-2.f, data[0]);
		EXPECT_EQ(-1.f, data[1]);
		EXPECT_EQ(2.f, data[2]);
		EXPECT_EQ(3.f, data[3]);
	}
	{
		SIMD_ALIGN
		static float data[] = { 2.f, -2.f, 3.f, 0.f, -2.f, -1.f };
		SIMD_ALIGN
		static bool const is_valid[] = { false, false, false, false, false,
		false };
		STATIC_ASSERT(ELEMENTSOF(data) == ELEMENTSOF(is_valid));
		size_t new_elements = size_t(-1);
		result = LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(ELEMENTSOF(data),
				is_valid, data, &new_elements);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
		EXPECT_EQ(0, new_elements);
	}
	{
		SIMD_ALIGN
		static float data[] = { 2.f, -2.f, 3.f, 0.f, -2.f, -1.f };
		size_t new_elements = size_t(-1);
		result = LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(ELEMENTSOF(data),
				nullptr, data, &new_elements);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
		EXPECT_EQ(size_t(-1), new_elements);
	}
	{
		SIMD_ALIGN
		static float data[] = { 2.f, -2.f, 3.f, 0.f, -2.f, -1.f };
		SIMD_ALIGN
		static bool const is_valid[] = { true, false, true, false, true, true };
		size_t new_elements = size_t(-1);
		result = LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(
		ELEMENTSOF(data) - 1, &is_valid[1], &data[1], &new_elements);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
		EXPECT_EQ(-2.f, data[1]);
		EXPECT_EQ(-1.f, data[2]);
		EXPECT_EQ(3.f, data[3]);
	}
	LIBSAKURA_SYMBOL(CleanUp)();
}

TEST(Statistics, ComputeStatistics) {
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
			nullptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

	{
		SIMD_ALIGN
		static std::array<float, 256> data;
		SIMD_ALIGN
		static std::array<bool, data.size()> is_valid;
		{
			is_valid.fill(true);
			std::iota(data.begin(), data.end(), 0);
			assert(data.size() % 2 == 0);
			LIBSAKURA_SYMBOL(StatisticsResultFloat) ref;
			ref.count = data.size();
			ref.index_of_min = 0;
			ref.min = 0;
			ref.index_of_max = data.size() - 1;
			ref.max = data.size() - 1;
			ref.sum = (data.size() - 1) * (data.size() / 2);
			ref.mean = (data.size() - 1) / 2.;
			auto rms2 = Rms2(data.size(), data.data(), is_valid.data())
					/ data.size();
			ref.rms = std::sqrt(rms2);
			ref.stddev = std::sqrt(std::abs(rms2 - ref.mean * ref.mean));

			CallAndTestResult(ref, data.size(), data.data(), is_valid.data(),
					TestResult, TestResult);
		}
		size_t min_max_index[16][2];
		for (size_t i = 0; i < ELEMENTSOF(min_max_index); ++i) {
			min_max_index[i][0] = 64 + i;
			min_max_index[i][1] = 127 - i;
		}
		std::fill(std::begin(data), std::end(data), 0.f);
		assert(data.size() % 2 == 0);
		std::for_each(std::begin(min_max_index), std::end(min_max_index),
				[](size_t const (&idx)[2]) {
					size_t max_index = idx[0];
					size_t min_index = idx[1];
					data[max_index] = 20.f;
					data[min_index] = -120.f;

					LIBSAKURA_SYMBOL(StatisticsResultFloat) ref;
					ref.count = data.size();
					ref.index_of_min = min_index;
					ref.min = data[min_index];
					ref.index_of_max = max_index;
					ref.max = data[max_index];
					ref.sum = data[max_index] + data[min_index];
					ref.mean = double(data[max_index] + data[min_index]) / data.size();

					double v = data[max_index];
					auto rms2 = v * v;
					v = data[min_index];
					rms2 += v * v;
					rms2 /= data.size();
					ref.rms = std::sqrt(rms2);
					ref.stddev = std::sqrt(std::abs(rms2 - ref.mean * ref.mean));

					CallAndTestResult(ref, data.size(), data.data(), is_valid.data(), TestResult, TestResult);

					data[max_index] = 0;
					data[min_index] = 0;
				});
	}
	{
		SIMD_ALIGN
		static std::array<float, 4> data;
		SIMD_ALIGN
		static std::array<bool, data.size()> is_valid;
		is_valid.fill(true);
		std::iota(data.begin(), data.end(), 0);

		LIBSAKURA_SYMBOL(StatisticsResultFloat) ref;
		ref.count = data.size();
		ref.index_of_min = 0;
		ref.min = 0;
		ref.index_of_max = data.size() - 1;
		ref.max = data.size() - 1;
		ref.sum = (data.size() - 1) * (data.size() / 2);
		ref.mean = (data.size() - 1) / 2.;
		auto rms2 = Rms2(data.size(), data.data(), is_valid.data())
				/ data.size();
		ref.rms = std::sqrt(rms2);
		ref.stddev = std::sqrt(std::abs(rms2 - ref.mean * ref.mean));

		CallAndTestResult(ref, data.size(), data.data(), is_valid.data(),
				TestResult, TestResult);
	}
	{
		SIMD_ALIGN
		static std::array<float, 1> data = { { 3.f } };
		SIMD_ALIGN
		static std::array<bool, data.size()> is_valid = { { true } };

		LIBSAKURA_SYMBOL(StatisticsResultFloat) ref;
		ref.count = data.size();
		ref.index_of_min = 0;
		ref.min = data[0];
		ref.index_of_max = data.size() - 1;
		ref.max = data[0];
		ref.sum = data[0];
		ref.mean = data[0];
		ref.rms = data[0];
		ref.stddev = 0.;

		CallAndTestResult(ref, data.size(), data.data(), is_valid.data(),
				TestResult, TestResult);
	}
	{
		SIMD_ALIGN
		static float data[0];
		SIMD_ALIGN
		static bool is_valid[ELEMENTSOF(data)];

		LIBSAKURA_SYMBOL(StatisticsResultFloat) ref;
		ref.count = ELEMENTSOF(data);
		ref.index_of_min = -1;
		ref.min = NAN;
		ref.index_of_max = -1;
		ref.max = NAN;
		ref.sum = 0;
		ref.mean = NAN;
		ref.rms = NAN;
		ref.stddev = NAN;

		CallAndTestResult(ref, ELEMENTSOF(data), data, is_valid, TestResult,
				TestResult);
	}
	{
		SIMD_ALIGN
		static std::array<float, 1024> data;
		SIMD_ALIGN
		static std::array<bool, data.size()> is_valid;
		is_valid.fill(false);
		std::iota(data.begin(), data.end(), 0);

		LIBSAKURA_SYMBOL(StatisticsResultFloat) ref;
		ref.count = 0;
		ref.index_of_min = -1;
		ref.min = NAN;
		ref.index_of_max = -1;
		ref.max = NAN;
		ref.sum = 0;
		ref.mean = NAN;
		ref.rms = NAN;
		ref.stddev = NAN;

		CallAndTestResult(ref, data.size(), data.data(), is_valid.data(),
				TestResult, TestResult);

		for (size_t end = data.size() - 128; end < data.size(); ++end) {
			for (size_t pos = end - 128; pos < end; ++pos) {
				is_valid[pos] = true;

				ref.count = 1;
				ref.index_of_min = pos;
				ref.min = data[pos];
				ref.index_of_max = pos;
				ref.max = data[pos];
				ref.sum = data[pos];
				ref.mean = data[pos];
				ref.rms = data[pos];
				ref.stddev = 0.;

				CallAndTestResult(ref, end, data.data(), is_valid.data(),
						TestResult, TestResult);

				is_valid[pos] = false;
			}
		}
	}
	{
		SIMD_ALIGN
		static std::array<float, 1024> data;
		SIMD_ALIGN
		static std::array<bool, data.size()> is_valid;
		size_t count = 0;
		double sum = 0;
		float max = -1;
		float min = data.size() + 1;
		size_t max_idx = -1;
		size_t min_idx = -1;
		double sqsum = 0.;
		std::iota(data.begin(), data.end(), 0);
		for (size_t i = 0; i < is_valid.size(); ++i) {
			is_valid[i] = i % 7 == 0;
			if (is_valid[i]) {
				++count;
				double v = data[i];
				sum += v;
				sqsum += v * v;
				if (max < data[i]) {
					max = data[i];
					max_idx = i;
				}
				if (min > data[i]) {
					min = data[i];
					min_idx = i;
				}
			}
		}
		auto rms2 = sqsum / count;
		auto mean = sum / count;

		LIBSAKURA_SYMBOL(StatisticsResultFloat) ref;
		ref.count = count;
		ref.index_of_min = min_idx;
		ref.min = min;
		ref.index_of_max = max_idx;
		ref.max = max;
		ref.sum = sum;
		ref.mean = mean;
		ref.rms = std::sqrt(rms2);
		ref.stddev = std::sqrt(std::abs(rms2 - mean * mean));

		CallAndTestResult(ref, data.size(), data.data(), is_valid.data(),
				TestResult, TestResult);

		is_valid.fill(false);
		sum = 0.;
		sqsum = 0.;
		count = 0;
		size_t base = data.size() - 16;
		data[base + 5] = data[base + 2] - 2;
		auto put = [&](size_t offset) {
			is_valid[base + offset] = true;
			++count;
			sum += data[base + offset];
			double v = data[base + offset];
			sqsum += v * v;
		};
		put(2);
		put(3);
		put(5);

		rms2 = sqsum / count;
		mean = sum / count;

		ref.count = count;
		ref.index_of_min = base + 5;
		ref.min = data[base + 5];
		ref.index_of_max = base + 3;
		ref.max = data[base + 3];
		ref.sum = sum;
		ref.mean = mean;
		ref.rms = std::sqrt(rms2);
		ref.stddev = std::sqrt(std::abs(rms2 - mean * mean));

		CallAndTestResult(ref, base + 7, data.data(), is_valid.data(),
				TestResult, TestResult);
	}
	{
		SIMD_ALIGN
		static std::array<float, 256> data;
		SIMD_ALIGN
		static std::array<bool, data.size()> is_valid;
		LIBSAKURA_SYMBOL(StatisticsResultFloat) ref;
		ref.min = 1.;
		ref.max = -1.;

		auto compare_min =
				[](LIBSAKURA_SYMBOL(StatisticsResultFloat) const &ref,
						LIBSAKURA_SYMBOL(StatisticsResultFloat) const &result) {
					EXPECT_EQ(ref.index_of_min, result.index_of_min);
					EXPECT_EQ(ref.min, result.min);
				};
		auto compare_max =
				[](LIBSAKURA_SYMBOL(StatisticsResultFloat) const &ref,
						LIBSAKURA_SYMBOL(StatisticsResultFloat) const &result) {
					EXPECT_EQ(ref.index_of_max, result.index_of_max);
					EXPECT_EQ(ref.max, result.max);
				};
		is_valid.fill(true);
		for (size_t end = data.size() / 2; end < data.size(); ++end) {
			for (size_t peak = 0; peak < end; ++peak) {
				for (size_t i = 0; i < end; ++i) {
					float v = i;
					data[i] = (v - peak) * (v - peak) + 1;
				}
				ref.index_of_min = peak;
				CallAndTestResult(ref, end, data.data(), is_valid.data(),
						compare_min, compare_min);

				for (size_t i = 0; i < end; ++i) {
					float v = i;
					data[i] = -(v - peak) * (v - peak) - 1;
				}
				ref.index_of_max = peak;
				CallAndTestResult(ref, end, data.data(), is_valid.data(),
						compare_max, compare_max);
			}
		}
	}
	LIBSAKURA_SYMBOL(CleanUp)();
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeStddevFloat)(
		size_t count, double mean, size_t elements, float const data[],
		bool const is_valid[], double *result);

TEST(Statistics, ComputeStatistics_Accuracy) {
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
			nullptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

	SIMD_ALIGN
	static std::array<float, 8192LU * 10000> data;
	SIMD_ALIGN
	static std::array<bool, data.size()> is_valid;
	{
		for (size_t i = 0; i < data.size(); ++i) {
			data[i] = i - float(data.size() + 8);
			is_valid[i] = i % 2 != 0;
		}
		LIBSAKURA_SYMBOL(StatisticsResultFloat) ref;
		assert(data.size() % 4 == 0);
		ref.count = data.size() / 2;
		ref.index_of_min = 1;
		ref.min = data[ref.index_of_min];
		ref.index_of_max = data.size() - 1 - 2;
		ref.max = data[ref.index_of_max];
		ref.sum = -1677721927680000.; // double(ref.max) + double(ref.min) * ref.count / 2.;
		ref.mean = -40960008; // (double(ref.max) + double(ref.min)) / 2.;
		ref.rms = 47296540.9802176008221133754;
		ref.stddev = 23648267.02600718885921762;

		auto compare = [](LIBSAKURA_SYMBOL(StatisticsResultFloat) const &ref,
				LIBSAKURA_SYMBOL(StatisticsResultFloat) const &result) {
			EXPECT_EQ(ref.count, result.count);
			EXPECT_LE(ref.index_of_max, result.index_of_max);
			EXPECT_EQ(ref.index_of_min, result.index_of_min);
			ExpectEQ(ref.max, result.max);
			ExpectEQ(ref.min, result.min);
			ExpectRealEQ(ref.sum, result.sum);
			ExpectRealEQ(ref.mean, result.mean);
			EXPECT_NEAR(ref.rms, result.rms, 7.283e-05);
			EXPECT_NEAR(ref.stddev, result.stddev, 0.00014565);
		};
		auto compare_accurate =
				[](LIBSAKURA_SYMBOL(StatisticsResultFloat) const &ref,
						LIBSAKURA_SYMBOL(StatisticsResultFloat) const &result) {
					EXPECT_EQ(ref.count, result.count);
					EXPECT_LE(ref.index_of_max, result.index_of_max);
					EXPECT_EQ(ref.index_of_min, result.index_of_min);
					ExpectEQ(ref.max, result.max);
					ExpectEQ(ref.min, result.min);
					ExpectRealEQ(ref.sum, result.sum);
					ExpectRealEQ(ref.mean, result.mean);
					ExpectRealEQ(ref.rms, result.rms);
					ExpectRealEQ(ref.stddev, result.stddev);
				};

		CallAndTestResult(ref, data.size(), data.data(), is_valid.data(),
				compare, compare_accurate);
	}
	{
		STATIC_ASSERT(data.size() % 4 == 0);
		constexpr float base = 9.0e+5f;
		for (size_t i = 0; i < data.size(); i += 4) {
			data[i] = base + 4;
			data[i + 1] = base + 7;
			data[i + 2] = base + 13;
			data[i + 3] = base + 16;
			is_valid[i] = true;
			is_valid[i + 1] = true;
			is_valid[i + 2] = true;
			is_valid[i + 3] = true;
		}

		LIBSAKURA_SYMBOL(StatisticsResultFloat) ref;
		assert(data.size() % 4 == 0);
		ref.count = data.size();
		ref.index_of_min = 0;
		ref.min = data[ref.index_of_min];
		ref.index_of_max = 3;
		ref.max = data[ref.index_of_max];
		ref.mean = base + 10;
		ref.sum = ref.mean * data.size();
		constexpr double rms2 = (double(base + 4) * double(base + 4)
				+ double(base + 7) * double(base + 7)
				+ double(base + 13) * double(base + 13)
				+ double(base + 16) * double(base + 16)) / 4.;
		ref.rms = std::sqrt(rms2);
		constexpr double variance = ((4. - 10.) * (4. - 10.)
				+ (7. - 10.) * (7. - 10.) + (13. - 10.) * (13. - 10.)
				+ (16. - 10.) * (16. - 10.)) / 4.;
		ref.stddev = std::sqrt(variance);

		auto compare = [](LIBSAKURA_SYMBOL(StatisticsResultFloat) const &ref,
				LIBSAKURA_SYMBOL(StatisticsResultFloat) const &result) {
			EXPECT_EQ(ref.count, result.count);
			EXPECT_TRUE(ref.index_of_max == result.index_of_max % 4);
			EXPECT_TRUE(ref.index_of_min == result.index_of_min % 4);
			ExpectEQ(ref.max, result.max);
			ExpectEQ(ref.min, result.min);
			ExpectRealEQ(ref.sum, result.sum);
			ExpectRealEQ(ref.mean, result.mean);
			EXPECT_NEAR(ref.rms, result.rms, 0.0001332);
			EXPECT_NEAR(ref.stddev, result.stddev, 11.4505);
		};
		auto compare_accurate =
				[](LIBSAKURA_SYMBOL(StatisticsResultFloat) const &ref,
						LIBSAKURA_SYMBOL(StatisticsResultFloat) const &result) {
					EXPECT_EQ(ref.count, result.count);
					EXPECT_TRUE(ref.index_of_max == result.index_of_max % 4);
					EXPECT_TRUE(ref.index_of_min == result.index_of_min % 4);
					ExpectEQ(ref.max, result.max);
					ExpectEQ(ref.min, result.min);
					ExpectRealEQ(ref.sum, result.sum);
					ExpectRealEQ(ref.mean, result.mean);
					ExpectRealEQ(ref.rms, result.rms);
					ExpectRealEQ(ref.stddev, result.stddev);
				};

		CallAndTestResult(ref, data.size(), data.data(), is_valid.data(),
				compare, compare_accurate);

		constexpr size_t start = 10000;
		for (auto i = start; i < data.size(); ++i) {
			if (CallAndIsDifferent(i, data.data(), is_valid.data())) {
				//std::cout << "ERROR index: " << i << std::endl;
				EXPECT_LT(start, i); // if this fails, make 'start' more smaller value.
				EXPECT_LE(11000, i);
				break;
			}
		}

		LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
		{
			LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL (ComputeAccurateStatisticsFloat)(data.size(),
					data.data(), is_valid.data(), &result);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		}

		double stddev = 0;
		LIBSAKURA_SYMBOL(ComputeStddevFloat)(result.count, result.mean,
				data.size(), data.data(), is_valid.data(), &stddev);
		//std::cout << std::setprecision(16) << std::sqrt(variance) << std::endl << stddev << std::endl << result.stddev << std::endl;
		EXPECT_DOUBLE_EQ(std::sqrt(variance), stddev);
	}
	{
		size_t elements = data.size() - 4 + 1;
		size_t spike_pos = elements - 1;
		for (size_t i = 0; i < spike_pos; i += 4) {
			data[i] = -3;
			data[i + 1] = -7;
			data[i + 2] = 3;
			data[i + 3] = 7;
			is_valid[i] = true;
			is_valid[i + 1] = true;
			is_valid[i + 2] = true;
			is_valid[i + 3] = true;
		}
		constexpr float spike = 9.0e+8f;
		data[spike_pos] = spike;

		LIBSAKURA_SYMBOL(StatisticsResultFloat) ref;
		ref.count = elements;
		ref.index_of_min = 1;
		ref.min = data[ref.index_of_min];
		ref.index_of_max = spike_pos;
		ref.max = spike;
		ref.sum = spike;
		ref.mean = ref.sum / ref.count;
		double rms2 = (double(data[0]) * double(data[0])
				+ double(data[1]) * double(data[1])
				+ double(data[2]) * double(data[2])
				+ double(data[3]) * double(data[3])) * ((elements - 1) / 4);
		rms2 += double(spike) * double(spike);
		rms2 /= elements;
		ref.rms = std::sqrt(rms2);
		double variance = (double(data[0] - ref.mean)
				* double(data[0] - ref.mean)
				+ double(data[1] - ref.mean) * double(data[1] - ref.mean)
				+ double(data[2] - ref.mean) * double(data[2] - ref.mean)
				+ double(data[3] - ref.mean) * double(data[3] - ref.mean))
				* ((elements - 1) / 4);
		variance += (spike - ref.mean) * (spike - ref.mean);
		variance /= elements;
		ref.stddev = std::sqrt(variance);

		auto compare = [](LIBSAKURA_SYMBOL(StatisticsResultFloat) const &ref,
				LIBSAKURA_SYMBOL(StatisticsResultFloat) const &result) {
			EXPECT_EQ(ref.count, result.count);
			EXPECT_EQ(ref.index_of_max, result.index_of_max);
			EXPECT_TRUE(ref.index_of_min == result.index_of_min % 4);
			ExpectEQ(ref.max, result.max);
			ExpectEQ(ref.min, result.min);
			ExpectRealEQ(ref.sum, result.sum);
			ExpectRealEQ(ref.mean, result.mean);
			EXPECT_NEAR(ref.rms, result.rms, 0.0001459);
			EXPECT_NEAR(ref.stddev, result.stddev, 0.0001459);
		};
		auto compare_accurate =
				[](LIBSAKURA_SYMBOL(StatisticsResultFloat) const &ref,
						LIBSAKURA_SYMBOL(StatisticsResultFloat) const &result) {
					EXPECT_EQ(ref.count, result.count);
					EXPECT_EQ(ref.index_of_max, result.index_of_max);
					EXPECT_TRUE(ref.index_of_min == result.index_of_min % 4);
					ExpectEQ(ref.max, result.max);
					ExpectEQ(ref.min, result.min);
					ExpectRealEQ(ref.sum, result.sum);
					ExpectRealEQ(ref.mean, result.mean);
					EXPECT_NEAR(ref.rms, result.rms, 1.6008e-10);
					EXPECT_NEAR(ref.stddev, result.stddev, 1.6008e-10);
				};

		CallAndTestResult(ref, elements, data.data(), is_valid.data(), compare,
				compare_accurate);

		std::swap(data[0], data[spike_pos]);
		spike_pos = 0;
		ref.index_of_max = spike_pos;
		CallAndTestResult(ref, elements, data.data(), is_valid.data(), compare,
				compare_accurate);

		size_t pos = elements / 2;
		std::swap(data[pos], data[spike_pos]);
		spike_pos = pos;
		ref.index_of_max = spike_pos;
		CallAndTestResult(ref, elements, data.data(), is_valid.data(), compare,
				compare_accurate);

		pos = elements / 2 + 13;
		std::swap(data[pos], data[spike_pos]);
		spike_pos = pos;
		ref.index_of_max = spike_pos;
		CallAndTestResult(ref, elements, data.data(), is_valid.data(), compare,
				compare_accurate);

		pos = elements / 4 - 13;
		std::swap(data[pos], data[spike_pos]);
		spike_pos = pos;
		ref.index_of_max = spike_pos;
		CallAndTestResult(ref, elements, data.data(), is_valid.data(), compare,
				compare_accurate);

		LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
		{
			LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL (ComputeAccurateStatisticsFloat)(elements,
					data.data(), is_valid.data(), &result);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		}

		double stddev = 0;
		LIBSAKURA_SYMBOL(ComputeStddevFloat)(result.count, result.mean,
				elements, data.data(), is_valid.data(), &stddev);
		//std::cout << std::setprecision(16) << std::sqrt(variance) << std::endl << stddev << std::endl << result.stddev << std::endl;
		//EXPECT_DOUBLE_EQ(std::sqrt(variance), stddev);
		EXPECT_NEAR(std::sqrt(variance), stddev, 3.89e-05);
	}
	LIBSAKURA_SYMBOL(CleanUp)();
}

namespace {
template<typename T, typename FuncType>
void TestComputeMedianAbsoluteDeviation(FuncType target) {
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
			nullptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

	{
		SIMD_ALIGN
		static T data[] = { T(-3), T(-2), T(-1), T(0), T(2), T(3) };
		SIMD_ALIGN
		T new_data[ELEMENTSOF(data)];

		{ // even
			result = target(ELEMENTSOF(data), data, new_data);
			constexpr auto median = (T(-1) + T(0)) / 2;
			auto const test = [&]() {
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
				EXPECT_EQ(std::abs(T(0) - median), new_data[0]);
				EXPECT_EQ(std::abs(T(-1) - median), new_data[1]);
				EXPECT_EQ(std::abs(T(0) - median), new_data[1]);
				EXPECT_EQ(std::abs(T(-1) - median), new_data[0]);
				EXPECT_EQ(std::abs(T(-2) - median), new_data[2]);
				EXPECT_EQ(std::abs(T(2) - median), new_data[3]);
				EXPECT_EQ(std::abs(T(-3) - median), new_data[4]);
				EXPECT_EQ(std::abs(T(2) - median), new_data[4]);
				EXPECT_EQ(std::abs(T(-3) - median), new_data[3]);
				EXPECT_EQ(std::abs(T(3) - median), new_data[5]);
			};
			test();

			std::copy_n(data, ELEMENTSOF(data), new_data);
			result = target(ELEMENTSOF(data), new_data, new_data);
			test();
		}

		{ // odd
			result = target(ELEMENTSOF(data) - 1, data, new_data);

			constexpr auto median = T(-1);
			auto const test = [&]() {
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
				EXPECT_EQ(std::abs(T(-1) - median), new_data[0]);
				EXPECT_EQ(std::abs(T(-2) - median), new_data[1]);
				EXPECT_EQ(std::abs(T(0) - median), new_data[2]);
				EXPECT_EQ(std::abs(T(-2) - median), new_data[2]);
				EXPECT_EQ(std::abs(T(0) - median), new_data[1]);
				EXPECT_EQ(std::abs(T(-3) - median), new_data[3]);
				EXPECT_EQ(std::abs(T(2) - median), new_data[4]);
			};
			test();

			std::copy_n(data, ELEMENTSOF(data), new_data);
			result = target(ELEMENTSOF(data) - 1, new_data, new_data);
			test();
		}
		{ // 0
			result = target(0, data, new_data);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

			std::copy_n(data, ELEMENTSOF(data), new_data);
			result = target(0, new_data, new_data);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
		}
		{ // 1
			result = target(1, data, new_data);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
			EXPECT_EQ(T(0), new_data[0]);

			std::copy_n(data, ELEMENTSOF(data), new_data);
			result = target(1, new_data, new_data);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
			EXPECT_EQ(T(0), new_data[0]);
		}
		{ // 2
			result = target(2, data, new_data);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
			constexpr auto median = (T(-3) + T(-2)) / 2;
			auto const test = [&]() {
				EXPECT_EQ(std::abs(data[0] - median), new_data[0]);
				EXPECT_EQ(std::abs(data[1] - median), new_data[1]);
				EXPECT_EQ(std::abs(data[0] - median), new_data[1]);
				EXPECT_EQ(std::abs(data[1] - median), new_data[0]);
			};
			test();

			std::copy_n(data, ELEMENTSOF(data), new_data);
			result = target(2, new_data, new_data);
			test();
		}
		{ // unaligned
			result = target(ELEMENTSOF(data), &data[1], new_data);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = target(ELEMENTSOF(data), data, &new_data[1]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = target(ELEMENTSOF(data), &new_data[1], &new_data[1]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
		}
		{ // nullptr
			result = target(ELEMENTSOF(data), nullptr, new_data);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = target(ELEMENTSOF(data), data, nullptr);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = target(ELEMENTSOF(data), nullptr, nullptr);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
		}
	}
	LIBSAKURA_SYMBOL(CleanUp)();}

}

TEST(Statistics, ComputeMedianAbsoluteDeviation) {
	TestComputeMedianAbsoluteDeviation<float,
			decltype(LIBSAKURA_SYMBOL(ComputeMedianAbsoluteDeviationFloat))>(
			LIBSAKURA_SYMBOL(ComputeMedianAbsoluteDeviationFloat));
}

namespace {

template<typename MessageType>
void ReportBenchmark(MessageType const &key, double sec) {
	std::cout << std::setprecision(5) << "#x# benchmark Stat_" << key << " " << sec
			<< std::endl;
}

template<typename T>
void Timing(char const str[], T stats_func) {
	double start = LIBSAKURA_SYMBOL(GetCurrentTime)();

	for (size_t i = 0; i < 15; ++i) {
		LIBSAKURA_SYMBOL(Status) status = stats_func();
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}
	double end = LIBSAKURA_SYMBOL(GetCurrentTime)();
	ReportBenchmark(str, end - start);
}

}

TEST(Statistics, ComputeStatistics_Performance) {
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
			nullptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

	{
		SIMD_ALIGN
		static std::array<float, 8192LU * 10000> data;
		SIMD_ALIGN
		static std::array<bool, data.size()> is_valid;
		std::iota(data.begin(), data.end(), 0);
		for (size_t i = 0; i < is_valid.size(); ++i) {
			is_valid[i] = i % 2 == 0;
		}
		LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
		Timing("ComputeStatisticsFloat",
				[&]() {return LIBSAKURA_SYMBOL (ComputeStatisticsFloat)(data.size(), data.data(),
							is_valid.data(), &result);});
		Timing("ComputeAccurateStatisticsFloat",
				[&]() {return LIBSAKURA_SYMBOL (ComputeAccurateStatisticsFloat)(data.size(), data.data(),
							is_valid.data(), &result);});
		assert(data.size() % 2 == 0);
		EXPECT_EQ(data.size() / 2, result.count);

		size_t new_elements = 0;
		Timing("SortValidValuesDenselyFloat",
				[&]() {
					return LIBSAKURA_SYMBOL (SortValidValuesDenselyFloat)(data.size() / 24, is_valid.data(), data.data(),
							&new_elements);
				});
		Timing("ComputeMedianAbsoluteDeviationFloat",
				[&]() {return LIBSAKURA_SYMBOL (ComputeMedianAbsoluteDeviationFloat)(new_elements, data.data(),
							data.data());});
}
	LIBSAKURA_SYMBOL(CleanUp)();
}
