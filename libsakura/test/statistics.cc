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
#include <algorithm>

#include <libsakura/sakura.h>

#include <libsakura/localdef.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"

namespace {
float Rms2(size_t elements, float const data[], bool const is_valid[]) {
	float result = 0;
	for (size_t i = 0; i < elements; ++i) {
		if (is_valid[i]) {
			result += data[i] * data[i];
		}
	}
	return result;
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
		size_t new_elements = static_cast<size_t>(-1);
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

		new_elements = static_cast<size_t>(-1);
		result = LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(0, is_valid, data,
				&new_elements);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
		EXPECT_EQ(0, new_elements);
	}
	{
		SIMD_ALIGN
		static float data[] = { 2.f, -2.f, 3.f, 0.f, -2.f, -1.f };
		SIMD_ALIGN
		static bool const is_valid[] = { true, false, true, false, true, true };
		STATIC_ASSERT(ELEMENTSOF(data) == ELEMENTSOF(is_valid));
		size_t new_elements = static_cast<size_t>(-1);
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
		size_t new_elements = static_cast<size_t>(-1);
		result = LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(ELEMENTSOF(data),
				is_valid, data, &new_elements);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
		EXPECT_EQ(0, new_elements);
	}
	{
		SIMD_ALIGN
		static float data[] = { 2.f, -2.f, 3.f, 0.f, -2.f, -1.f };
		size_t new_elements = static_cast<size_t>(-1);
		result = LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(ELEMENTSOF(data),
				nullptr, data, &new_elements);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
		EXPECT_EQ(static_cast<size_t>(-1), new_elements);
	}
	LIBSAKURA_SYMBOL(CleanUp)();
}

TEST(Statistics, ComputeStatistics) {
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
			nullptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

	{
		SIMD_ALIGN
		static float data[256];
		SIMD_ALIGN
		static bool is_valid[ELEMENTSOF(data)];
		{
			for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
				data[i] = i;
				is_valid[i] = true;
			}
			LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
			LIBSAKURA_SYMBOL(Status) status =
					LIBSAKURA_SYMBOL (ComputeStatisticsFloat)(ELEMENTSOF(data),
							data, is_valid, &result);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
			EXPECT_EQ(ELEMENTSOF(data), result.count);
			EXPECT_EQ(ELEMENTSOF(data) - 1, result.index_of_max);
			EXPECT_EQ(0, result.index_of_min);
			EXPECT_EQ(ELEMENTSOF(data) - 1, result.max);
			EXPECT_EQ(0, result.min);
			assert(ELEMENTSOF(data) % 2 == 0);
			EXPECT_EQ((ELEMENTSOF(data) - 1) * (ELEMENTSOF(data) / 2),
					result.sum);
			float const mean = (ELEMENTSOF(data) - 1) / 2.f;
			EXPECT_EQ(mean, result.mean);
			float rms2 = Rms2(ELEMENTSOF(data), data,
					is_valid) / ELEMENTSOF(data);
			EXPECT_DOUBLE_EQ(sqrt(rms2), result.rms);
			EXPECT_DOUBLE_EQ(sqrt(rms2 - mean * mean), result.stddev);
		}
		size_t min_max_index[16][2];
		for (size_t i = 0; i < ELEMENTSOF(min_max_index); ++i) {
			min_max_index[i][0] = 64 + i;
			min_max_index[i][1] = 127 - i;
		}
		std::fill(std::begin(data), std::end(data), 0.f);
		std::for_each(std::begin(min_max_index),std::end(min_max_index), [](size_t const (&idx)[2]) {
			size_t max_index = idx[0];
			size_t min_index = idx[1];
			data[max_index] = 20.f;
			data[min_index] = -120.f;
			LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
			LIBSAKURA_SYMBOL(Status) status =
					LIBSAKURA_SYMBOL (ComputeStatisticsFloat)(ELEMENTSOF(data),
							data, is_valid, &result);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
			EXPECT_EQ(ELEMENTSOF(data), result.count);
			EXPECT_EQ(max_index, result.index_of_max);
			EXPECT_EQ(min_index, result.index_of_min);
			EXPECT_EQ(data[max_index], result.max);
			EXPECT_EQ(data[min_index], result.min);
			assert(ELEMENTSOF(data) % 2 == 0);
			EXPECT_EQ(data[max_index] + data[min_index],
					result.sum);
			float const mean = (data[max_index] + data[min_index]) / ELEMENTSOF(data);
			EXPECT_EQ(mean, result.mean);
			float rms2 = (data[max_index] * data[max_index] + data[min_index] * data[min_index]) / ELEMENTSOF(data);
			EXPECT_DOUBLE_EQ(sqrt(rms2), result.rms);
			EXPECT_DOUBLE_EQ(sqrt(rms2 - mean * mean), result.stddev);
			data[max_index] = 0;
			data[min_index] = 0;
		});
	}
	{
		SIMD_ALIGN
		static float data[4];
		SIMD_ALIGN
		static bool is_valid[ELEMENTSOF(data)];
		for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
			data[i] = i;
			is_valid[i] = true;
		}
		LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL (ComputeStatisticsFloat)(
		ELEMENTSOF(data), data, is_valid, &result);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		EXPECT_EQ(ELEMENTSOF(data), result.count);
		EXPECT_EQ(ELEMENTSOF(data) - 1, result.index_of_max);
		EXPECT_EQ(0, result.index_of_min);
		EXPECT_EQ(ELEMENTSOF(data) - 1, result.max);
		EXPECT_EQ(0, result.min);
		assert(ELEMENTSOF(data) % 2 == 0);
		EXPECT_EQ((ELEMENTSOF(data) - 1) * (ELEMENTSOF(data) / 2), result.sum);
		float const mean = (ELEMENTSOF(data) - 1) / 2.f;
		EXPECT_EQ(mean, result.mean);
		float rms2 = Rms2(ELEMENTSOF(data), data, is_valid) / ELEMENTSOF(data);
		EXPECT_DOUBLE_EQ(sqrt(rms2), result.rms);
		EXPECT_DOUBLE_EQ(sqrt(rms2 - mean * mean), result.stddev);
	}

	{
		SIMD_ALIGN
		static float data[1] = { 3.f };
		SIMD_ALIGN
		static bool is_valid[ELEMENTSOF(data)] = { true };
		LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL (ComputeStatisticsFloat)(
		ELEMENTSOF(data), data, is_valid, &result);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		EXPECT_EQ(ELEMENTSOF(data), result.count);
		EXPECT_EQ(ELEMENTSOF(data) - 1, result.index_of_max);
		EXPECT_EQ(0, result.index_of_min);
		EXPECT_EQ(data[0], result.max);
		EXPECT_EQ(data[0], result.min);
		EXPECT_EQ(data[0], result.sum);
		float const mean = data[0];
		EXPECT_EQ(mean, result.mean);
		float rms2 = (data[0] * data[0]) / 1.f;
		EXPECT_DOUBLE_EQ(sqrt(rms2), result.rms);
		EXPECT_DOUBLE_EQ(sqrt(rms2 - mean * mean), result.stddev);
	}
	{
		SIMD_ALIGN
		static float data[0];
		SIMD_ALIGN
		static bool is_valid[ELEMENTSOF(data)];
		LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL (ComputeStatisticsFloat)(
		ELEMENTSOF(data), data, is_valid, &result);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		EXPECT_EQ(ELEMENTSOF(data), result.count);
		EXPECT_EQ(-1, result.index_of_max);
		EXPECT_EQ(-1, result.index_of_min);
		EXPECT_TRUE(std::isnan(result.max));
		EXPECT_TRUE(std::isnan(result.min));
		EXPECT_EQ(0.f, result.sum);
		EXPECT_TRUE(std::isnan(result.mean));
		EXPECT_TRUE(std::isnan(result.rms));
		EXPECT_TRUE(std::isnan(result.stddev));
	}
	{
		SIMD_ALIGN
		static float data[1024];
		SIMD_ALIGN
		static bool is_valid[ELEMENTSOF(data)];
		for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
			data[i] = i;
			is_valid[i] = false;
		}
		LIBSAKURA_SYMBOL (StatisticsResultFloat) result;
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL (ComputeStatisticsFloat)(
		ELEMENTSOF(data), data, is_valid, &result);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		EXPECT_EQ(0, result.count);
		EXPECT_EQ(-1, result.index_of_max);
		EXPECT_EQ(-1, result.index_of_min);
		EXPECT_TRUE(std::isnan(result.max));
		EXPECT_TRUE(std::isnan(result.min));
		EXPECT_EQ(0.f, result.sum);
		EXPECT_TRUE(std::isnan(result.mean));
		EXPECT_TRUE(std::isnan(result.rms));
		EXPECT_TRUE(std::isnan(result.stddev));
	}
	{
		SIMD_ALIGN
		static float data[1024];
		SIMD_ALIGN
		static bool is_valid[ELEMENTSOF(data)];
		size_t count = 0;
		float sum = 0;
		float max = -1;
		float min = ELEMENTSOF(data) + 1;
		size_t max_idx = -1;
		size_t min_idx = -1;
		float sqsum = 0;
		for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
			data[i] = i;
			is_valid[i] = i % 7 == 0;
			if (is_valid[i]) {
				count++;
				sum += data[i];
				sqsum += data[i] * data[i];
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
		float const rms2 = sqsum / count;
		float const mean = sum / count;

		LIBSAKURA_SYMBOL (StatisticsResultFloat) result;
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL (ComputeStatisticsFloat)(
		ELEMENTSOF(data), data, is_valid, &result);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		EXPECT_EQ(count, result.count);
		EXPECT_EQ(max_idx, result.index_of_max);
		EXPECT_EQ(min_idx, result.index_of_min);
		EXPECT_EQ(max, result.max);
		EXPECT_EQ(min, result.min);
		EXPECT_DOUBLE_EQ(sum, result.sum);
		EXPECT_DOUBLE_EQ(mean, result.mean);
		EXPECT_EQ(static_cast<int>(1000 * std::sqrt(rms2)),
				static_cast<int>(1000 * result.rms));
		EXPECT_EQ(static_cast<int>(1000 * std::sqrt(rms2 - mean * mean)),
				static_cast<int>(1000 * result.stddev));
	}
	LIBSAKURA_SYMBOL(CleanUp)();
}

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeStddevFloat)(
		size_t count,
		double mean, size_t elements, float const data[], bool const is_valid[],
		double *result) {
	CHECK_ARGS(elements <= INT32_MAX);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(is_valid != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(is_valid));
	CHECK_ARGS(result != nullptr);

#if 1
	double sq_diff = 0.;
	for (size_t i = 0; i < elements; ++i) {
		if (is_valid[i]) {
			double diff = data[i] - mean;
			sq_diff += diff * diff;
		}
	}
	*result = sqrt(sq_diff / count);
#else
	{
		size_t n = 0;
		double mean = 0;
		double M2 = 0;

		for (size_t i = 0; i < elements; ++i) {
			if (is_valid[i]) {
				++n;
				double delta = data[i] - mean;
				mean += delta / n;
				M2 += delta * (data[i] - mean);
			}
		}

		if (n < 2) {
			*result = 0.;
		}
		*result = M2 / (n - 1);
	}
#endif
	return LIBSAKURA_SYMBOL(Status_kOK);
}

TEST(Statistics, ComputeStatistics_Accuracy) {
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
			nullptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

	{
		SIMD_ALIGN
		static float data[8192 * 10000];
		SIMD_ALIGN
		static bool is_valid[ELEMENTSOF(data)];
		for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
			data[i] = i;
			is_valid[i] = i % 2 == 0;
		}
		LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
		{
			double start = LIBSAKURA_SYMBOL(GetCurrentTime)();

			for (size_t i = 0; i < 1; ++i) {
				LIBSAKURA_SYMBOL(Status) status =
				LIBSAKURA_SYMBOL (ComputeStatisticsFloat)(ELEMENTSOF(data), data,
						is_valid, &result);
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
			}
			double end = LIBSAKURA_SYMBOL(GetCurrentTime)();
			printf("Statistics time: %lf\n", end - start);
		}
		assert(ELEMENTSOF(data) % 2 == 0);
		EXPECT_EQ(ELEMENTSOF(data) / 2, result.count);

		STATIC_ASSERT(ELEMENTSOF(data) % 4 == 0);
		constexpr float base = 9.0e+5f;
		for (size_t i = 0; i < ELEMENTSOF(data); i += 4) {
			data[i] = base + 4;
			data[i+1] = base + 7;
			data[i+2] = base + 13;
			data[i+3] = base + 16;
			is_valid[i] = true;
			is_valid[i+1] = true;
			is_valid[i+2] = true;
			is_valid[i+3] = true;
		}
		{
			LIBSAKURA_SYMBOL(Status) status =
			LIBSAKURA_SYMBOL (ComputeStatisticsFloat)(ELEMENTSOF(data), data,
					is_valid, &result);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		}
		assert(ELEMENTSOF(data) % 2 == 0);
		EXPECT_EQ(ELEMENTSOF(data), result.count);
		EXPECT_EQ(base + 4, result.min);
		EXPECT_EQ(base + 16, result.max);
		EXPECT_EQ(base + 10, result.mean);
		EXPECT_NEAR((base + 10) * ELEMENTSOF(data), result.sum, 1400000);
		constexpr double rms2 = (double(base + 4) * double(base + 4)
				+ double(base + 7) * double(base + 7)
				+ double(base + 13) * double(base + 13)
				+ double(base + 16) * double(base + 16)) / 4.;
		constexpr auto rms = float(sqrt(rms2));
		EXPECT_NEAR(rms, result.rms, 0.00006);
		constexpr double variance = ((4. - 10.) * (4. - 10.)
				+ (7. - 10.) * (7. - 10.) + (13. - 10.) * (13. - 10.)
				+ (16. - 10.) * (16. - 10.)) / 4.;

		double stddev = 0;
		LIBSAKURA_SYMBOL(ComputeStddevFloat)(result.count, result.mean, ELEMENTSOF(data), data,
					is_valid, &stddev);
		//std::cout << std::setprecision(16) << sqrt(variance) << std::endl << stddev << std::endl << result.stddev << std::endl;
		EXPECT_DOUBLE_EQ(sqrt(variance), stddev);
		EXPECT_NEAR(sqrt(variance), result.stddev, 5.5);
	}
	LIBSAKURA_SYMBOL(CleanUp)();
}

TEST(Statistics, ComputeStatistics_Performance) {
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
			nullptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

	{
		SIMD_ALIGN
		static float data[8192 * 10];
		SIMD_ALIGN
		static bool is_valid[ELEMENTSOF(data)];
		for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
			data[i] = i;
			is_valid[i] = i % 2 == 0;
		}
		LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
		{
			double start = LIBSAKURA_SYMBOL(GetCurrentTime)();

			for (size_t i = 0; i < 20000; ++i) {
				LIBSAKURA_SYMBOL(Status) status =
				LIBSAKURA_SYMBOL (ComputeStatisticsFloat)(ELEMENTSOF(data), data,
						is_valid, &result);
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
			}
			double end = LIBSAKURA_SYMBOL(GetCurrentTime)();
			printf("Statistics time: %lf\n", end - start);
		}
		assert(ELEMENTSOF(data) % 2 == 0);
		EXPECT_EQ(ELEMENTSOF(data) / 2, result.count);
	}
	LIBSAKURA_SYMBOL(CleanUp)();
}
