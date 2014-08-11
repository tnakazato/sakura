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
		result = LIBSAKURA_SYMBOL(SortValidValuesDensely)(ELEMENTSOF(data),
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
		result = LIBSAKURA_SYMBOL(SortValidValuesDensely)(0, is_valid, data,
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
		result = LIBSAKURA_SYMBOL(SortValidValuesDensely)(ELEMENTSOF(data),
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
		result = LIBSAKURA_SYMBOL(SortValidValuesDensely)(ELEMENTSOF(data),
				is_valid, data, &new_elements);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
		EXPECT_EQ(0, new_elements);
	}
	{
		SIMD_ALIGN
		static float data[] = { 2.f, -2.f, 3.f, 0.f, -2.f, -1.f };
		size_t new_elements = static_cast<size_t>(-1);
		result = LIBSAKURA_SYMBOL(SortValidValuesDensely)(ELEMENTSOF(data),
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
		for (size_t i = 0; i < ELEMENTSOF(data); ++i) {
			data[i] = i;
			is_valid[i] = true;
		}
		LIBSAKURA_SYMBOL(StatisticsResult) result;
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL (ComputeStatistics)(
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
		EXPECT_FLOAT_EQ(sqrt(rms2), result.rms);
		EXPECT_FLOAT_EQ(sqrt(rms2 - mean * mean), result.stddev);
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
		LIBSAKURA_SYMBOL(StatisticsResult) result;
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL (ComputeStatistics)(
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
		EXPECT_FLOAT_EQ(sqrt(rms2), result.rms);
		EXPECT_FLOAT_EQ(sqrt(rms2 - mean * mean), result.stddev);
	}

	{
		SIMD_ALIGN
		static float data[1] = { 3.f };
		SIMD_ALIGN
		static bool is_valid[ELEMENTSOF(data)] = { true };
		LIBSAKURA_SYMBOL(StatisticsResult) result;
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL (ComputeStatistics)(
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
		EXPECT_FLOAT_EQ(sqrt(rms2), result.rms);
		EXPECT_FLOAT_EQ(sqrt(rms2 - mean * mean), result.stddev);
	}
	{
		SIMD_ALIGN
		static float data[0];
		SIMD_ALIGN
		static bool is_valid[ELEMENTSOF(data)];
		LIBSAKURA_SYMBOL(StatisticsResult) result;
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL (ComputeStatistics)(
		ELEMENTSOF(data), data, is_valid, &result);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		EXPECT_EQ(ELEMENTSOF(data), result.count);
		EXPECT_EQ(-1, result.index_of_max);
		EXPECT_EQ(-1, result.index_of_min);
		EXPECT_TRUE(isnanf(result.max));
		EXPECT_TRUE(isnanf(result.min));
		EXPECT_EQ(0.f, result.sum);
		EXPECT_TRUE(isnanf(result.mean));
		EXPECT_TRUE(isnanf(result.rms));
		EXPECT_TRUE(isnanf(result.stddev));
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
		LIBSAKURA_SYMBOL (StatisticsResult) result;
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL (ComputeStatistics)(
		ELEMENTSOF(data), data, is_valid, &result);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		EXPECT_EQ(0, result.count);
		EXPECT_EQ(-1, result.index_of_max);
		EXPECT_EQ(-1, result.index_of_min);
		EXPECT_TRUE(isnanf(result.max));
		EXPECT_TRUE(isnanf(result.min));
		EXPECT_EQ(0.f, result.sum);
		EXPECT_TRUE(isnanf(result.mean));
		EXPECT_TRUE(isnanf(result.rms));
		EXPECT_TRUE(isnanf(result.stddev));
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

		LIBSAKURA_SYMBOL (StatisticsResult) result;
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL (ComputeStatistics)(
		ELEMENTSOF(data), data, is_valid, &result);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
		EXPECT_EQ(count, result.count);
		EXPECT_EQ(max_idx, result.index_of_max);
		EXPECT_EQ(min_idx, result.index_of_min);
		EXPECT_EQ(max, result.max);
		EXPECT_EQ(min, result.min);
		EXPECT_FLOAT_EQ(sum, result.sum);
		EXPECT_FLOAT_EQ(mean, result.mean);
		EXPECT_EQ(static_cast<int>(1000 * std::sqrt(rms2)),
				static_cast<int>(1000 * result.rms));
		EXPECT_EQ(static_cast<int>(1000 * std::sqrt(rms2 - mean * mean)),
				static_cast<int>(1000 * result.stddev));
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
		LIBSAKURA_SYMBOL(StatisticsResult) result;
		{
			double start = LIBSAKURA_SYMBOL(GetCurrentTime)();

			for (size_t i = 0; i < 20000; ++i) {
				LIBSAKURA_SYMBOL(Status) status =
				LIBSAKURA_SYMBOL (ComputeStatistics)(ELEMENTSOF(data), data,
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
