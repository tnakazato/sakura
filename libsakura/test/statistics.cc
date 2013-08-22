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
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)();
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

	{
		SIMD_ALIGN
		static float data[] = { 2.f, -2.f, 3.f, 0.f, -2.f, -1.f };
		SIMD_ALIGN
		static bool const is_valid[] = { true, true, true, true, true, true };
		ASSERT_EQ(ELEMENTSOF(data), ELEMENTSOF(is_valid));
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
		assert(ELEMENTSOF(data) == ELEMENTSOF(is_valid));
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
		assert(ELEMENTSOF(data) == ELEMENTSOF(is_valid));
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
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)();
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
	LIBSAKURA_SYMBOL(CleanUp)();
}

TEST(Statistics, ComputeStatistics_Speed) {
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)();
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

			for (size_t i = 0; i < 30000; ++i) {
				LIBSAKURA_SYMBOL(Status) status =
						LIBSAKURA_SYMBOL (ComputeStatistics)(ELEMENTSOF(data),
								data, is_valid, &result);
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
