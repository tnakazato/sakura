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
#include <cassert>
#include <cstddef>
#include <cmath>
#include <climits>
#include <algorithm>
#include <sys/types.h>

#include "libsakura/sakura.h"
#include "libsakura/localdef.h"

#define FORCE_EIGEN 0

namespace {

template<size_t kUnitSize, size_t kBlockSize, typename Reducer>
typename Reducer::MiddleLevelAccumulator BlockWiseTraverse(size_t offset,
		int levels, size_t stride, size_t num_blocks, size_t full_blocks,
		typename Reducer::MiddleLevelAccumulator const &residual,
		Reducer *reducer) {
	//std::cout << "stride = " << stride << std::endl;
	typename Reducer::MiddleLevelAccumulator acc_stack[levels + 1][kUnitSize];
	acc_stack[0][0].Clear();
	typename Reducer::MiddleLevelAccumulator (*acc)[kUnitSize] = &acc_stack[1];

	size_t indices_stack[levels + 1];
	indices_stack[0] = 0;
	size_t *indices = &indices_stack[1];
	size_t start = offset;

	int leaf = levels - 1;
	int level = 0;
	indices[level] = 0;
	auto const upper_limit = offset + num_blocks * kBlockSize;
	while (level >= 0) {
		//std::cout << "level: " << level << ", i: " << indices[level] << " full_blocks: " << full_blocks << " limit: " << upper_limit << std::endl;
		assert(0 <= level && level < levels);
		if (level >= leaf) {
			auto &accumlator = acc[level - 1][indices[level - 1]];
			accumlator.Clear();
			if (start == offset) {
				reducer->MiddleLevelReduce(residual, accumlator);
			}
			auto delta = start;
			for (size_t i = 0; i < full_blocks; ++i) {
				for (size_t j = 0; j < kUnitSize; ++j) {
					reducer->LastLevelReduce(delta + j * stride, kBlockSize,
							accumlator);
				}
				delta += stride * kUnitSize;
			}
			for (size_t j = 0; j < kUnitSize; ++j) {
				auto index = delta + j * stride;
				if (index >= upper_limit) {
					break;
				}
				reducer->LastLevelReduce(index, kBlockSize, accumlator);
			}
			start += kBlockSize;
			--level;
			++indices[level];
		} else {
			if (indices[level] < kUnitSize) {
				++level;
				indices[level] = 0;
			} else {
				auto &accumlator = acc[level - 1][indices[level - 1]];
				accumlator = acc[level][0];
				auto acc_local = acc[level];
				for (size_t j = 1; j < kUnitSize; ++j) {
					reducer->MiddleLevelReduce(acc_local[j], accumlator);
				}
				--level;
				++indices[level];
			}
		}
	}
	return acc_stack[0][0];
}

/**
 * @brief Traverses data in a range [@a offset, @a offset + @a num_data)
 * to reduce them to @a Reducer::TopLevelAccumulator using @a reducer.
 *
 * This function reduces blocks to @a Reducer::MiddleLevelAccumulator
 * using @a Reducer::LastLevelReduce() or remaining consecutive data less than @a kBlockSize
 * to @a Reducer::MiddleLevelAccumulator using @a Reducer::SequentialReduce() .
 * Every @a kUnitSize @a Reducer::MiddleLevelAccumulators are repeatedly reduced
 * to a @a Reducer::MiddleLevelAccumulator
 * using @a Reducer::MiddleLevelReduce() in the depth first order until all of them are reduced to a single one.
 * Finally, the @a Reducer::MiddleLevelAccumulator is reduced to @a Reducer::TopLevelAccumulator
 * using @a Reducer::TopLevelReduce() .
 *
 * @see ../doc/misc/SakuraStatistics-algorithm.xlsx for detail of Traverse algorithm.
 *
 * @tparam kUnitSize	The number of blocks in a unit. Unit is a processing unit of middle level reducer.
 * @tparam kBlockSize	Size of a block. Block is a processing unit of low level reducer
 * and consecutive @a kBlockSize data.
 * @tparam Reducer	@a Reducer defines how to reduce data.
 * @param offset	@a offset indicates a starting point of the range to be traversed.
 * @param num_data	@a num_data indicates the number of data in the range to be traversed.
 * @param reducer	@a reducer defines how to reduce data in the range.
 * @return a reduced data
 */
template<size_t kUnitSize, size_t kBlockSize, typename Reducer>
typename Reducer::TopLevelAccumulator Traverse(size_t offset,
		size_t const num_data, Reducer *reducer) {
	STATIC_ASSERT(kUnitSize > 1 && kBlockSize > 0);
	typename Reducer::TopLevelAccumulator top_level_result;
	top_level_result.Clear();

	auto num_blocks = num_data / kBlockSize;
	int levels = 0;
	size_t n = 1;
	auto stride = kBlockSize;
	while (n * kUnitSize <= num_blocks) {
		++levels;
		n *= kUnitSize;
		stride *= kUnitSize;
	}
	auto full_blocks = num_blocks / n;

	typename Reducer::MiddleLevelAccumulator residual;
	residual.Clear();
	{
		if (levels == 0) {
			num_blocks = 0;
		}
		auto rest = num_data - num_blocks * kBlockSize;
		reducer->SequentialReduce(offset + num_blocks * kBlockSize, rest,
				residual);
	}

	if (levels > 0) {
		stride /= kUnitSize;
		reducer->TopLevelReduce(
				reducer->TopLevelReduce(
						BlockWiseTraverse<kUnitSize, kBlockSize, Reducer>(
								offset, levels, stride, num_blocks, full_blocks,
								residual, reducer)), top_level_result);
	} else {
		reducer->TopLevelReduce(reducer->TopLevelReduce(residual),
				top_level_result);
	}
	return top_level_result;
}

template<typename Scalar, typename Accumulator>
struct ScalarStats {
	Scalar const *data;bool const *mask;
	ScalarStats(Scalar const *data_arg, bool const *mask_arg) :
			data(data_arg), mask(mask_arg) {
	}
	struct SIMD_ALIGN TopLevelAccumulator {
		size_t count;
		Accumulator sum;
		Accumulator square_sum;
		Scalar min, max;
		int32_t index_of_min, index_of_max;

		void Clear() {
			count = 0;
			sum = 0.;
			square_sum = 0.;
			min = max = NAN;
			index_of_min = index_of_max = -1;
		}
	};

	typedef TopLevelAccumulator MiddleLevelAccumulator;

	static void StaticTopLevelReduce(TopLevelAccumulator const &increment,
			TopLevelAccumulator&accumulator) {
		accumulator.count += increment.count;
		accumulator.sum += increment.sum;
		accumulator.square_sum += increment.square_sum;
		if (accumulator.index_of_min < 0) {
			accumulator.index_of_min = increment.index_of_min;
			accumulator.index_of_max = increment.index_of_max;
			accumulator.min = increment.min;
			accumulator.max = increment.max;
		} else {
			if (increment.min < accumulator.min) {
				assert(increment.index_of_min >= 0);
				accumulator.index_of_min = increment.index_of_min;
				accumulator.min = increment.min;
			}
			if (increment.max > accumulator.max) {
				assert(increment.index_of_max >= 0);
				accumulator.index_of_max = increment.index_of_max;
				accumulator.max = increment.max;
			}
		}
	}

	void TopLevelReduce(TopLevelAccumulator const &increment,
			TopLevelAccumulator&accumulator) {
		StaticTopLevelReduce(increment, accumulator);
	}

	TopLevelAccumulator TopLevelReduce(
			MiddleLevelAccumulator const &accumulator) {
		return accumulator;
	}

	void MiddleLevelReduce(MiddleLevelAccumulator const &increment,
			MiddleLevelAccumulator &accumulator) {
		TopLevelReduce(increment, accumulator);
	}

	template<typename T>
	static void Accumulate(size_t pos, size_t negative_offset, Scalar value,
			T &count, int32_t &index_of_min, int32_t &index_of_max,
			Scalar &min, Scalar &max, Accumulator &sum,
			Accumulator &square_sum) {
		++count;
		Accumulator v = value;
		sum += v;
		square_sum += v * v;
		if (index_of_min < 0) {
			index_of_min = index_of_max = pos - negative_offset;
			min = max = value;
		} else {
			if (value < min) {
				index_of_min = pos - negative_offset;
				min = value;
			}
			if (value > max) {
				index_of_max = pos - negative_offset;
				max = value;
			}
		}
	}

	void Reduce(size_t pos, MiddleLevelAccumulator &accumulator) {
		if (mask[pos]) {
			Accumulate(pos, 0u, data[pos], accumulator.count,
					accumulator.index_of_min, accumulator.index_of_max,
					accumulator.min, accumulator.max, accumulator.sum,
					accumulator.square_sum);
		}
	}

	void LastLevelReduce(size_t position, size_t block_size,
			MiddleLevelAccumulator &accumulator) {
		for (size_t i = 0; i < block_size; ++i) {
			Reduce(position + i, accumulator);
		}
	}

	void SequentialReduce(size_t offset, size_t elements,
			MiddleLevelAccumulator &accumulator) {
		for (size_t i = 0; i < elements; ++i) {
			Reduce(offset + i, accumulator);
		}
	}
};

}

#if defined(__AVX__) && !defined(ARCH_SCALAR) && (! FORCE_EIGEN)
#include <immintrin.h>
#include <cstdint>

namespace {

#include "libsakura/packed_operation.h"

union m256 {
	__m256 m256;
	__m256i m256i;
	__m256d m256d;
	float floatv[8];
	double doublev[4];
	int32_t intv[8];
	char charv[32];
};

union m128 {
	__m128 m128;
	__m128i m128i;
	__m128d m128d;
	float floatv[4];
	double doublev[2];
	int32_t intv[4];
	char charv[16];
};

#if defined(__AVX2__)
struct SIMDWordForInt {
	m256 value;
	auto IntValue() const -> decltype(value.m256i) {
		return value.m256i;
	}
	void Clear() {
		value.m256 = _mm256_setzero_ps();
	}
	void Add(decltype(value.m256i) const &v) {
		value.m256i += v;
	}
};
#else
struct SIMDWordForInt {
	m128 value;
	auto IntValue() const -> decltype(value.m128i) {
		return value.m128i;
	}
	void Clear() {
		value.m128 = _mm_setzero_ps();
	}
	void Add(decltype(value.m128i) const &v) {
		value.m128i += v;
	}
};
#endif

inline void StatsBlock(size_t i, __m256 const *data_arg, double const *mask_arg,
		SIMDWordForInt &count, __m256d &sum, __m256d &square_sum, __m256 &min,
		__m256 &max, __m256i &index_of_min, __m256i &index_of_max) {
	auto const data = AssumeAligned(data_arg);

	auto const zero = _mm256_setzero_ps();
	auto const nan = _mm256_set1_ps(NAN);

#if defined(__AVX2__)
	auto const mask = AssumeAligned(mask_arg);
	auto const zero256i = _mm256_setzero_si256();
	auto mask8 = _mm256_cvtepi8_epi32(_mm_castpd_si128(_mm_load1_pd(mask)));
	count.Add(mask8);
	mask8 = _mm256_cmpeq_epi32(mask8, zero256i);
#else
	auto const mask = AssumeAligned<float const *, sizeof(decltype(count.IntValue()))>(reinterpret_cast<float const *>(mask_arg));
	auto const zero128i = _mm_setzero_si128();
	auto mask0 = _mm_castps_si128(_mm_load_ss(&mask[0]));
	auto mask1 = _mm_castps_si128(_mm_load_ss(&mask[1]));
	mask0 = _mm_cvtepu8_epi32(mask0);
	mask1 = _mm_cvtepu8_epi32(mask1);
	count.Add(_mm_add_epi32(mask0, mask1));

	mask0 = _mm_cmpeq_epi32(mask0, zero128i);
	mask1 = _mm_cmpeq_epi32(mask1, zero128i);
	auto mask8 = _mm256_insertf128_si256(_mm256_castsi128_si256(mask0),
			mask1, 1);
#endif
	auto maskf = _mm256_cvtepi32_ps(mask8);
	/* maskf: 0xffffffff means invalid data, 0 means valid data. */

	auto value = _mm256_blendv_ps(data[0], nan, maskf);
	auto const index = _mm256_castps_si256(
			_mm256_blendv_ps(_mm256_castsi256_ps(_mm256_set1_epi32(i)), maskf,
					maskf));
	{
		auto take_increment = _mm256_or_ps(
				_mm256_cmp_ps(value, min, _CMP_LT_OQ),
				_mm256_castsi256_ps(index_of_min));
		min = _mm256_blendv_ps(min, value, take_increment);
		index_of_min = _mm256_castps_si256(
				_mm256_blendv_ps(_mm256_castsi256_ps(index_of_min),
						_mm256_castsi256_ps(index), take_increment));
	}
	{
		auto take_increment = _mm256_or_ps(
				_mm256_cmp_ps(value, max, _CMP_GT_OQ),
				_mm256_castsi256_ps(index_of_max));
		max = _mm256_blendv_ps(max, value, take_increment);
		index_of_max = _mm256_castps_si256(
				_mm256_blendv_ps(_mm256_castsi256_ps(index_of_max),
						_mm256_castsi256_ps(index), take_increment));
	}

	value = _mm256_blendv_ps(value, zero, maskf);
	auto v = _mm256_cvtps_pd(_mm256_castps256_ps128(value));
	sum = _mm256_add_pd(sum, v);
	square_sum = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
	LIBSAKURA_SYMBOL(SimdPacketAVX), double>(v, v, square_sum);

	v = _mm256_cvtps_pd(_mm256_extractf128_ps(value, 1));
	sum = _mm256_add_pd(sum, v);
	square_sum = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
	LIBSAKURA_SYMBOL(SimdPacketAVX), double>(v, v, square_sum);
}

inline double AddHorizontally(__m256d packed_values) {
	packed_values = _mm256_hadd_pd(packed_values, packed_values);
	__m128d sum2 = _mm256_extractf128_pd(packed_values, 1);
	return _mm_cvtsd_f64(
			_mm_add_sd(_mm256_castpd256_pd128(packed_values), sum2));
}

#if defined(__AVX2__)
inline int32_t AddHorizontally(__m256i packed_values) {
	packed_values = _mm256_hadd_epi32(packed_values, packed_values);
	packed_values = _mm256_hadd_epi32(packed_values, packed_values);
	__m128i sum2 = _mm256_extractf128_si256(packed_values, 1);
	return _mm_cvtsi128_si32(
			_mm_add_epi32(_mm256_castsi256_si128(packed_values), sum2));
}
#else
inline int32_t AddHorizontally128(__m128i packed_values) {
	packed_values = _mm_hadd_epi32(packed_values, packed_values);
	packed_values = _mm_hadd_epi32(packed_values, packed_values);
	int32_t total = _mm_extract_epi32(packed_values, 0);
	return total;
}
#endif

template<typename Scalar, typename Accumulator>
struct SIMDStats {
};

template<>
struct SIMDStats<float, double> {
	typedef float Scalar;
	typedef double Accumulator;
	typedef ScalarStats<Scalar, Accumulator> ScalarVersion;

	Scalar const *data;bool const *is_valid;
	SIMDStats(Scalar const *data_arg, bool const *mask_arg) :
			data(data_arg), is_valid(mask_arg) {
	}

	typedef ScalarVersion::TopLevelAccumulator TopLevelAccumulator;

	struct SIMD_ALIGN MiddleLevelAccumulator {
		SIMDWordForInt count;
		m256 sum;
		m256 square_sum;
		m256 min, max;
		m256 index_of_min, index_of_max;

		void Clear() {
			auto const zero = _mm256_setzero_ps();
			count.Clear();
			sum.m256 = zero;
			square_sum.m256 = zero;
			auto const nan = _mm256_set1_ps(NAN);
			min.m256 = nan;
			max.m256 = nan;
			auto const neg1 = _mm256_set1_epi32(-1);
			index_of_min.m256i = neg1;
			index_of_max.m256i = neg1;
		}
	};

	void TopLevelReduce(TopLevelAccumulator const &increment,
			TopLevelAccumulator&accumulator) {
		ScalarVersion::StaticTopLevelReduce(increment, accumulator);
	}

	TopLevelAccumulator TopLevelReduce(
			MiddleLevelAccumulator const &accumulator) {
		TopLevelAccumulator result;
#if defined(__AVX2__)
		result.count = static_cast<size_t>(AddHorizontally(
				accumulator.count.IntValue()));
#else
		result.count = static_cast<size_t>(AddHorizontally128(accumulator.count.IntValue()));
#endif

		result.sum = AddHorizontally(accumulator.sum.m256d);
		result.square_sum = AddHorizontally(accumulator.square_sum.m256d);

		{
			m256 tmp;
			tmp.m256 = accumulator.min.m256;
			m256 tmp_index;
			tmp_index.m256i = accumulator.index_of_min.m256i;

			float r = tmp.floatv[0];
			auto result_index = tmp_index.intv[0] + 0;
			for (unsigned i = 1; i < ELEMENTSOF(tmp.floatv); ++i) {
				if (!std::isnan(tmp.floatv[i])) {
					assert(tmp_index.intv[i] >= 0);
					if (std::isnan(r) || tmp.floatv[i] < r) {
						r = tmp.floatv[i];
						result_index = tmp_index.intv[i] + i;
					}
				}
			}
			result.min = r;
			result_index = std::max(static_cast<ssize_t>(-1), static_cast<ssize_t>(result_index));
			result.index_of_min = result_index;
		}
		{
			m256 tmp;
			tmp.m256 = accumulator.max.m256;
			m256 tmp_index;
			tmp_index.m256i = accumulator.index_of_max.m256i;

			float r = tmp.floatv[0];
			auto result_index = tmp_index.intv[0] + 0;
			for (unsigned i = 1; i < ELEMENTSOF(tmp.floatv); ++i) {
				if (!std::isnan(tmp.floatv[i])) {
					assert(tmp_index.intv[i] >= 0);
					if (std::isnan(r) || tmp.floatv[i] > r) {
						r = tmp.floatv[i];
						result_index = tmp_index.intv[i] + i;
					}
				}
			}
			result.max = r;
			result_index = std::max(static_cast<ssize_t>(-1), static_cast<ssize_t>(result_index));
			result.index_of_max = result_index;
		}
		return result;
	}

	void MiddleLevelReduce(MiddleLevelAccumulator const &increment,
			MiddleLevelAccumulator &accumulator) {
		accumulator.count.Add(increment.count.IntValue());
		accumulator.sum.m256d += increment.sum.m256d;
		accumulator.square_sum.m256d += increment.square_sum.m256d;

		{
			auto take_increment = _mm256_or_ps(
					_mm256_cmp_ps(increment.min.m256, accumulator.min.m256,
							_CMP_LT_OQ),
					_mm256_castsi256_ps(accumulator.index_of_min.m256i));
			accumulator.min.m256 = _mm256_blendv_ps(accumulator.min.m256,
					increment.min.m256, take_increment);
			accumulator.index_of_min.m256i = _mm256_castps_si256(
					_mm256_blendv_ps(
							_mm256_castsi256_ps(accumulator.index_of_min.m256i),
							_mm256_castsi256_ps(increment.index_of_min.m256i),
							take_increment));
		}
		{
			auto take_increment = _mm256_or_ps(
					_mm256_cmp_ps(increment.max.m256, accumulator.max.m256,
							_CMP_GT_OQ),
					_mm256_castsi256_ps(accumulator.index_of_max.m256i));
			accumulator.max.m256 = _mm256_blendv_ps(accumulator.max.m256,
					increment.max.m256, take_increment);
			accumulator.index_of_max.m256i = _mm256_castps_si256(
					_mm256_blendv_ps(
							_mm256_castsi256_ps(accumulator.index_of_max.m256i),
							_mm256_castsi256_ps(increment.index_of_max.m256i),
							take_increment));
		}
	}

	void LastLevelReduce(size_t position, size_t block_size,
			MiddleLevelAccumulator &accumulator) {
		auto mask_ptr = reinterpret_cast<double const *>(&is_valid[position]);
		StatsBlock(position, reinterpret_cast<__m256 const *>(&data[position]),
				mask_ptr, accumulator.count, accumulator.sum.m256d,
				accumulator.square_sum.m256d, accumulator.min.m256,
				accumulator.max.m256, accumulator.index_of_min.m256i,
				accumulator.index_of_max.m256i);
		//std::cout << position << " - " << position + 7 << std::endl;
	}

	void SequentialReduce(size_t offset, size_t elements,
			MiddleLevelAccumulator &accumulator) {
		for (size_t j = 0; j < elements; ++j) {
			auto pos = offset + j;
			if (is_valid[pos]) {
				auto i_float = j % ELEMENTSOF(accumulator.min.floatv);
				auto i_double = j % ELEMENTSOF(accumulator.min.doublev);
				ScalarStats<Scalar, Accumulator>::Accumulate(pos, j, data[pos],
						accumulator.count.value.intv[i_float
								% ELEMENTSOF(accumulator.count.value.intv)],
						accumulator.index_of_min.intv[i_float],
						accumulator.index_of_max.intv[i_float],
						accumulator.min.floatv[i_float],
						accumulator.max.floatv[i_float],
						accumulator.sum.doublev[i_double],
						accumulator.square_sum.doublev[i_double]);
			}
		}
	}
};

void ComputeStatisticsSimdFloat(size_t num_data, float const data[],
		bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result_arg) {
	STATIC_ASSERT(sizeof(m256) == sizeof(__m256 ));
	auto const zero_d = _mm256_setzero_pd();
	auto sum = zero_d;
	auto square_sum = zero_d;
	auto const nan = _mm256_set1_ps(NAN);
	auto min = nan;
	auto max = nan;
	auto index_of_min = _mm256_set1_epi32(-1);
	auto index_of_max = _mm256_set1_epi32(-1);
	SIMDWordForInt count;
	count.Clear();
	STATIC_ASSERT(
			true == 1 && false == 0 && sizeof(bool) == 1 && sizeof(bool) * 8 == sizeof(double));
	double const *mask_ptr = reinterpret_cast<double const *>(is_valid);
	auto const *data_ptr = AssumeAligned(reinterpret_cast<__m256 const *>(data));
	for (size_t i = 0; i < num_data / (sizeof(__m256 ) / sizeof(float)); ++i) {
		StatsBlock(i, &data_ptr[i], &mask_ptr[i], count, sum, square_sum, min, max,
				index_of_min, index_of_max);
	}
#if defined(__AVX2__)
	size_t counted = static_cast<size_t>(AddHorizontally(count.IntValue()));
#else
	size_t counted = static_cast<size_t>(AddHorizontally128(count.IntValue()));
#endif

	double total = AddHorizontally(sum);
	double square_total = AddHorizontally(square_sum);

	LIBSAKURA_SYMBOL(StatisticsResultFloat) result;
	{
		m256 tmp;
		tmp.m256 = min;
		m256 tmp_index;
		tmp_index.m256i = index_of_min;

		float r = tmp.floatv[0];
		ssize_t result_index =
				tmp_index.intv[0] == -1 ?
						-1 : tmp_index.intv[0] * ELEMENTSOF(tmp.intv);
		for (unsigned i = 1; i < ELEMENTSOF(tmp.floatv); ++i) {
			if (!std::isnan(tmp.floatv[i])) {
				assert(tmp_index.intv[i] >= 0);
				if (std::isnan(r) || tmp.floatv[i] < r) {
					r = tmp.floatv[i];
					result_index = tmp_index.intv[i] * ELEMENTSOF(tmp.intv) + i;
				}
			}
		}
		result.min = r;
		result_index = std::max(static_cast<ssize_t>(-1), result_index);
		result.index_of_min = result_index;
	}

	{
		m256 tmp;
		tmp.m256 = max;
		m256 tmp_index;
		tmp_index.m256i = index_of_max;

		float r = tmp.floatv[0];
		ssize_t result_index =
				tmp_index.intv[0] == -1 ?
						-1 : tmp_index.intv[0] * ELEMENTSOF(tmp.intv);
		for (unsigned i = 1; i < ELEMENTSOF(tmp.floatv); ++i) {
			if (!std::isnan(tmp.floatv[i])) {
				assert(tmp_index.intv[i] >= 0);
				if (std::isnan(r) || tmp.floatv[i] > r) {
					r = tmp.floatv[i];
					result_index = tmp_index.intv[i] * ELEMENTSOF(tmp.intv) + i;
				}
			}
		}
		result.max = r;
		result_index = std::max(static_cast<ssize_t>(-1), result_index);
		result.index_of_max = result_index;
	}
	int32_t start = (num_data / (sizeof(__m256 ) / sizeof(float)))
			* (sizeof(__m256 ) / sizeof(float));
	assert(num_data <= INT32_MAX);
	for (int32_t i = start; i < static_cast<int32_t>(num_data); ++i) {
		if (is_valid[i]) {
			auto const float_data = data[i];
			double const double_data = float_data;
			++counted;
			total += double_data;
			square_total += double_data * double_data;
			if (!std::isnan(float_data)) {
				if (std::isnan(result.min) || float_data < result.min) {
					result.min = float_data;
					result.index_of_min = i;
				}
				if (std::isnan(result.max) || float_data > result.max) {
					result.max = float_data;
					result.index_of_max = i;
				}
			}
		}
	}
	result.count = counted;
	result.sum = total;
	result.square_sum = square_total;
	*result_arg = result;
}

} /* anonymous namespace */

#else /* defined(__AVX__) && !defined(ARCH_SCALAR) && (! FORCE_EIGEN) */

#define EIGEN_DENSEBASE_PLUGIN "eigen_binary_visitor_plugin.h"
#include <Eigen/Core>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

template<typename Accumulator, typename Scalar, typename ScararOther,
		typename Index>
struct StatVisitor {
	StatVisitor() :
			count(0), sum(0), square_sum(0), min(NAN), max(NAN), index_of_min(
					-1), index_of_max(-1) {
	}
// called for the first coefficient
	inline bool Init(const Scalar& value_f, const ScararOther &is_valid, int i,
			int j) {
		assert(j == 0); // support only 1 dimension array
		if (!is_valid) {
			return false;
		}
		Accumulator const value = Accumulator(value_f);
		count = 1;
		sum = value;
		min = max = value_f;
		square_sum = value * value;
		index_of_min = index_of_max = i;
		return true;
	}
// called for all other coefficients
	inline void operator()(const Scalar& value_f, const ScararOther & is_valid,
			int i, int j) {
		assert(j == 0); // support only 1 dimension array
		if (is_valid) {
			Accumulator const value = Accumulator(value_f);
			++count;
			assert(!std::isnan(value));
			sum += value;
			square_sum += value * value;
			if (value_f < min) {
				min = value_f;
				index_of_min = i;
			}
			if (value_f > max) {
				max = value_f;
				index_of_max = i;
			}
		}
	}

	size_t count;
	Accumulator sum;
	Accumulator square_sum;
	Scalar min, max;
	int index_of_min, index_of_max;
};

template<typename InternalDataType, typename DataType>
inline void ComputeStatisticsEigen(size_t num_data, DataType const *data,
bool const *is_valid,
LIBSAKURA_SYMBOL(StatisticsResultFloat) *result_) {
	LIBSAKURA_SYMBOL(StatisticsResultFloat) &result = *result_;

	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	Map<Array<DataType, Dynamic, 1>, Aligned> data_(const_cast<float *>(data),
			num_data);

	assert(LIBSAKURA_SYMBOL(IsAligned)(is_valid));
	Map<Array<bool, Dynamic, 1>, Aligned> is_valid_(
			const_cast<bool *>(is_valid), num_data);

	StatVisitor<InternalDataType, DataType, bool,
			typename Map<Array<DataType, Dynamic, 1>, Aligned>::Index> visitor;
	data_.VisitWith(is_valid_, visitor);
	result.count = visitor.count;
	result.sum = visitor.sum;
	result.square_sum = visitor.square_sum;
	result.min = visitor.min;
	result.index_of_min = visitor.index_of_min;
	result.max = visitor.max;
	result.index_of_max = visitor.index_of_max;
}

} /* anonymous namespace */

#endif /* defined(__AVX__) && !defined(ARCH_SCALAR) && (! FORCE_EIGEN) */

namespace {

template<typename T, typename Result>
void ComputeStatistics(size_t num_data, T const data[],
bool const is_valid[], Result *result) {
	assert(((void)"Not yet implemented", false));
}

template<>
void ComputeStatistics<float, LIBSAKURA_SYMBOL(StatisticsResultFloat)>(
		size_t num_data, float const data[],
		bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) {
#if defined(__AVX__) && !defined(ARCH_SCALAR) && (! FORCE_EIGEN)
	ComputeStatisticsSimdFloat(num_data, data, is_valid, result);
#else
	ComputeStatisticsEigen<double, float>(num_data, data, is_valid, result);
#endif
}

template<typename T, typename Result>
void ComputeAccurateStatistics(size_t num_data, T const data[],
bool const is_valid[], Result *result) {
	assert(((void)"Not yet implemented", false));
}

template<>
void ComputeAccurateStatistics<float, LIBSAKURA_SYMBOL(StatisticsResultFloat)>(
		size_t num_data, float const data[],
		bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) {
#if defined(__AVX__) && !defined(ARCH_SCALAR) && (! FORCE_EIGEN)
	SIMDStats<float, double> reducer(data, is_valid);
#else
	ScalarStats<float, double> reducer(data, is_valid);
#endif
	auto stats = Traverse<4u, 8u, decltype(reducer)>(0u, num_data, &reducer);

	result->count = stats.count;
	result->sum = stats.sum;
	result->square_sum = stats.square_sum;
	result->min = stats.min;
	result->index_of_min = stats.index_of_min;
	result->max = stats.max;
	result->index_of_max = stats.index_of_max;
}

} /* anonymous namespace */

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

template<typename Func>
LIBSAKURA_SYMBOL(Status) ComputeStatisticsFloatGateKeeper(Func func,
		size_t num_data, float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) {
	CHECK_ARGS(num_data <= INT32_MAX);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(is_valid != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(is_valid));
	CHECK_ARGS(result != nullptr);

	try {
		func(/*num_data, data, is_valid, result*/);
	} catch (...) {
		assert(false); // No exception should be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeStatisticsFloat)(
		size_t num_data, float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) noexcept {
	return ComputeStatisticsFloatGateKeeper(
			[=] {ComputeStatistics(num_data, data, is_valid, result);},
			num_data, data, is_valid, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeAccurateStatisticsFloat)(
		size_t num_data, float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) noexcept {
	return ComputeStatisticsFloatGateKeeper(
			[=] {ComputeAccurateStatistics(num_data, data, is_valid, result);},
			num_data, data, is_valid, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeStddevFloat)(
		size_t degree_of_freedom, double mean, size_t num_data, float const data[],
		bool const is_valid[], double *result) noexcept {
	CHECK_ARGS(num_data <= INT32_MAX);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(is_valid != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(is_valid));
	CHECK_ARGS(degree_of_freedom > 0);
	CHECK_ARGS(result != nullptr);

#if 1
	double sq_diff = 0.;
	for (size_t i = 0; i < num_data; ++i) {
		if (is_valid[i]) {
			double diff = data[i] - mean;
			sq_diff += diff * diff;
		}
	}
	*result = std::sqrt(sq_diff / degree_of_freedom);
#else
	{
		size_t n = 0;
		double mean = 0;
		double M2 = 0;

		for (size_t i = 0; i < num_data; ++i) {
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
		*result = std::sqrt(M2 / (n - 1));
	}
#endif
	return LIBSAKURA_SYMBOL(Status_kOK);
}

namespace {

template<typename T, typename COMPARATOR>
void QuickSort(size_t num_data, T data[]) {
	// TODO implement quick sort using expression template to avoid overhead of calling COMPARATOR.
	qsort(data, num_data, sizeof(T),
			reinterpret_cast<int (*)(void const*,
					void const*)>(COMPARATOR::Compare));}

template<typename T>
struct AscendingOrder {
	static int Compare(T const*a, T const*b) {
		if (false) { // no branch version
			auto tmp = *a - *b;
			int sign = static_cast<int>(std::abs(tmp) / tmp);
			if (static_cast<int>(0. / 0.) == INT_MIN) {
				sign <<= 1; // make INT_MIN 0
			} else if (static_cast<int>(0. / 0.) == 0) {
			} else {
				assert(false);
			}
			return sign;
		} else {
			if (*a < *b) {
				return -1;
			} else if (*a > *b) {
				return 1;
			} else {
				return 0;
			}
		}
	}
};

}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(
		size_t num_data, bool const is_valid[], float data[],
		size_t *new_num_data) noexcept {
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(is_valid != nullptr);
	CHECK_ARGS(new_num_data != nullptr);

	try {
		size_t valid_count = 0;
		for (size_t i = 0; i < num_data; ++i) {
			if (is_valid[i]) {
				auto const the_data = data[i];
				assert(!std::isnan(the_data) && !std::isinf(the_data));
				data[valid_count] = the_data;
				++valid_count;
			}
		}

		QuickSort<float, AscendingOrder<float> >(valid_count, data);
		*new_num_data = valid_count;
	} catch (...) {
		assert(false); // No exception should be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

namespace {

template<typename T>
LIBSAKURA_SYMBOL(Status) ComputeMedianAbsoluteDeviation(size_t num_data,
		T const data[], T new_data[]) noexcept {
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(new_data != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(new_data));

	try {
		if (num_data < 1) {
			return LIBSAKURA_SYMBOL(Status_kOK);
		}
		auto data_aligned = AssumeAligned(data);
		auto const new_data_aligned = AssumeAligned(new_data);
		auto median = data_aligned[num_data / 2];
		if (num_data % 2 == 0) {
			median = (median + data_aligned[num_data / 2 - 1]) / 2;
		}
		for (size_t i = 0; i < num_data; ++i) {
			new_data_aligned[i] = std::abs(data_aligned[i] - median);
		}

		QuickSort<T, AscendingOrder<T> >(num_data, new_data);
	} catch (...) {
		assert(false); // No exception should be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeMedianAbsoluteDeviationFloat)(
		size_t num_data, float const data[], float new_data[]) noexcept {
	return ComputeMedianAbsoluteDeviation(num_data, data, new_data);
}
