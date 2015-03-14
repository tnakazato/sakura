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
#include <cassert>
#include <cstddef>
#include <cmath>
#include <climits>
#include <algorithm>

#include "libsakura/sakura.h"
#include "libsakura/localdef.h"

#define FORCE_EIGEN 0

namespace {

template<size_t kUnitSize, size_t kBlockSize, typename Reducer>
typename Reducer::MiddleLevelAccumulator BlockWiseTraverse(size_t offset,
		int levels, size_t stride, size_t num_blocks, size_t full_blocks,
		typename Reducer::MiddleLevelAccumulator const &residual,
		Reducer *reducer) {
	//cout << "stride = " << stride << endl;
	typename Reducer::MiddleLevelAccumulator acc_stack[levels + 1][kUnitSize];
	acc_stack[0][0].clear();
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
		//cout << "level: " << level << ", i: " << indices[level] << " full_blocks: " << full_blocks << " limit: " << upper_limit << endl;
		assert(0 <= level && level < levels);
		if (level >= leaf) {
			auto &accumlator = acc[level - 1][indices[level - 1]];
			accumlator.clear();
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

template<size_t kUnitSize, size_t kBlockSize, typename Reducer>
typename Reducer::TopLevelAccumulator Traverse(size_t offset,
		size_t const num_data, Reducer *reducer) {
	STATIC_ASSERT(kUnitSize > 1 && kBlockSize > 0);
	typename Reducer::TopLevelAccumulator top_level_result;
	top_level_result.clear();

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
	residual.clear();
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
		int index_of_min, index_of_max;

		void clear() {
			count = 0;
			sum = 0.;
			square_sum = 0.;
			min = max = NAN;
			index_of_min = index_of_max = -1;
		}
	};

	typedef TopLevelAccumulator MiddleLevelAccumulator;

	void TopLevelReduce(TopLevelAccumulator const &increment,
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
				accumulator.index_of_min = increment.index_of_min;
				accumulator.min = increment.min;
			}
			if (increment.max > accumulator.max) {
				accumulator.index_of_max = increment.index_of_max;
				accumulator.max = increment.max;
			}
		}
	}

	TopLevelAccumulator TopLevelReduce(
			MiddleLevelAccumulator const &accumulator) {
		return accumulator;
	}

	void MiddleLevelReduce(MiddleLevelAccumulator const &increment,
			MiddleLevelAccumulator &accumulator) {
		TopLevelReduce(increment, accumulator);
	}

	void Reduce(size_t pos, MiddleLevelAccumulator &accumulator) {
#if 1
		if (mask[pos]) {
			++accumulator.count;
			Accumulator v = data[pos];
			accumulator.sum += v;
			accumulator.square_sum += v * v;
			if (accumulator.index_of_min < 0) {
				accumulator.index_of_min = accumulator.index_of_max = pos;
				accumulator.min = accumulator.max = data[pos];
			} else {
				if (data[pos] < accumulator.min) {
					accumulator.index_of_min = pos;
					accumulator.min = data[pos];
				}
				if (data[pos] > accumulator.max) {
					accumulator.index_of_max = pos;
					accumulator.max = data[pos];
				}
			}
		}
#else
		accumulator.sum += pos;
#endif
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

template<typename Scalar, typename Accumulator>
struct SIMDStats {
	Scalar const *data;
	bool const *is_valid;
	SIMDStats(Scalar const *data_arg, bool const *mask_arg) :
			data(data_arg), is_valid(mask_arg) {
	}
	struct SIMD_ALIGN TopLevelAccumulator {
		size_t count;
		Accumulator sum;
		Accumulator square_sum;
		Scalar min, max;
		int index_of_min, index_of_max;

		void clear() {
			count = 0;
			sum = 0.;
			square_sum = 0.;
			min = max = NAN;
			index_of_min = index_of_max = -1;
		}
	};

	struct SIMD_ALIGN MiddleLevelAccumulator {
		m256 count;
		m256 sum;
		m256 square_sum;
		m256 min, max;
		m256 index_of_min, index_of_max;

		void clear() {
			auto const zero = _mm256_setzero_ps();
			count.m256 = zero;
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
				accumulator.index_of_min = increment.index_of_min;
				accumulator.min = increment.min;
			}
			if (increment.max > accumulator.max) {
				accumulator.index_of_max = increment.index_of_max;
				accumulator.max = increment.max;
			}
		}
	}

	TopLevelAccumulator TopLevelReduce(
			MiddleLevelAccumulator const &accumulator) {
		return accumulator;
	}

	void MiddleLevelReduce(MiddleLevelAccumulator const &increment,
			MiddleLevelAccumulator &accumulator) {
		TopLevelReduce(increment, accumulator);
	}

	void Reduce(size_t pos, MiddleLevelAccumulator &accumulator) {
#if 1
		if (is_valid[pos]) {
			++accumulator.count;
			Accumulator v = data[pos];
			accumulator.sum += v;
			accumulator.square_sum += v * v;
			if (accumulator.index_of_min < 0) {
				accumulator.index_of_min = accumulator.index_of_max = pos;
				accumulator.min = accumulator.max = data[pos];
			} else {
				if (data[pos] < accumulator.min) {
					accumulator.index_of_min = pos;
					accumulator.min = data[pos];
				}
				if (data[pos] > accumulator.max) {
					accumulator.index_of_max = pos;
					accumulator.max = data[pos];
				}
			}
		}
#else
		accumulator.sum += pos;
#endif
	}

	void LastLevelReduce(size_t position, size_t block_size,
			MiddleLevelAccumulator &accumulator) {
		__m256 const *data_ = AssumeAligned(
				reinterpret_cast<__m256 const *>(data));
		__m256 const zero = _mm256_setzero_ps();
#if defined(__AVX2__)
		__m256i const zero256i = _mm256_setzero_si256();
		__m256 const nan = _mm256_set1_ps(NAN);
		double const *mask_ = AssumeAligned(
				reinterpret_cast<double const *>(is_valid));
#else
		__m128i const zero128i = _mm_setzero_si128();
		float const *mask_ = AssumeAligned(reinterpret_cast<float const *>(is_valid));
#endif
		for (size_t j = 0; j < block_size; ++j) {
			auto i = position + j;
			STATIC_ASSERT(sizeof(double) == 8);
#if defined(__AVX2__)
			__m256i mask8 = _mm256_cvtepi8_epi32(
					_mm_castpd_si128(_mm_load1_pd(&mask_[i])));
			accumulator.count = _mm256_add_epi32(accumulator.count, mask8);
			mask8 = _mm256_cmpeq_epi32(mask8, zero256i);
#else
			__m128i mask0 = _mm_castps_si128(_mm_load_ss(&mask_[i * 2]));
			__m128i mask1 = _mm_castps_si128(_mm_load_ss(&mask_[i * 2 + 1]));
			mask0 = _mm_cvtepu8_epi32(mask0);
			mask1 = _mm_cvtepu8_epi32(mask1);
			accumulator.count = _mm_add_epi32(_mm_add_epi32(accumulator.count, mask0), mask1);

			mask0 = _mm_cmpeq_epi32(mask0, zero128i);
			mask1 = _mm_cmpeq_epi32(mask1, zero128i);
			__m256i mask8 = _mm256_insertf128_si256(_mm256_castsi128_si256(mask0),
					mask1, 1);
#endif
			__m256 maskf = _mm256_cvtepi32_ps(mask8);
			/* maskf: 0xffffffff means invalid data, 0 means valid data. */

			__m256 value = data_[i];
			accumulator.min = _mm256_min_ps(accumulator.min,
					_mm256_blendv_ps(value, accumulator.min, maskf));
			accumulator.max = _mm256_max_ps(accumulator.max,
					_mm256_blendv_ps(value, accumulator.max, maskf));

			{
				__m256 value_nan = _mm256_blendv_ps(value, nan, maskf);
				__m256i   const index = _mm256_set1_epi32(i);
				accumulator.index_of_min = _mm256_castps_si256(
						_mm256_blendv_ps(_mm256_castsi256_ps(index),
								_mm256_castsi256_ps(accumulator.index_of_min),
								_mm256_cmp_ps(accumulator.min, value_nan,
										_CMP_NEQ_UQ)));
				accumulator.index_of_max = _mm256_castps_si256(
						_mm256_blendv_ps(_mm256_castsi256_ps(index),
								_mm256_castsi256_ps(accumulator.index_of_max),
								_mm256_cmp_ps(accumulator.max, value_nan,
										_CMP_NEQ_UQ)));
			}

			value = _mm256_blendv_ps(value, zero, maskf);
			__m256d v = _mm256_cvtps_pd(_mm256_castps256_ps128(value));
			accumulator.sum = _mm256_add_pd(accumulator.sum, v);
			accumulator.square_sum = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
			LIBSAKURA_SYMBOL(SimdPacketAVX), double>(v, v,
					accumulator.square_sum);

			v = _mm256_cvtps_pd(_mm256_extractf128_ps(value, 1));
			accumulator.sum = _mm256_add_pd(accumulator.sum, v);
			accumulator.square_sum = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
			LIBSAKURA_SYMBOL(SimdPacketAVX), double>(v, v,
					accumulator.square_sum);
		}
	}

	void SequentialReduce(size_t offset, size_t elements,
			MiddleLevelAccumulator &accumulator) {
		for (size_t i = 0; i < elements; ++i) {
			Reduce(offset + i, accumulator);
		}
	}
};

inline float AddHorizontally(__m256 packed_values) {
	packed_values = _mm256_hadd_ps(packed_values, packed_values);
	packed_values = _mm256_hadd_ps(packed_values, packed_values);
	__m128 sum2 = _mm256_extractf128_ps(packed_values, 1);
	return _mm_cvtss_f32(
			_mm_add_ss(_mm256_castps256_ps128(packed_values), sum2));
}

inline double AddHorizontally(__m256d packed_values) {
	packed_values = _mm256_hadd_pd(packed_values, packed_values);
	__m128d sum2 = _mm256_extractf128_pd(packed_values, 1);
	return _mm_cvtsd_f64(
			_mm_add_sd(_mm256_castpd256_pd128(packed_values), sum2));
}

inline int32_t AddHorizontally(__m256i packed_values) {
#if defined(__AVX2__)
	packed_values = _mm256_hadd_epi32(packed_values, packed_values);
	packed_values = _mm256_hadd_epi32(packed_values, packed_values);
	__m128i sum2 = _mm256_extractf128_si256(packed_values, 1);
	return _mm_cvtsi128_si32(
			_mm_add_epi32(_mm256_castsi256_si128(packed_values), sum2));
#else
	__m128i count2 = _mm256_castsi256_si128(
			_mm256_permute2f128_si256(packed_values, packed_values, 1));
	__m128i count4w = _mm_hadd_epi32(_mm256_castsi256_si128(packed_values),
			count2);
	count4w = _mm_hadd_epi32(count4w, count4w);
	count4w = _mm_hadd_epi32(count4w, count4w);
	int32_t total = _mm_extract_epi32(count4w, 0);
	return total;
#endif
}

inline int32_t AddHorizontally128(__m128i packed_values) {
	packed_values = _mm_hadd_epi32(packed_values, packed_values);
	packed_values = _mm_hadd_epi32(packed_values, packed_values);
	int32_t total = _mm_extract_epi32(packed_values, 0);
	return total;
}

void ComputeStatisticsSimdFloat(float const data[], bool const is_valid[],
		size_t elements, LIBSAKURA_SYMBOL(StatisticsResultFloat) *result_arg) {
	STATIC_ASSERT(sizeof(m256) == sizeof(__m256 ));
	__m256d const zero_d = _mm256_setzero_pd();
	__m256d sum = zero_d;
	__m256d square_sum = zero_d;
	__m256 const zero = _mm256_setzero_ps();
	__m256 const nan = _mm256_set1_ps(NAN);
	__m256 min = nan;
	__m256 max = nan;
	__m256i index_of_min = _mm256_set1_epi32(-1);
	__m256i index_of_max = _mm256_set1_epi32(-1);
#if defined(__AVX2__)
	__m256i const zero256i = _mm256_setzero_si256();
	__m256i count = zero256i;
	double const *mask_ = AssumeAligned(
			reinterpret_cast<double const *>(is_valid));
#else
	__m128i const zero128i = _mm_setzero_si128();
	__m128i count = zero128i;
	float const *mask_ = AssumeAligned(reinterpret_cast<float const *>(is_valid));
#endif
	__m256 const *data_ = AssumeAligned(reinterpret_cast<__m256  const *>(data));
	for (size_t i = 0; i < elements / (sizeof(__m256 ) / sizeof(float)); ++i) {
#if defined(__AVX2__)
		__m256i mask8 = _mm256_cvtepi8_epi32(
				_mm_castpd_si128(_mm_load1_pd(&mask_[i])));
		count = _mm256_add_epi32(count, mask8);
		mask8 = _mm256_cmpeq_epi32(mask8, zero256i);
#else
		__m128i mask0 = _mm_castps_si128(_mm_load_ss(&mask_[i * 2]));
		__m128i mask1 = _mm_castps_si128(_mm_load_ss(&mask_[i * 2 + 1]));
		mask0 = _mm_cvtepu8_epi32(mask0);
		mask1 = _mm_cvtepu8_epi32(mask1);
		count = _mm_add_epi32(_mm_add_epi32(count, mask0), mask1);

		mask0 = _mm_cmpeq_epi32(mask0, zero128i);
		mask1 = _mm_cmpeq_epi32(mask1, zero128i);
		__m256i mask8 = _mm256_insertf128_si256(_mm256_castsi128_si256(mask0),
				mask1, 1);
#endif
		__m256 maskf = _mm256_cvtepi32_ps(mask8);
		/* maskf: 0xffffffff means invalid data, 0 means valid data. */

		__m256 value = data_[i];
		min = _mm256_min_ps(min, _mm256_blendv_ps(value, min, maskf));
		max = _mm256_max_ps(max, _mm256_blendv_ps(value, max, maskf));

		{
			__m256 value_nan = _mm256_blendv_ps(value, nan, maskf);
			__m256i  const index = _mm256_set1_epi32(i);
			index_of_min = _mm256_castps_si256(
					_mm256_blendv_ps(_mm256_castsi256_ps(index),
							_mm256_castsi256_ps(index_of_min),
							_mm256_cmp_ps(min, value_nan, _CMP_NEQ_UQ)));
			index_of_max = _mm256_castps_si256(
					_mm256_blendv_ps(_mm256_castsi256_ps(index),
							_mm256_castsi256_ps(index_of_max),
							_mm256_cmp_ps(max, value_nan, _CMP_NEQ_UQ)));
		}

		value = _mm256_blendv_ps(value, zero, maskf);
		__m256d v = _mm256_cvtps_pd(_mm256_castps256_ps128(value));
		sum = _mm256_add_pd(sum, v);
		square_sum = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
				LIBSAKURA_SYMBOL(SimdPacketAVX), double>(v, v, square_sum);

		v = _mm256_cvtps_pd(_mm256_extractf128_ps(value, 1));
		sum = _mm256_add_pd(sum, v);
		square_sum = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
				LIBSAKURA_SYMBOL(SimdPacketAVX), double>(v, v, square_sum);
	}
#if defined(__AVX2__)
	size_t counted = static_cast<size_t>(AddHorizontally(count));
#else
	size_t counted = static_cast<size_t>(AddHorizontally128(count));
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
		int result_index = tmp_index.intv[0] * ELEMENTSOF(tmp.intv);
		for (unsigned i = 1; i < ELEMENTSOF(tmp.floatv); ++i) {
			if (std::isnan(r)
					|| (!std::isnan(tmp.floatv[i]) && tmp.floatv[i] < r)) {
				r = tmp.floatv[i];
				result_index = tmp_index.intv[i] * ELEMENTSOF(tmp.intv) + i;
			}
		}
		result.min = r;
		result_index = std::max(-1, result_index);
		result.index_of_min = result_index;
	}

	{
		m256 tmp;
		tmp.m256 = max;
		m256 tmp_index;
		tmp_index.m256i = index_of_max;

		float r = tmp.floatv[0];
		int result_index = tmp_index.intv[0] * ELEMENTSOF(tmp.intv);
		for (unsigned i = 1; i < ELEMENTSOF(tmp.floatv); ++i) {
			if (std::isnan(r)
					|| (!std::isnan(tmp.floatv[i]) && tmp.floatv[i] > r)) {
				r = tmp.floatv[i];
				result_index = tmp_index.intv[i] * ELEMENTSOF(tmp.intv) + i;
			}
		}
		result.max = r;
		result_index = std::max(-1, result_index);
		result.index_of_max = result_index;
	}
	int start = (elements / (sizeof(__m256 ) / sizeof(float)))
			* (sizeof(__m256 ) / sizeof(float));
	assert(elements <= INT32_MAX);
	for (int32_t i = start; i < static_cast<int32_t>(elements); ++i) {
		if (is_valid[i]) {
			++counted;
			total += data[i];
			square_total += double(data[i]) * double(data[i]);
			if (std::isnan(result.min)
					|| (!std::isnan(data[i]) && data[i] < result.min)) {
				result.min = data[i];
				result.index_of_min = i;
			}
			if (std::isnan(result.max)
					|| (!std::isnan(data[i]) && data[i] > result.max)) {
				result.max = data[i];
				result.index_of_max = i;
			}
		}
	}
	result.count = counted;
	result.sum = total;
	double double_count = static_cast<double>(result.count);
	double mean = NAN;
	double rms2 = NAN;
	if (result.count != 0) {
		mean = total / double_count;
		rms2 = square_total / double_count;
	}
	result.mean = mean;
	result.rms = std::sqrt(rms2);
	result.stddev = std::sqrt(std::abs(rms2 - mean * mean));
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
	class StatVisitor {
	public:
		StatVisitor() :
		count(0), sum(0), square_sum(0), min(NAN), max(NAN), index_of_min(
				-1), index_of_max(-1) {
		}
// called for the first coefficient
		inline bool Init(const Scalar& value_f, const ScararOther &is_valid,
				int i, int j) {
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
		inline void operator()(const Scalar& value_f,
				const ScararOther & is_valid, int i, int j) {
			assert(j == 0); // support only 1 dimension array
			if (is_valid) {
				Accumulator const value = Accumulator(value_f);
				++count;
				assert(!std::isnan(value));
				sum += value;
				square_sum += Accumulator(value) * Accumulator(value);
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
	inline void ComputeStatisticsEigen(DataType const *data, bool const *is_valid,
			size_t elements,
			LIBSAKURA_SYMBOL(StatisticsResultFloat) *result_) {
		LIBSAKURA_SYMBOL(StatisticsResultFloat) &result = *result_;

		assert(LIBSAKURA_SYMBOL(IsAligned)(data));
		Map<Array<DataType, Dynamic, 1>, Aligned> data_(const_cast<float *>(data),
				elements);

		assert(LIBSAKURA_SYMBOL(IsAligned)(is_valid));
		Map<Array<bool, Dynamic, 1>, Aligned> is_valid_(
				const_cast<bool *>(is_valid), elements);

		StatVisitor<InternalDataType, DataType, bool,
		typename Map<Array<DataType, Dynamic, 1>, Aligned>::Index> visitor;
		data_.VisitWith(is_valid_, visitor);
		result.count = visitor.count;
		result.sum = visitor.sum;
		result.min = visitor.min;
		result.index_of_min = visitor.index_of_min;
		result.max = visitor.max;
		result.index_of_max = visitor.index_of_max;
		InternalDataType mean = NAN;
		InternalDataType rms2 = NAN;
		if (visitor.count != 0) {
			mean = visitor.sum / result.count;
			rms2 = visitor.square_sum / result.count;
		}
		result.mean = mean;
		result.rms = std::sqrt(rms2);
		result.stddev = std::sqrt(std::abs(rms2 - mean * mean));
	}

} /* anonymous namespace */

#endif /* defined(__AVX__) && !defined(ARCH_SCALAR) && (! FORCE_EIGEN) */

namespace {

template<typename T, typename Result>
void ComputeStatistics(T const data[],
bool const is_valid[], size_t elements, Result *result) {
	assert(("Not yet implemented", false));
}

template<>
void ComputeStatistics<float, LIBSAKURA_SYMBOL(StatisticsResultFloat)>(
		float const data[],
		bool const is_valid[], size_t elements,
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) {
#if defined( __AVX__) && !defined(ARCH_SCALAR) && (! FORCE_EIGEN)
	ComputeStatisticsSimdFloat(data, is_valid, elements, result);
#else
	ComputeStatisticsEigen<double, float>(data, is_valid, elements, result);
#endif
}

template<typename T, typename Result>
void ComputeAccurateStatistics(T const data[],
bool const is_valid[], size_t elements, Result *result) {
	assert(("Not yet implemented", false));
}

template<>
void ComputeAccurateStatistics<float, LIBSAKURA_SYMBOL(StatisticsResultFloat)>(
		float const data[],
		bool const is_valid[], size_t elements,
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) {
	ScalarStats<float, double> reducer(data, is_valid);
	auto stats = Traverse<4u, 8u, decltype(reducer)>(0u, elements, &reducer);

	result->count = stats.count;
	result->sum = stats.sum;
	result->min = stats.min;
	result->index_of_min = stats.index_of_min;
	result->max = stats.max;
	result->index_of_max = stats.index_of_max;
	double mean = NAN;
	double rms2 = NAN;
	if (stats.count != 0) {
		mean = stats.sum / stats.count;
		rms2 = stats.square_sum / stats.count;
	}
	result->mean = mean;
	result->rms = std::sqrt(rms2);
	result->stddev = std::sqrt(std::abs(rms2 - mean * mean));

}
} /* anonymous namespace */

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

template<typename Func>
LIBSAKURA_SYMBOL(Status) ComputeStatisticsFloatGateKeeper(
		Func func,
		size_t elements, float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) {
	CHECK_ARGS(elements <= INT32_MAX);
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(is_valid != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(is_valid));
	CHECK_ARGS(result != nullptr);

	try {
		func(/*data, is_valid, elements, result*/);
	} catch (...) {
		assert(false); // no exception should not be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeStatisticsFloat)(
		size_t elements, float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) {
	return ComputeStatisticsFloatGateKeeper(
			[=] {ComputeStatistics(data, is_valid, elements, result);}, elements,
			data, is_valid, result);
}

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeAccurateStatisticsFloat)(
		size_t elements, float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) {
	return ComputeStatisticsFloatGateKeeper(
			[=] {ComputeAccurateStatistics(data, is_valid, elements, result);}, elements,
			data, is_valid, result);
}

namespace {

template<typename T, typename COMPARATOR>
void QuickSort(T data[], size_t elements) {
	// TODO implement quick sort using expression template to avoid overhead of calling COMPARATOR.
	qsort(data, elements, sizeof(T),
			reinterpret_cast<int (*)(void const*,
					void const*)>(COMPARATOR::Compare));}

template<typename T>
class AscendingOrder {
public:
	static int Compare(T const*a, T const*b) {
		if (false) {
			auto tmp = *a - *b;
			int sign = static_cast<int>(std::abs(tmp) / tmp);
			if (static_cast<int>(0. / 0.) == INT_MIN) {
				sign <<= 1;
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
		size_t elements, bool const is_valid[], float data[],
		size_t *new_elements) {
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(is_valid != nullptr);
	CHECK_ARGS(new_elements != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(is_valid));

	try {
		size_t valid_count = 0;
		for (size_t i = 0; i < elements; ++i) {
			if (is_valid[i]) {
				assert(!std::isnan(data[i]));
				data[valid_count] = data[i];
				++valid_count;
			}
		}

		QuickSort<float, AscendingOrder<float> >(data, valid_count);
		*new_elements = valid_count;
	} catch (...) {
		assert(false); // no exception should not be raised for the current implementation.
		return LIBSAKURA_SYMBOL(Status_kUnknownError);
	}
	return LIBSAKURA_SYMBOL(Status_kOK);
}
