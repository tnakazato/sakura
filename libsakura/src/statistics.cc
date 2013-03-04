#include <cassert>
#include <cmath>
#include <algorithm>
#include <libsakura/sakura.h>
#include <libsakura/OptimizedImplementationFactoryImpl.h>
#define EIGEN_DENSEBASE_PLUGIN "EigenBinaryVisitorPlugin.h"
#include <Eigen/Core>
#include <libsakura/localdef.h>

using namespace Eigen;

namespace {

template<typename Scalar, typename Index>
struct UnconditionalVisitor {
	size_t count;
	Scalar sum;
	Scalar min, max;
	Scalar sqSum;
	UnconditionalVisitor() :
			count(0), sum(0), sqSum(0) {
	}
// called for the first coefficient
	inline void init(const Scalar& value, int i, int j) {
		count++;
		sum += value;
		min = max = value;
		sqSum += value * value;
	}
// called for all other coefficients
	inline void operator()(const Scalar& value, int i, int j) {
		count++;
		sum += value;
		min = std::min(min, value);
		max = std::max(max, value);
		sqSum += value * value;
	}
};

template<typename Scalar, typename ScararOther, typename Index>
struct StatVisitor {
	size_t count;
	Scalar sum;
	Scalar min, max;
	Scalar sqSum;
	StatVisitor() :
			count(0), sum(0), min(NAN), max(NAN), sqSum(0) {
	}
// called for the first coefficient
	inline void init(const Scalar& value, const ScararOther &valueOther, int i,
			int j) {
		count += valueOther ? 1 : 0;
		sum += valueOther ? value : 0;
		min = valueOther ? (min == NAN ? value : std::min(min, value)) : min;
		max = valueOther ? (max == NAN ? value : std::max(max, value)) : max;
		sqSum += valueOther ? value * value : 0;
	}
// called for all other coefficients
	inline void operator()(const Scalar& value, const ScararOther & valueOther,
			int i, int j) {
		count += valueOther ? 1 : 0;
		sum += valueOther ? value : 0;
		min = valueOther ? (min == NAN ? value : std::min(min, value)) : min;
		max = valueOther ? (max == NAN ? value : std::max(max, value)) : max;
		sqSum += valueOther ? value * value : 0;
	}
};

template<typename Scalar, typename ScararOther, typename Index>
struct StatVisitor2 {
	size_t count;
	Scalar sum;
	Scalar min, max;
	Scalar sqSum;
	StatVisitor2() :
			count(0), sum(0), min(NAN), max(NAN), sqSum(0) {
	}
// called for the first coefficient
	inline void init(const Scalar& value, const ScararOther &valueOther, int i,
			int j) {
		(*this)(value, valueOther, i, j);
	}
// called for all other coefficients
	inline void operator()(const Scalar& value, const ScararOther & valueOther,
			int i, int j) {
		int tmp = static_cast<int>(valueOther);
		count += tmp;
		float tmpf = static_cast<float>(tmp);
		sum += value * tmpf;
		min = valueOther ? (min == NAN ? value : std::min(min, value)) : min;
		max = valueOther ? (max == NAN ? value : std::max(max, value)) : max;
		sqSum += value * value * tmpf;
	}
};

template<typename DataType>
inline void stat(libsakura_symbol(statistics_result) &result,
		DataType const *data, bool const *mask, size_t elements) {
	assert(false);

	assert(libsakura_symbol(is_aligned)(data));
	Map<Array<DataType, Dynamic, 1>, Aligned> data_(const_cast<float *>(data),
			elements);

	assert(libsakura_symbol(is_aligned)(mask));
	Map<Array<bool, Dynamic, 1>, Aligned> mask_(const_cast<bool *>(mask),
			elements);

#if 0 // SimpleMean
#warning Simple Stats
	UnconditionalVisitor<DataType,
	typename Map<Array<DataType, Dynamic, 1> >::Index> visitor;
	data_.visit(visitor);
#else // masked mean
#warning Masked Stats
	StatVisitor<DataType, bool,
	typename Map<Array<DataType, Dynamic, 1> >::Index> visitor;
	data_.visitWith(mask_, visitor);
#endif
	result.count = visitor.count;
	result.sum = visitor.sum;
	result.mean = result.count == 0 ? NAN : result.sum / result.count;
	result.min = visitor.min;
	result.max = visitor.max;
	float rms2 = result.count == 0 ? NAN : visitor.sqSum / result.count;
	result.rms = std::sqrt(rms2);
	result.stddev = std::sqrt(rms2 - result.mean * result.mean);
	//printf("%f\n", result.mean);
}

/*
 template<>
 void sdstat<float>(float const *data_, char const *mask, size_t elements) {
 Map<ArrayXf> data(const_cast<float *>(data_), elements);
 }
 */

#if 0
// actual statistic calculations
statDictOut["max"] = get_stats(scantableIn, channelMask, "max")
statDictOut["min"] = get_stats(scantableIn, channelMask, "min")
statDictOut["sum"] = get_stats(scantableIn, channelMask, "sum")
statDictOut["mean"] = get_stats(scantableIn, channelMask, "mean")
statDictOut["median"] = get_stats(scantableIn, channelMask, "median")
statDictOut["rms"] = get_stats(scantableIn, channelMask, "rms")
statDictOut["stddev"] = get_stats(scantableIn, channelMask, "stddev")
statDictOut["max_abc"] = get_stats_pos(scantableIn, channelMask, "max_abc")
statDictOut["min_abc"] = get_stats_pos(scantableIn, channelMask, "min_abc")
#endif

#if defined( __AVX__)
#include <immintrin.h>
#include <cstdint>
#warning Optimized for SIMD

union m256 {
	__m256 m256;
	__m256i m256i;
	__m256d m256d;
	float floatv[8];
	double doublev[4];
	int32_t intv[8];
	char charv[32];
};

inline float add_horizontally(__m256 packedValues) {
	packedValues = _mm256_hadd_ps(packedValues, packedValues);
	packedValues = _mm256_hadd_ps(packedValues, packedValues);
	__m256 sum2 = _mm256_permute2f128_ps(packedValues, packedValues, 1);
	float total;
	_mm_store_ss(&total, _mm_add_ss(_mm256_castps256_ps128(packedValues), _mm256_castps256_ps128(sum2)));
	return total;
}

inline int32_t add_horizontally(__m256i packedValues) {
	__m128i count2 = _mm256_castsi256_si128(_mm256_permute2f128_si256(packedValues, packedValues, 1));
	__m128i count4w = _mm_hadd_epi32(_mm256_castsi256_si128(packedValues), count2);
	count4w = _mm_hadd_epi32(count4w, count4w);
	count4w = _mm_hadd_epi32(count4w, count4w);
	int32_t total;
	_mm_store_ss((float*)&total, (__m128)count4w);
	return total;
}

inline int32_t add_horizontally_128(__m128i packedValues) {
	packedValues = _mm_hadd_epi32(packedValues, packedValues);
	packedValues = _mm_hadd_epi32(packedValues, packedValues);
	int32_t total;
	_mm_store_ss((float*)&total, (__m128)packedValues);
	return total;
}

void simd(libsakura_symbol(statistics_result) &result,float const *data, bool const *mask, size_t elements) {
	assert(elements % sizeof(__m256) == 0); // TODO support arbitrary elements
	assert(sizeof(m256) == sizeof(__m256));

	__m256 sum = _mm256_setzero_ps();
	__m256 zero = _mm256_setzero_ps();
	__m256i zeroi = _mm256_setzero_si256();
	__m256i onei = _mm256_set1_epi32(1);
	__m128i count = _mm_setzero_si128();
	__m256 nan = _mm256_set1_ps(NAN);
	__m256 min = nan;
	__m256 max = nan;
	__m256 sqSum = zero;
	float const *mask_ = (float const *)mask;
	__m256 const *data_ = (__m256 const *)data;
	for (int i = 0; i < elements/(sizeof(__m256)/sizeof(float)); i++) {
		__m128i mask0 = _mm_castps_si128(_mm_load_ss(&mask_[i*2]));
		__m128i mask1 = _mm_castps_si128(_mm_load_ss(&mask_[i*2+1]));
		__m128i first4 = _mm_cvtepu8_epi32(mask0); // first 4bytes
		__m128i second4 = _mm_cvtepu8_epi32(mask1);// second 4bytes
		count = _mm_add_epi32(_mm_add_epi32(first4, second4), count);

		__m256i mask8 = _mm256_insertf128_si256(_mm256_castsi128_si256(first4),
				second4, 1);
		__m256 maskf = _mm256_cmp_ps(_mm256_cvtepi32_ps(mask8),
				zero, _CMP_NEQ_UQ);
		__m256 value = data_[i];
		sum += _mm256_blendv_ps(zero, value, maskf);
		min = _mm256_min_ps(min, _mm256_blendv_ps(min, value, maskf));
		max = _mm256_max_ps(max, _mm256_blendv_ps(max, value, maskf));

		sqSum += _mm256_blendv_ps(zero, _mm256_mul_ps(value, value), maskf);
	}
	{

		size_t numTotal = static_cast<size_t>(add_horizontally_128(count));
		float total = add_horizontally(sum);
		result.sum = total;
		result.count = numTotal;
		if (result.count == 0) {
			result.mean = NAN;
		} else {
			result.mean = result.sum / result.count;
		}
	}

	{
		m256 tmp;
		tmp.m256 = min;
		float r = tmp.floatv[0];
		for (int i = 1; i < elementsof(tmp.floatv); i++) {
			if (r == NAN) {
				r = tmp.floatv[i];
			} else if (tmp.floatv[i] != NAN) {
				r = std::min(r, tmp.floatv[i]);
			}
		}
		result.min = r;
	}

	{
		m256 tmp;
		tmp.m256 = max;
		float r = tmp.floatv[0];
		for (int i = 1; i < elementsof(tmp.floatv); i++) {
			if (r == NAN) {
				r = tmp.floatv[i];
			} else if (tmp.floatv[i] != NAN) {
				r = std::max(r, tmp.floatv[i]);
			}
		}
		result.max = r;
	}

	float total = add_horizontally(sqSum);
	float rms2 = result.count == 0 ? NAN : total / result.count;
	result.rms = std::sqrt(rms2);
	result.stddev = std::sqrt(rms2 - result.mean * result.mean);
}
#endif
}

#include <cstdio>
namespace libsakura_PREFIX {
void ADDSUFFIX(Statistics, ARCH_SUFFIX)::reduce(
	libsakura_symbol(statistics_result) &result, float const *data,
	bool const *mask, size_t elements) const {
#if 0
stat(result, data, mask, elements);
//printf("%f\n", result);
#else
#if defined( __AVX__)
simd(result, data, mask, elements);
#endif
#endif
}
}
