#include <cassert>
#include <cstddef>
#include <cmath>
#include <algorithm>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory_impl.h"
#include "libsakura/localdef.h"

#if defined(__AVX__)
#include <immintrin.h>
#include <cstdint>

namespace {

union m256 {
	__m256 m256;
	__m256i m256i;
	__m256d m256d;
	float floatv[8];
	double doublev[4];
	int32_t intv[8];
	char charv[32];
};

inline float AddHorizontally(__m256 packed_values) {
	packed_values = _mm256_hadd_ps(packed_values, packed_values);
	packed_values = _mm256_hadd_ps(packed_values, packed_values);
	__m256 sum2 = _mm256_permute2f128_ps(packed_values, packed_values, 1);
	float total;
	_mm_store_ss(&total,
			_mm_add_ss(_mm256_castps256_ps128(packed_values),
					_mm256_castps256_ps128(sum2)));
	return total;
}

inline int32_t AddHorizontally(__m256i packed_values) {
	__m128i count2 = _mm256_castsi256_si128(
			_mm256_permute2f128_si256(packed_values, packed_values, 1) );
	__m128i count4w = _mm_hadd_epi32(_mm256_castsi256_si128(packed_values),
			count2);
	count4w = _mm_hadd_epi32(count4w, count4w);
	count4w = _mm_hadd_epi32(count4w, count4w);
	int32_t total = _mm_extract_epi32(count4w, 0);
	return total;
}

inline int32_t AddHorizontally128(__m128i packed_values) {
	packed_values = _mm_hadd_epi32(packed_values, packed_values);
	packed_values = _mm_hadd_epi32(packed_values, packed_values);
	int32_t total = _mm_extract_epi32(packed_values, 0);
	return total;
}

void ComputeStatisticsSimd(float const data[], bool const is_valid[], size_t elements,
		LIBSAKURA_SYMBOL(StatisticsResult) &result) {
	assert(sizeof(m256) == sizeof(__m256 ));
	__m256 sum = _mm256_setzero_ps();
	__m256   const zero = _mm256_setzero_ps();
	__m256i   const zeroi = _mm256_setzero_si256();
	__m128i   const zero128i = _mm_setzero_si128();
	__m256i   const onei = _mm256_set1_epi32(1);
	__m128i   const one128i = _mm_set1_epi32(1);
	__m128i count = _mm_setzero_si128();
	__m256   const nan = _mm256_set1_ps(NAN);
	__m256 min = nan;
	__m256 max = nan;
	__m256 square_sum = zero;
	float const *mask_ = reinterpret_cast<float const *>(is_valid);
	__m256   const *data_ = reinterpret_cast<__m256   const *>(data);
	for (int i = 0; i < elements / (sizeof(__m256 ) / sizeof(float)); ++i) {
		__m128i mask0 = _mm_castps_si128(_mm_load_ss(&mask_[i * 2]));
		__m128i mask1 = _mm_castps_si128(_mm_load_ss(&mask_[i * 2 + 1]));
		mask0 = _mm_cvtepu8_epi32(mask0);
		mask1 = _mm_cvtepu8_epi32(mask1);
		count = _mm_add_epi32(_mm_add_epi32(count, mask0), mask1);

		mask0 = _mm_cmpeq_epi32(mask0, zero128i);
		mask1 = _mm_cmpeq_epi32(mask1, zero128i);
		/* mask[01] == 0xffffffff means invalid data, 0 means valid data. */

		__m256i mask8 = _mm256_insertf128_si256(_mm256_castsi128_si256(mask0),
				mask1, 1);
		__m256 maskf = _mm256_cvtepi32_ps(mask8);

		__m256 value = data_[i];
		sum += _mm256_blendv_ps(value, zero, maskf);
		min = _mm256_min_ps(min, _mm256_blendv_ps(value, min, maskf));
		max = _mm256_max_ps(max, _mm256_blendv_ps(value, max, maskf));

		square_sum += _mm256_blendv_ps(_mm256_mul_ps(value, value), zero, maskf);
	}
	{

		size_t numTotal = static_cast<size_t>(AddHorizontally128(count));
		float total = AddHorizontally(sum);
		result.sum = total;
		result.count = numTotal;
	}

	{
		m256 tmp;
		tmp.m256 = min;
		float r = tmp.floatv[0];
		for (int i = 1; i < ELEMENTSOF(tmp.floatv); ++i) {
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
		for (int i = 1; i < ELEMENTSOF(tmp.floatv); ++i) {
			if (r == NAN) {
				r = tmp.floatv[i];
			} else if (tmp.floatv[i] != NAN) {
				r = std::max(r, tmp.floatv[i]);
			}
		}
		result.max = r;
	}
	float total = AddHorizontally(square_sum);
	int start = (elements / (sizeof(__m256 ) / sizeof(float)))
			* (sizeof(__m256 ) / sizeof(float));
	for (int i = start; i < elements; ++i) {
		if (is_valid[i]) {
			++result.count;
			result.sum += data[i];
			result.min = std::min(result.min, data[i]);
			result.max = std::max(result.max, data[i]);
			total += data[i] * data[i];
		}
	}
	float float_count = static_cast<float>(result.count);
	result.mean = result.count == 0 ? NAN : result.sum / float_count;

	float rms2 = result.count == 0 ? NAN : total / float_count;
	result.rms = std::sqrt(rms2);
	result.stddev = std::sqrt(rms2 - result.mean * result.mean);
}

} /* anonymous namespace */

#else /* defined(__AVX__) */

#define EIGEN_DENSEBASE_PLUGIN "eigen_binary_visitor_plugin.h"
#include <Eigen/Core>

using ::Eigen::Map;
using ::Eigen::Array;
using ::Eigen::Dynamic;
using ::Eigen::Aligned;

namespace {

template<typename Scalar, typename ScararOther, typename Index>
class StatVisitor {
public:
	StatVisitor() :
			count(0), sum(0), min(NAN), max(NAN), square_sum(0) {
	}
// called for the first coefficient
	inline bool Init(const Scalar& value, const ScararOther &is_valid, int i,
			int j) {
		if (!is_valid) {
			return false;
		}
		count = 1;
		sum = value;
		min = max = value;
		square_sum = value * value;
		return true;
	}
// called for all other coefficients
	inline void operator()(const Scalar& value, const ScararOther & is_valid,
			int i, int j) {
		if (is_valid) {
			++count;
			assert(!isnanf(value));
			sum += value;
			min = std::min(min, value);
			max = std::max(max, value);
			square_sum += value * value;
		}
	}

	size_t count;
	Scalar sum;
	Scalar min, max;
	Scalar square_sum;
};

template<typename DataType>
inline void ComputeStatisticsEigen(DataType const *data, bool const *is_valid, size_t elements,
		LIBSAKURA_SYMBOL(StatisticsResult) &result) {
	assert(false);

	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	Map<Array<DataType, Dynamic, 1>, Aligned> data_(const_cast<float *>(data),
			elements);

	assert(LIBSAKURA_SYMBOL(IsAligned)(is_valid));
	Map<Array<bool, Dynamic, 1>, Aligned> is_valid_(
			const_cast<bool *>(is_valid), elements);

	StatVisitor<DataType, bool,
			typename Map<Array<DataType, Dynamic, 1> >::Index> visitor;
	data_.VisitWith(is_valid_, visitor);
	result.count = visitor.count;
	result.sum = visitor.sum;
	result.mean = result.sum / result.count;
	result.min = visitor.min;
	result.max = visitor.max;
	float rms2 = visitor.square_sum / result.count;
	result.rms = std::sqrt(rms2);
	result.stddev = std::sqrt(rms2 - result.mean * result.mean);
	//printf("%f\n", result.mean);
}

#if 0
/* actual statistic calculations
 statDictOut["median"] = get_stats(scantableIn, channelMask, "median")
 statDictOut["max_abc"] = get_stats_pos(scantableIn, channelMask, "max_abc")
 statDictOut["min_abc"] = get_stats_pos(scantableIn, channelMask, "min_abc")
 */
#endif

} /* anonymous namespace */

#endif /* defined(__AVX__) */

namespace LIBSAKURA_PREFIX {
void ADDSUFFIX(Statistics, ARCH_SUFFIX)::ComputeStatistics(float const data[],
		bool const is_valid[], size_t elements,
		LIBSAKURA_SYMBOL(StatisticsResult) &result) const {
#if defined( __AVX__)
	ComputeStatisticsSimd(data, is_valid, elements, result);
#else
	ComputeStatisticsEigen(data, is_valid, elements, result);
#endif
}
}
