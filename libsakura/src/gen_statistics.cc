#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

#define CHECK_ARGS(x) do { \
	if (!(x)) { \
		return LIBSAKURA_SYMBOL(Status_kInvalidArgument); \
	} \
} while (false)

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeStatistics)(
		size_t elements, float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResult) *result) {
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(is_valid != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(is_valid));
	CHECK_ARGS(result != nullptr);

	auto stat =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetStatisticsImpl();
	stat->ComputeStatistics(data, is_valid, elements, result);
	return LIBSAKURA_SYMBOL(Status_kOK);
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

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SortValidValuesDensely)(
		size_t elements, bool const is_valid[], float data[],
		size_t *new_elements) {
	CHECK_ARGS(data != nullptr);
	CHECK_ARGS(is_valid != nullptr);
	CHECK_ARGS(new_elements != nullptr);
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(data));
	CHECK_ARGS(LIBSAKURA_SYMBOL(IsAligned)(is_valid));

	size_t valid_count = 0;
	for (size_t i = 0; i < elements; ++i) {
		if (is_valid[i]) {
			assert(!isnanf(data[i]));
			data[valid_count] = data[i];
			++valid_count;
		}
	}

	QuickSort<float, AscendingOrder<float> >(data, valid_count);
	*new_elements = valid_count;
	return LIBSAKURA_SYMBOL(Status_kOK);
}
