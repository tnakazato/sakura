#include <cassert>
#include <cstdlib>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" void LIBSAKURA_SYMBOL(ComputeStatistics)(float const data[],
		bool const is_valid[], size_t elements,
		LIBSAKURA_SYMBOL(StatisticsResult) *result) {
	assert(data != nullptr);
	assert(is_valid != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(is_valid));
	assert(result != nullptr);

	auto stat =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetStatisticsImpl();
	stat->ComputeStatistics(data, is_valid, elements, result);
}

namespace {

template <typename T>
class AscendingOrder {
public:
	static int Compare(T const*a, T const*b) {
		// TODO try to optimize below to eliminate branch
		if (*a < *b) {
			return -1;
		} else		if (*a > *b) {
			return 1;
		} else {
			return 0;
		}
	}
};

template <typename T, typename COMPARATOR >
void QuickSort(T data[], size_t elements) {
	// TODO implement so that calling COMPARATOR inline for better performance.
	qsort(data, elements, sizeof(T), reinterpret_cast<int (*)(void const*, void const*)>(COMPARATOR::Compare));
}

}

extern "C" size_t LIBSAKURA_SYMBOL(SortValidValuesDensely)(
		bool const is_valid[], size_t elements, float data[]) {
	assert(data != nullptr);
	assert(is_valid != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(is_valid));

	size_t valid_count = 0;
	for (size_t i = 0; i < elements; ++i) {
		if (is_valid[i]) {
			data[valid_count] = data[i];
			++valid_count;
		}
	}

	QuickSort<float, AscendingOrder<float> >(data, valid_count);
	return valid_count;
}
