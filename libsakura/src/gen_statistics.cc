#include <cassert>
#include <cstdlib>
#include <climits>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" void LIBSAKURA_SYMBOL(ComputeStatistics)(size_t elements,
		float const data[], bool const is_valid[],
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

extern "C" size_t LIBSAKURA_SYMBOL(SortValidValuesDensely)(size_t elements,
		bool const is_valid[], float data[]) {
	assert(data != nullptr);
	assert(is_valid != nullptr);
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(is_valid));

	size_t valid_count = 0;
	for (size_t i = 0; i < elements; ++i) {
		if (is_valid[i]) {
			assert(!isnanf(data[i]));
			data[valid_count] = data[i];
			++valid_count;
		}
	}

	QuickSort<float, AscendingOrder<float> >(data, valid_count);
	return valid_count;
}
