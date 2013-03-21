#include <cassert>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

namespace {
void test() {
	static float data[128] __attribute__((aligned(32)));
	static bool is_invalid[128] __attribute__((aligned(32)));
	LIBSAKURA_SYMBOL(StatisticsResult) result;
	return LIBSAKURA_SYMBOL(ComputeStatistics)(data, is_invalid, ELEMENTSOF(data), &result);
}
}

extern "C" void LIBSAKURA_SYMBOL(ComputeStatistics)(float const data[],
		bool const is_invalid[], size_t elements,
		LIBSAKURA_SYMBOL(StatisticsResult) *result) {
	assert(result != nullptr);
	auto stat =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetStatisticsImpl();
	stat->Reduce(data, is_invalid, elements, *result);
}
