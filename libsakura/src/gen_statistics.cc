#include <cassert>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"
#include "libsakura/localdef.h"

extern "C" void LIBSAKURA_SYMBOL(ComputeStatistics)(float const data[],
		bool const is_valid[], size_t elements,
		LIBSAKURA_SYMBOL(StatisticsResult) *result) {
	assert(result != nullptr);
	auto stat =
			::LIBSAKURA_PREFIX::OptimizedImplementationFactory::GetFactory()->GetStatisticsImpl();
	stat->ComputeStatistics(data, is_valid, elements, *result);
}
