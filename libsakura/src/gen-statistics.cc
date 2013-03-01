#include <cassert>
#include <libsakura/OptimizedImplementationFactory.h>
#include <libsakura/sakura.h>
#include <libsakura/localdef.h>

using namespace std;
using namespace libsakura_PREFIX;

namespace {
void test() {
	static float data[128] __attribute__((aligned(32)));
	static bool mask[128] __attribute__((aligned(32)));
	libsakura_symbol(statistics_result) result;
	return libsakura_symbol(statistics)(&result, data, mask, elementsof(data));
}
}

extern "C"
void libsakura_symbol(statistics)(libsakura_symbol(statistics_result) *result,
		float const data[], bool const mask[],
		size_t elements) {
	assert(result != nullptr);
	Statistics const *stat =
			OptimizedImplementationFactory::getFactory()->getStatisticsImpl();
	stat->reduce(*result, data, mask, elements);
}

