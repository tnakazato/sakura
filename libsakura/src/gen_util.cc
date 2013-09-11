/*
 * util.c
 *
 *  Created on: 2013/02/27
 *      Author: kohji
 */

#include <string>
#include <algorithm>
#include <sys/time.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "libsakura/sakura.h"
#include "libsakura/optimized_implementation_factory.h"

extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Initialize)() {
	std::string env_val_name(LIBSAKURA_PREFIX_STRING "_SIMD");
	std::transform(env_val_name.begin(), env_val_name.end(),
			env_val_name.begin(), ::toupper);
	char const *simd_spec = getenv(env_val_name.c_str());
	if (simd_spec == nullptr) {
		simd_spec = "adaptive";
	}
	LIBSAKURA_PREFIX::OptimizedImplementationFactory::InitializeFactory(
			simd_spec);

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" void LIBSAKURA_SYMBOL(CleanUp)() {
	LIBSAKURA_PREFIX::OptimizedImplementationFactory::CleanUpFactory();
}

extern "C" double LIBSAKURA_SYMBOL(GetCurrentTime)() {
	struct timeval tv;
	int result = gettimeofday(&tv, NULL);
	(void) (0 && result); // to avoid warning
	return tv.tv_sec + ((double) tv.tv_usec) / 1000000.;
}

/*
 * ALIGNMENT must be a power of 2.
 */
#define ALIGNMENT (256/* avx 256bits */ / 8)

extern "C" size_t LIBSAKURA_SYMBOL (GetAlignment)() {
	return ALIGNMENT;
}

extern "C" bool LIBSAKURA_SYMBOL(IsAligned)(void const *ptr) {
	uint64_t addr = (uint64_t) ptr;
	return addr % ALIGNMENT == 0;
}

extern "C" void *LIBSAKURA_SYMBOL(AlignAny)(size_t size_of_arena, void *vp,
		size_t size_required) {
	if (vp == nullptr) {
		return nullptr;
	}
	uint64_t addr = (uint64_t) vp;
	uint64_t new_addr = (addr + (ALIGNMENT - 1)) & ~(ALIGNMENT - 1);
	if (size_of_arena - (size_t) (new_addr - addr) < size_required) {
		return nullptr;
	}
	return (void *) new_addr;
}

extern "C" float *LIBSAKURA_SYMBOL(AlignFloat)(size_t elements_in_arena,
		float *fp, size_t elements_required) {
	return static_cast<float *>(LIBSAKURA_SYMBOL(AlignAny)(
			elements_in_arena * sizeof(float), fp,
			elements_required * sizeof(float)));
}

extern "C" double *LIBSAKURA_SYMBOL(AlignDouble)(size_t elements_in_arena,
		double *dp, size_t elements_required) {
	return static_cast<double *>(LIBSAKURA_SYMBOL(AlignAny)(
			elements_in_arena * sizeof(double), dp,
			elements_required * sizeof(double)));
}
