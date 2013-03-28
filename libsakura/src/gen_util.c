/*
 * util.c
 *
 *  Created on: 2013/02/27
 *      Author: kohji
 */

#include <sys/time.h>
#include <stdint.h>

#include "libsakura/sakura.h"

double LIBSAKURA_SYMBOL(GetCurrentTime)() {
	struct timeval tv;
	int result = gettimeofday(&tv, NULL );
	return tv.tv_sec + ((double) tv.tv_usec) / 1000000.;
}

/*
 * ALIGNMENT must be a power of 2.
 */
#define ALIGNMENT (256/* avx 256bits */ / 8)

size_t LIBSAKURA_SYMBOL (GetAlignment)() {
	return ALIGNMENT;
}

bool LIBSAKURA_SYMBOL(IsAligned)(void const *ptr) {
	uint64_t addr = (uint64_t) ptr;
	return addr % ALIGNMENT == 0;
}

void const *LIBSAKURA_SYMBOL(AlignAny)(void const *vp, size_t size_of_arena,
		size_t size_required) {
	assert(vp != NULL );
	if (vp == NULL ) {
		return NULL ;
	}
	uint64_t addr = (uint64_t) vp;
	uint64_t new_addr = (addr + (ALIGNMENT - 1)) & ~(ALIGNMENT - 1);
	if (size_of_arena - (size_t) (new_addr - addr) < size_required) {
		return NULL ;
	}
	return (void const *) new_addr;
}

float const *LIBSAKURA_SYMBOL(AlignFloat)(float const *fp,
		size_t elements_of_arena, size_t elements_required) {
	return LIBSAKURA_SYMBOL(AlignAny)(fp, elements_of_arena * sizeof(float),
			elements_required * sizeof(float));
}

double const *LIBSAKURA_SYMBOL(AlignDouble)(double const *dp,
		size_t elements_of_arena, size_t elements_required) {
	return LIBSAKURA_SYMBOL(AlignAny)(dp, elements_of_arena * sizeof(double),
			elements_required * sizeof(double));
}
