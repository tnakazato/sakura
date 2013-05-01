/*
 * util.c
 *
 *  Created on: 2013/02/27
 *      Author: kohji
 */

#include <sys/time.h>
#include <assert.h>
#include <stdint.h>

#include "libsakura/sakura.h"

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Initialize)() {
	/* any code to initialize libsakura */

	return LIBSAKURA_SYMBOL(Status_kOK);
}

void LIBSAKURA_SYMBOL(CleanUp)() {
	/* any code to cleanup libsakura */
}

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

void const *LIBSAKURA_SYMBOL(AlignAny)(size_t size_of_arena, void const *vp,
		size_t size_required) {
	assert(vp != NULL);
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

float const *LIBSAKURA_SYMBOL(AlignFloat)(size_t elements_in_arena,
		float const *fp, size_t elements_required) {
	return LIBSAKURA_SYMBOL(AlignAny)(elements_in_arena * sizeof(float), fp,
			elements_required * sizeof(float));
}

double const *LIBSAKURA_SYMBOL(AlignDouble)(size_t elements_in_arena,
		double const *dp, size_t elements_required) {
	return LIBSAKURA_SYMBOL(AlignAny)(elements_in_arena * sizeof(double), dp,
			elements_required * sizeof(double));
}
