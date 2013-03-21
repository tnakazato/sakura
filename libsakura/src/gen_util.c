/*
 * util.c
 *
 *  Created on: 2013/02/27
 *      Author: kohji
 */

#include <sys/time.h>

#include "libsakura/sakura.h"

double LIBSAKURA_SYMBOL(GetCurrentTime)() {
	struct timeval tv;
	int result = gettimeofday(&tv, NULL );
	return tv.tv_sec + ((double) tv.tv_usec) / 1000000.;
}
