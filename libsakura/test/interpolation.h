/*
 * interpolation.h
 *
 *  Created on: Sep 18, 2013
 *      Author: nakazato
 */

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include <stdarg.h>

void InitializeDoubleArray(size_t num_array, double array[], ...) {
	va_list arguments_list;
	va_start(arguments_list, array);
	for (size_t i = 0; i < num_array; ++i) {
		array[i] = va_arg(arguments_list, double);
	}
}

void InitializeFloatArray(size_t num_array, float array[], ...) {
	va_list arguments_list;
	va_start(arguments_list, array);
	for (size_t i = 0; i < num_array; ++i) {
		array[i] = static_cast<float>(va_arg(arguments_list, double));
	}
}

#endif /* INTERPOLATION_H_ */
