/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2016
 * National Astronomical Observatory of Japan
 * 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
 * 
 * This file is part of Sakura.
 * 
 * Sakura is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * Sakura is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
 * @SAKURA_LICENSE_HEADER_END@
 */
/*
 * c_interface.c
 *
 *  Created on: 2014/12/22
 *      Author: kohji
 */

#include <libsakura/sakura.h>
#include <stdlib.h>
#include <stdio.h>

/* Just a header file conformance check to C language */
int main() {
	unsigned const test_cases = 1;
	unsigned const tests = 1;
	unsigned tests_tried = 0;
	unsigned tests_passed = 0;
	printf("[==========] Running %d tests from %d test case.\n", tests,
			test_cases);
	{
		++tests_tried;
		sakura_Status status = sakura_Initialize(NULL, NULL);
		if (status == sakura_Status_kOK) {
			++tests_passed;
		}
		sakura_CleanUp();
	}
	printf("[==========] %d tests from %d test case ran.\n", tests_tried,
			test_cases);
	printf("[  PASSED  ] %d tests.\n", tests_passed);
	if (tests_tried - tests_passed > 0) {
		printf("[  FAILED  ] %d tests.\n", tests_tried - tests_passed);
		return 1;
	}

	return 0;
}
