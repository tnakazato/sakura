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
 * init.cc
 *
 *  Created on: 2013/05/01
 *      Author: kohji
 */

#include <libsakura/sakura.h>
#include <unistd.h>
#include <iostream>
#include <cstdlib>
#include "loginit.h"
#include "gtest/gtest.h"

using namespace std;

namespace {

size_t alloc_count = 0;
size_t free_count = 0;

void *MyAllocate(size_t size) {
	++alloc_count;
	return malloc(size);
}

void MyFree(void *ptr) {
	++free_count;
	free(ptr);
}

}

TEST(Global, Init) {
	sakura_Status result = LIBSAKURA_SYMBOL(Initialize)(nullptr, nullptr);
	EXPECT_EQ(result, sakura_Status_kOK);
	LIBSAKURA_SYMBOL(CleanUp)();

	result = LIBSAKURA_SYMBOL(Initialize)(nullptr, nullptr);
	EXPECT_EQ(result, sakura_Status_kOK);
	LIBSAKURA_SYMBOL(CleanUp)();

	result = LIBSAKURA_SYMBOL(Initialize)(MyAllocate, MyFree);
	EXPECT_EQ(result, sakura_Status_kOK);
	{
		alloc_count = 0;
		free_count = 0;
		//uint16_t const nw = {1};
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context = nullptr;
		auto status = LIBSAKURA_SYMBOL(CreateLSQFitContextPolynomialFloat)(
				LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), 5, 100, &context);
		EXPECT_EQ(status, sakura_Status_kOK);
		EXPECT_NE(nullptr, context);
		status = LIBSAKURA_SYMBOL(DestroyLSQFitContextFloat)(context);
		EXPECT_EQ(status, sakura_Status_kOK);
		EXPECT_LT(0, alloc_count);
		EXPECT_LT(0, free_count);
		EXPECT_EQ(alloc_count, free_count);
	}
	LIBSAKURA_SYMBOL(CleanUp)();
}

TEST(Global, Align) {
	sakura_Status result = sakura_Initialize(nullptr, nullptr);
	EXPECT_EQ(result, sakura_Status_kOK);

	cout << "Alignment is " << sakura_GetAlignment() << endl;

	static_assert(sizeof(uintptr_t) >= sizeof(void *), "");
	{
		uintptr_t offset = 0;
		uintptr_t base = 0x10000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 50);
		auto aaddr = (uintptr_t) aptr;
		EXPECT_EQ(base, aaddr);
		EXPECT_TRUE(sakura_IsAligned(ptr));

		auto fptr = sakura_AlignFloat(100, static_cast<float *>(ptr), 50);
		aaddr = (uintptr_t) fptr;
		EXPECT_EQ(base, aaddr);

		auto dptr = sakura_AlignDouble(100, static_cast<double *>(ptr), 50);
		aaddr = (uintptr_t) dptr;
		EXPECT_EQ(base, aaddr);
	}
	{
		uintptr_t offset = 1;
		uintptr_t base = 0x10000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 50);
		auto aaddr = (uintptr_t) aptr;
		EXPECT_EQ(base, aaddr);
		EXPECT_FALSE(sakura_IsAligned(ptr));
		EXPECT_TRUE(sakura_IsAligned(aptr));

		auto fptr = sakura_AlignFloat(100, static_cast<float *>(ptr), 50);
		aaddr = (uintptr_t) fptr;
		EXPECT_EQ(base, aaddr);

		auto dptr = sakura_AlignDouble(100, static_cast<double *>(ptr), 50);
		aaddr = (uintptr_t) dptr;
		EXPECT_EQ(base, aaddr);
	}
	{
		uintptr_t offset = sakura_GetAlignment() - 1;
		uintptr_t base = 0x10000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 50);
		auto aaddr = (uintptr_t) aptr;
		EXPECT_EQ(base, aaddr);
	}

	{
		uintptr_t offset = 0;
		uintptr_t base = 0xf000000010000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 50);
		auto aaddr = (uintptr_t) aptr;
		EXPECT_EQ(base, aaddr);
	}
	{
		uintptr_t offset = 1;
		uintptr_t base = 0xf000000010000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 50);
		auto aaddr = (uintptr_t) aptr;
		EXPECT_EQ(base, aaddr);
	}
	{
		uintptr_t offset = sakura_GetAlignment() - 1;
		uintptr_t base = 0xf000000010000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 50);
		auto aaddr = (uintptr_t) aptr;
		EXPECT_EQ(base, aaddr);
	}

	{
		void *aptr = sakura_AlignAny(100, nullptr, 50);
		EXPECT_TRUE(nullptr == aptr);
	}

	sakura_CleanUp();
}

TEST(Global, AlignShort) {
	sakura_Status result = sakura_Initialize(nullptr, nullptr);
	EXPECT_EQ(result, sakura_Status_kOK);

	static_assert(sizeof(uintptr_t) >= sizeof(void *), "");
	{
		uintptr_t offset = 0;
		uintptr_t base = 0x10000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(10, ptr, 100);
		EXPECT_EQ(nullptr, aptr);

		auto fptr = sakura_AlignFloat(10, static_cast<float *>(ptr), 100);
		EXPECT_EQ(nullptr, fptr);

		auto dptr = sakura_AlignDouble(10, static_cast<double *>(ptr), 100);
		EXPECT_EQ(nullptr, dptr);
	}
	{
		uintptr_t offset = 0;
		uintptr_t base = 0x10000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 100);
		auto aaddr = (uintptr_t) aptr;
		EXPECT_EQ(base, aaddr);

		auto fptr = sakura_AlignFloat(100, static_cast<float *>(ptr), 100);
		aaddr = (uintptr_t) fptr;
		EXPECT_EQ(base, aaddr);

		auto dptr = sakura_AlignDouble(100, static_cast<double *>(ptr), 100);
		aaddr = (uintptr_t) dptr;
		EXPECT_EQ(base, aaddr);
	}
	{
		uintptr_t offset = 1;
		uintptr_t base = 0x10000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 100);
		EXPECT_EQ(nullptr, aptr);

		auto fptr = sakura_AlignFloat(100, static_cast<float *>(ptr), 100);
		EXPECT_EQ(nullptr, fptr);

		auto dptr = sakura_AlignDouble(100, static_cast<double *>(ptr), 100);
		EXPECT_EQ(nullptr, dptr);
	}
	{
		uintptr_t offset = sakura_GetAlignment() - 1;
		uintptr_t base = 0x10000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 100);
		EXPECT_EQ(nullptr, aptr);
	}
	{
		uintptr_t offset = sakura_GetAlignment() - 1;
		uintptr_t base = 0x10000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(1, ptr, 100);
		EXPECT_EQ(nullptr, aptr);
	}

	{
		uintptr_t offset = 0;
		uintptr_t base = 0xf000000010000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 100);
		auto aaddr = (uintptr_t) aptr;
		EXPECT_EQ(base, aaddr);
	}
	{
		uintptr_t offset = 1;
		uintptr_t base = 0xf000000010000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 100);
		EXPECT_EQ(nullptr, aptr);
	}
	{
		uintptr_t offset = sakura_GetAlignment() - 1;
		uintptr_t base = 0xf000000010000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 100);
		EXPECT_EQ(nullptr, aptr);
	}
	{
		uintptr_t offset = sakura_GetAlignment() - 1;
		uintptr_t base = 0xf000000010000000ULL;
		uintptr_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(1, ptr, 100);
		EXPECT_EQ(nullptr, aptr);
	}

	{
		void *aptr = sakura_AlignAny(100, nullptr, 100);
		EXPECT_TRUE(nullptr == aptr);
	}

	sakura_CleanUp();
}

