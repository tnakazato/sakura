/*
 * init.cc
 *
 *  Created on: 2013/05/01
 *      Author: kohji
 */

#include <libsakura/sakura.h>
#include "loginit.h"
#include "gtest/gtest.h"

TEST(Global, Init) {
	sakura_Status result = sakura_Initialize(nullptr, nullptr);
	EXPECT_EQ(result, sakura_Status_kOK);

	sakura_CleanUp();

	result = sakura_Initialize(nullptr, nullptr);
	EXPECT_EQ(result, sakura_Status_kOK);

	sakura_CleanUp();
}

TEST(Global, Align) {
	sakura_Status result = sakura_Initialize(nullptr, nullptr);
	EXPECT_EQ(result, sakura_Status_kOK);

	static_assert(sizeof(uint64_t) >= sizeof(void *), "");
	{
		uint64_t offset = 0;
		uint64_t base = 0x10000000ull;
		uint64_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 50);
		auto aaddr = (uint64_t) aptr;
		EXPECT_EQ(base, aaddr);
	}
	{
		uint64_t offset = 1;
		uint64_t base = 0x10000000ull;
		uint64_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 50);
		auto aaddr = (uint64_t) aptr;
		EXPECT_EQ(base, aaddr);
	}
	{
		uint64_t offset = sakura_GetAlignment() - 1;
		uint64_t base = 0x10000000ull;
		uint64_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 50);
		auto aaddr = (uint64_t) aptr;
		EXPECT_EQ(base, aaddr);
	}

	{
		uint64_t offset = 0;
		uint64_t base = 0xf000000010000000ull;
		uint64_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 50);
		auto aaddr = (uint64_t) aptr;
		EXPECT_EQ(base, aaddr);
	}
	{
		uint64_t offset = 1;
		uint64_t base = 0xf000000010000000ull;
		uint64_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 50);
		auto aaddr = (uint64_t) aptr;
		EXPECT_EQ(base, aaddr);
	}
	{
		uint64_t offset = sakura_GetAlignment() - 1;
		uint64_t base = 0xf000000010000000ull;
		uint64_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 50);
		auto aaddr = (uint64_t) aptr;
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

	static_assert(sizeof(uint64_t) >= sizeof(void *), "");
	{
		uint64_t offset = 0;
		uint64_t base = 0x10000000ull;
		uint64_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 100);
		auto aaddr = (uint64_t) aptr;
		EXPECT_EQ(base, aaddr);
	}
	{
		uint64_t offset = 1;
		uint64_t base = 0x10000000ull;
		uint64_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 100);
		EXPECT_EQ(nullptr, aptr);
	}
	{
		uint64_t offset = sakura_GetAlignment() - 1;
		uint64_t base = 0x10000000ull;
		uint64_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 100);
		EXPECT_EQ(nullptr, aptr);
	}

	{
		uint64_t offset = 0;
		uint64_t base = 0xf000000010000000ull;
		uint64_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 100);
		auto aaddr = (uint64_t) aptr;
		EXPECT_EQ(base, aaddr);
	}
	{
		uint64_t offset = 1;
		uint64_t base = 0xf000000010000000ull;
		uint64_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 100);
		EXPECT_EQ(nullptr, aptr);
	}
	{
		uint64_t offset = sakura_GetAlignment() - 1;
		uint64_t base = 0xf000000010000000ull;
		uint64_t addr = base - offset;
		auto ptr = (void *) addr;
		auto aptr = sakura_AlignAny(100, ptr, 100);
		EXPECT_EQ(nullptr, aptr);
	}

	{
		void *aptr = sakura_AlignAny(100, nullptr, 100);
		EXPECT_TRUE(nullptr == aptr);
	}

	sakura_CleanUp();
}
