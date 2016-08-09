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
 * util.c
 *
 *  Created on: 2013/02/27
 *      Author: kohji
 */

#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <memory>
#include <algorithm>
#include <chrono>

#if LIBSAKURA_HAS_LOG4CXX
#include <log4cxx/propertyconfigurator.h>
#endif

#include "libsakura/sakura.h"
#include "libsakura/localdef.h"
#include "libsakura/logger.h"
#include "libsakura/memory_manager.h"

namespace {

auto memory_logger = LIBSAKURA_PREFIX::Logger::GetLogger("memory");

void *DefaultAllocator(size_t size) {
	void *mem_ptr = nullptr;
	auto result = posix_memalign(&mem_ptr, LIBSAKURA_ALIGNMENT,
			std::max(static_cast<size_t>(1), size));
	return result == 0 ? mem_ptr : nullptr;
}

void DefaultFree(void *ptr) {
	free(ptr);
}

} /* namespace */

namespace LIBSAKURA_PREFIX {

#if LIBSAKURA_HAS_LOG4CXX
::log4cxx::LoggerPtr Logger::GetLogger(char const *suffix) {
	std::string logger = LIBSAKURA_PREFIX_STRING;
	if (*suffix != '\0') {
		logger.append(".");
		logger.append(suffix);
	}
	return log4cxx::Logger::getLogger(logger);
}
#endif

LIBSAKURA_SYMBOL(UserAllocator) Memory::allocator_;
LIBSAKURA_SYMBOL(UserDeallocator) Memory::deallocator_;

void * Memory::Allocate(size_t size) noexcept {
	void *ptr = allocator_(size);
	LOG4CXX_DEBUG(memory_logger, "Memory::Allocate: " << ptr);
	return ptr;
}

void Memory::Free(void *ptr) noexcept {
	LOG4CXX_DEBUG(memory_logger, "Memory::Free: " << ptr);
	deallocator_(ptr);
}

} /* namespace LIBSAKURA_PREFIX */


extern "C" LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Initialize)(
LIBSAKURA_SYMBOL(UserAllocator) allocator,
LIBSAKURA_SYMBOL(UserDeallocator) deallocator) noexcept {
	LIBSAKURA_PREFIX::Memory::allocator_ =
			allocator == nullptr ? DefaultAllocator : allocator;
	LIBSAKURA_PREFIX::Memory::deallocator_ =
			deallocator == nullptr ? DefaultFree : deallocator;

	return LIBSAKURA_SYMBOL(Status_kOK);
}

extern "C" void LIBSAKURA_SYMBOL(CleanUp)() noexcept {
}

extern "C" size_t LIBSAKURA_SYMBOL (GetAlignment)() noexcept {
	return LIBSAKURA_ALIGNMENT;
}

extern "C" bool LIBSAKURA_SYMBOL(IsAligned)(void const *ptr) noexcept {
	uintptr_t addr = (uintptr_t) ptr;
	return addr % LIBSAKURA_ALIGNMENT == 0;
}

extern "C" void *LIBSAKURA_SYMBOL(AlignAny)(size_t size_of_arena, void *vp,
		size_t size_required) noexcept {
	if (vp == nullptr) {
		return nullptr;
	}

#if GCC_SUPPORTS_STD_ALIGN // If gcc starts supporting std::align,
	return std::align(LIBSAKURA_ALIGNMENT, size_required, vp, size_of_arena);
#else
	uintptr_t addr = (uintptr_t) vp;
	constexpr uintptr_t kMaxPadding = LIBSAKURA_ALIGNMENT - 1u;
	uintptr_t new_addr = (addr + kMaxPadding) & ~kMaxPadding;
	if (size_of_arena < size_required + size_t(new_addr - addr)) {
		return nullptr;
	}
	return (void *) new_addr;
#endif
}

extern "C" float *LIBSAKURA_SYMBOL(AlignFloat)(size_t elements_in_arena,
		float *fp, size_t elements_required) noexcept {
	return static_cast<float *>(LIBSAKURA_SYMBOL(AlignAny)(
			elements_in_arena * sizeof(float), fp,
			elements_required * sizeof(float)));
}

extern "C" double *LIBSAKURA_SYMBOL(AlignDouble)(size_t elements_in_arena,
		double *dp, size_t elements_required) noexcept {
	return static_cast<double *>(LIBSAKURA_SYMBOL(AlignAny)(
			elements_in_arena * sizeof(double), dp,
			elements_required * sizeof(double)));
}
