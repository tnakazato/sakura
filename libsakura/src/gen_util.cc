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

#if LIBSAKURA_HAS_LOG4CXX
#include <log4cxx/propertyconfigurator.h>
#endif

#include "libsakura/sakura.h"
#include "libsakura/localdef.h"
#include "libsakura/logger.h"
#include "libsakura/memory_manager.h"
#include "libsakura/optimized_implementation_factory.h"

namespace {

auto memory_logger = LIBSAKURA_PREFIX::Logger::GetLogger("memory");

void *DefaultAllocator(size_t size) {
	return malloc(std::max(static_cast<size_t>(1), size));
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
LIBSAKURA_SYMBOL(UserDeallocator) deallocator) {
	LIBSAKURA_PREFIX::Memory::allocator_ =
			allocator == nullptr ? DefaultAllocator : allocator;
	LIBSAKURA_PREFIX::Memory::deallocator_ =
			deallocator == nullptr ? DefaultFree : deallocator;

	{
		std::string env_val_name(LIBSAKURA_PREFIX_STRING "_SIMD");
		std::transform(env_val_name.begin(), env_val_name.end(),
				env_val_name.begin(), ::toupper);
		char const *simd_spec = getenv(env_val_name.c_str());
		if (simd_spec == nullptr) {
			simd_spec = "adaptive";
		}
		LIBSAKURA_PREFIX::OptimizedImplementationFactory::InitializeFactory(
				simd_spec);
	}

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

extern "C" size_t LIBSAKURA_SYMBOL (GetAlignment)() {
	return LIBSAKURA_ALIGNMENT;
}

extern "C" bool LIBSAKURA_SYMBOL(IsAligned)(void const *ptr) {
	STATIC_ASSERT(sizeof(uint64_t) >= sizeof(void *));
	uint64_t addr = (uint64_t) ptr;
	return addr % LIBSAKURA_ALIGNMENT == 0;
}

extern "C" void *LIBSAKURA_SYMBOL(AlignAny)(size_t size_of_arena, void *vp,
		size_t size_required) {
	if (vp == nullptr) {
		return nullptr;
	}

	STATIC_ASSERT(sizeof(uint64_t) >= sizeof(void *));
	STATIC_ASSERT(sizeof(uint64_t) >= sizeof(size_t));

	uint64_t addr = (uint64_t) vp;
	uint64_t max_padding = LIBSAKURA_ALIGNMENT - 1u;
	uint64_t new_addr = (addr + max_padding) & ~max_padding;
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
