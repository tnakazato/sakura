/*
 * @SAKURA_LICENSE_HEADER_START@
 * @SAKURA_LICENSE_HEADER_END@
 */
/*
 * aligned_memory.h
 *
 *  Created on: 2013/08/09
 *      Author: kohji
 */

#ifndef ALIGNED_MEMORY_H_
#define ALIGNED_MEMORY_H_

#include <cstdlib>
#include <cstdint>
#include <cassert>

template<size_t ALIGN>
class AlignedMemory {
	char buf[ALIGN - 1];
public:
	AlignedMemory() = default;
	static void *operator new(size_t size, size_t realSize) {
		void *p = new char[size + realSize];
		return p;
	}
	static void operator delete(void *p) {
		delete[] reinterpret_cast<char *>(p);
	}

	static AlignedMemory *newAlignedMemory(size_t size) {
		AlignedMemory *p = new (size) AlignedMemory();
		return p;
	}

	inline void operator()(void *ptr) const noexcept {
		delete[] reinterpret_cast<char *>(ptr);
	}

	template<typename T>
	static inline void *AlignedAllocateOrException(size_t size_in_bytes,
			T **aligned_address) throw (std::bad_alloc) {
		assert(aligned_address != nullptr);
		auto ptr = newAlignedMemory(size_in_bytes);
		*aligned_address = ptr->template memory<T>();
		assert(*aligned_address != nullptr);
		return reinterpret_cast<void *>(ptr);
	}

	template<typename T>
	T *memory() {
		uint64_t addr = reinterpret_cast<uint64_t>(&buf[ALIGN - 1]);
		addr &= ~(ALIGN - 1);
		return reinterpret_cast<T *>(addr);
	}
};

typedef AlignedMemory<32> DefaultAlignedMemory; // for AVX

#ifndef SIMD_ALIGN
# if defined(__INTEL_COMPILER) && __INTEL_COMPILER < 1600
#  define SIMD_ALIGN /* nothing */
# else
#  define SIMD_ALIGN alignas(32) // for AVX
# endif
#endif

#endif /* ALIGNED_MEMORY_H_ */
