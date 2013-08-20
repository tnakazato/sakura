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

template<size_t ALIGN>
class AlignedMemory {
	char buf[ALIGN - 1];
	AlignedMemory() {
	}
public:
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

	template<typename T>
	T *memory() {
		uint64_t addr = reinterpret_cast<uint64_t>(&buf[ALIGN - 1]);
		addr &= ~(ALIGN - 1);
		return reinterpret_cast<T *>(addr);
	}
};

typedef AlignedMemory<32> DefaultAlignedMemory; // for AVX

#define SIMD_ALIGN alignas(32) // for AVX
#endif /* ALIGNED_MEMORY_H_ */
