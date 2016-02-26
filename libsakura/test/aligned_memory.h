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
#ifndef ALIGNED_MEMORY_H_
#define ALIGNED_MEMORY_H_

#include <cstdlib>
#include <cstdint>
#include <cassert>
#include <new>

template<size_t kAlign>
class AlignedMemory {
	char buf[kAlign - 1];
public:
	AlignedMemory() = default;
	static void *operator new(size_t size, size_t real_size) {
		void *p = new char[size + real_size];
		return p;
	}
	static void operator delete(void *p) {
		delete[] reinterpret_cast<char *>(p);
	}

	static void operator delete(void *p, size_t real_size) {
		delete[] reinterpret_cast<char *>(p);
	}

	static AlignedMemory *NewAlignedMemory(size_t size) {
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
		auto ptr = NewAlignedMemory(size_in_bytes);
		*aligned_address = ptr->template Memory<T>();
		assert(*aligned_address != nullptr);
		return reinterpret_cast<void *>(ptr);
	}

	template<typename T>
	T *Memory() {
		uintptr_t addr = reinterpret_cast<uintptr_t>(&buf[kAlign - 1]);
		addr &= ~(kAlign - 1);
		return reinterpret_cast<T *>(addr);
	}
};

typedef AlignedMemory<32> DefaultAlignedMemory; // for AVX

#ifndef SIMD_ALIGN
# if defined(DISABLE_ALIGNAS) || defined(__INTEL_COMPILER) && __INTEL_COMPILER < 1600
#  define SIMD_ALIGN /* nothing */
# else
#  define SIMD_ALIGN alignas(32) // for AVX
# endif
#endif

#endif /* ALIGNED_MEMORY_H_ */
