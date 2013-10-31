/*
 * memory_manager.h
 *
 *  Created on: 2013/09/13
 *      Author: kohji
 */

/**
 * @file
 * Contains a memory manager class
 *
 * This file is an internal header file for libsakura.
 * This is not a part of libsakura API.
 */

#ifndef LIBSAKURA_LIBSAKURA_MEMORY_MANAGER_H_
#define LIBSAKURA_LIBSAKURA_MEMORY_MANAGER_H_

#include <new>
#include <libsakura/sakura.h>

namespace LIBSAKURA_PREFIX {
class Memory {
public:
	typedef void (*TypeOfFree)(void*);

	/**
	 * @~
	 * @brief
	 * Allocates a memory
	 *
	 * MT-safe
	 *
	 * @param[in] size 0 is acceptable and a valid address for the area with size 0 would be returned.
	 */
	static inline void * Allocate(size_t size) noexcept {
		void *ptr = allocator_(size);
		// printf("Memory::Allocate: %p\n", ptr);
		return ptr;
	}
	/**
	 * @~
	 * @brief
	 * Frees the memory
	 *
	 * MT-safe
	 *
	 * @param[in] ptr nullptr is acceptable.
	 */
	static inline void Free(void *ptr) noexcept {
		// printf("Memory::Free: %p\n", ptr);
		deallocator_(ptr);
	}

	/**
	 * @~
	 * @brief
	 * Allocate an aligned memory with @a size_in_bytes size and the aligned address is assigned to @a aligned_address
	 *
	 * @param[in] size_in_bytes size in bytes you requires, not a number of T.
	 * @param[out] aligned_address @a *aligned_address points aligned memory on successful call, otherwise @a *aligned_address is set to nullptr.
	 * @return the address of the top of allocated area which may be unaligned. This address is required when freeing the area by @ref Free(void *ptr).
	 *
	 * MT-safe
	 */
	template<typename T>
	static inline void *AlignedAllocateOrException(size_t size_in_bytes,
			T **aligned_address) throw (std::bad_alloc) {
		assert(aligned_address != nullptr);
		size_t alignment = LIBSAKURA_SYMBOL(GetAlignment)();

		size_t const size_of_arena = size_in_bytes + alignment - 1;
		void *ptr = Allocate(size_of_arena);
		if (ptr == nullptr) {
			*aligned_address = nullptr;
			throw std::bad_alloc();
		}
		*aligned_address = reinterpret_cast<T *>(LIBSAKURA_SYMBOL(AlignAny)(
				size_of_arena, ptr, size_in_bytes));
		assert(*aligned_address != nullptr);
		return ptr;
	}

private:
	friend ::LIBSAKURA_SYMBOL(Status) (::LIBSAKURA_SYMBOL(Initialize))(
			::LIBSAKURA_SYMBOL(UserAllocator) allocator,
			::LIBSAKURA_SYMBOL(UserDeallocator) deallocator);
	Memory() = delete;
	static ::LIBSAKURA_SYMBOL(UserAllocator) allocator_;
	static ::LIBSAKURA_SYMBOL(UserDeallocator) deallocator_;
};

} /* namespace LIBSAKURA_PREFIX */

#endif /* LIBSAKURA_LIBSAKURA_MEMORY_MANAGER_H_ */
