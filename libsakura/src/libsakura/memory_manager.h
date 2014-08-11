/*
 * @SAKURA_LICENSE_HEADER_START@
 * @SAKURA_LICENSE_HEADER_END@
 */
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
#include <libsakura/localdef.h>

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
	static void * Allocate(size_t size) noexcept;

	/**
	 * @~
	 * @brief
	 * Frees the memory
	 *
	 * MT-safe
	 *
	 * @param[in] ptr nullptr is acceptable.
	 */
	static void Free(void *ptr) noexcept;

	/**
	 * @~
	 * Frees the memory by calling @ref Free(void *ptr)
	 *
	 * This method is intended to use with std::unique_ptr.
	 *
	 * MT-safe
	 *
	 * @param[in] ptr nullptr is acceptable.
	 */
	inline void operator()(void *ptr) const noexcept {
		Free(ptr);
	}

	/**
	 * @~
	 * @brief
	 * Allocate an aligned memory with @a size_in_bytes size and the aligned address is assigned to @a aligned_address
	 *
	 * @param[in] size_in_bytes size in bytes you requires, not a number of objects of type T.
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
		*aligned_address = AssumeAligned(
				reinterpret_cast<T *>(LIBSAKURA_SYMBOL(AlignAny)(size_of_arena,
						ptr, size_in_bytes)));
		assert(*aligned_address != nullptr);
		return ptr;
	}

	Memory() = default;
private:
	friend ::LIBSAKURA_SYMBOL(Status) (::LIBSAKURA_SYMBOL(Initialize))(
			::LIBSAKURA_SYMBOL(UserAllocator) allocator,
			::LIBSAKURA_SYMBOL(UserDeallocator) deallocator);
	static ::LIBSAKURA_SYMBOL(UserAllocator) allocator_;
	static ::LIBSAKURA_SYMBOL(UserDeallocator) deallocator_;
};

} /* namespace LIBSAKURA_PREFIX */

#endif /* LIBSAKURA_LIBSAKURA_MEMORY_MANAGER_H_ */
