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

#include <libsakura/sakura.h>

namespace LIBSAKURA_PREFIX {
class Memory {
public:
	/**
	 * @~
	 * @brief
	 * Allocates a memory
	 *
	 * MT-safe
	 *
	 * @param size 0 is acceptable and a valid address for the area with size 0 would be returned.
	 */
	static inline void * Allocate(size_t size) noexcept {
		return allocator_(size);
	}
	/**
	 * @~
	 * @brief
	 * Frees the memory
	 *
	 * MT-safe
	 *
	 * @param ptr nullptr is acceptable.
	 */
	static inline void Free(void *ptr) noexcept {
		deallocator_(ptr);
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
