/*
 * memory_manager.h
 *
 *  Created on: 2013/09/13
 *      Author: kohji
 */

#ifndef LIBSAKURA_LIBSAKURA_MEMORY_MANAGER_H_
#define LIBSAKURA_LIBSAKURA_MEMORY_MANAGER_H_

#include <libsakura/sakura.h>

namespace LIBSAKURA_PREFIX {
class Memory {
public:
	static inline void * Allocate(size_t size) throw () {
		return allocator_(size);
	}
	static inline void Free(void *ptr) throw () {
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
