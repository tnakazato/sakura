/*
 * @SAKURA_LICENSE_HEADER_START@
 * @SAKURA_LICENSE_HEADER_END@
 */
#ifndef _LIBSAKURA_INTERPOLATION_UTILS_TCC_
#define _LIBSAKURA_INTERPOLATION_UTILS_TCC_

#include <libsakura/memory_manager.h>

namespace {

template<class DataType>
struct StorageAndAlignedPointer {
	StorageAndAlignedPointer() :
			storage(nullptr), pointer(nullptr) {
	}
	~StorageAndAlignedPointer() {
		if (storage != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(storage);
			storage = nullptr;
		}
	}
	void *storage;
	DataType *pointer;
};

template<class DataType>
inline void AllocateAndAlign(size_t num_array,
		StorageAndAlignedPointer<DataType> *holder) {
	holder->storage = LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
			num_array * sizeof(DataType), &(holder->pointer));
}

}

#endif /* _LIBSAKURA_INTERPOLATION_UTILS_TCC_ */
