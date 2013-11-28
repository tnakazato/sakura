#ifndef _LIBSAKURA_INTERPOLATION_UTILS_H_
#define _LIBSAKURA_INTERPOLATION_UTILS_H_

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

#endif /* _LIBSAKURA_INTERPOLATION_UTILS_H_ */
