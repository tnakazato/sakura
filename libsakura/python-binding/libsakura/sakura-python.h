#ifndef LIBSAKURA_LIBSAKURA_SAKURA_PYTHON_H_
#define LIBSAKURA_LIBSAKURA_SAKURA_PYTHON_H_

#include <libsakura/sakura.h>
#include <libsakura/config-python.h>

#ifdef __cplusplus
extern "C" {
#endif

struct LIBSAKURA_SYMBOL(PyAlignedBuffer);

typedef enum {
	LIBSAKURA_SYMBOL(TypeId_kBool),
	LIBSAKURA_SYMBOL(TypeId_kInt8),
	LIBSAKURA_SYMBOL(TypeId_kInt32),
	LIBSAKURA_SYMBOL(TypeId_kInt64),
	LIBSAKURA_SYMBOL(TypeId_kFloat),
	LIBSAKURA_SYMBOL(TypeId_kDouble),
	LIBSAKURA_SYMBOL(TypeId_kLongDouble),
	LIBSAKURA_SYMBOL(TypeId_kEnd)
}LIBSAKURA_SYMBOL(PyTypeId);

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
LIBSAKURA_SYMBOL(PyTypeId) type, void *original_addr, void *aligned_addr,
		size_t dimensions, size_t elements[], void (*destructor)(void *),
		LIBSAKURA_SYMBOL(PyAlignedBuffer) **buffer);

void LIBSAKURA_SYMBOL(PyAlignedBufferDestroy)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer);

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferDimensions)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, size_t *dimensions);

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferElements)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, size_t dimensions,
		size_t elements[]);

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferAlignedAddr)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, void **aligned_addr);

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferType)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, LIBSAKURA_SYMBOL(PyTypeId) *type);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LIBSAKURA_LIBSAKURA_SAKURA_PYTHON_H_ */
