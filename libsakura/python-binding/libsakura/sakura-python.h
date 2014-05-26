/**
 * @file
 *  Sakura python binding header file.
 *
 * sakura-python.h
 *
 *  Created on: 2014/05/26
 *      Author: kohji
 */

#ifndef LIBSAKURA_LIBSAKURA_SAKURA_PYTHON_H_
#define LIBSAKURA_LIBSAKURA_SAKURA_PYTHON_H_

#include <Python.h>
#include <libsakura/sakura.h>
#include <libsakura/config-python.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @~japanese
 * @brief アラインされたバッファー
 *
 * 内部に次の属性を持つ。
 *
 * @li @c 型 @ref sakura_PyTypeId
 * @li @c サイズ(次元(最大5次元)、各次元の要素数)
 * @li @c 領域の先頭ポインター
 * @li @c アラインされたアドレス
 * @li @c 領域の解放関数
 */
struct LIBSAKURA_SYMBOL(PyAlignedBuffer);

/**
 * @~japanese
 * @brief @ref sakura_PyAlignedBuffer に格納される要素の型
 *
 */
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

/**
 * @~japanese
 * @brief @ref PyAlignedBuffer を  PyObject(PyCapsule) にカプセル化する
 *
 * @param[in] buffer カプセル化する @ref sakura_PyAlignedBuffer。カプセル化に成功した後は、 @a buffer の解放は PyObjectに委ねられる。
 * @param[out] capsule 成功した場合は、カプセル化した @a PyObject のアドレス。それ以外の場合は不定。
 * @return 終了ステータス
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferEncapsulate)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, PyObject **capsule);

/**
 * @~japanese
 * @brief PyObject(PyCapsule) を @ref PyAlignedBuffer にデカプセル化する
 *
 * @param[in] capsule デカプセル化する @a PyObject。
 * @param[in] buffer 成功した場合は、デカプセル化された@ref sakura_PyAlignedBuffer。それ以外の場合は不定。
 * デカプセル化に成功した場合は、capsuleが存続する間だけ@ref PyAlignedBufferは有効である。
 * @return 終了ステータス
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferDecapsulate)(
		PyObject *capsule,
		LIBSAKURA_SYMBOL(PyAlignedBuffer) **buffer);

/**
 * @~japanese
 * @brief @ref PyAlignedBuffer を作成する
 *
 * @param[in] elements 各次元の要素数。どういう順序で格納するかは、利用者の自由。
 * @return 終了ステータス
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
LIBSAKURA_SYMBOL(PyTypeId) type, void *original_addr, void *aligned_addr,
		size_t dimensions, size_t elements[], void (*destructor)(void *),
		LIBSAKURA_SYMBOL(PyAlignedBuffer) **buffer);

void LIBSAKURA_SYMBOL(PyAlignedBufferDestroy)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer);

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferDimensions)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, size_t *dimensions);

/**
 * @~japanese
 * @brief @ref PyAlignedBuffer から各次元の要素数を取得する
 *
 * @param[in] dimensions	各次元の要素数を先頭から何個取り出したいか(@a elementsに用意されている要素の数)。
 * @a sakura_PyAlignedBufferCreate を呼び出した際の@a dimensions 以下であること。
 * @return 終了ステータス
 */
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
