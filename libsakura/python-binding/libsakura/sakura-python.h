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
/** @file
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
	LIBSAKURA_SYMBOL(PyTypeId_kBool),
	LIBSAKURA_SYMBOL(PyTypeId_kInt8),
	LIBSAKURA_SYMBOL(PyTypeId_kInt32),
	LIBSAKURA_SYMBOL(PyTypeId_kInt64),
	LIBSAKURA_SYMBOL(PyTypeId_kFloat),
	LIBSAKURA_SYMBOL(PyTypeId_kDouble),
	LIBSAKURA_SYMBOL(PyTypeId_kLongDouble),
	LIBSAKURA_SYMBOL(PyTypeId_kEnd)
}LIBSAKURA_SYMBOL(PyTypeId);

/**
 * @~japanese
 * @brief @ref sakura_PyAlignedBuffer を  PyObject(PyCapsule) にカプセル化する
 *
 * @param[in] buffer カプセル化する @ref sakura_PyAlignedBuffer。カプセル化に成功した後は、 @a buffer の解放は PyObjectに委ねられる。
 * @param[out] capsule 成功した場合は、カプセル化した @a PyObject のアドレス(リファレンスカウントは1)。それ以外の場合は不定。
 * @return 終了ステータス
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferEncapsulate)(
		LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, PyObject **capsule) LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief @a PyObject(@a PyCapsule) から、 @a sakura_PyAlignedBuffer を取り出す(デカプセルする)
 *
 * @param[in] capsule デカプセル化する @a PyObject。
 * @param[out] buffer 成功した場合は、デカプセル化された@a sakura_PyAlignedBuffer。それ以外の場合は不定。
 * デカプセル化に成功した場合は、@a capsule が存続する間だけ@a sakura_PyAlignedBuffer は有効である。
 * @a buffer は@a capsule が所有したままなので、勝手に@a buffer を解放してはならない。
 * @return 終了ステータス
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferDecapsulate)(
		PyObject *capsule,
		LIBSAKURA_SYMBOL(PyAlignedBuffer) **buffer) LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief @ref sakura_PyAlignedBuffer を作成する
 *
 * @param[in] type 要素の型。signed/unsignedは、無視。
 * @param[in] original_addr 領域の先頭アドレス。
 * @param[in] aligned_addr 領域中のアラインされたアドレス。
 * このアドレスの利用者間で合意されていれば、必ずしもアラインされていなくてもよい(単に@a original_addr と同じアドレスでもよい)。
 * @param[in] dimensions 次元数(@a elementsの要素数)。現在の実装では、最大次元数は5である。
 * @param[in] elements 各次元の要素数。どういう順序で格納するかは、利用者の自由。
 * @param[in] destructor @a sakura_PyAlignedBuffer が破棄されるときに @a original_addr を解放するために呼び出されるデストラクター。
 * @param[out] buffer @a sakura_PyAlignedBuffer の作成に成功した場合は、そのアドレスが格納される。その他の場合は不定。
 * @return 終了ステータス
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferCreate)(
LIBSAKURA_SYMBOL(PyTypeId) type, void *original_addr, void *aligned_addr,
		size_t dimensions, size_t elements[], void (*destructor)(void *),
		LIBSAKURA_SYMBOL(PyAlignedBuffer) **buffer)
LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief @ref sakura_PyAlignedBuffer を破棄する
 * @ref sakura_PyAlignedBufferCreate で指定した@a original_addrは、@a destructorにより、共に破棄される。
 * @a buffer がNULLの場合は、何もしない。
 *
 * @param[in] buffer @a sakura_PyAlignedBuffer のアドレス。
 */
void LIBSAKURA_SYMBOL(PyAlignedBufferDestroy)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer);

/**
 * @~japanese
 * @brief @ref sakura_PyAlignedBuffer から次元数を取得する
 *
 * @param[in] buffer @a sakura_PyAlignedBuffer のアドレス。
 * @param[out] dimensions 成功した場合は、次元数が格納される。その他の場合は不定。
 * @return 終了ステータス
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferDimensions)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, size_t *dimensions);

/**
 * @~japanese
 * @brief @ref sakura_PyAlignedBuffer から各次元の要素数を取得する
 *
 * @param[in] buffer @a sakura_PyAlignedBuffer のアドレス。
 * @param[in] dimensions	各次元の要素数を先頭から何個取り出したいか(@a elementsに用意されている要素の数)。
 * @ref sakura_PyAlignedBufferCreate を呼び出した際の@a dimensions 以下であること。
 * @param[out] elements 成功した場合は、各次元の要素数が格納される。その他の場合は不定。
 * @return 終了ステータス
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferElements)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, size_t dimensions,
		size_t elements[]);

/**
 * @~japanese
 * @brief @ref sakura_PyAlignedBuffer からアラインされたアドレスを取得する
 *
 * @param[in] buffer @a sakura_PyAlignedBuffer のアドレス。
 * @param[out] aligned_addr 成功した場合は、アラインされたアドレスが格納される。その他の場合は不定。
 * @return 終了ステータス
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferAlignedAddr)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, void **aligned_addr);

/**
 * @~japanese
 * @brief @ref sakura_PyAlignedBuffer から要素の型を取得する
 *
 * @param[in] buffer @a sakura_PyAlignedBuffer のアドレス。
 * @param[out] type 成功した場合は、要素の型が格納される。その他の場合は不定。
 * @return 終了ステータス
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(PyAlignedBufferType)(
LIBSAKURA_SYMBOL(PyAlignedBuffer) *buffer, LIBSAKURA_SYMBOL(PyTypeId) *type);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LIBSAKURA_LIBSAKURA_SAKURA_PYTHON_H_ */
