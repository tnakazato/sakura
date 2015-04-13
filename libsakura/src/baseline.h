/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2014
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
/**
 * baseline.h
 *
 *  Created on: 2013/11/14
 *      Author: wataru
 */

#ifndef BASELINE_H_
#define BASELINE_H_

#include <libsakura/sakura.h>

extern "C" {
struct LIBSAKURA_SYMBOL(BaselineContext) {
	LIBSAKURA_SYMBOL(BaselineType) baseline_type; /**< ベースラインの関数形 */
	size_t num_bases; /**< ベースラインモデルを構成する基底関数の数 */
	size_t num_basis_data; /**< 個々の基底関数を表現するデータ点の数 */
	void *basis_data_storage; /**< 基底関数データを格納する領域を確保した時に返されるアドレス。アラインされておらず、解放時に用いる目的のみのために保持される。 */
	double *basis_data; /**< 基底関数データを指すポインタ。アドレスはアラインされている。要素数は @a num_bases * @a num_basis_data 。i 番目の点の j 番目の基底関数のデータは [ @a num_basis_data * i + j ] 番目の要素に格納される。 */
};
}

#endif /* BASELINE_H_ */
