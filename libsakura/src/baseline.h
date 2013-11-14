/*
 * baseline.h
 *
 *  Created on: 2013/11/14
 *      Author: wataru
 */

#ifndef BASELINE_H_
#define BASELINE_H_

#include <libsakura/sakura.h>

struct LIBSAKURA_SYMBOL(BaselineContext) {
	LIBSAKURA_SYMBOL(BaselineType) baseline_type; /**< ベースラインの関数形 */
	size_t num_bases; /**< ベースラインモデルを構成する基底関数の数 */
	size_t num_basis_data; /**< 個々の基底関数を表現するデータ点の数 */
	void *basis_data_storage; /**< 基底関数データを格納する領域を確保した時に返されるアドレス。アラインされておらず、解放時に用いる目的のみのために保持される。 */
	double *basis_data; /**< 基底関数データを指すポインタ。アドレスはアラインされている。要素数は @a num_bases * @a num_basis_data 。i 番目の基底関数の j 番目の点のデータは [ @a num_basis_data * i + j ] 番目の要素に格納される。 */
};

#endif /* BASELINE_H_ */
