/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2015
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
struct LIBSAKURA_SYMBOL(BaselineContextFloat) {
	size_t num_bases; /**< ベースラインモデルを構成する基底関数の数 */
	size_t num_basis_data; /**< 個々の基底関数を表現するデータ点の数 */
	size_t num_lsq_bases_max; /** Maximum number of simultaneous equations for least-square fitting: @a baseline_param + 1 @a for Polynomial and Chebyshev, @a baseilne_param + 3 @a for CubicSpline, and @a 2 * baseline_param + 1 @a for Sinusoid. */
	void *basis_data_storage; /**< 基底関数データを格納する領域を確保した時に返されるアドレス。アラインされておらず、解放時に用いる目的のみのために保持される。 */
	double *basis_data; /**< 基底関数データを指すポインタ。アドレスはアラインされている。要素数は @a num_bases * @a num_basis_data 。i 番目の点の j 番目の基底関数のデータは [ @a num_basis_data * i + j ] 番目の要素に格納される。 */
	void *lsq_matrix_storage; /** An address returned when allocating @a lsq_matrix @a . This one itself is not aligned, and is hold just for deallocating. */
	double *lsq_matrix; /** A pointer to aligned 1D array for matrix data used in least-square fitting. Its size is @a num_lsq_bases_max * num_lsq_bases_max @a . */
	void *lsq_vector_storage; /** An address returned when allocating @a lsq_vector @a . This one itself is not aligned, and is hold just for deallocating. */
	double *lsq_vector; /** A pointer to aligned 1D array for vector data used in least-square fitting. Its size is @a num_lsq_bases_max @a . */
	void *clipped_indices_storage; /** An address returned when allocating @a clipped_indices @a . This one itself is not aligned, and is hold just for deallocating. */
	size_t *clipped_indices; /** A pointer to aligned 1D array for clipped indices. Its size is @a num_basis_data @a . */
	void *best_fit_model_storage; /** An address returned when allocating @a best_fit_model @a . This one itself is not aligned, and is hold just for deallocating. */
	float *best_fit_model; /** A pointer to aligned 1D array for best-fit model. Its size is @a num_basis_data @a . */
	void *residual_data_storage; /** An address returned when allocating @a residual_data @a . This one itself is not aligned, and is hold just for deallocating. */
	float *residual_data; /** A pointer to aligned 1D array for baseline-subtracted data. Its size is @a num_basis_data @a . */
	void *use_bases_idx_storage; /** An address returned when allocating @a use_bases_idx @a . This one itself is not aligned, and is hold just for deallocating. */
	size_t *use_bases_idx; /** A pointer to aligned 1D array for indices of bases that are used for fitting. Its size is @a num_lsq_bases_max @a . For baseline type other than Sinusoid, the element values should be equal to their index. */
	void *coeff_full_storage; /** An address returned when allocating @a coeff_full @a . This one itself is not aligned, and is hold just for deallocating. */
	double *coeff_full; /** A pointer to aligned 1D array for storing full coefficients. Its size is @a baseline_param * 4 @a for CubicSpline, and @a num_lsq_bases_max @a for other baseline types. */
	void *cspline_basis_storage; /** An address returned when allocating @a cspline_basis @a . This one itself is not aligned, and is hold just for deallocating. */
	double *cspline_basis; /** A pointer to aligned 1D array for basis functions used in cubic spline fitting. Its size is @a num_basis_data * num_lsq_bases_max @a . */
	void *cspline_lsq_coeff_storage; /** An address returned when allocating @a cspline_lsq_coeff @a . This one itself is not aligned, and is hold just for deallocating. */
	double *cspline_lsq_coeff; /** A pointer to aligned 1D array for least-square fitting results. Its size is @a num_lsq_bases_max @a . */
	LIBSAKURA_SYMBOL(BaselineType) baseline_type; /**< ベースラインの関数形 */
	uint16_t baseline_param; /** order for Polynomial or Chebyshev, npiece for CubicSpline, or maximum wave number for Sinusoid.*/
};
}

#endif /* BASELINE_H_ */
