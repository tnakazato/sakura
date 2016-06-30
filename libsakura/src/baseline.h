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
/**
 * baseline.h
 *
 *  Created on: 2013/11/14
 *      Author: wataru
 */

#ifndef BASELINE_H_
#define BASELINE_H_

#include <libsakura/sakura.h>

typedef enum {
	/**
	 * @brief Polynomial
	 */LSQFitTypeInternal_kPolynomial,

	/**
	 * @brief Chebyshev Polynomial
	 */LSQFitTypeInternal_kChebyshev,

	/**
	 * @brief Natural Cubic Spline
	 */LSQFitTypeInternal_kCubicSpline,

	/**
	 * @brief Sinusoids
	 */LSQFitTypeInternal_kSinusoid,

	/**
	 * @brief Number of baseline functions implemented
	 */LSQFitTypeInternal_kNumElements
} LSQFitTypeInternal;

extern "C" {
struct LIBSAKURA_SYMBOL(LSQFitContextFloat) {
	size_t num_bases; /**< Number of basis functions. */
	size_t num_basis_data; /**< Number of discrete data points of each basis function. */
	size_t num_lsq_bases_max; /**< Maximum number of simultaneous equations for least-square fitting: @a lsqfit_param + 1 @a for Polynomial and Chebyshev, @a lsqfit_param + 3 @a for CubicSpline, and @a 2 * lsqfit_param + 1 @a for Sinusoid. */
	void *basis_data_storage; /**< An address returned when allocating @a basis_data @a . This one itself is not aligned, and is hold just for deallocating. */
	double *basis_data; /**< A pointer to aligned 1D array for basis data. Its size is @a num_bases * num_basis_data @a . The basis data of the @a j @a th basis at the @a i @a th point is stored at the [ @a num_basis_data * i + j @a ]-th element. */
	void *lsq_matrix_storage; /**< An address returned when allocating @a lsq_matrix @a . This one itself is not aligned, and is hold just for deallocating. */
	double *lsq_matrix; /**< A pointer to aligned 1D array for matrix data used in least-square fitting. Its size is @a num_lsq_bases_max * num_lsq_bases_max @a . */
	void *lsq_vector_storage; /**< An address returned when allocating @a lsq_vector @a . This one itself is not aligned, and is hold just for deallocating. */
	double *lsq_vector; /**< A pointer to aligned 1D array for vector data used in least-square fitting. Its size is @a num_lsq_bases_max @a . */
	void *clipped_indices_storage; /**< An address returned when allocating @a clipped_indices @a . This one itself is not aligned, and is hold just for deallocating. */
	size_t *clipped_indices; /**< A pointer to aligned 1D array for clipped indices. Its size is @a num_basis_data @a . */
	void *best_fit_model_storage; /**< An address returned when allocating @a best_fit_model @a . This one itself is not aligned, and is hold just for deallocating. */
	float *best_fit_model; /**< A pointer to aligned 1D array for best-fit model. Its size is @a num_basis_data @a . */
	void *residual_data_storage; /**< An address returned when allocating @a residual_data @a . This one itself is not aligned, and is hold just for deallocating. */
	float *residual_data; /**< A pointer to aligned 1D array for baseline-subtracted data. Its size is @a num_basis_data @a . */
	void *use_bases_idx_storage; /**< An address returned when allocating @a use_bases_idx @a . This one itself is not aligned, and is hold just for deallocating. */
	size_t *use_bases_idx; /**< A pointer to aligned 1D array for indices of bases that are used for fitting. Its size is @a num_lsq_bases_max @a . For baseline type other than Sinusoid, the element values should be equal to their index. */
	void *coeff_full_storage; /**< An address returned when allocating @a coeff_full @a . This one itself is not aligned, and is hold just for deallocating. */
	double *coeff_full; /**< A pointer to aligned 1D array for storing full coefficients. Its size is @a lsqfit_param * 4 @a for CubicSpline, and @a num_lsq_bases_max @a for other baseline types. */
	void *cspline_basis_storage; /**< An address returned when allocating @a cspline_basis @a . This one itself is not aligned, and is hold just for deallocating. */
	double *cspline_basis; /**< A pointer to aligned 1D array for basis functions used in cubic spline fitting. Its size is @a num_basis_data * num_lsq_bases_max @a . */
	void *cspline_lsq_coeff_storage; /**< An address returned when allocating @a cspline_lsq_coeff @a . This one itself is not aligned, and is hold just for deallocating. */
	double *cspline_lsq_coeff; /**< A pointer to aligned 1D array for least-square fitting results. Its size is @a num_lsq_bases_max @a . */
	LSQFitTypeInternal lsqfit_type; /**< Functional type of lsq fitting. */
	uint16_t lsqfit_param; /**< order for Polynomial or Chebyshev, npiece for CubicSpline, or maximum wave number for Sinusoid.*/
};
}

#endif /* BASELINE_H_ */
