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
#ifndef _LIBSAKURA_INTERPOLATION_POLYNOMIAL_TCC_
#define _LIBSAKURA_INTERPOLATION_POLYNOMIAL_TCC_

#include "interpolation_utils.tcc"

namespace {

// TODO: documentation must be added
template<class XDataType, class YDataType>
class PolynomialXInterpolatorImpl {
public:
	PolynomialXInterpolatorImpl() :
			polynomial_order_(0), num_elements_(0), holder_(2) {
	}
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
		polynomial_order_ =
				(polynomial_order + 1u >= num_base) ?
						static_cast<uint8_t>(num_base - 1) : polynomial_order;
		num_elements_ = polynomial_order_ + 1;
		AllocateAndAlign < XDataType > (num_elements_, &holder_[0]);
		AllocateAndAlign < XDataType > (num_elements_, &holder_[1]);
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset) {
		for (size_t k = 1; k < num_location; ++k) {
			int left_edge1 = offset + k - 1 - polynomial_order_ / 2;
			size_t left_edge2 = num_base - num_elements_;
			size_t left_edge = static_cast<size_t>(
					(left_edge1 > 0) ? left_edge1 : 0);
			left_edge = (left_edge > left_edge2) ? left_edge2 : left_edge;
			for (size_t j = 0; j < num_array; ++j) {
				for (size_t i = location[k - 1]; i < location[k]; ++i) {
					PerformNevilleAlgorithm(num_base, base_position, base_data,
							left_edge, j, interpolated_position[i],
							&interpolated_data[j * num_interpolated + i]);
				}
			}
		}
	}
private:
	void PerformNevilleAlgorithm(size_t num_base,
			XDataType const base_position[], YDataType const base_data[],
			size_t left_index, size_t array_index,
			XDataType interpolated_position, YDataType *interpolated_data) {

		// working pointers
		XDataType const *x_ptr = &(base_position[left_index]);
		YDataType const *y_ptr =
				&(base_data[num_base * array_index + left_index]);

		XDataType *c = holder_[0].pointer;
		XDataType *d = holder_[1].pointer;

		for (size_t i = 0; i < num_elements_; ++i) {
			c[i] = static_cast<XDataType>(y_ptr[i]);
			d[i] = c[i];
		}

		// Neville's algorithm
		XDataType work = c[0];
		for (size_t m = 1; m < num_elements_; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = m; i < num_elements_; ++i) {
				XDataType cd = c[i + 1 - m] - d[i - m];
				XDataType dx = x_ptr[i - m] - x_ptr[i];
				assert(dx != 0);
				cd /= dx;
				c[i - m] = (x_ptr[i - m] - interpolated_position) * cd;
				d[i - m] = (x_ptr[i] - interpolated_position) * cd;
			}

			// In each step, c[0] holds Cm1 which is a correction between
			// P12...m and P12...[m+1]. Thus, the following repeated update
			// corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
			work += c[0];
		}

		*interpolated_data = static_cast<YDataType>(work);
	}

	uint8_t polynomial_order_;
	size_t num_elements_;
	std::vector<StorageAndAlignedPointer<XDataType> > holder_;
};

template<class XDataType, class YDataType>
class PolynomialYInterpolatorImpl {
public:
	PolynomialYInterpolatorImpl() :
			polynomial_order_(0), num_elements_(0), holder_(3) {
	}
	void PrepareForInterpolation(uint8_t polynomial_order, size_t num_base,
			size_t num_array, XDataType const base_position[],
			YDataType const base_data[]) {
		polynomial_order_ =
				(polynomial_order + 1u >= num_base) ?
						static_cast<uint8_t>(num_base - 1) : polynomial_order;
		num_elements_ = polynomial_order_ + 1;
		AllocateAndAlign(num_elements_ * num_array, &holder_[0]);
		AllocateAndAlign(num_elements_ * num_array, &holder_[1]);
		AllocateAndAlign(num_array, &holder_[2]);
	}
	void Interpolate1D(size_t num_base, XDataType const base_position[],
			size_t num_array, YDataType const base_data[],
			size_t num_interpolated, XDataType const interpolated_position[],
			YDataType interpolated_data[], size_t num_location,
			size_t const location[], size_t offset) {
		for (size_t k = 1; k < num_location; ++k) {
			int left_edge1 = offset + k - 1 - polynomial_order_ / 2;
			size_t left_edge2 = num_base - num_elements_;
			size_t left_edge = static_cast<size_t>(
					(left_edge1 > 0) ? left_edge1 : 0);
			left_edge = (left_edge > left_edge2) ? left_edge2 : left_edge;
			for (size_t i = location[k - 1]; i < location[k]; ++i) {
				PerformNevilleAlgorithm(num_array, base_position, base_data,
						left_edge, interpolated_position[i],
						&interpolated_data[num_array * i]);
			}
		}
	}
private:
	void PerformNevilleAlgorithm(size_t num_array,
			XDataType const base_position[], YDataType const base_data[],
			size_t left_index, XDataType interpolated_position,
			YDataType interpolated_data[]) {
		// working pointers
		XDataType const *x_ptr = &(base_position[left_index]);
		XDataType *c = holder_[0].pointer;
		XDataType *d = holder_[1].pointer;
		XDataType *work = holder_[2].pointer;

		size_t start = left_index * num_array;
		size_t num_elements = num_elements_ * num_array;
		for (size_t i = 0; i < num_elements; ++i) {
			c[i] = static_cast<XDataType>(base_data[start + i]);
			d[i] = c[i];
		}

		// Neville's algorithm
		for (size_t i = 0; i < num_array; ++i) {
			work[i] = c[i];
		}
		for (size_t m = 1; m < num_elements_; ++m) {
			// Evaluate Cm1, Cm2, Cm3, ... Cm[n-m] and Dm1, Dm2, Dm3, ... Dm[n-m].
			// Those are stored to c[0], c[1], ..., c[n-m-1] and d[0], d[1], ...,
			// d[n-m-1].
			for (size_t i = m; i < num_elements_; ++i) {
				XDataType dx = x_ptr[i - m] - x_ptr[i];
				assert(dx != 0);
				size_t offset = (i - m) * num_array;
				for (size_t j = 0; j < num_array; ++j) {
					XDataType cd = (c[offset + num_array + j] - d[offset + j])
							/ dx;
					c[offset + j] = (x_ptr[i - m] - interpolated_position) * cd;
					d[offset + j] = (x_ptr[i] - interpolated_position) * cd;
				}
			}

			// In each step, c[0] holds Cm1 which is a correction between
			// P12...m and P12...[m+1]. Thus, the following repeated update
			// corresponds to the route P1 -> P12 -> P123 ->...-> P123...n.
			for (size_t i = 0; i < num_array; ++i) {
				work[i] += c[i];
			}
		}

		for (size_t i = 0; i < num_array; ++i) {
			interpolated_data[i] = static_cast<YDataType>(work[i]);
		}
	}

	uint8_t polynomial_order_;
	size_t num_elements_;
	std::vector<StorageAndAlignedPointer<XDataType> > holder_;
};

}

#endif /* _LIBSAKURA_INTERPOLATION_POLYNOMIAL_TCC_ */
