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
/*
 * EigenBinaryVisitorPlugin.h
 *
 *  Created on: 2013/02/22
 *      Author: kohji
 */

#ifndef LIBSAKURA_LIBSAKURA_EIGENBINARYVISITORPLUGIN_H_
#define LIBSAKURA_LIBSAKURA_EIGENBINARYVISITORPLUGIN_H_

template<typename Visitor_, typename Derived_, typename OtherDerived,
		int kUnrollCount>
struct BinaryVisitorImpl {
	static constexpr decltype(Derived_::RowsAtCompileTime) col = (kUnrollCount
			- 1) / Derived_::RowsAtCompileTime;
	static constexpr decltype(col) row = (kUnrollCount - 1)
			% Derived_::RowsAtCompileTime;

	static inline void Run(const Derived_ &mat, const OtherDerived &mat2,
			Visitor_& visitor) {
		BinaryVisitorImpl<Visitor_, Derived_, OtherDerived, kUnrollCount - 1>::Run(
				mat, mat2, visitor);
		visitor(mat.coeff(row, col), mat2.coeff(row, col), row, col);
	}
};

template<typename Visitor_, typename Derived_, typename OtherDerived>
struct BinaryVisitorImpl<Visitor_, Derived_, OtherDerived, 1> {
	static inline void Run(const Derived_ &mat, const OtherDerived &mat2,
			Visitor_& visitor) {
		visitor.Init(mat.coeff(0, 0), mat2.coeff(0, 0), 0, 0);
	}
};

template<typename Visitor_, typename Derived_, typename OtherDerived>
struct BinaryVisitorImpl<Visitor_, Derived_, OtherDerived, Dynamic> {
	typedef typename Derived_::Index Index;
	static inline void Run(const Derived_& mat, const OtherDerived &mat2,
			Visitor_& visitor) {
		Index i, j;
		for (j = 0; j < mat.cols(); ++j) {
			for (i = 0; i < mat.rows(); ++i) {
				bool init_phase_finished = visitor.Init(mat.coeff(i, j), mat2.coeff(i, j), i, j);
				if (init_phase_finished) {
					goto doRest;
				}
			}
		}
doRest:
		for (++i; i < mat.rows(); ++i) {
			visitor(mat.coeff(i, j), mat2.coeff(i, j), i, j);
		}
		for (++j; j < mat.cols(); ++j) {
			for (Index i = 0; i < mat.rows(); ++i) {
				visitor(mat.coeff(i, j), mat2.coeff(i, j), i, j);
			}
		}
	}
};

template<typename OtherDerived, typename Visitor>
void VisitWith(const DenseBase<OtherDerived> &other, Visitor &visitor) {
	eigen_assert(
			other.rows() == this->rows() && other.cols() == this->cols()
					&& "DenseBase::visitWith(): inconsistent size.");
	enum {
		unroll = SizeAtCompileTime != Dynamic && CoeffReadCost != Dynamic
				&& (SizeAtCompileTime == 1
						|| internal::functor_traits < Visitor > ::Cost
								!= Dynamic)
				&& SizeAtCompileTime * CoeffReadCost
						+ (SizeAtCompileTime - 1) * internal::functor_traits
						< Visitor > ::Cost <= EIGEN_UNROLLING_LIMIT
	};
	return BinaryVisitorImpl<Visitor, Derived, OtherDerived,
			unroll ? int(SizeAtCompileTime) : Dynamic>::Run(derived(), other.derived(),
			visitor);
}

#endif /* LIBSAKURA_LIBSAKURA_EIGENBINARYVISITORPLUGIN_H_ */
