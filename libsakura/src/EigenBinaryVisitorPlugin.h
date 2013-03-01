/*
 * EigenBinaryVisitorPlugin.h
 *
 *  Created on: 2013/02/22
 *      Author: kohji
 */

#ifndef EIGENBINARYVISITORPLUGIN_H_
#define EIGENBINARYVISITORPLUGIN_H_

template<typename Visitor_, typename Derived_, typename OtherDerived,
		int UnrollCount>
struct binary_visitor_impl {
	enum {
		col = (UnrollCount - 1) / Derived_::RowsAtCompileTime,
		row = (UnrollCount - 1) % Derived::RowsAtCompileTime
	};

	static inline void run(const Derived_ &mat, const OtherDerived &mat2,
			Visitor_& visitor) {
		binary_visitor_impl<Visitor_, Derived_, OtherDerived, UnrollCount - 1>::run(
				mat, mat2, visitor);
		visitor(mat.coeff(row, col), mat2.coeff(row, col), row, col);
	}
};

template<typename Visitor_, typename Derived_, typename OtherDerived>
struct binary_visitor_impl<Visitor_, Derived_, OtherDerived, 1> {
	static inline void run(const Derived_ &mat, const OtherDerived &mat2,
			Visitor_& visitor) {
		visitor.init(mat.coeff(0, 0), mat2.coeff(0, 0), 0, 0);
	}
};

template<typename Visitor_, typename Derived_, typename OtherDerived>
struct binary_visitor_impl<Visitor_, Derived_, OtherDerived, Dynamic> {
	typedef typename Derived_::Index Index;
	static inline void run(const Derived_& mat, const OtherDerived &mat2,
			Visitor_& visitor) {
		visitor.init(mat.coeff(0, 0), mat2.coeff(0, 0), 0, 0);
		for (Index i = 1; i < mat.rows(); ++i)
			visitor(mat.coeff(i, 0), mat2.coeff(i, 0), i, 0);
		for (Index j = 1; j < mat.cols(); ++j)
			for (Index i = 0; i < mat.rows(); ++i)
				visitor(mat.coeff(i, j), mat2.coeff(i, j), i, j);
	}
};

template<typename OtherDerived, typename Visitor>
void visitWith(const DenseBase<OtherDerived> &other, Visitor &visitor) {
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
	return binary_visitor_impl<Visitor, Derived, OtherDerived,
			unroll ? int(SizeAtCompileTime) : Dynamic>::run(derived(), other.derived(),
			visitor);
}

#endif /* EIGENBINARYVISITORPLUGIN_H_ */
