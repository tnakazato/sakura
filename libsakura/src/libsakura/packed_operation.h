/**
 * @~japanese
 * @file
 * これはSakuraのAPIの一部ではない。
 * Eigenやintrinsicsを使わずに、ベクトル演算をポータブルに記述するための
 * 型やテンプレートを提供する。
 *
 *  Created on: 2013/04/04
 *      Author: kohji
 */

#ifndef LIBSAKURA_LIBSAKURA_PACKED_OPERATION_H_
#define LIBSAKURA_LIBSAKURA_PACKED_OPERATION_H_

#include <cassert>
#include <libsakura/sakura.h>
#include <libsakura/packed_type.h>

// namespace {

template<typename Arch>
class LIBSAKURA_SYMBOL(SimdConvert) {
public:
	static inline typename Arch::PacketType
	ByteToInt32SignExtend(typename Arch::PriorArch::PacketType bytes) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}

	static inline typename Arch::PacketType
	ByteToFloat(typename Arch::PriorArch::PacketType bytes) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}

	static inline typename Arch::PacketType
	FloatToDouble(typename Arch::PriorArch::PacketType floats) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
};


/**
 * @~japanese
 * ベクトル演算する要素の型の特性を表現する型
 * @tparam Arch SIMDアーキテクチャーを識別する型。通常は、@ref sakura_SimdArchNative を指定すれば良い。
 * @tparam ElementType ベクトル演算する要素の型。float or double or int32_t or int64_t
 */
template<typename Arch, typename ElementType>
struct LIBSAKURA_SYMBOL(SimdScalarType) {
	typedef ElementType Type;
	enum LIBSAKURA_SYMBOL(SimdElementsInPacket) {
		kElementsInPacket = Arch::PacketType::kSize / sizeof(Type)
	};
};

/**
 * @~japanese
 *
 * @tparam Arch SIMDアーキテクチャーを識別する型。通常は、@ref sakura_SimdArchNative を指定すれば良い。
 * @tparam ElementType ベクトル演算する要素の型。float or double or int32_t or int64_t
 */
template<typename Arch, typename ElementType>
class LIBSAKURA_SYMBOL(SimdBlend) {
public:
	/**
	 * @~japanese
	 * @a condition の該当要素の最上位bitが1なら@a true_value を、それ以外なら@ false_value を返す。
	 * @a condition と @a false_value, @a false_value の要素のサイズは揃っていなければならない。
	 * 型は違っても良い。
	 */
	static inline typename Arch::PacketType Blend(
			typename Arch::PacketType const &condition,
			typename Arch::PacketType const &false_value,
			typename Arch::PacketType const &true_value) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, ElementType>::kElementsInPacket;
				++i) {
			switch (LIBSAKURA_SYMBOL(SimdScalarType)<Arch, ElementType>::kElementsInPacket) {
			case Arch::PacketType::kNumInt32:
				result.v_int32.v[i] =
						(condition.v_int32.v[i]
								& (static_cast<uint32_t>(1) << 31)) ?
								true_value.v_int32.v[i] :
								false_value.v_int32.v[i];
				break;
			case Arch::PacketType::kNumInt64:
				result.v_int64.v[i] =
						(condition.v_int64.v[i]
								& (static_cast<uint64_t>(1) << 63)) ?
								true_value.v_int64.v[i] :
								false_value.v_int64.v[i];
				break;
			default:
				assert(false);
				break;
			}
		}
		return result;
	}
};

template<typename Arch, typename ElementType>
class LIBSAKURA_SYMBOL(SimdMath) {
public:
	static inline typename Arch::PacketType Add(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType Sub(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType Mul(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType Div(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType Not(
			typename Arch::PacketType operand) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType And(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType Or(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType Xor(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
};

template<typename Arch>
class LIBSAKURA_SYMBOL(SimdMath)<Arch, float> {
	typedef float Type;
public:
	static inline typename Arch::PacketType Add(
			typename Arch::PacketType const &lhs,
			typename Arch::PacketType const &rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_float.v[i] = lhs.v_float.v[i] + rhs.v_float.v[i];
		}
		return result;
	}
	static inline typename Arch::PacketType Sub(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_float.v[i] = lhs.v_float.v[i] - rhs.v_float.v[i];
		}
		return result;
	}
	static inline typename Arch::PacketType Mul(
			typename Arch::PacketType const &lhs,
			typename Arch::PacketType const &rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_float.v[i] = lhs.v_float.v[i] * rhs.v_float.v[i];
		}
		return result;
	}
	static inline typename Arch::PacketType Div(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_float.v[i] = lhs.v_float.v[i] / rhs.v_float.v[i];
		}
		return result;
	}
};

template<typename Arch, typename ElementType>
class LIBSAKURA_SYMBOL(SimdCompare) {
public:
	static inline typename Arch::PacketType Not(
			typename Arch::PacketType operand) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType Equal(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType NotEqual(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType LessThan(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType LessOrEqual(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
};

template<typename Arch>
class LIBSAKURA_SYMBOL(SimdCompare)<Arch, float> {
	typedef float Type;
public:
	static inline typename Arch::PacketType Not(
			typename Arch::PacketType operand) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int32.v[i] = ~ operand.v_int32.v[i];
		}
		return result;
	}
	static inline typename Arch::PacketType Equal(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int32.v[i] =
					lhs.v_float.v[i] == rhs.v_float.v[i] ? -1LL : 0;
		}
		return result;
	}
	static inline typename Arch::PacketType NotEqual(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int32.v[i] =
					lhs.v_float.v[i] != rhs.v_float.v[i] ? -1LL : 0;
		}
		return result;
	}
	static inline typename Arch::PacketType LessThan(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int32.v[i] =
					lhs.v_float.v[i] < rhs.v_float.v[i] ? -1LL : 0;
		}
		return result;
	}
	static inline typename Arch::PacketType LessOrEqual(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int32.v[i] =
					lhs.v_float.v[i] <= rhs.v_float.v[i] ? -1LL : 0;
		}
		return result;
	}
};

template<typename Arch>
class LIBSAKURA_SYMBOL(SimdCompare)<Arch, double> {
	typedef double Type;
public:
	static inline typename Arch::PacketType Not(
			typename Arch::PacketType operand) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int64.v[i] = ~ operand.v_int64.v[i];
		}
		return result;
	}
	static inline typename Arch::PacketType Equal(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int64.v[i] =
					lhs.v_double.v[i] == rhs.v_double.v[i] ? -1LL : 0;
		}
		return result;
	}
	static inline typename Arch::PacketType NotEqual(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int64.v[i] =
					lhs.v_double.v[i] != rhs.v_double.v[i] ? -1LL : 0;
		}
		return result;
	}
	static inline typename Arch::PacketType LessThan(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int64.v[i] =
					lhs.v_double.v[i] < rhs.v_double.v[i] ? -1LL : 0;
		}
		return result;
	}
	static inline typename Arch::PacketType LessOrEqual(
			typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i < LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int64.v[i] =
					lhs.v_double.v[i] <= rhs.v_double.v[i] ? -1LL : 0;
		}
		return result;
	}
};

#if defined(__SSE4_2__)

template<>
class LIBSAKURA_SYMBOL(SimdConvert)<LIBSAKURA_SYMBOL(SimdArchSSE)> {
public:
	static inline LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType
	ByteToInt32(LIBSAKURA_SYMBOL(SimdArchSSE)::PriorArch::PacketType bytes) {
		LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType result;
		result.v_prior.v[0] = bytes;
		result.raw_int32 = _mm_cvtepi8_epi32(result.raw_int32);
		return result;
	}

	static inline LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType
	ByteToFloat(LIBSAKURA_SYMBOL(SimdArchSSE)::PriorArch::PacketType bytes) {
		LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType result;
		result.raw_float = _mm_cvtpi8_ps(bytes.raw_int64);
		return result;
	}

	static inline LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType
	FloatToDouble(LIBSAKURA_SYMBOL(SimdArchSSE)::PriorArch::PacketType floats) {
		LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType result;
		result.raw_double = _mm_cvtps_pd((__m128)_mm_set_epi64(floats.raw_int64, floats.raw_int64));
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdBlend)<LIBSAKURA_SYMBOL(SimdArchSSE), float> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
public:
	static inline Arch::PacketType Blend(
			Arch::PacketType condition,
			Arch::PacketType false_value,
			Arch::PacketType true_value) {
		Arch::PacketType result;
		result.raw_float = _mm_blendv_ps(false_value.raw_float,
				true_value.raw_float, condition.raw_float);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdBlend)<LIBSAKURA_SYMBOL(SimdArchSSE), double> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
public:
	static inline Arch::PacketType Blend(
			Arch::PacketType condition,
			Arch::PacketType false_value,
			Arch::PacketType true_value) {
		Arch::PacketType result;
		result.raw_double = _mm_blendv_pd(false_value.raw_double,
				true_value.raw_double, condition.raw_double);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdBlend)<LIBSAKURA_SYMBOL(SimdArchSSE), int32_t> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
public:
	static inline Arch::PacketType Blend(
			Arch::PacketType condition,
			Arch::PacketType false_value,
			Arch::PacketType true_value) {
		Arch::PacketType result;
		result.raw_float = _mm_blendv_ps(false_value.raw_float,
				true_value.raw_float, condition.raw_float);
		return result;
	}
};

/*----------------------------------------------------*/

template<>
class LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchSSE), int32_t> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
	typedef int32_t Type;
public:
	static inline Arch::PacketType Add(
			Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_add_epi32(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Sub(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_sub_epi32(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Mul(
			Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_mul_epi32(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Div(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		//result.raw_int32 = _mm_div_epi32(lhs.raw_int32, rhs.raw_int32);
		result.v_int32.v[0] = lhs.v_int32.v[0] / rhs.v_int32.v[0];
		result.v_int32.v[1] = lhs.v_int32.v[1] / rhs.v_int32.v[1];
		result.v_int32.v[2] = lhs.v_int32.v[2] / rhs.v_int32.v[2];
		result.v_int32.v[3] = lhs.v_int32.v[3] / rhs.v_int32.v[3];
		return result;
	}
	static inline Arch::PacketType Not(
			Arch::PacketType operand) {
		Arch::PacketType all_one;
		all_one.set1(~ static_cast<int32_t>(0));
		return Xor(operand, all_one);
	}
	static inline Arch::PacketType And(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_and_si128(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Or(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_or_si128(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Xor(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_xor_si128(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchSSE), float> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
	typedef float Type;
public:
	static inline Arch::PacketType Add(
			Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_float = _mm_add_ps(lhs.raw_float,
				rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Sub(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm_sub_ps(lhs.raw_float,
				rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Mul(
			Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_float = _mm_mul_ps(lhs.raw_float,
				rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Div(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm_div_ps(lhs.raw_float,
				rhs.raw_float);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchSSE), double> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
	typedef double Type;
public:
	static inline Arch::PacketType Add(
			Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_double = _mm_add_pd(lhs.raw_double,
				rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType Sub(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double = _mm_sub_pd(lhs.raw_double,
				rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType Mul(
			Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_double = _mm_mul_pd(lhs.raw_double,
				rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType Div(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double = _mm_div_pd(lhs.raw_double,
				rhs.raw_double);
		return result;
	}
};

/*----------------------------------------------------*/

template<>
class LIBSAKURA_SYMBOL(SimdCompare)<LIBSAKURA_SYMBOL(SimdArchSSE), int32_t> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
public:
	static inline Arch::PacketType Not(
			Arch::PacketType operand) {
		return LIBSAKURA_SYMBOL(SimdMath)<Arch, int32_t>::Not(operand);
	}
	static inline Arch::PacketType Equal(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 =
				_mm_cmpeq_epi32(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType NotEqual(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		return Not(Equal(lhs, rhs));
	}
	static inline Arch::PacketType LessThan(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_cmplt_epi32(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType LessOrEqual(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		return Not(LessThan(rhs, lhs));
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdCompare)<LIBSAKURA_SYMBOL(SimdArchSSE), float> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
public:
	static inline Arch::PacketType Not(
			Arch::PacketType operand) {
		return LIBSAKURA_SYMBOL(SimdMath)<Arch, int32_t>::Not(operand);
	}
	static inline Arch::PacketType Equal(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float =
				_mm_cmpeq_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType NotEqual(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float =
				_mm_cmpneq_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType LessThan(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float =
				_mm_cmplt_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType LessOrEqual(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float =
				_mm_cmple_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdCompare)<LIBSAKURA_SYMBOL(SimdArchSSE), double> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
public:
	static inline Arch::PacketType Not(
			Arch::PacketType operand) {
		return LIBSAKURA_SYMBOL(SimdMath)<Arch, int32_t>::Not(operand);
	}
	static inline Arch::PacketType Equal(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double =
				_mm_cmpeq_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType NotEqual(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double =
				_mm_cmpneq_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType LessThan(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double =
				_mm_cmplt_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType LessOrEqual(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double =
				_mm_cmple_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
};

#endif /* defined(__SSE__) */

#if defined(__AVX__)

template<>
class LIBSAKURA_SYMBOL(SimdConvert)<LIBSAKURA_SYMBOL(SimdArchAVX)> {
public:
	static inline LIBSAKURA_SYMBOL(SimdArchAVX)::PacketType
	ByteToInt32SignExtend(LIBSAKURA_SYMBOL(SimdArchAVX)::PriorArch::PacketType bytes) {
		LIBSAKURA_SYMBOL(SimdArchAVX)::PacketType result;
#if defined(__AVX2__)
		result.raw_int32 =
				_mm256_cvtepi8_epi32(bytes.raw_int32);
#else
		result.v_prior.v[0].raw_int32 = _mm_cvtepi8_epi32(bytes.raw_int32);
		result.v_prior.v[1].raw_int32 = _mm_cvtepi8_epi32(
				_mm_shuffle_epi32(bytes.raw_int32, _MM_SHUFFLE(0,0,0,1)));
#endif
		return result;
	}

	static inline LIBSAKURA_SYMBOL(SimdArchAVX)::PacketType
	ByteToFloat(LIBSAKURA_SYMBOL(SimdArchAVX)::PriorArch::PacketType bytes) {
		LIBSAKURA_SYMBOL(SimdArchAVX)::PacketType result;
#if defined(__AVX2__)
		result.raw_float = _mm256_cvtepi32_ps(
				_mm256_cvtepi8_epi32(bytes.raw_int32));

#else
		result.v_prior.v[0].raw_int32 = _mm_cvtepi8_epi32(bytes.raw_int32);
		result.v_prior.v[1].raw_int32 = _mm_cvtepi8_epi32(
				_mm_shuffle_epi32(bytes.raw_int32, _MM_SHUFFLE(0,0,0,1)));
		result.raw_float = _mm256_cvtepi32_ps(result.raw_int32);
#endif
		return result;
	}

	static inline LIBSAKURA_SYMBOL(SimdArchAVX)::PacketType
	FloatToDouble(LIBSAKURA_SYMBOL(SimdArchAVX)::PriorArch::PacketType floats) {
		LIBSAKURA_SYMBOL(SimdArchAVX)::PacketType result;
		result.raw_double = _mm256_cvtps_pd(floats.raw_float);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdBlend)<LIBSAKURA_SYMBOL(SimdArchAVX), float> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
public:
	static inline Arch::PacketType Blend(
			Arch::PacketType condition,
			Arch::PacketType false_value,
			Arch::PacketType true_value) {
		Arch::PacketType result;
		result.raw_float = _mm256_blendv_ps(false_value.raw_float,
				true_value.raw_float, condition.raw_float);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdBlend)<LIBSAKURA_SYMBOL(SimdArchAVX), double> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
public:
	static inline Arch::PacketType Blend(
			Arch::PacketType condition,
			Arch::PacketType false_value,
			Arch::PacketType true_value) {
		Arch::PacketType result;
		result.raw_double = _mm256_blendv_pd(false_value.raw_double,
				true_value.raw_double, condition.raw_double);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdBlend)<LIBSAKURA_SYMBOL(SimdArchAVX), int32_t> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
public:
	static inline Arch::PacketType Blend(
			Arch::PacketType condition,
			Arch::PacketType false_value,
			Arch::PacketType true_value) {
		Arch::PacketType result;
		result.raw_float = _mm256_blendv_ps(false_value.raw_float,
				true_value.raw_float, condition.raw_float);
		return result;
	}
};

/*----------------------------------------------------*/

template<>
class LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchAVX), int32_t> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
	typedef int32_t Type;
public:
	static inline Arch::PacketType Add(
			Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_int32 = _mm256_add_epi32(lhs.raw_int32,
				rhs.raw_int32);
#else
		result.v_prior.v[0].raw_int32 = _mm_add_epi32(lhs.v_prior.v[0].raw_int32, rhs.v_prior.v[0].raw_int32);
		result.v_prior.v[1].raw_int32 = _mm_add_epi32(lhs.v_prior.v[1].raw_int32, rhs.v_prior.v[1].raw_int32);
#endif
		return result;
	}
	static inline Arch::PacketType Sub(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_int32 = _mm256_sub_epi32(lhs.raw_int32,
				rhs.raw_int32);
#else
		result.v_prior.v[0].raw_int32 = _mm_sub_epi32(lhs.v_prior.v[0].raw_int32, rhs.v_prior.v[0].raw_int32);
		result.v_prior.v[1].raw_int32 = _mm_sub_epi32(lhs.v_prior.v[1].raw_int32, rhs.v_prior.v[1].raw_int32);
#endif
		return result;
	}
	static inline Arch::PacketType Mul(
			Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_int32 = _mm256_mul_epi32(lhs.raw_int32,
				rhs.raw_int32);
#else
		result.v_prior.v[0].raw_int32 = _mm_mul_epi32(lhs.v_prior.v[0].raw_int32, rhs.v_prior.v[0].raw_int32);
		result.v_prior.v[1].raw_int32 = _mm_mul_epi32(lhs.v_prior.v[1].raw_int32, rhs.v_prior.v[1].raw_int32);
#endif
		return result;
	}
	static inline Arch::PacketType Div(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		// result.raw_int32 = _mm256_div_epi32(lhs.raw_int32, rhs.raw_int32);
		result.v_prior.v[0] = LIBSAKURA_SYMBOL(SimdMath)<Arch::PriorArch, int32_t>::Div(lhs.v_prior.v[0], rhs.v_prior.v[0]);
		result.v_prior.v[1] = LIBSAKURA_SYMBOL(SimdMath)<Arch::PriorArch, int32_t>::Div(lhs.v_prior.v[1], rhs.v_prior.v[1]);
		return result;
	}
	static inline Arch::PacketType Not(
			Arch::PacketType operand) {
		Arch::PacketType all_one;
		all_one.set1(~ static_cast<int32_t>(0));
		return Xor(operand, all_one);
	}
	static inline Arch::PacketType And(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_and_ps(lhs.raw_float,
				rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Or(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_or_ps(lhs.raw_float,
				rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Xor(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_xor_ps(lhs.raw_float,
				rhs.raw_float);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchAVX), float> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
	typedef float Type;
public:
	static inline Arch::PacketType Add(
			Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_add_ps(lhs.raw_float,
				rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Sub(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_sub_ps(lhs.raw_float,
				rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Mul(
			Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_mul_ps(lhs.raw_float,
				rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Div(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_div_ps(lhs.raw_float,
				rhs.raw_float);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchAVX), double> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
	typedef double Type;
public:
	static inline Arch::PacketType Add(
			Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_double = _mm256_add_pd(lhs.raw_double,
				rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType Sub(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double = _mm256_sub_pd(lhs.raw_double,
				rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType Mul(
			Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_double = _mm256_mul_pd(lhs.raw_double,
				rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType Div(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double = _mm256_div_pd(lhs.raw_double,
				rhs.raw_double);
		return result;
	}
};

/*----------------------------------------------------*/

template<>
class LIBSAKURA_SYMBOL(SimdCompare)<LIBSAKURA_SYMBOL(SimdArchAVX), int32_t> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
public:
	static inline Arch::PacketType Not(
			Arch::PacketType operand) {
		return LIBSAKURA_SYMBOL(SimdMath)<Arch, int32_t>::Not(operand);
	}
	static inline Arch::PacketType Equal(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_float =
				 _mm256_cmpeq_epi32(lhs.raw_int32, rhs.raw_int32);
#else
		result.v_prior.v[0].raw_int32 = _mm_cmpeq_epi32(lhs.v_prior.v[0].raw_int32, rhs.v_prior.v[0].raw_int32);
		result.v_prior.v[1].raw_int32 = _mm_cmpeq_epi32(lhs.v_prior.v[1].raw_int32, rhs.v_prior.v[1].raw_int32);
#endif
		return result;
	}
	static inline Arch::PacketType NotEqual(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		return Not(Equal(lhs, rhs));
	}
	static inline Arch::PacketType LessThan(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_float =
				 _mm256_cmpgt_epi32(rhs.raw_int32, lhs.raw_int32);
#else
		result.v_prior.v[0].raw_int32 = _mm_cmplt_epi32(lhs.v_prior.v[0].raw_int32, rhs.v_prior.v[0].raw_int32);
		result.v_prior.v[1].raw_int32 = _mm_cmplt_epi32(lhs.v_prior.v[1].raw_int32, rhs.v_prior.v[1].raw_int32);
#endif
		return result;
	}
	static inline Arch::PacketType LessOrEqual(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		return Not(LessThan(rhs, lhs));
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdCompare)<LIBSAKURA_SYMBOL(SimdArchAVX), float> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
public:
	static inline Arch::PacketType Not(
			Arch::PacketType operand) {
		return LIBSAKURA_SYMBOL(SimdMath)<Arch, int32_t>::Not(operand);
	}
	static inline Arch::PacketType Equal(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float =
				_mm256_cmp_ps(lhs.raw_float, rhs.raw_float, _CMP_EQ_UQ);
		return result;
	}
	static inline Arch::PacketType NotEqual(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float =
				_mm256_cmp_ps(lhs.raw_float, rhs.raw_float, _CMP_NEQ_UQ);
		return result;
	}
	static inline Arch::PacketType LessThan(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float =
				_mm256_cmp_ps(lhs.raw_float, rhs.raw_float, _CMP_NGE_UQ);
		return result;
	}
	static inline Arch::PacketType LessOrEqual(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float =
				_mm256_cmp_ps(lhs.raw_float, rhs.raw_float, _CMP_NGT_UQ);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdCompare)<LIBSAKURA_SYMBOL(SimdArchAVX), double> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
public:
	static inline Arch::PacketType Not(
			Arch::PacketType operand) {
		return LIBSAKURA_SYMBOL(SimdMath)<Arch, int32_t>::Not(operand);
	}
	static inline Arch::PacketType Equal(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double =
				_mm256_cmp_pd(lhs.raw_double, rhs.raw_double, _CMP_EQ_UQ);
		return result;
	}
	static inline Arch::PacketType NotEqual(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double =
				_mm256_cmp_pd(lhs.raw_double, rhs.raw_double, _CMP_NEQ_UQ);
		return result;
	}
	static inline Arch::PacketType LessThan(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double =
				_mm256_cmp_pd(lhs.raw_double, rhs.raw_double, _CMP_NGE_UQ);
		return result;
	}
	static inline Arch::PacketType LessOrEqual(
			Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double =
				_mm256_cmp_pd(lhs.raw_double, rhs.raw_double, _CMP_NGT_UQ);
		return result;
	}
};

#endif /* defined(__AVX__) */

/**
 * @~japanese
 * @brief アラインされたデータ(要素数はSIMDパケット長の倍数でなくても良い)を、ベクトル処理とスカラー処理を組み合わせて、繰り返し処理する。
 *
 * 次のように処理することで、@a data 配列を処理する。
 *
 * @a PacketAction::prologueを1回、@a PacketAction::action をn回(elementsが一定数(実装依存)以上あり、SIMDによる速度改善が見込める場合は n &gt; 0、そうでない場合は n = 0)、@a PacketAction::epilogueを一回呼ぶ。
 *
 * 次に、@a ScalarAction::prologueを1回、@a ScalarAction::action をn回(パケット単位で処理できない端数がある場合は n &gt; 0、そうでない場合は n = 0)、@a ScalarAction::epilogueを一回呼ぶ。
 *
 * @tparam Arch	SIMDアーキテクチャーを識別する型。通常は、@ref sakura_SimdArchNative を指定すれば良い。
 *  @ref sakura_SimdArchAVX or @ref sakura_SimdArchSSE or @ref sakura_SimdArchMMX
 * @tparam ScalarType	@a data の要素(スカラー)の型
 * @tparam PacketAction	SIMD処理の内容。次のメソッドを備えた型であること。
 * <pre>
 * template &lt;typename Arch, typename Context&gt;
 * struct PacketAction {
 * 	static inline void prologue(Context *context) {
 * 	}
 * 	static inline void action(size_t idx, typename Arch::PacketType *data, Context *context) {
 * 	}
 * 	static inline void epilogue(Context *context) {
 * 	}
 * };
 * </pre>
 * @tparam ScalarAction	スカラー処理の内容。次のメソッドを備えた型であること。
 * <pre>
 * template &lt;typename Arch, typename Context&gt;
 * struct ScalarAction {
 * 	static inline void prologue(Context *context) {
 * 	}
 * 	static inline void action(size_t idx, ScalarType *data, Context *context) {
 * 	}
 * 	static inline void epilogue(Context *context) {
 * 	}
 * };
 * </pre>
 * @tparam Context	処理のコンテキスト情報の型
 * @param data	アライメントされたデータ
 * @param elements	@a data の要素数
 * @param context		処理のコンテキスト情報へのポインタ
 */
template<typename Arch, typename ScalarType, typename PacketAction, typename ScalarAction,
		typename Context>
void LIBSAKURA_SYMBOL(SimdIterate)(size_t elements, ScalarType data[],
		Context *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	size_t const kUnit =
			LIBSAKURA_SYMBOL(SimdScalarType)<Arch, ScalarType>::kElementsInPacket;
	size_t const packet_count = elements >= kUnit * 4 ? elements / kUnit : 0;
	auto ptr = reinterpret_cast<typename Arch::PacketType *>(data);
	PacketAction::prologue(context);
	for (size_t i = 0; i < packet_count; ++i) {
		PacketAction::action(i, &ptr[i], context);
	}
	PacketAction::epilogue(context);
	ScalarAction::prologue(context);
	for (size_t i = packet_count * kUnit; i < elements; ++i) {
		ScalarAction::action(i, &data[i], context);
	}
	ScalarAction::epilogue(context);
}

//} /* namespace */

#endif /* LIBSAKURA_LIBSAKURA_PACKED_OPERATION_H_ */
