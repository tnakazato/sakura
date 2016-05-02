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
	static inline typename Arch::PacketType ByteToInt32SignExtend(
			typename Arch::PriorArch::PacketType bytes) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}

	static inline typename Arch::PriorArch::PacketType Int32ToLSByte(
			typename Arch::PacketType ints) {
		assert(false);
		typename Arch::PriorArch::PacketType result = { 0 };
		return result;
	}

	static inline typename Arch::PacketType ByteToFloat(
			typename Arch::PriorArch::PacketType bytes) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}

	static inline typename Arch::PacketType FloatToDouble(
			typename Arch::PriorArch::PacketType floats) {
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
	static constexpr unsigned kElementsInPacket = Arch::PacketType::kSize
			/ sizeof(Type);
};

/**
 * @~japanese
 * ポータブルな FMAdd/FMSub 機能を提供する
 */
struct LIBSAKURA_SYMBOL(FMA) {
	template<typename Packet, typename T>
	struct GetType {
		typedef void type;
	};
	template<typename Packet>
	struct GetType<Packet, float> {
		typedef typename Packet::RawFloat type;
	};
	template<typename Packet>
	struct GetType<Packet, double> {
		typedef typename Packet::RawDouble type;
	};
	template<typename Packet, typename T>
	static typename GetType<Packet, T>::type MultiplyAdd(
			typename GetType<Packet, T>::type const &a,
			typename GetType<Packet, T>::type const &b,
			typename GetType<Packet, T>::type const &c) {
		return (a * b) + c;	// compiler select fmadd instruction if possible
	}
	template<typename Packet, typename T>
	static typename GetType<Packet, T>::type MultiplySub(
			typename GetType<Packet, T>::type const &a,
			typename GetType<Packet, T>::type const &b,
			typename GetType<Packet, T>::type const &c) {
		return (a * b) - c;	// compiler select fmsub instruction if possible
	}
};

#if defined(__AVX2__)
template<>
typename LIBSAKURA_SYMBOL(FMA)::GetType<LIBSAKURA_SYMBOL(SimdPacketAVX), double>::type LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
		LIBSAKURA_SYMBOL(SimdPacketAVX), double>(
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketAVX), double>::type const &a,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketAVX), double>::type const &b,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketAVX), double>::type const &c) {
	return _mm256_fmadd_pd(a, b, c);
}

template<>
typename LIBSAKURA_SYMBOL(FMA)::GetType<LIBSAKURA_SYMBOL(SimdPacketAVX), double>::type LIBSAKURA_SYMBOL(FMA)::MultiplySub<
		LIBSAKURA_SYMBOL(SimdPacketAVX), double>(
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketAVX), double>::type const &a,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketAVX), double>::type const &b,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketAVX), double>::type const &c) {
	return _mm256_fmsub_pd(a, b, c);
}
template<>
typename LIBSAKURA_SYMBOL(FMA)::GetType<LIBSAKURA_SYMBOL(SimdPacketAVX), float>::type LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
		LIBSAKURA_SYMBOL(SimdPacketAVX), float>(
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketAVX), float>::type const &a,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketAVX), float>::type const &b,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketAVX), float>::type const &c) {
	return _mm256_fmadd_ps(a, b, c);
}

template<>
typename LIBSAKURA_SYMBOL(FMA)::GetType<LIBSAKURA_SYMBOL(SimdPacketAVX), float>::type LIBSAKURA_SYMBOL(FMA)::MultiplySub<
		LIBSAKURA_SYMBOL(SimdPacketAVX), float>(
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketAVX), float>::type const &a,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketAVX), float>::type const &b,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketAVX), float>::type const &c) {
	return _mm256_fmsub_ps(a, b, c);
}

template<>
typename LIBSAKURA_SYMBOL(FMA)::GetType<LIBSAKURA_SYMBOL(SimdPacketSSE), double>::type LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
		LIBSAKURA_SYMBOL(SimdPacketSSE), double>(
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketSSE), double>::type const &a,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketSSE), double>::type const &b,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketSSE), double>::type const &c) {
	return _mm_fmadd_pd(a, b, c);
}

template<>
typename LIBSAKURA_SYMBOL(FMA)::GetType<LIBSAKURA_SYMBOL(SimdPacketSSE), double>::type LIBSAKURA_SYMBOL(FMA)::MultiplySub<
		LIBSAKURA_SYMBOL(SimdPacketSSE), double>(
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketSSE), double>::type const &a,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketSSE), double>::type const &b,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketSSE), double>::type const &c) {
	return _mm_fmsub_pd(a, b, c);
}

template<>
typename LIBSAKURA_SYMBOL(FMA)::GetType<LIBSAKURA_SYMBOL(SimdPacketSSE), float>::type LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
		LIBSAKURA_SYMBOL(SimdPacketSSE), float>(
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketSSE), float>::type const &a,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketSSE), float>::type const &b,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketSSE), float>::type const &c) {
	return _mm_fmadd_ps(a, b, c);
}

template<>
typename LIBSAKURA_SYMBOL(FMA)::GetType<LIBSAKURA_SYMBOL(SimdPacketSSE), float>::type LIBSAKURA_SYMBOL(FMA)::MultiplySub<
		LIBSAKURA_SYMBOL(SimdPacketSSE), float>(
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketSSE), float>::type const &a,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketSSE), float>::type const &b,
		typename LIBSAKURA_SYMBOL(FMA)::GetType<
		LIBSAKURA_SYMBOL(SimdPacketSSE), float>::type const &c) {
	return _mm_fmsub_ps(a, b, c);
}
#endif

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
	static inline typename Arch::PacketType Reciprocal(
			typename Arch::PacketType value) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType Add(typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType Sub(typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType Mul(typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType MulAdd(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		return Add(Mul(lhs, rhs), diff);
	}
	static inline typename Arch::PacketType MulSub(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		return Sub(Mul(lhs, rhs), diff);
	}
	static inline typename Arch::PacketType Div(typename Arch::PacketType lhs,
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
	static inline typename Arch::PacketType And(typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType Or(typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType Xor(typename Arch::PacketType lhs,
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
	static inline typename Arch::PacketType Reciprocal(
			typename Arch::PacketType value) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_float.v[i] = float(1) / value.v_float.v[i];
		}
		return result;
	}
	static inline typename Arch::PacketType Add(
			typename Arch::PacketType const &lhs,
			typename Arch::PacketType const &rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_float.v[i] = lhs.v_float.v[i] + rhs.v_float.v[i];
		}
		return result;
	}
	static inline typename Arch::PacketType Sub(typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
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
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_float.v[i] = lhs.v_float.v[i] * rhs.v_float.v[i];
		}
		return result;
	}
	static inline typename Arch::PacketType MulAdd(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		return Add(Mul(lhs, rhs), diff);
	}
	static inline typename Arch::PacketType MulSub(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		return Sub(Mul(lhs, rhs), diff);
	}
	static inline typename Arch::PacketType Div(typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
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
	static inline typename Arch::PacketType Equal(typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType NotEqual(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType LessThan(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs) {
		assert(false);
		typename Arch::PacketType result = { 0 };
		return result;
	}
	static inline typename Arch::PacketType LessOrEqual(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs) {
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
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int32.v[i] = ~operand.v_int32.v[i];
		}
		return result;
	}
	static inline typename Arch::PacketType Equal(typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int32.v[i] =
					lhs.v_float.v[i] == rhs.v_float.v[i] ? -1LL : 0;
		}
		return result;
	}
	static inline typename Arch::PacketType NotEqual(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int32.v[i] =
					lhs.v_float.v[i] != rhs.v_float.v[i] ? -1LL : 0;
		}
		return result;
	}
	static inline typename Arch::PacketType LessThan(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int32.v[i] =
					lhs.v_float.v[i] < rhs.v_float.v[i] ? -1LL : 0;
		}
		return result;
	}
	static inline typename Arch::PacketType LessOrEqual(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
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
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int64.v[i] = ~operand.v_int64.v[i];
		}
		return result;
	}
	static inline typename Arch::PacketType Equal(typename Arch::PacketType lhs,
			typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int64.v[i] =
					lhs.v_double.v[i] == rhs.v_double.v[i] ? -1LL : 0;
		}
		return result;
	}
	static inline typename Arch::PacketType NotEqual(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int64.v[i] =
					lhs.v_double.v[i] != rhs.v_double.v[i] ? -1LL : 0;
		}
		return result;
	}
	static inline typename Arch::PacketType LessThan(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
				++i) {
			result.v_int64.v[i] =
					lhs.v_double.v[i] < rhs.v_double.v[i] ? -1LL : 0;
		}
		return result;
	}
	static inline typename Arch::PacketType LessOrEqual(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs) {
		typename Arch::PacketType result;
		for (size_t i = 0;
				i
						< LIBSAKURA_SYMBOL(SimdScalarType)<Arch, Type>::kElementsInPacket;
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
	static inline LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType ByteToInt32(
	LIBSAKURA_SYMBOL(SimdArchSSE)::PriorArch::PacketType bytes) {
		LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType result;
		result.v_prior.v[0] = bytes;
		result.raw_int32 = _mm_cvtepi8_epi32(result.raw_int32);
		return result;
	}

	static inline LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType ByteToFloat(
	LIBSAKURA_SYMBOL(SimdArchSSE)::PriorArch::PacketType bytes) {
		LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType result;
		result.raw_float = _mm_cvtpi8_ps(bytes.raw_int64);
		return result;
	}

	static inline LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType FloatToDouble(
	LIBSAKURA_SYMBOL(SimdArchSSE)::PriorArch::PacketType floats) {
		LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType result;
		result.raw_double = _mm_cvtps_pd(
				(__m128 ) _mm_set_epi64(floats.raw_int64, floats.raw_int64));
		return result;
	}
	static inline typename LIBSAKURA_SYMBOL(SimdArchSSE)::PriorArch::PacketType Int32ToLSByte(
			LIBSAKURA_SYMBOL(SimdArchSSE)::PacketType ints) {
		LIBSAKURA_SYMBOL(SimdArchSSE)::PriorArch::PacketType result;
		constexpr int32_t kZero = 0x80808080;
		constexpr int32_t kLSBytes = 0x0c080400;
		result.raw_int64 = _mm_cvtsi64_m64(_mm_extract_epi64(_mm_shuffle_epi8(ints.raw_int32,
				_mm_set_epi32(kZero, kLSBytes, kZero, kZero)), 1));
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdBlend)<LIBSAKURA_SYMBOL(SimdArchSSE), float> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
public:
	static inline Arch::PacketType Blend(Arch::PacketType condition,
			Arch::PacketType false_value, Arch::PacketType true_value) {
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
	static inline Arch::PacketType Blend(Arch::PacketType condition,
			Arch::PacketType false_value, Arch::PacketType true_value) {
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
	static inline Arch::PacketType Blend(Arch::PacketType condition,
			Arch::PacketType false_value, Arch::PacketType true_value) {
		Arch::PacketType result;
		result.raw_float = _mm_blendv_ps(false_value.raw_float,
				true_value.raw_float, condition.raw_float);
		return result;
	}
};

/*----------------------------------------------------*/

template<>
class LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchSSE), int8_t> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
	typedef int8_t Type;
public:
	static inline Arch::PacketType Reciprocal(Arch::PacketType lhs) {
		Arch::PacketType one;
		one.set1(Type(1));
		return Div(one, lhs);
	}
	static inline Arch::PacketType Add(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_add_epi8(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Sub(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_sub_epi8(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Mul(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		int8_t const *bytes_lhs = reinterpret_cast<int8_t const*>(lhs.v_int32.v);
		int8_t const *bytes_rhs = reinterpret_cast<int8_t const*>(rhs.v_int32.v);
		int8_t *bytes_result = reinterpret_cast<int8_t *>(result.v_int32.v);
		for (int i = 0; i < Arch::PacketType::kSize; ++i) {
			bytes_result[i] = bytes_lhs[i] * bytes_rhs[i];
		}
		return result;
	}
	static inline typename Arch::PacketType MulAdd(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		return Add(Mul(lhs, rhs), diff);
	}
	static inline typename Arch::PacketType MulSub(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		return Sub(Mul(lhs, rhs), diff);
	}
	static inline Arch::PacketType Div(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		int8_t const *bytes_lhs = reinterpret_cast<int8_t const*>(lhs.v_int32.v);
		int8_t const *bytes_rhs = reinterpret_cast<int8_t const*>(rhs.v_int32.v);
		int8_t *bytes_result = reinterpret_cast<int8_t *>(result.v_int32.v);
		for (int i = 0; i < Arch::PacketType::kSize; ++i) {
			bytes_result[i] = bytes_lhs[i] / bytes_rhs[i];
		}
		return result;
	}
	static inline Arch::PacketType Not(Arch::PacketType operand) {
		Arch::PacketType all_one;
		all_one.set1(~static_cast<int8_t>(0));
		return Xor(operand, all_one);
	}
	static inline Arch::PacketType And(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_and_si128(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Or(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_or_si128(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Xor(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_xor_si128(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchSSE), int32_t> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
	typedef int32_t Type;
public:
	static inline Arch::PacketType Reciprocal(Arch::PacketType lhs) {
		Arch::PacketType one;
		one.set1(Type(1));
		return Div(one, lhs);
	}
	static inline Arch::PacketType Add(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_add_epi32(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Sub(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_sub_epi32(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Mul(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_mul_epi32(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline typename Arch::PacketType MulAdd(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		return Add(Mul(lhs, rhs), diff);
	}
	static inline typename Arch::PacketType MulSub(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		return Sub(Mul(lhs, rhs), diff);
	}
	static inline Arch::PacketType Div(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		//result.raw_int32 = _mm_div_epi32(lhs.raw_int32, rhs.raw_int32);
		result.v_int32.v[0] = lhs.v_int32.v[0] / rhs.v_int32.v[0];
		result.v_int32.v[1] = lhs.v_int32.v[1] / rhs.v_int32.v[1];
		result.v_int32.v[2] = lhs.v_int32.v[2] / rhs.v_int32.v[2];
		result.v_int32.v[3] = lhs.v_int32.v[3] / rhs.v_int32.v[3];
		return result;
	}
	static inline Arch::PacketType Not(Arch::PacketType operand) {
		Arch::PacketType all_one;
		all_one.set1(~static_cast<int32_t>(0));
		return Xor(operand, all_one);
	}
	static inline Arch::PacketType And(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_and_si128(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Or(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_or_si128(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType Xor(Arch::PacketType lhs,
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
	static inline Arch::PacketType Reciprocal(Arch::PacketType lhs) {
		Arch::PacketType result;
		result.raw_float = _mm_rcp_ps(lhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Add(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_float = _mm_add_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Sub(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm_sub_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Mul(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_float = _mm_mul_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline typename Arch::PacketType MulAdd(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		Arch::PacketType result;
		result.raw_float = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
				typename Arch::PacketType, Type>(lhs.raw_float, rhs.raw_float,
				diff.raw_float);
		return result;
	}
	static inline typename Arch::PacketType MulSub(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		Arch::PacketType result;
		result.raw_float = LIBSAKURA_SYMBOL(FMA)::MultiplySub<
				typename Arch::PacketType, Type>(lhs.raw_float, rhs.raw_float,
				diff.raw_float);
		return result;
	}
	static inline Arch::PacketType Div(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm_div_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchSSE), double> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
	typedef double Type;
public:
	static inline Arch::PacketType Reciprocal(Arch::PacketType lhs) {
		Arch::PacketType one;
		one.set1(Type(1.));
		return Div(one, lhs);
	}
	static inline Arch::PacketType Add(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_double = _mm_add_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType Sub(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double = _mm_sub_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType Mul(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_double = _mm_mul_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
	static inline typename Arch::PacketType MulAdd(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		Arch::PacketType result;
		result.raw_double = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
				typename Arch::PacketType, Type>(lhs.raw_double, rhs.raw_double,
				diff.raw_double);
		return result;
	}
	static inline typename Arch::PacketType MulSub(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		Arch::PacketType result;
		result.raw_double = LIBSAKURA_SYMBOL(FMA)::MultiplySub<
				typename Arch::PacketType, Type>(lhs.raw_double, rhs.raw_double,
				diff.raw_double);
		return result;
	}
	static inline Arch::PacketType Div(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double = _mm_div_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
};

/*----------------------------------------------------*/

template<>
class LIBSAKURA_SYMBOL(SimdCompare)<LIBSAKURA_SYMBOL(SimdArchSSE), int32_t> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
public:
	static inline Arch::PacketType Not(Arch::PacketType operand) {
		return LIBSAKURA_SYMBOL(SimdMath)<Arch, int32_t>::Not(operand);
	}
	static inline Arch::PacketType Equal(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_cmpeq_epi32(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType NotEqual(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		return Not(Equal(lhs, rhs));
	}
	static inline Arch::PacketType LessThan(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_int32 = _mm_cmplt_epi32(lhs.raw_int32, rhs.raw_int32);
		return result;
	}
	static inline Arch::PacketType LessOrEqual(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		return Not(LessThan(rhs, lhs));
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdCompare)<LIBSAKURA_SYMBOL(SimdArchSSE), float> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
public:
	static inline Arch::PacketType Not(Arch::PacketType operand) {
		return LIBSAKURA_SYMBOL(SimdMath)<Arch, int32_t>::Not(operand);
	}
	static inline Arch::PacketType Equal(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm_cmpeq_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType NotEqual(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm_cmpneq_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType LessThan(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm_cmplt_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType LessOrEqual(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm_cmple_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdCompare)<LIBSAKURA_SYMBOL(SimdArchSSE), double> {
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) Arch;
public:
	static inline Arch::PacketType Not(Arch::PacketType operand) {
		return LIBSAKURA_SYMBOL(SimdMath)<Arch, int32_t>::Not(operand);
	}
	static inline Arch::PacketType Equal(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double = _mm_cmpeq_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType NotEqual(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double = _mm_cmpneq_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType LessThan(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double = _mm_cmplt_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType LessOrEqual(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double = _mm_cmple_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
};

#endif /* defined(__SSE4_2__) */

#if defined(__AVX__)

template<>
class LIBSAKURA_SYMBOL(SimdConvert)<LIBSAKURA_SYMBOL(SimdArchAVX)> {
public:
	static inline LIBSAKURA_SYMBOL(SimdArchAVX)::PacketType ByteToInt32SignExtend(
	LIBSAKURA_SYMBOL(SimdArchAVX)::PriorArch::PacketType bytes) {
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

	static inline typename LIBSAKURA_SYMBOL(SimdArchAVX)::PriorArch::PacketType Int32ToLSByte(
			LIBSAKURA_SYMBOL(SimdArchAVX)::PacketType ints) {
		LIBSAKURA_SYMBOL(SimdArchAVX)::PriorArch::PacketType result;
		constexpr int32_t kZero = 0x80808080;
		constexpr int32_t kLSBytes = 0x0c080400;
#if defined(__AVX2__)
		const auto idx = _mm256_set_epi32(kZero, kZero, kLSBytes, kZero,
				kZero, kZero, kZero, kLSBytes);
		auto lsbytes = _mm256_shuffle_epi8(ints.raw_int32, idx);
		result.raw_int32 = _mm_or_si128(_mm256_castsi256_si128(lsbytes), _mm256_extracti128_si256(lsbytes, 1));
#else
		result.raw_int32 = _mm_or_si128(
				_mm_shuffle_epi8(ints.v_prior.v[0].raw_int32,
						_mm_set_epi32(kZero, kZero, kZero, kLSBytes)),
				_mm_shuffle_epi8(ints.v_prior.v[1].raw_int32,
						_mm_set_epi32(kZero, kZero, kLSBytes, kZero)));
#endif
		return result;
	}

	static inline LIBSAKURA_SYMBOL(SimdArchAVX)::PacketType ByteToFloat(
	LIBSAKURA_SYMBOL(SimdArchAVX)::PriorArch::PacketType bytes) {
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

	static inline LIBSAKURA_SYMBOL(SimdArchAVX)::PacketType FloatToDouble(
	LIBSAKURA_SYMBOL(SimdArchAVX)::PriorArch::PacketType floats) {
		LIBSAKURA_SYMBOL(SimdArchAVX)::PacketType result;
		result.raw_double = _mm256_cvtps_pd(floats.raw_float);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdBlend)<LIBSAKURA_SYMBOL(SimdArchAVX), float> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
public:
	static inline Arch::PacketType Blend(Arch::PacketType condition,
			Arch::PacketType false_value, Arch::PacketType true_value) {
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
	static inline Arch::PacketType Blend(Arch::PacketType condition,
			Arch::PacketType false_value, Arch::PacketType true_value) {
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
	static inline Arch::PacketType Blend(Arch::PacketType condition,
			Arch::PacketType false_value, Arch::PacketType true_value) {
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
	static inline Arch::PacketType Reciprocal(Arch::PacketType lhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_rcp_ps(lhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Add(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_int32 = _mm256_add_epi32(lhs.raw_int32,
				rhs.raw_int32);
#else
		result.v_prior.v[0].raw_int32 = _mm_add_epi32(
				lhs.v_prior.v[0].raw_int32, rhs.v_prior.v[0].raw_int32);
		result.v_prior.v[1].raw_int32 = _mm_add_epi32(
				lhs.v_prior.v[1].raw_int32, rhs.v_prior.v[1].raw_int32);
#endif
		return result;
	}
	static inline Arch::PacketType Sub(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_int32 = _mm256_sub_epi32(lhs.raw_int32,
				rhs.raw_int32);
#else
		result.v_prior.v[0].raw_int32 = _mm_sub_epi32(
				lhs.v_prior.v[0].raw_int32, rhs.v_prior.v[0].raw_int32);
		result.v_prior.v[1].raw_int32 = _mm_sub_epi32(
				lhs.v_prior.v[1].raw_int32, rhs.v_prior.v[1].raw_int32);
#endif
		return result;
	}
	static inline Arch::PacketType Mul(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_int32 = _mm256_mul_epi32(lhs.raw_int32,
				rhs.raw_int32);
#else
		result.v_prior.v[0].raw_int32 = _mm_mul_epi32(
				lhs.v_prior.v[0].raw_int32, rhs.v_prior.v[0].raw_int32);
		result.v_prior.v[1].raw_int32 = _mm_mul_epi32(
				lhs.v_prior.v[1].raw_int32, rhs.v_prior.v[1].raw_int32);
#endif
		return result;
	}
	static inline typename Arch::PacketType MulAdd(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		return Add(Mul(lhs, rhs), diff);
	}
	static inline typename Arch::PacketType MulSub(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		return Sub(Mul(lhs, rhs), diff);
	}
	static inline Arch::PacketType Div(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		// result.raw_int32 = _mm256_div_epi32(lhs.raw_int32, rhs.raw_int32);
		result.v_prior.v[0] = LIBSAKURA_SYMBOL(SimdMath)<Arch::PriorArch,
				int32_t>::Div(lhs.v_prior.v[0], rhs.v_prior.v[0]);
		result.v_prior.v[1] = LIBSAKURA_SYMBOL(SimdMath)<Arch::PriorArch,
				int32_t>::Div(lhs.v_prior.v[1], rhs.v_prior.v[1]);
		return result;
	}
	static inline Arch::PacketType Not(Arch::PacketType operand) {
		Arch::PacketType all_one;
		all_one.set1(~static_cast<int32_t>(0));
		return Xor(operand, all_one);
	}
	static inline Arch::PacketType And(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_int32 = _mm256_and_si256(lhs.raw_int32, rhs.raw_int32);
#else
		result.raw_float = _mm256_and_ps(lhs.raw_float, rhs.raw_float);
#endif
		return result;
	}
	static inline Arch::PacketType Or(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_int32 = _mm256_or_si256(lhs.raw_int32, rhs.raw_int32);
#else
		result.raw_float = _mm256_or_ps(lhs.raw_float, rhs.raw_float);
#endif
		return result;
	}
	static inline Arch::PacketType Xor(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_int32 = _mm256_xor_si256(lhs.raw_int32, rhs.raw_int32);
#else
		result.raw_float = _mm256_xor_ps(lhs.raw_float, rhs.raw_float);
#endif
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchAVX), float> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
	typedef float Type;
public:
	static inline Arch::PacketType Reciprocal(Arch::PacketType lhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_rcp_ps(lhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Add(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_add_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Sub(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_sub_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Mul(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_mul_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline typename Arch::PacketType MulAdd(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		Arch::PacketType result;
		result.raw_float = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
				typename Arch::PacketType, Type>(lhs.raw_float, rhs.raw_float,
				diff.raw_float);
		return result;
	}
	static inline typename Arch::PacketType MulSub(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		Arch::PacketType result;
		result.raw_float = LIBSAKURA_SYMBOL(FMA)::MultiplySub<
				typename Arch::PacketType, Type>(lhs.raw_float, rhs.raw_float,
				diff.raw_float);
		return result;
	}
	static inline Arch::PacketType Div(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_div_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Not(Arch::PacketType operand) {
		Arch::PacketType all_one;
		all_one.set1(~static_cast<int32_t>(0));
		return Xor(operand, all_one);
	}
	static inline Arch::PacketType And(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_and_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Or(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_or_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
	static inline Arch::PacketType Xor(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float = _mm256_xor_ps(lhs.raw_float, rhs.raw_float);
		return result;
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdMath)<LIBSAKURA_SYMBOL(SimdArchAVX), double> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
	typedef double Type;
public:
	static inline Arch::PacketType Reciprocal(Arch::PacketType lhs) {
		Arch::PacketType one;
		one.set1(Type(1.));
		return Div(one, lhs);
	}
	static inline Arch::PacketType Add(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_double = _mm256_add_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType Sub(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double = _mm256_sub_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
	static inline Arch::PacketType Mul(Arch::PacketType const &lhs,
			Arch::PacketType const &rhs) {
		Arch::PacketType result;
		result.raw_double = _mm256_mul_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
	static inline typename Arch::PacketType MulAdd(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		Arch::PacketType result;
		result.raw_double = LIBSAKURA_SYMBOL(FMA)::MultiplyAdd<
				typename Arch::PacketType, Type>(lhs.raw_double, rhs.raw_double,
				diff.raw_double);
		return result;
	}
	static inline typename Arch::PacketType MulSub(
			typename Arch::PacketType lhs, typename Arch::PacketType rhs,
			typename Arch::PacketType diff) {
		Arch::PacketType result;
		result.raw_double = LIBSAKURA_SYMBOL(FMA)::MultiplySub<
				typename Arch::PacketType, Type>(lhs.raw_double, rhs.raw_double,
				diff.raw_double);
		return result;
	}
	static inline Arch::PacketType Div(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double = _mm256_div_pd(lhs.raw_double, rhs.raw_double);
		return result;
	}
};

/*----------------------------------------------------*/

template<>
class LIBSAKURA_SYMBOL(SimdCompare)<LIBSAKURA_SYMBOL(SimdArchAVX), int32_t> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
public:
	static inline Arch::PacketType Not(Arch::PacketType operand) {
		return LIBSAKURA_SYMBOL(SimdMath)<Arch, int32_t>::Not(operand);
	}
	static inline Arch::PacketType Equal(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_int32 =
		_mm256_cmpeq_epi32(lhs.raw_int32, rhs.raw_int32);
#else
		result.v_prior.v[0].raw_int32 = _mm_cmpeq_epi32(
				lhs.v_prior.v[0].raw_int32, rhs.v_prior.v[0].raw_int32);
		result.v_prior.v[1].raw_int32 = _mm_cmpeq_epi32(
				lhs.v_prior.v[1].raw_int32, rhs.v_prior.v[1].raw_int32);
#endif
		return result;
	}
	static inline Arch::PacketType NotEqual(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		return Not(Equal(lhs, rhs));
	}
	static inline Arch::PacketType LessThan(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
#if defined(__AVX2__)
		result.raw_int32 =
		_mm256_cmpgt_epi32(rhs.raw_int32, lhs.raw_int32);
#else
		result.v_prior.v[0].raw_int32 = _mm_cmplt_epi32(
				lhs.v_prior.v[0].raw_int32, rhs.v_prior.v[0].raw_int32);
		result.v_prior.v[1].raw_int32 = _mm_cmplt_epi32(
				lhs.v_prior.v[1].raw_int32, rhs.v_prior.v[1].raw_int32);
#endif
		return result;
	}
	static inline Arch::PacketType LessOrEqual(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		return Not(LessThan(rhs, lhs));
	}
};

template<>
class LIBSAKURA_SYMBOL(SimdCompare)<LIBSAKURA_SYMBOL(SimdArchAVX), float> {
	typedef LIBSAKURA_SYMBOL(SimdArchAVX) Arch;
public:
	static inline Arch::PacketType Not(Arch::PacketType operand) {
		return LIBSAKURA_SYMBOL(SimdMath)<Arch, int32_t>::Not(operand);
	}
	static inline Arch::PacketType Equal(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float =
		_mm256_cmp_ps(lhs.raw_float, rhs.raw_float, _CMP_EQ_UQ);
		return result;
	}
	static inline Arch::PacketType NotEqual(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float =
		_mm256_cmp_ps(lhs.raw_float, rhs.raw_float, _CMP_NEQ_UQ);
		return result;
	}
	static inline Arch::PacketType LessThan(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_float =
		_mm256_cmp_ps(lhs.raw_float, rhs.raw_float, _CMP_NGE_UQ);
		return result;
	}
	static inline Arch::PacketType LessOrEqual(Arch::PacketType lhs,
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
	static inline Arch::PacketType Not(Arch::PacketType operand) {
		return LIBSAKURA_SYMBOL(SimdMath)<Arch, int32_t>::Not(operand);
	}
	static inline Arch::PacketType Equal(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double =
		_mm256_cmp_pd(lhs.raw_double, rhs.raw_double, _CMP_EQ_UQ);
		return result;
	}
	static inline Arch::PacketType NotEqual(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double =
		_mm256_cmp_pd(lhs.raw_double, rhs.raw_double, _CMP_NEQ_UQ);
		return result;
	}
	static inline Arch::PacketType LessThan(Arch::PacketType lhs,
			Arch::PacketType rhs) {
		Arch::PacketType result;
		result.raw_double =
		_mm256_cmp_pd(lhs.raw_double, rhs.raw_double, _CMP_NGE_UQ);
		return result;
	}
	static inline Arch::PacketType LessOrEqual(Arch::PacketType lhs,
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
 * @a PacketAction::Prologueを1回、@a PacketAction::Action をn回(elementsが一定数(実装依存)以上あり、SIMDによる速度改善が見込める場合は n &gt; 0、そうでない場合は n = 0)、@a PacketAction::Epilogueを一回呼ぶ。
 *
 * 次に、@a ScalarAction::Prologueを1回、@a ScalarAction::Action をn回(パケット単位で処理できない端数がある場合は n &gt; 0、そうでない場合は n = 0)、@a ScalarAction::Epilogueを一回呼ぶ。
 *
 * @tparam Arch	SIMDアーキテクチャーを識別する型。通常は、@ref sakura_SimdArchNative を指定すれば良い。
 *  @ref sakura_SimdArchAVX or @ref sakura_SimdArchSSE or @ref sakura_SimdArchMMX
 * @tparam ScalarType	@a data の要素(スカラー)の型
 * @tparam PacketAction	SIMD処理の内容。次のメソッドを備えた型であること。
 * <pre>
 * template &lt;typename Arch, typename Context&gt;
 * struct PacketAction {
 * 	static inline void Prologue(Context *context) {
 * 	}
 * 	static inline void Action(size_t idx, typename Arch::PacketType *data, Context *context) {
 * 	}
 * 	static inline void Epilogue(Context *context) {
 * 	}
 * };
 * </pre>
 * @tparam ScalarAction	スカラー処理の内容。次のメソッドを備えた型であること。
 * <pre>
 * template &lt;typename ScalarType, typename Context&gt;
 * struct ScalarAction {
 * 	static inline void Prologue(Context *context) {
 * 	}
 * 	static inline void Action(size_t idx, ScalarType *data, Context *context) {
 * 	}
 * 	static inline void Epilogue(Context *context) {
 * 	}
 * };
 * </pre>
 * @tparam Context	処理のコンテキスト情報の型
 * @param data	アライメントされたデータ
 * @param elements	@a data の要素数
 * @param context		処理のコンテキスト情報へのポインタ
 */
template<typename Arch, typename ScalarType, typename PacketAction,
		typename ScalarAction, typename Context>
void LIBSAKURA_SYMBOL(SimdIterate)(size_t elements, ScalarType data[],
		Context *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data));
	constexpr size_t kUnit =
	LIBSAKURA_SYMBOL(SimdScalarType)<Arch, ScalarType>::kElementsInPacket;
	size_t const packet_count = elements >= kUnit * 1 ? elements / kUnit : 0;
	auto ptr =
			const_cast<typename Arch::PacketType *>(reinterpret_cast<typename Arch::PacketType const *>(data));
	PacketAction::Prologue(context);
	for (size_t i = 0; i < packet_count; ++i) {
		PacketAction::Action(i, &ptr[i], context);
	}
	PacketAction::Epilogue(context);
	ScalarAction::Prologue(context);
	for (size_t i = packet_count * kUnit; i < elements; ++i) {
		ScalarAction::Action(i, &data[i], context);
	}
	ScalarAction::Epilogue(context);
}

/**
 * @~japanese
 * @brief アラインされたデータ(要素数はSIMDパケット長の倍数でなくても良い)を、ベクトル処理とスカラー処理を組み合わせて、繰り返し処理する。
 *
 * 次のように処理することで、@a data 配列を処理する。
 *
 * @a PacketAction::Prologueを1回、@a PacketAction::Action をn回(elementsが一定数(実装依存)以上あり、SIMDによる速度改善が見込める場合は n &gt; 0、そうでない場合は n = 0)、@a PacketAction::Epilogueを一回呼ぶ。
 *
 * 次に、@a ScalarAction::Prologueを1回、@a ScalarAction::Action をn回(パケット単位で処理できない端数がある場合は n &gt; 0、そうでない場合は n = 0)、@a ScalarAction::Epilogueを一回呼ぶ。
 *
 * @tparam Arch1	SIMDアーキテクチャーを識別する型。通常は、@ref sakura_SimdArchNative を指定すれば良い。
 *  @ref sakura_SimdArchAVX or @ref sakura_SimdArchSSE or @ref sakura_SimdArchMMX。
 *  PacketActionによる処理単位はArch1によって決定される。
 * @tparam Arch2	SIMDアーキテクチャーを識別する型。通常は、@ref sakura_SimdArchNative を指定すれば良い。
 * @tparam ScalarType1	@a data1 の要素(スカラー)の型
 * @tparam ScalarType2	@a data2 の要素(スカラー)の型
 * @tparam PacketAction	SIMD処理の内容。次のメソッドを備えた型であること。
 * <pre>
 * template &lt;typename Arch1, typename Arch2, typename Context&gt;
 * struct PacketAction {
 * 	static inline void Prologue(Context *context) {
 * 	}
 * 	static inline void Action(size_t idx, typename Arch1::PacketType *data1, typename Arch2::PacketType *data2, Context *context) {
 * 	}
 * 	static inline void Epilogue(Context *context) {
 * 	}
 * };
 * </pre>
 * @tparam ScalarAction	スカラー処理の内容。次のメソッドを備えた型であること。
 * <pre>
 * template &lt;typename ScalarType1, typename ScalarType2, typename Context&gt;
 * struct ScalarAction {
 * 	static inline void Prologue(Context *context) {
 * 	}
 * 	static inline void Action(size_t idx, ScalarType1 *data1, ScalarType2 *data2, Context *context) {
 * 	}
 * 	static inline void Epilogue(Context *context) {
 * 	}
 * };
 * </pre>
 * @tparam Context	処理のコンテキスト情報の型
 * @param data1	アライメントされたデータ
 * @param data2	アライメントされたデータ
 * @param elements	@a data1, @a data2 の要素数
 * @param context		処理のコンテキスト情報へのポインタ
 */
template<typename Arch1, typename ScalarType1, typename Arch2,
		typename ScalarType2, typename PacketAction, typename ScalarAction,
		typename Context>
void LIBSAKURA_SYMBOL(SimdIterate)(size_t elements, ScalarType1 data1[],
		ScalarType2 data2[], Context *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data2));
	constexpr size_t kUnit =
	LIBSAKURA_SYMBOL(SimdScalarType)<Arch1, ScalarType1>::kElementsInPacket;
	size_t const packet_count = elements >= kUnit * 1 ? elements / kUnit : 0;
	auto ptr1 =
			const_cast<typename Arch1::PacketType *>(reinterpret_cast<typename Arch1::PacketType const *>(data1));
	auto ptr2 =
			const_cast<typename Arch2::PacketType *>(reinterpret_cast<typename Arch2::PacketType const *>(data2));
	PacketAction::Prologue(context);
	for (size_t i = 0; i < packet_count; ++i) {
		PacketAction::Action(i, &ptr1[i], &ptr2[i], context);
	}
	PacketAction::Epilogue(context);
	ScalarAction::Prologue(context);
	for (size_t i = packet_count * kUnit; i < elements; ++i) {
		ScalarAction::Action(i, &data1[i], &data2[i], context);
	}
	ScalarAction::Epilogue(context);
}

/**
 * @~japanese
 * @brief アラインされたデータ(要素数はSIMDパケット長の倍数でなくても良い)を、ベクトル処理とスカラー処理を組み合わせて、繰り返し処理する。
 *
 * 次のように処理することで、@a data 配列を処理する。
 *
 * @a PacketAction::Prologueを1回、@a PacketAction::Action をn回(elementsが一定数(実装依存)以上あり、SIMDによる速度改善が見込める場合は n &gt; 0、そうでない場合は n = 0)、@a PacketAction::Epilogueを一回呼ぶ。
 *
 * 次に、@a ScalarAction::Prologueを1回、@a ScalarAction::Action をn回(パケット単位で処理できない端数がある場合は n &gt; 0、そうでない場合は n = 0)、@a ScalarAction::Epilogueを一回呼ぶ。
 *
 * @tparam Arch1	SIMDアーキテクチャーを識別する型。通常は、@ref sakura_SimdArchNative を指定すれば良い。
 *  @ref sakura_SimdArchAVX or @ref sakura_SimdArchSSE or @ref sakura_SimdArchMMX。
 *  PacketActionによる処理単位はArch1によって決定される。
 * @tparam Arch2	SIMDアーキテクチャーを識別する型。通常は、@ref sakura_SimdArchNative を指定すれば良い。
 * @tparam Arch3	SIMDアーキテクチャーを識別する型。通常は、@ref sakura_SimdArchNative を指定すれば良い。
 * @tparam ScalarType1	@a data1 の要素(スカラー)の型
 * @tparam ScalarType2	@a data2 の要素(スカラー)の型
 * @tparam ScalarType3	@a data3 の要素(スカラー)の型
 * @tparam PacketAction	SIMD処理の内容。次のメソッドを備えた型であること。
 * <pre>
 * template &lt;typename Arch1, typename Arch2, typename Arch3, typename Context&gt;
 * struct PacketAction {
 * 	static inline void Prologue(Context *context) {
 * 	}
 * 	static inline void Action(size_t idx, typename Arch1::PacketType *data1, typename Arch2::PacketType *data2, typename Arch3::PacketType *data3, Context *context) {
 * 	}
 * 	static inline void Epilogue(Context *context) {
 * 	}
 * };
 * </pre>
 * @tparam ScalarAction	スカラー処理の内容。次のメソッドを備えた型であること。
 * <pre>
 * template &lt;typename ScalarType1, typename ScalarType2, typename ScalarType3, typename Context&gt;
 * struct ScalarAction {
 * 	static inline void Prologue(Context *context) {
 * 	}
 * 	static inline void Action(size_t idx, ScalarType1 *data1, ScalarType2 *data2, ScalarType3 *data3, Context *context) {
 * 	}
 * 	static inline void Epilogue(Context *context) {
 * 	}
 * };
 * </pre>
 * @tparam Context	処理のコンテキスト情報の型
 * @param data1	アライメントされたデータ
 * @param data2	アライメントされたデータ
 * @param data3	アライメントされたデータ
 * @param elements	@a data1, @a data2, @a data3 の要素数
 * @param context		処理のコンテキスト情報へのポインタ
 */
template<typename Arch1, typename ScalarType1, typename Arch2,
		typename ScalarType2, typename Arch3, typename ScalarType3,
		typename PacketAction, typename ScalarAction, typename Context>
void LIBSAKURA_SYMBOL(SimdIterate)(size_t elements, ScalarType1 data1[],
		ScalarType2 data2[], ScalarType3 data3[], Context *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data2));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data3));
	constexpr size_t kUnit =
	LIBSAKURA_SYMBOL(SimdScalarType)<Arch1, ScalarType1>::kElementsInPacket;
	size_t const packet_count = elements >= kUnit * 1 ? elements / kUnit : 0;
	auto ptr1 =
			const_cast<typename Arch1::PacketType *>(reinterpret_cast<typename Arch1::PacketType const *>(data1));
	auto ptr2 =
			const_cast<typename Arch2::PacketType *>(reinterpret_cast<typename Arch2::PacketType const *>(data2));
	auto ptr3 =
			const_cast<typename Arch3::PacketType *>(reinterpret_cast<typename Arch3::PacketType const *>(data3));
	PacketAction::Prologue(context);
	for (size_t i = 0; i < packet_count; ++i) {
		PacketAction::Action(i, &ptr1[i], &ptr2[i], &ptr3[i], context);
	}
	PacketAction::Epilogue(context);
	ScalarAction::Prologue(context);
	for (size_t i = packet_count * kUnit; i < elements; ++i) {
		ScalarAction::Action(i, &data1[i], &data2[i], &data3[i], context);
	}
	ScalarAction::Epilogue(context);
}

template<typename Arch1, typename ScalarType1, typename Arch2,
		typename ScalarType2, typename Arch3, typename ScalarType3,
		typename Arch4, typename ScalarType4, typename PacketAction,
		typename ScalarAction, typename Context>
void LIBSAKURA_SYMBOL(SimdIterate)(size_t elements, ScalarType1 data1[],
		ScalarType2 data2[], ScalarType3 data3[], ScalarType4 data4[],
		Context *context) {
	assert(LIBSAKURA_SYMBOL(IsAligned)(data1));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data2));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data3));
	assert(LIBSAKURA_SYMBOL(IsAligned)(data4));
	constexpr size_t kUnit =
	LIBSAKURA_SYMBOL(SimdScalarType)<Arch1, ScalarType1>::kElementsInPacket;
	size_t const packet_count = elements >= kUnit * 1 ? elements / kUnit : 0;
	auto ptr1 =
			const_cast<typename Arch1::PacketType *>(reinterpret_cast<typename Arch1::PacketType const *>(data1));
	auto ptr2 =
			const_cast<typename Arch2::PacketType *>(reinterpret_cast<typename Arch2::PacketType const *>(data2));
	auto ptr3 =
			const_cast<typename Arch3::PacketType *>(reinterpret_cast<typename Arch3::PacketType const *>(data3));
	auto ptr4 =
			const_cast<typename Arch4::PacketType *>(reinterpret_cast<typename Arch4::PacketType const *>(data4));
	PacketAction::Prologue(context);
	for (size_t i = 0; i < packet_count; ++i) {
		PacketAction::Action(i, &ptr1[i], &ptr2[i], &ptr3[i], &ptr4[i],
				context);
	}
	PacketAction::Epilogue(context);
	ScalarAction::Prologue(context);
	for (size_t i = packet_count * kUnit; i < elements; ++i) {
		ScalarAction::Action(i, &data1[i], &data2[i], &data3[i], &data4[i],
				context);
	}
	ScalarAction::Epilogue(context);
}

//} /* namespace */

#endif /* LIBSAKURA_LIBSAKURA_PACKED_OPERATION_H_ */
