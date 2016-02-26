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

#ifndef LIBSAKURA_LIBSAKURA_PACKED_TYPE_H_
#define LIBSAKURA_LIBSAKURA_PACKED_TYPE_H_

#include <stdint.h>
#include <libsakura/config.h>

#if defined(__MMX__) || defined(__x86_64__)
/* for 64 bit CPU */
#include <mmintrin.h>

/**
 * @~japanese
 * @brief 	一度のベクトル演算で処理するデータの固まり(Packet)を格納する型(MMX用)である。
 */
union LIBSAKURA_SYMBOL(SimdPacketMMX) {
	enum {
		/**
		 * @~japanese
		 *  パケットのサイズ
		 */
		kSize = sizeof(__m64 ),
		/**
		 * @~japanese
		 *  パケットに格納できるfloatの要素数
		 */
		kNumFloat = kSize / sizeof(float),
		/**
		 * @~japanese
		 *  パケットに格納できるdoubleの要素数
		 */
		kNumDouble = kSize / sizeof(double),
		/**
		 * @~japanese
		 * パケットに格納できるint32_tの要素数
		 */
		kNumInt32 = kSize / sizeof(int32_t),
		/**
		 * @~japanese
		 * パケットに格納できるint64_tの要素数
		 */
		kNumInt64 = kSize / sizeof(int64_t),
	};

	typedef struct {
		float dummy[kNumFloat];
	} RawFloat;
	typedef struct {
		double dummy[kNumDouble];
	} RawDouble;
	typedef struct {
		int32_t dummy[kNumInt32];
	} RawInt32;
	typedef __m64 RawInt64;
	RawInt64 raw_int64;
	RawFloat raw_float;
	RawDouble raw_double;
	RawInt32 raw_int32;

	/**
	 * @~japanese
	 * @brief @a value を各要素に設定して初期化する。
	 */
	inline void set1(double value) {
		v_double.v[0] = value;
	}
	/**
	 * @~japanese
	 * @brief @a value を各要素に設定して初期化する。
	 */
	inline void set1(float value) {
		v_float.v[0] = value;
		v_float.v[1] = value;
	}
	/**
	 * @~japanese
	 * @brief @a value を各要素に設定して初期化する。
	 */
	inline void set1(int8_t value) {
		raw_int64 = _mm_set_pi8(value, value, value, value, value, value, value,
				value);
	}
	/**
	 * @~japanese
	 * @brief @a value を各要素に設定して初期化する。
	 */
	inline void set1(int16_t value) {
		raw_int64 = _mm_set_pi16(value, value, value, value);
	}
	/**
	 * @~japanese
	 * @brief @a value を各要素に設定して初期化する。
	 */
	inline void set1(int32_t value) {
		raw_int64 = _mm_set_pi32(value, value);
	}
	/**
	 * @~japanese
	 * @brief @a value を各要素に設定して初期化する。
	 */
	inline void set1(int64_t value) {
		raw_int64 = _mm_cvtsi64_m64(value);
	}

	/**
	 * @~japanese
	 * パケットに格納されているfloat型要素にアクセスするための構造体
	 */
	typedef struct {
		float v[kNumFloat];
	} VFloat;
	/**
	 * @~japanese
	 * パケットに格納されているdouble型要素にアクセスするための構造体
	 */
	typedef struct {
		double v[kNumDouble];
	} VDouble;
	/**
	 * @~japanese
	 * パケットに格納されているint32_t型要素にアクセスするための構造体
	 */
	typedef struct {
		int32_t v[kNumInt32];
	} VInt32;
	/**
	 * @~japanese
	 * パケットに格納されているint64_t型要素にアクセスするための構造体
	 */
	typedef struct {
		int64_t v[kNumInt64];
	} VInt64;
	/**
	 * @~japanese
	 *  パケットに格納されているfloat型データ
	 */
	VFloat v_float;
	/**
	 * @~japanese
	 *  パケットに格納されているdouble型データ
	 */
	VDouble v_double;
	/**
	 * @~japanese
	 * パケットに格納されているint32_t型データ
	 */
	VInt32 v_int32;
	/**
	 * @~japanese
	 * パケットに格納されているにnt64_t型データ
	 */
	VInt64 v_int64;
};

/**
 * SIMDアーキテクチャーを識別する型(MMX)
 */
typedef struct {
	// 64bit
	/**
	 * @~japanese
	 * @brief
	 * このアーキテクチャーにおけるPacket型
	 */
	typedef LIBSAKURA_SYMBOL(SimdPacketMMX) PacketType;
}LIBSAKURA_SYMBOL(SimdArchMMX);

#if defined(__SSE4_2__)

#include <smmintrin.h>

/**
 * @~japanese
 * @brief 	一度のベクトル演算で処理するデータの固まり(Packet)を格納する型(SSE用)である。
 *
 * 各メンバーの詳細は@ref sakura_SimdPacketMMX を参照。
 */
union LIBSAKURA_SYMBOL(SimdPacketSSE) {
	enum {
		kSize = sizeof(__m128 ),
		kNumFloat = kSize / sizeof(float),
		kNumDouble = kSize / sizeof(double),
		kNumInt32 = kSize / sizeof(int32_t),
		kNumInt64 = kSize / sizeof(int64_t),
	};
	typedef __m128 RawFloat;
	typedef __m128d RawDouble;
	typedef __m128i RawInt32;
	typedef __m128i RawInt64;
	RawInt64 raw_int64;
	RawFloat raw_float;
	RawDouble raw_double;
	RawInt32 raw_int32;

	inline void set1(double value) {
		raw_double = _mm_set1_pd(value);
	}
	inline void set1(float value) {
		raw_float = _mm_set1_ps(value);
	}
	inline void set1(int8_t value) {
		raw_int32 = _mm_set1_epi8(value);
	}
	inline void set1(int16_t value) {
		raw_int32 = _mm_set1_epi16(value);
	}
	inline void set1(int32_t value) {
		raw_int32 = _mm_set1_epi32(value);
	}
	inline void set1(int64_t value) {
		raw_int64 = _mm_set1_epi64x(value);
	}

	typedef struct {
		float v[kNumFloat];
	} VFloat;
	typedef struct {
		double v[kNumDouble];
	} VDouble;
	typedef struct {
		int32_t v[kNumInt32];
	} VInt32;
	typedef struct {
		int64_t v[kNumInt64];
	} VInt64;
	VFloat v_float;
	VDouble v_double;
	VInt32 v_int32;
	VInt64 v_int64;

	/**
	 * @~japanese
	 * @brief 	前の世代のSIMDアーキテクチャーのPacket型(@ref sakura_SimdPacketMMX)
	 */
	typedef struct {
		LIBSAKURA_SYMBOL(SimdPacketMMX) v[kSize
				/ LIBSAKURA_SYMBOL(SimdPacketMMX)::kSize];
	} VPrior;
	/**
	 * @~japanese
	 * @brief 	前の世代のSIMDアーキテクチャーのPacket型(@ref sakura_SimdPacketMMX)の配列としてアクセスするためのメンバー
	 */
	VPrior v_prior;
};

/**
 * SIMDアーキテクチャーを識別する型(SSE)
 */
typedef struct {
	// 128bit
	/**
	 * @~japanese
	 * @brief
	 * このアーキテクチャーにおけるPacket型
	 */
	typedef LIBSAKURA_SYMBOL(SimdPacketSSE) PacketType;
	/**
	 * @~japanese
	 * @brief
	 * 前の世代のSIMDアーキテクチャーの型
	 */
	typedef LIBSAKURA_SYMBOL(SimdArchMMX) PriorArch;
}LIBSAKURA_SYMBOL(SimdArchSSE);

#if defined(__AVX__)

#include <immintrin.h>

/**
 * @~japanese
 * @brief 	一度のベクトル演算で処理するデータの固まり(Packet)を格納する型(AVX用)である。
 *
 * 各メンバーの詳細は@ref sakura_SimdPacketMMX を参照。
 */
union LIBSAKURA_SYMBOL(SimdPacketAVX) {
	enum {
		kSize = sizeof(__m256 ),
		kNumFloat = kSize / sizeof(float),
		kNumDouble = kSize / sizeof(double),
		kNumInt32 = kSize / sizeof(int32_t),
		kNumInt64 = kSize / sizeof(int64_t),
	};

	typedef __m256 RawFloat;
	typedef __m256d RawDouble;
	typedef __m256i RawInt32;
	typedef __m256i RawInt64;
	RawInt64 raw_int64;
	RawFloat raw_float;
	RawDouble raw_double;
	RawInt32 raw_int32;

	inline void set1(double value) {
		raw_double = _mm256_set1_pd(value);
	}
	inline void set1(float value) {
		raw_float = _mm256_set1_ps(value);
	}
	inline void set1(int8_t value) {
		raw_int32 = _mm256_set1_epi8(value);
	}
	inline void set1(int16_t value) {
		raw_int32 = _mm256_set1_epi16(value);
	}
	inline void set1(int32_t value) {
		raw_int32 = _mm256_set1_epi32(value);
	}
	inline void set1(int64_t value) {
		raw_int64 = _mm256_set1_epi64x(value);
	}

	typedef struct {
		float v[kNumFloat];
	} VFloat;
	typedef struct {
		double v[kNumDouble];
	} VDouble;
	typedef struct {
		int32_t v[kNumInt32];
	} VInt32;
	typedef struct {
		int64_t v[kNumInt64];
	} VInt64;
	VFloat v_float;
	VDouble v_double;
	VInt32 v_int32;
	VInt64 v_int64;

	/**
	 * @~japanese
	 * @brief 	前の世代のSIMDアーキテクチャーのPacket型(@ref sakura_SimdPacketSSE)
	 */
	typedef struct {
		LIBSAKURA_SYMBOL(SimdPacketSSE) v[kSize
				/ LIBSAKURA_SYMBOL(SimdPacketSSE)::kSize];
	} VPrior;
	/**
	 * @~japanese
	 * @brief 	前の世代のSIMDアーキテクチャーのPacket型(@ref sakura_SimdPacketSSE)の配列としてアクセスするためのメンバー
	 */
	VPrior v_prior;
};

/**
 * SIMDアーキテクチャーを識別する型(AVX)
 */
typedef struct {
	// 256bit
	/**
	 * @~japanese
	 * @brief
	 * このアーキテクチャーにおけるPacket型
	 */
	typedef LIBSAKURA_SYMBOL(SimdPacketAVX) PacketType;
	/**
	 * @~japanese
	 * @brief
	 * 前の世代のSIMDアーキテクチャーの型
	 */
	typedef LIBSAKURA_SYMBOL(SimdArchSSE) PriorArch;
}LIBSAKURA_SYMBOL(SimdArchAVX);
#endif /* defined(__AVX__) */

#if defined(__AVX__)
/**
 * @~japanese
 * @brief サポートされている最新のSIMDアーキテクチャー
 */
typedef LIBSAKURA_SYMBOL(SimdArchAVX) LIBSAKURA_SYMBOL(SimdArchNative);
/**
 * @~japanese
 * @brief サポートされている最新のSIMDアーキテクチャーのパケット型
 */
typedef LIBSAKURA_SYMBOL(SimdPacketAVX) LIBSAKURA_SYMBOL(SimdPacketNative);
#define LIBSAKURA_PACKET_NATIVE_DEFINED 1

#endif /* defined(__AVX__) */

#if ! defined(LIBSAKURA_PACKET_NATIVE_DEFINED)
/**
 * @~japanese
 * @brief サポートされている最新のSIMDアーキテクチャー
 */
typedef LIBSAKURA_SYMBOL(SimdArchSSE) LIBSAKURA_SYMBOL(SimdArchNative);
/**
 * @~japanese
 * @brief サポートされている最新のSIMDアーキテクチャーのパケット型
 */
typedef LIBSAKURA_SYMBOL(SimdPacketSSE) LIBSAKURA_SYMBOL(SimdPacketNative);
#define LIBSAKURA_PACKET_NATIVE_DEFINED 1
#endif

#endif /* defined(__SSE4_2__) */

#if ! defined(LIBSAKURA_PACKET_NATIVE_DEFINED)
/**
 * @~japanese
 * @brief サポートされている最新のSIMDアーキテクチャー
 */
typedef LIBSAKURA_SYMBOL(SimdArchMMX) LIBSAKURA_SYMBOL(SimdArchNative);
/**
 * @~japanese
 * @brief サポートされている最新のSIMDアーキテクチャーのパケット型
 */
typedef LIBSAKURA_SYMBOL(SimdPacketMMX) LIBSAKURA_SYMBOL(SimdPacketNative);
#define LIBSAKURA_PACKET_NATIVE_DEFINED 1
#endif

#endif /* defined(__MMX__) || defined(__x86_64__) */

#endif /* LIBSAKURA_LIBSAKURA_PACKED_TYPE_H_ */
