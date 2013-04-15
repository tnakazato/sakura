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

#if defined(__AVX__)

#include <immintrin.h>

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
	RawFloat raw_float;
	RawDouble raw_double;
	RawInt32 raw_int32;
	RawInt32 raw_int64;

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
};

typedef struct {
	// 256bit
	typedef LIBSAKURA_SYMBOL(SimdPacketAVX) PacketType;
}LIBSAKURA_SYMBOL(SimdArchAVX);

typedef LIBSAKURA_SYMBOL(SimdArchAVX) LIBSAKURA_SYMBOL(SimdArchNative);
typedef LIBSAKURA_SYMBOL(SimdPacketAVX) LIBSAKURA_SYMBOL(SimdPacketNative);
#define LIBSAKURA_PACKET_NATIVE_DEFINED 1

#endif /* defined(__AVX__) */

#if defined(__SSE__)

#include <xmmintrin.h>

union LIBSAKURA_SYMBOL(SimdPacketSSE) {
	enum {
		kSize = sizeof(__m128),
		kNumFloat = kSize / sizeof(float),
		kNumDouble = kSize / sizeof(double),
		kNumInt32 = kSize / sizeof(int32_t),
		kNumInt64 = kSize / sizeof(int64_t),
	};
	typedef __m128 RawFloat;
	typedef __m128d RawDouble;
	typedef __m128i RawInt32;
	typedef __m128i RawInt64;
	RawFloat raw_float;
	RawDouble raw_double;
	RawInt32 raw_int32;
	RawInt32 raw_int64;

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
};

typedef struct {
	// 128bit
	typedef LIBSAKURA_SYMBOL(SimdPacketSSE) PacketType;
}LIBSAKURA_SYMBOL(SimdArchSSE);

#if ! defined(LIBSAKURA_PACKET_NATIVE_DEFINED)
typedef LIBSAKURA_SYMBOL(SimdArchSSE) LIBSAKURA_SYMBOL(SimdArchNative);
typedef LIBSAKURA_SYMBOL(SimdPacketSSE) LIBSAKURA_SYMBOL(SimdPacketNative);
#define LIBSAKURA_PACKET_NATIVE_DEFINED 1
#endif

#endif /* defined(__SSE__) */

#if defined(__MMX__)
/* for 64 bit CPU */

/**
 * @brief
 * @~japanese
 * 	一度のベクトル演算で処理するデータの固まり(Packet)を格納する型である。
 */
union LIBSAKURA_SYMBOL(SimdPacketMMX) {
	enum {
		/**
		 * @~japanese
		 *  パケットのサイズ
		 */
		kSize = sizeof(double),
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

	typedef struct {float dummy[kNumFloat];}RawFloat;
	typedef struct {double dummy[kNumDouble];}RawDouble;
	typedef struct {int32_t dummy[kNumInt32];}RawInt32;
	typedef struct {int64_t dummy[kNumInt64];}RawInt64;
	RawFloat raw_float;
	RawDouble raw_double;
	RawInt32 raw_int32;
	RawInt32 raw_int64;

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

/*
 * SIMDアーキテクチャーを識別する型
 */
typedef struct {
	// 64bit
	typedef LIBSAKURA_SYMBOL(SimdPacketMMX) PacketType;
}LIBSAKURA_SYMBOL(SimdArchMMX);

#if ! defined(LIBSAKURA_PACKET_NATIVE_DEFINED)
/**
 * @brief
 * @~japanese
 * 	現在のコンパイル条件におけるネイティブなSIMDアーキテクチャの型を表す。
 */
typedef LIBSAKURA_SYMBOL(SimdArchMMX) LIBSAKURA_SYMBOL(SimdArchNative);
/**
 * @brief
 * @~japanese
 * 	現在のコンパイル条件における最大のパケット型を表す。
 */
	typedef LIBSAKURA_SYMBOL(SimdPacketMMX) LIBSAKURA_SYMBOL(SimdPacketNative);
#define LIBSAKURA_PACKET_NATIVE_DEFINED 1
#endif

#endif /* defined(__MMX__) */

#endif /* LIBSAKURA_LIBSAKURA_PACKED_TYPE_H_ */
