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
	typedef __m64 RawInt64;
	RawFloat raw_float;
	RawDouble raw_double;
	RawInt32 raw_int32;
	RawInt64 raw_int64;

	inline void set1(double value) {
		v_double.v[0] = value;
	}
	inline void set1(float value) {
		v_float.v[0] = value;
		v_float.v[1] = value;
	}
	inline void set1(int8_t value) {
		raw_int64 = _mm_set_pi8(value, value, value, value,
				value, value, value, value);
	}
	inline void set1(int16_t value) {
		raw_int64 =  _mm_set_pi16(value, value, value, value);
	}
	inline void set1(int32_t value) {
		raw_int64 =  _mm_set_pi32(value, value);
	}
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

/*
 * SIMDアーキテクチャーを識別する型
 */
typedef struct {
	// 64bit
	typedef LIBSAKURA_SYMBOL(SimdPacketMMX) PacketType;
}LIBSAKURA_SYMBOL(SimdArchMMX);

#endif /* defined(__MMX__) || defined(__x86_64__) */

#if defined(__SSE4_1__)

#include <smmintrin.h>

/**
 * @~japanese
 * @brief 	一度のベクトル演算で処理するデータの固まり(Packet)を格納する型(SSE用)である。
 *
 * 各メンバーの詳細は@ref sakura_SimdPacketMMX を参照。
 */
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

	LIBSAKURA_SYMBOL(SimdPacketMMX) v_prior[kSize / LIBSAKURA_SYMBOL(SimdPacketMMX)::kSize];
};

typedef struct {
	// 128bit
	typedef LIBSAKURA_SYMBOL(SimdPacketSSE) PacketType;
}LIBSAKURA_SYMBOL(SimdArchSSE);

#endif /* defined(__SSE__) */

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
	RawFloat raw_float;
	RawDouble raw_double;
	RawInt32 raw_int32;
	RawInt32 raw_int64;


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

	LIBSAKURA_SYMBOL(SimdPacketSSE) v_prior[kSize / LIBSAKURA_SYMBOL(SimdPacketSSE)::kSize];
};

typedef struct {
	// 256bit
	typedef LIBSAKURA_SYMBOL(SimdPacketAVX) PacketType;
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


#if defined(__SSE4_1__)
#if ! defined(LIBSAKURA_PACKET_NATIVE_DEFINED)
typedef LIBSAKURA_SYMBOL(SimdArchSSE) LIBSAKURA_SYMBOL(SimdArchNative);
typedef LIBSAKURA_SYMBOL(SimdPacketSSE) LIBSAKURA_SYMBOL(SimdPacketNative);
#define LIBSAKURA_PACKET_NATIVE_DEFINED 1
#endif

#endif /* defined(__SSE__) */


#if defined(__MMX__) || defined(__x86_64__)
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

#endif /* defined(__MMX__) || defined(__x86_64__) */

#endif /* LIBSAKURA_LIBSAKURA_PACKED_TYPE_H_ */
