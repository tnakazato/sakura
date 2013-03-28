/**
 * @file
 *  Sakura main header file.
 *
 * sakura.h
 *
 *  Created on: 2013/02/20
 *      Author: kohji
 */

#ifndef LIBSAKURA_LIBSAKURA_SAKURA_H_
#define LIBSAKURA_LIBSAKURA_SAKURA_H_

#include <stddef.h>
#include <stdbool.h>

#include <libsakura/config.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
	LIBSAKURA_SYMBOL(Status_kOK) = 0, LIBSAKURA_SYMBOL(Status_kNG) = 1
}LIBSAKURA_SYMBOL(Status);
/**
 * @~english
 * @brief Initializes Sakura Library
 * @return Only when sakura_Status_kOK is returned, you can use Sakura Library.
 * @~japanese
 * @brief Sakuraライブラリを初期化する。
 * @~
 * MT-unsafe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Initialize)();

/**
 * @~english
 * @brief Returns the current time.
 * @~japanese
 * @brief 現在時刻(単位は秒)を返す。精度はgettimeofday(2)依存。
 * @~
 * MT-safe
 */
double LIBSAKURA_SYMBOL(GetCurrentTime)();

/*
 * memory alignment(for SIMD)
 */
/**
 * @~japanese
 * @brief SAKURAが想定するアライメントに、@a ptr が合っているか調べる
 * @param ptr アラインされているか調べたいアドレス
 * @return アラインされているなら true , そうでないなら false
 * @~
 * MT-safe
 */

bool LIBSAKURA_SYMBOL(IsAligned)(void const *ptr);

/**
 * @~japanese
 * @brief SAKURAがベクトル演算を行う配列に期待するアライメントを返す
 * @return 戻り値の倍数アドレスにアラインされることを期待する
 * @~
 * MT-safe
 */
size_t LIBSAKURA_SYMBOL (GetAlignment)();

/**
 * @~japanese
 * @brief @a vp がアラインされていないならば、
 * アドレスを必要最小限増加させ、アラインされたアドレスを返す。
 * @a vp がアラインされていれば@a vp を返す。
 * @param vp 確保されている領域の先頭アドレス
 * @param size_of_arena 確保されている領域のサイズ
 * @param size_required アライン後も利用可能でなければならないサイズ
 * @return アラインされたアドレス。もし、 @a size_required を格納するのに
 * 十分な大きさの@a size_of_arena が無いならば、 nullptr を返す。
 * @~
 * MT-safe
 */
void const *LIBSAKURA_SYMBOL(AlignAny)(void const *vp, size_t size_of_arena,
		size_t size_required);
float const *LIBSAKURA_SYMBOL(AlignFloat)(float const *fp,
		size_t elements_of_arena, size_t elements_required);
double const *LIBSAKURA_SYMBOL(AlignDouble)(double const *dp,
		size_t elements_of_arena, size_t elements_required);

/**
 * @ref sakura_ComputeStatistics の結果を格納する構造体。
 */
typedef struct {
	size_t count; /**< 個数 */
	float sum; /**< 合計 */
	float mean; /**< 平均 */
	float rms; /**< 二乗平均平方根 */
	float stddev; /**< 分散 */
	float min; /**< 最小 */
	float max; /**< 最大 */
	int index_of_min; /**< 最小値のインデックス(有効な値がなければ-1) */
	int index_of_max; /**< 最大値のインデックス(有効な値がなければ-1) */
}LIBSAKURA_SYMBOL(StatisticsResult);

/**
 * @~japanese
 * @brief 統計値を計算する。どのような統計値を算出するかは
 * @ref sakura_StatisticsResult を参照。
 * @param data 対象となるデータ。対応する@a is_valid がtrueの場合、NaNであってはならない。
 * @param is_valid データのマスク。この値が false だと、
 * 対応する@a data の要素が無視される
 * @param elements @a data 及び@a is_valid の要素の数
 * @param result 結果の格納先。計算不可能な場合は、構造体内のメンバーの値にNaNが設定される。
 * 同じ値が複数あった場合、どの値のインデックスが@a index_of_min, @a index_of_maxに格納されるかは不定である。
 *
 * @~
 * MT-safe
 */
void LIBSAKURA_SYMBOL(ComputeStatistics)(float const data[],
		bool const is_valid[], size_t elements,
		LIBSAKURA_SYMBOL(StatisticsResult) *result);

/**
 * @~japanese
 * @brief validな値のみを先頭に詰めて昇順にソートする。
 * @param is_valid データのマスク。この値が false だと、
 * 対応する@a data の要素が無視される
 * @param elements @a data 及び@a is_valid の要素の数
 * @param data ソート対象のデータ。In placeでソートするので、この配列内の順序は変更される。
 * 対応する@a is_valid がtrueの場合、NaNであってはならない。
 * @return (validでないデータを除いた)ソートされた要素数( <= elements)
 * @~
 * MT-safe
 */
size_t LIBSAKURA_SYMBOL(SortValidValuesDensely)(bool const is_valid[],
		size_t elements, float data[]);

#ifdef __cplusplus
}
/* extern "C" */
#endif

#endif /* LIBSAKURA_LIBSAKURA_SAKURA_H_ */
