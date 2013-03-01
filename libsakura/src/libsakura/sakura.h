/**
 * @file
 *  Sakura main header file.
 *
 * sakura.h
 *
 *  Created on: 2013/02/20
 *      Author: kohji
 */

#ifndef SAKURA_H_
#define SAKURA_H_

#include <stddef.h>
#include <stdbool.h>
#include <libsakura/config.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @~english
 * @brief Returns the current time.
 * @~japanese
 * @brief 現在時刻(単位は秒)を返す。精度はgettimeofday(2)依存。
 * @~
 * @author Kohji Nakamura
 */
double libsakura_symbol(currenttime)();

/*
 * memory alignment(for SIMD)
 */
/**
 * @~japanese
 * @brief SAKURAが想定するアライメントに、@a ptr が合っているか調べる
 * @param ptr アラインされているか調べたいアドレス
 * @return アラインされているなら true , そうでないなら false
 */bool libsakura_symbol(is_aligned)(void const *ptr);

/**
 * @~japanese
 * @brief SAKURAがベクトル演算を行う配列に期待するアライメントを返す
 * @return 戻り値の倍数アドレスにアラインされることを期待する
 */
size_t libsakura_symbol (get_alignment)();

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
 */
void const *libsakura_symbol(align_any)(void const *vp, size_t size_of_arena,
		size_t size_required);
float const *libsakura_symbol(align_float)(float const *fp,
		size_t elements_of_arena, size_t elements_required);
double const *libsakura_symbol(align_double)(double const *dp,
		size_t elements_of_arena, size_t elements_required);

/**
 * @ref sakura_statistics の結果を格納する構造体。
 */
typedef struct {
	size_t count; /**< 個数 */
	float sum; /**< 合計 */
	float mean; /**< 平均 */
	float min; /**< 最小 */
	float max; /**< 最大 */
	float rms; /**< 二乗平均平方根 */
	float stddev; /**< 分散 */
} libsakura_symbol(statistics_result);

/**
 * @~japanese
 * @brief 統計値を計算する。どのような統計値を算出するかは
 * @ref sakura_statistics_result を参照。
 * @param result 結果の格納先
 * @param data 対象となるデータ
 * @param mask データのマスク。この値が true だと、
 * 対応する@a data の要素が無視される
 * @param elements @a data 及び@a mask の要素の数
 */
void libsakura_symbol(statistics)
(libsakura_symbol(statistics_result) *result,
		float const data[], bool const mask[], size_t elements);

#ifdef __cplusplus
}
/* extern "C" */
#endif

#endif /* SAKURA_H_ */
