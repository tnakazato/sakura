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
#include <stdint.h>

#include <libsakura/config.h>

#if defined(__GNUC__) || defined(__GNUG__)
#	define LIBSAKURA_WARN_UNUSED_RESULT __attribute__((warn_unused_result))
#else
#	define LIBSAKURA_WARN_UNUSED_RESULT /* Don't ignore result value */
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @~japanese
 * @brief 関数の呼び出し結果を表す
 *
 */
typedef enum {
	/**
	 * @~japanese
	 * @brief 成功または正常
	 */LIBSAKURA_SYMBOL(Status_kOK) = 0, /**
	 * @~japanese
	 * @brief 失敗または異常
	 */LIBSAKURA_SYMBOL(Status_kNG) = 1, /**
	 * @~japanese
	 * @brief 引数が不正だった
	 * @a must-be-aligned 制約に違反している場合も含む。
	 */LIBSAKURA_SYMBOL(Status_kInvalidArgument) = 2, /**
	 * @~japanese
	 * @brief メモリーが足りない
	 */LIBSAKURA_SYMBOL(Status_kNoMemory) = 3, /**
	 * @~japanese
	 * @brief 原因不明のエラー
	 */LIBSAKURA_SYMBOL(Status_kUnknownError) = 99
}LIBSAKURA_SYMBOL(Status);

/**
 * @~english
 * @brief A type of the allocator function called by Sakura Library.
 *
 * Implementation of function of this type must be reentrant.
 * And it must returns a valid pointer for memory region of size 0 if 0 is passed as @a size parameter.
 *
 * @param[in] size	Size of required memory in bytes
 * @return Allocated memory. NULL if failed to allocate.
 *
 * @~japanese
 * @brief Sakuraライブラリが動的にメモリーを確保するときに呼び出す関数の型
 *
 * 関数はリエントラントな実装でなければならない。
 * sizeが0の場合も、(メモリーが確保できるなら)長さ0の領域のポインタを返すこと(NULLを返さないこと)。
 * @param[in] size	確保するサイズ(バイト)
 * @return 確保した領域のアドレス。確保できなかった場合は、NULL。
 * @~
 * MT-safe
 */
typedef void *(*LIBSAKURA_SYMBOL(UserAllocator))(size_t size);

/**
 * @~english
 * @brief A type of the deallocator function called by Sakura Library.
 *
 * Implementation of function of this type must be reentrant.
 * And it must do nothing if NULL is passed as @a pointer parameter.
 *
 * @param[in] pointer	NULL or an address to be released which was allocated by the allocator of type @ref sakura_UserAllocator .
 *
 * @~japanese
 * @brief Sakuraライブラリが動的に確保したメモリーを開放するときに呼び出す関数の型
 *
 * 関数はリエントラントな実装でなければならない。
 * @a pointer にNULLが渡された場合、何も行わないこと。
 * @param[in] pointer	@ref sakura_UserAllocator が返したアドレスまたはNULL。
 * @~
 * MT-safe
 */
typedef void (*LIBSAKURA_SYMBOL(UserDeallocator))(void *pointer);

/**
 * @~english
 * @brief Initializes Sakura Library.
 *
 * Initialize libsakura by calling this function before calling any other function of Sakura Library.
 *
 * Without calling @ref sakura_CleanUp() , don't call this function again.
 * @param[in]	allocator	An allocator which is used when Sakura Library needs to allocate memory dynamically. malloc(3) is used if NULL is provided. See @ref sakura_UserAllocator .
 * @param[in]	deallocator	A deallocator which is used when Sakura Library needs to free dynamically allocated memory. free(3) is used if NULL is provided. See @ref sakura_UserDeallocator .
 * @return Only when sakura_Status_kOK is returned, you can use Sakura Library.
 *
 * @~japanese
 * @brief Sakuraライブラリを初期化する
 *
 * 他の全てのSakuraライブラリAPIの呼び出しに先立って、呼び出すこと。
 * マルチスレッドセーフではないので、単一のスレッドから呼び出すこと。
 * @ref sakura_CleanUp() の呼び出しを挟まず、複数回この関数を呼び出してはならない。
 * @param[in]	allocator	Sakuraライブラリ内で、メモリーを確保するときに呼び出されるアロケーター。NULLの場合はmalloc(3)が使用される。 @ref sakura_UserAllocator 参照
 * @param[in]	deallocator	Sakuraライブラリ内で、メモリーを開放するときに呼び出されるデアロケーター。NULLの場合はfree(3)が使用される。 @ref sakura_UserDeallocator 参照
 * @return @a sakura_Status_kOK が返されたときのみ、Sakuraライブラリを使用できる。
 * @~
 * MT-unsafe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Initialize)(
LIBSAKURA_SYMBOL(UserAllocator) allocator,
LIBSAKURA_SYMBOL(UserDeallocator) deallocator) LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~english
 * @brief Cleans up Sakura Library.
 *
 * When you call this function, no function of Sakura Library must not be running.
 * @~japanese
 * @brief Sakuraライブラリをクリーンアップする
 *
 * マルチスレッドセーフではないので、単一のスレッドから呼び出すこと。
 * また、この関数を呼び出す時点で、他のSakuraライブラリの呼び出しは全て完了しており、以後の呼び出しも発生しないこと。
 * @~
 * MT-unsafe
 */
void LIBSAKURA_SYMBOL(CleanUp)();

/**
 * @~english
 * @brief Returns a current time.
 * Precision of the time depends on gettimeofday(2).
 * @return Current time in seconds since the Epoch.
 *
 * @~japanese
 * @brief 現在時刻(単位は秒)を返す
 *
 * 精度はgettimeofday(2)依存。
 * @return 現在時刻(単位は秒)
 * @~
 * MT-safe
 */
double LIBSAKURA_SYMBOL(GetCurrentTime)();

/*
 * memory alignment(for SIMD)
 */
/**
 * @~japanese
 * @brief Sakuraライブラリが想定するアライメントに、@a ptr が合っているか調べる
 *
 * @param[in] ptr アラインされているか調べたいアドレス。NULL も受け付ける。
 * @return アラインされているなら true , そうでないなら false
 * @~
 * MT-safe
 */

bool LIBSAKURA_SYMBOL(IsAligned)(void const *ptr);

/**
 * @~english
 * @brief Returns an alignment that Sakura Library expects for arrays on which vector operations are performed.
 *
 * @return An expected alignment for arrays marked as must-be-aligned.
 *
 * @~japanese
 * @brief Sakuraライブラリがベクトル演算を行う配列に期待するアライメントを返す
 *
 * @return must-be-alignedとマークされた配列に期待されるアライメント
 * @~
 * MT-safe
 */
size_t LIBSAKURA_SYMBOL (GetAlignment)();

/**
 * @~japanese
 * @brief @a arena がアラインされていないならば、
 * アドレスを必要最小限増加させ、アラインされたアドレスを返す
 *
 * @a arena がアラインされていれば@a arena を返す。
 *
 * @param[in] arena 確保されている領域の先頭アドレス
 * @param[in] size_of_arena 確保されている領域のサイズ
 * @param[in] size_required アライン後も利用可能でなければならないサイズ
 * @return アラインされたアドレス。もし、 @a size_required を格納するのに
 * 十分な大きさの@a size_of_arena が無いならば、 NULL を返す。
 * @~
 * MT-safe
 */
void *LIBSAKURA_SYMBOL(AlignAny)(size_t size_of_arena, void *arena,
		size_t size_required);
/**
 * @~japanese
 * @brief @a arena がアラインされていないならば、
 * アドレスを必要最小限増加させ、アラインされたアドレスを返す
 *
 * @a arena がアラインされていれば@a arena を返す。
 *
 * @param[in] arena 確保されている領域の先頭アドレス
 * @param[in] elements_in_arena 確保されている要素数
 * @param[in] elements_required アライン後も利用可能でなければならない要素数
 * @return アラインされたアドレス。もし、 @a elements_required を格納するのに
 * 十分な大きさの@a elements_in_arena が無いならば、 NULL を返す。
 * @~
 * MT-safe
 */
float *LIBSAKURA_SYMBOL(AlignFloat)(size_t elements_in_arena, float *arena,
		size_t elements_required);
/**
 * @~japanese
 * @brief @a arena がアラインされていないならば、
 * アドレスを必要最小限増加させ、アラインされたアドレスを返す
 *
 * @a arena がアラインされていれば@a arena を返す。
 *
 * @param[in] arena 確保されている領域の先頭アドレス
 * @param[in] elements_in_arena 確保されている要素数
 * @param[in] elements_required アライン後も利用可能でなければならない要素数
 * @return アラインされたアドレス。もし、 @a elements_required を格納するのに
 * 十分な大きさの@a elements_in_arena が無いならば、 NULL を返す。
 * @~
 * MT-safe
 */
double *LIBSAKURA_SYMBOL(AlignDouble)(size_t elements_in_arena, double *arena,
		size_t elements_required);

/**
 * @~japanese
 * @brief @ref sakura_ComputeStatisticsFloat の結果を格納する構造体
 */
typedef struct {
	/**
	 * @~
	 * a number of valid data
	 */
	size_t count;
	/**
	 * @~
	 * a sum of valid data
	 */
	float sum;
	/**
	 * @~
	 * a mean of valid data
	 */
	float mean;
	/**
	 * @~
	 * a root-mean-square of valid data
	 */
	float rms;
	/**
	 * @~
	 * an stddev of valid data
	 */
	float stddev;
	/**
	 * @~
	 * a min value of valid data
	 */
	float min;
	/**
	 * @~
	 * a max value of valid data
	 */
	float max;
	/**
	 * @~
	 * one of index for min value. -1 if there is no valid data.
	 */
	int index_of_min;
	/**
	 * @~
	 * one of index for max value. -1 if there is no valid data.
	 */
	int index_of_max;
}LIBSAKURA_SYMBOL(StatisticsResultFloat);
/**
 * @~japanese
 * @brief 統計値を計算する。どのような統計値を算出するかは
 * @ref sakura_StatisticsResultFloat を参照
 *
 * @param[in] num_data @a data 及び@a is_valid の要素の数。@a num_data <= INT32_MAX
 * @param[in] data 対象となるデータ。対応する@a is_valid がtrueの場合、Inf/NaNであってはならない。
 * <br/>must-be-aligned
 * @param[in] is_valid データのマスク。この値が false だと、
 * 対応する@a data の要素が無視される。
 * <br/>must-be-aligned
 * @param[out] result 結果の格納先。計算不可能な場合は、構造体内のメンバーの値にNaNが設定される。
 * 同じ値が複数あった場合、どの値のインデックスが@a index_of_min, @a index_of_maxに格納されるかは不定である。
 * @return 終了ステータス
 *
 * @~english
 * @brief Computes statistics. Refer to
 * @ref sakura_StatisticsResultFloat to see what kind of statistics are computed.
 *
 * @param[in] num_data A number of elements in @a data and @a is_valid . @a num_data <= INT32_MAX
 * @param[in] data Data. If corresponding element in @a is_valid is true, the element in @a data must not be Inf nor NaN.
 * <br/>must-be-aligned
 * @param[in] is_valid Masks of @a data. If a value of element is false,
 * the corresponding element in @a data is ignored.
 * <br/>must-be-aligned
 * @param[out] result An address where the result should be stored. Some fields may be set to NaN if it is impossible to figure out.
 * If there is more than one occurrence of min or max value, it is undefined which index of the occurrence is selected for @a index_of_min or @a index_of_max.
 * @return status code
 *
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeStatisticsFloat)(
		size_t num_data, float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result)
				LIBSAKURA_WARN_UNUSED_RESULT;
/**
 * @~japanese
 * @brief validな値のみを先頭に詰めて昇順にソートする
 *
 * @param[in] is_valid データのマスク。この値が false だと、
 * 対応する@a data の要素が無視される
 * @param[in] num_data @a data 及び@a is_valid の要素の数
 * @param[in,out] data ソート対象のデータ。In placeでソートするので、この配列内の順序は変更される。
 * 対応する@a is_valid がtrueの場合、NaNであってはならない。
 * @param[out] new_num_data (validでないデータを除いた)ソートされた要素数( <= @a num_data ) の格納先
 * @return 終了ステータス
 *
 * @~english
 * @brief Sort only valid data in ascending order.
 *
 * @param[in] num_data A number of elements in @a data and @a is_valid .
 * @param[in] is_valid Masks of @a data. If a value of element is false,
 * the corresponding element in @a data is ignored.
 * @param[in,out] data Data to be sorted. Since data is sorted in place, contents of this array are not preserved.
 * If corresponding element in @a is_valid is true, the element in @a data must not be Inf nor NaN.
 * @param[out] new_num_data A number of sorted elements that don't include invalid data( <= @a num_data ) is stored here.
 * @return status code
 *
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(
		size_t num_data, bool const is_valid[], float data[],
		size_t *new_num_data) LIBSAKURA_WARN_UNUSED_RESULT;
/**
 * @~japanese
 * @brief 畳み込みしながらグリッドする
 *
 * @a start_spectrum から @a end_spectrum までの、@a x , @a y 座標上の @a value の値を
 * @a convolution_table で表される広がりを持った点として、@a grid 上にプロットする。
 * エッジから2*supportピクセル分の値は信頼できない(理由:処理の高速化のため)。
 * @a grid にプロットする際は、@a polarization_map , @a channel_map によって偏波とチャネルのマッピングが行われる。
 *
 * 各浮動小数点の数値はNaN/+-Infであってはならない。
 * @param[in] num_spectra 次の関係でなければならない。 0 <= start_spectrum <= end_spectrum <= num_spectra
 * @param[in] start_spectrum 開始spectrumの添字
 * @param[in] end_spectrum 終了spectrumの添字+1
 * @param[in] spectrum_mask	要素数はnum_spectra。falseのスペクトルは無視される。<br/>must-be-aligned
 * @param[in] x 要素数は@a num_spectra 。2次元平面に投射済みのx座標。範囲は、INT32_MIN < x[i] < INT32_MAX<br/>must-be-aligned
 * @param[in] y 要素数は@a num_spectra 。2次元平面に投射済みのy座標。範囲は、INT32_MIN < y[i] < INT32_MAX<br/>must-be-aligned
 * @param[in] support @a width x @a height 平面における畳み込みカーネルの広がり(中心か らのpixel数)。範囲は、0 < support <= 256<br/>
 * ただし、@a support * @a sampling に比例するサイズの領域をスタック上に確保するので、@a support * @a sampling が大きな値となる場合、スタックオーバーフローを起こす。
 * @param[in] sampling 畳み込みカーネルの精度(/pixel)。範囲は、0 < sampling <= INT32_MAX
 * @param[in] num_polarizations 範囲は、0 < num_polarizations <= INT32_MAX
 * @param[in] polarization_map	要素数は、num_polarizations。各要素の値は、[0,num_polarizations_for_grid)でなければならない。<br/>must-be-aligned
 * @param[in] num_channels 範囲は、0 < num_channels <= INT32_MAX
 * @param[in] channel_map	要素数は、num_channels。各要素の値は、[0,num_channels_for_grid)でなければならない。<br/>must-be-aligned
 * @param[in] mask	要素のレイアウトは、[num_spectra][num_polarizations][num_channels]。falseの場合は、該当するスペクトル、偏波、チ ャネルの組み合わせのデータは無視される。<br/>must-be-aligned
 * @param[in] value	要素のレイアウトは、[num_spectra][num_polarizations][num_channels]。グリッディングする値。<br/>must-be-aligned
 * @param[in] weight 要素のレイアウトは、[num_spectra][num_channels]。重み。<br/>must-be-aligned
 * @param[in] weight_only @a value に重みを掛けたものではなく、重み自体をグリッディングする場合は、true。
 * @param[in] num_convolution_table @a convolution_table の要素数。 範囲は、ceil(sqrt(2.)*(support+1)*sampling) <= num_convolution_table <= INT32_MAX / 32
 * @param[in] convolution_table	要素数は、@a num_convolution_table 。畳み込みに使用する重みカーブ。各要素の値は、NaN,Infであってはならない。要素0が中心を表す。<br/>must-be-aligned
 * @param[in] num_polarizations_for_grid 範囲は、0 < num_polarizations_for_grid <= INT32_MAX
 * @param[in] num_channels_for_grid 範囲は、0 < num_channels_for_grid <= INT32_MAX
 * @param[in] width 範囲は、0 < width <= INT32_MAX
 * @param[in] height 範囲は、0 < height <= INT32_MAX
 * @param[out] weight_sum	要素のレイアウトは、[num_polarizations_for_grid][num_channels_for_grid]。重みの合計。<br/>must-be-aligned
 * @param[out] weight_of_grid	要素のレイアウトは、[height][width][num_polarizations_for_grid][num_channels_for_grid]。グリッドの重み。<br/>must-be-aligned
 * @param[out] grid	要素のレイアウトは、[height][width][num_polarizations_for_grid][num_channels_for_grid]。グリッディング結果。<br/>must-be-aligned
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GridConvolvingFloat)(
		size_t num_spectra, size_t start_spectrum, size_t end_spectrum,
		bool const spectrum_mask[/*num_spectra*/],
		double const x[/*num_spectra*/], double const y[/*num_spectra*/],
		size_t support, size_t sampling, size_t num_polarizations,
		uint32_t const polarization_map[/*num_polarizations*/],
		size_t num_channels, uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_spectra][num_polarizations]*/[/*num_channels*/],
		float const value/*[num_spectra][num_polarizations]*/[/*num_channels*/],
		float const weight/*[num_spectra]*/[/*num_channels*/], bool weight_only,
		size_t num_convolution_table/*= ceil(sqrt(2.)*(support+1)*sampling)*/,
		float const convolution_table[/*num_convolution_table*/],
		size_t num_polarizations_for_grid, size_t num_channels_for_grid,
		size_t width, size_t height,
		double weight_sum/*[num_polarizations_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarizations_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarizations_for_grid]*/[/*num_channels_for_grid*/])
				LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~english
 * @brief Returns if the values in input array are in any of specified range (inclusive).
 * @details Returns true if the corresponding element in the input array
 * is in range of upper and lower boundary pairs,
 * @par
 * @a lower_bound[k] <= @a data[i] <= @a upper_bound[k].
 *
 * The function takes more than one upper and lower boundary pairs as arrays,
 * @a lower_bounds and @a upper_bounds.@n
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 * All elements in @a result are set to false when no condition is given, i.e., @a num_condition = 0.
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data An input array of size, @a num_data.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] num_condition The number of elements in the arrays, @a lower_bounds
 * and @a upper_bounds.
 * @param[in] lower_bounds The input array of size, @a num_condition.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] upper_bounds The input array of size, @a num_condition.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return status code
 * @~japanese
 * @brief 入力配列の値が、与えられた下限値と上限値の組の範囲に入っているかを検定する。(inclusive).
 * @details 複数の下限値( @a lower_bounds ) と上限値 ( @a upper_bounds ) の組を配列として取り、
 * 入力配列の要素の値がそれらのいずれかの範囲に含まれていれば真を返す。@n
 * すなわち、
 * @par
 * @a lower_bound[k] <= @a data[i] <= @a upper_bound[k].
 *
 * を検定する。@n
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。@n
 * 上限値と下限値が与えられなければ (@a num_condition = 0)、@a result は全ての要素がfalseとなる。
 *
 * @param[in] num_data 一次元配列@a data 及び@a result の要素の数。
 * @param[in] data 入力一次元配列。検定の対象となる値を格納する。要素数は@a num_data でなければならない。
 * 配列が浮動少数点型の場合、要素がInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[in] num_condition 一次元配列@a lower_bounds 及び@a upper_bounds の要素の数。
 * 下限値と上限値の組の数を表す。
 * @param[in] lower_bounds 入力一次元配列。検定条件の下限値を格納する。要素数は@a num_condition でなければならない。
 * 配列が浮動少数点型の場合、要素がInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[in] upper_bounds 入力一次元配列。検定条件の上限値を格納する。要素数は@a num_condition でなければならない。
 * 配列が浮動少数点型の場合、要素がInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatInRangesInclusive)(
		size_t num_data, float const data[/*num_data*/], size_t num_condition,
		float const lower_bounds[/*num_condition*/],
		float const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]);

/**
 * @copybrief sakura_SetTrueFloatInRangesInclusive
 * @copydetails sakura_SetTrueFloatInRangesInclusive
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntInRangesInclusive)(
		size_t num_data, int const data[/*num_data*/], size_t num_condition,
		int const lower_bounds[/*num_condition*/],
		int const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]);

/**
 * @~english
 * @brief Returns if the values in input array are in any of specified range (exclusive).
 * @details Returns true if the corresponding element in the input array
 * is in between upper and lower boundary pairs,
 * @par
 * @a lower_bound[k] < @a data[i] < @a upper_bound[k].
 *
 * The function takes more than one upper and lower boundary pairs as arrays,
 * @a lower_bounds and @a upper_bounds.@n
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 * All elements in @a result are set to false when no condition is given, i.e., @a num_condition = 0.
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data An input array of size, @a num_data.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] num_condition The number of elements in the arrays, @a lower_bounds
 * and @a upper_bounds.
 * @param[in] lower_bounds The input array of size, @a num_condition.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] upper_bounds The input array of size, @a num_condition.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return status code
 * @~japanese
 * @brief 入力配列の値が、与えられた下限値と上限値の組の範囲に入っているかを検定する。(exclusive).
 * @details 複数の下限値( @a lower_bounds ) と上限値 ( @a upper_bounds ) の組を配列として取り、
 * 入力配列の要素の値がそれらのいずれかの範囲に含まれていれば真を返す。@n
 * すなわち、
 * @par
 * @a lower_bound[k] < @a data[i] < @a upper_bound[k].
 *
 * を検定する。@n
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。@n
 * 上限値と下限値が与えられなければ (@a num_condition = 0)、@a result は全ての要素がfalseとなる。
 *
 * @param[in] num_data 一次元配列@a data 及び@a result の要素の数。
 * @param[in] data 入力一次元配列。検定の対象となる値を格納する。要素数は@a num_data でなければならない。
 * 配列が浮動少数点型の場合、要素がInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[in] num_condition 一次元配列@a lower_bounds 及び@a upper_bounds の要素の数。
 * 下限値と上限値の組の数を表す。
 * @param[in] lower_bounds 入力一次元配列。検定条件の下限値を格納する。要素数は@a num_condition でなければならない。
 * 配列が浮動少数点型の場合、要素がInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[in] upper_bounds 入力一次元配列。検定条件の上限値を格納する。要素数は@a num_condition でなければならない。
 * 配列が浮動少数点型の場合、要素がInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatInRangesExclusive)(
		size_t num_data, float const data[/*num_data*/], size_t num_condition,
		float const lower_bounds[/*num_condition*/],
		float const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]);

/**
 * @copybrief sakura_SetTrueFloatInRangesExclusive
 * @copydetails sakura_SetTrueFloatInRangesExclusive
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntInRangesExclusive)(
		size_t num_data, int const data[/*num_data*/], size_t num_condition,
		int const lower_bounds[/*num_condition*/],
		int const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]);

/**
 * @~english
 * @brief Returns if the values in input array are greater than a threshold.
 * @details Returns true if the corresponding element in the input array
 * is greater than a @a threshold,
 * @par
 * @a data[i] > @a threshold .
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data An input array of size, @a num_data.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] threshold the threshold of evaluation.
 * In case the parameter is floating-point type, the value should not be Inf nor NaN.
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return status code
 * @~japanese
 * @brief 入力配列の値が、与えられたしきい値より大きいかどうかを検定する。
 * @details 入力配列の要素の値がしきい値( @a threshold )より大きければ真を返す。@n
 * すなわち、
 * @par
 * @a data[i] > @a threshold .
 *
 * を検定する。@n
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。@n
 *
 * @param[in] num_data 一次元配列@a data 及び@a result の要素の数。
 * @param[in] data 入力一次元配列。検定の対象となる値を格納する。要素数は@a num_data でなければならない。
 * 配列が浮動少数点型の場合、要素がInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[in] threshold しきい値。引数が浮動少数点型の場合、値がInfやNaNであってはならない。InfやNaNの場合の動作は不定。
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatGreaterThan)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]);

/**
 * @copybrief sakura_SetTrueFloatGreaterThan
 * @copydetails sakura_SetTrueFloatGreaterThan
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntGreaterThan)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]);

/**
 * @~english
 * @brief Returns if the values in input array are greater than or equal to a threshold.
 * @details Returns true if the corresponding element in the input array
 * is greater than or equals to a @a threshold,
 * @par
 * @a data[i] >= @a threshold .
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data An input array of size, @a num_data.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] threshold the threshold of evaluation.
 * In case the parameter is floating-point type, the value should not be Inf nor NaN.
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return status code
 * @~japanese
 * @brief 入力配列の値が、与えられたしきい値以上かどうかを検定する。
 * @details 入力配列の要素の値がしきい値( @a threshold )以上であれば真を返す。@n
 * すなわち、
 * @par
 * @a data[i] >= @a threshold .
 *
 * を検定する。@n
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。@n
 *
 * @param[in] num_data 一次元配列@a data 及び@a result の要素の数。
 * @param[in] data 入力一次元配列。検定の対象となる値を格納する。要素数は@a num_data でなければならない。
 * 配列が浮動少数点型の場合、要素がInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[in] threshold しきい値。引数が浮動少数点型の場合、値がInfやNaNであってはならない。InfやNaNの場合の動作は不定。
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatGreaterThanOrEquals)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]);

/**
 * @copybrief sakura_SetTrueFloatGreaterThanOrEquals
 * @copydetails sakura_SetTrueFloatGreaterThanOrEquals
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntGreaterThanOrEquals)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]);

/**
 * @~english
 * @brief Returns if the values in input array are less than a threshold.
 * @details Returns true if the corresponding element in the input array
 * is less than a @a threshold,
 * @par
 * @a data[i] < @a threshold .
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data An input array of size, @a num_data.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] threshold the threshold of evaluation.
 * In case the parameter is floating-point type, the value should not be Inf nor NaN.
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return status code
 * @~japanese
 * @brief 入力配列の値が、与えられたしきい値より小さいかどうかを検定する。
 * @details 入力配列の要素の値がしきい値( @a threshold )より小さければ真を返す。@n
 * すなわち、
 * @par
 * @a data[i] < @a threshold .
 *
 * を検定する。@n
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。@n
 *
 * @param[in] num_data 一次元配列@a data 及び@a result の要素の数。
 * @param[in] data 入力一次元配列。検定の対象となる値を格納する。要素数は@a num_data でなければならない。
 * 配列が浮動少数点型の場合、要素がInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[in] threshold しきい値。引数が浮動少数点型の場合、値がInfやNaNであってはならない。InfやNaNの場合の動作は不定。
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatLessThan)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]);

/**
 * @copybrief sakura_SetTrueFloatLessThan
 * @copydetails sakura_SetTrueFloatLessThan
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntLessThan)(size_t num_data,
		int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]);

/**
 * @~english
 * @brief Returns if the values in input array are less than or equal to a threshold.
 * @details Returns true if the corresponding element in the input array
 * is less than or equals to a @a threshold,
 * @par
 * @a data[i] <= @a threshold .
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data An input array of size, @a num_data.
 * In case the array is floating-point type, the elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] threshold the threshold of evaluation.
 * In case the parameter is floating-point type, the value should not be Inf nor NaN.
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return status code
 * @~japanese
 * @brief 入力配列の値が、与えられたしきい値以下かどうかを検定する。
 * @details 入力配列の要素の値がしきい値( @a threshold )以下であれば真を返す。@n
 * すなわち、
 * @par
 * @a data[i] <= @a threshold .
 *
 * を検定する。@n
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。@n
 *
 * @param[in] num_data 一次元配列@a data 及び@a result の要素の数。
 * @param[in] data 入力一次元配列。検定の対象となる値を格納する。要素数は@a num_data でなければならない。
 * 配列が浮動少数点型の場合、要素がInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[in] threshold しきい値。引数が浮動少数点型の場合、値がInfやNaNであってはならない。InfやNaNの場合の動作は不定。
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueFloatLessThanOrEquals)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]);

/**
 * @copybrief sakura_SetTrueFloatLessThanOrEquals
 * @copydetails sakura_SetTrueFloatLessThanOrEquals
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIntLessThanOrEquals)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]);

/**
 * @~english
 * @brief Returns if the values in input array are finite numbers.
 * @details Returns false if the corresponding element in the input array
 * is not a number (NaN) or infinity.
 *
 * @note No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data The input array of of size, @a num_data.
 * @n must-be-aligned
 * @param[out] result The output array of of size, @a num_data.
 * @n must-be-aligned
 * @return status code
 * @~japanese
 * @brief 入力配列の値が、有限の値かどうかを検定する。
 * @details 入力配列の要素の値が非数（NaN）または無限大（Inf）ならば偽を返す。
 *
 * @note 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。
 *
 * @param[in] num_data @a data 及び@a result の要素の数。
 * @param[in] data 入力配列。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetFalseFloatIfNanOrInf)(
		size_t num_data, float const data[/*num_data*/],
		bool result[/*num_data*/]);

/**
 * @~english
 * @brief Convert an input array to a boolean array.
 * @details Returns true if the corresponding element in input array != 0.
 * @note No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data The input array of of size, @a num_data.
 * @n must-be-aligned
 * @param[out] result The output array of of size, @a num_data.
 * @n must-be-aligned
 * @return status code
 * @~japanese
 * @brief 入力配列を論理値の配列に変換する。
 * @details 入力配列の対応する要素に、値が1のビットがひとつでもあれば、trueを返す。
 * @note 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。
 *
 * @param[in] num_data @a data 及び@a result の要素の数。
 * @param[in] data 入力配列。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint8ToBool)(size_t num_data,
		uint8_t const data[/*num_data*/], bool result[/*num_data*/]);

/**
 * @copybrief sakura_Uint8ToBool
 * @copydetails sakura_Uint8ToBool
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint32ToBool)(size_t num_data,
		uint32_t const data[/*num_data*/], bool result[/*num_data*/]);

/**
 * @~english
 * @brief Invert a boolean array
 * @note No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data The input array of of size, @a num_data.
 * @n must-be-aligned
 * @param[out] result The output array of of size, @a num_data.
 * @n must-be-aligned
 * @return status code
 * @~japanese
 * @brief 入力配列を論理反転する。
 * @note 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。
 *
 * @param[in] num_data @a data 及び@a result の要素の数。
 * @param[in] data 入力配列。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InvertBool)(size_t num_data,
bool const data[/*num_data*/], bool result[/*num_data*/]);

/**
 * @~english
 * @brief Invoke bit operation AND between a a bit mask and an array.
 * @details Invokes the following bit operation to @a i- th element of @a result :
 * @code
 * result [i] = edit_mask[i] ? (data[i] & bit_maks) : data[i]
 * @endcode
 *
 * @note
 * No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @note
 * The function can also be used to invoke material nonimplication of @a data and @a bit_mask .
 * Input the complement of @a bit_mask (~@a bit_mask ) to invoke material nonimplication.
 * For details of mateial nonimplication, see, e.g.,@n
 * http://en.wikipedia.org/wiki/Truth_table
 *
 * @param[in] bit_mask A bit mask. The bit operation is invoked
 * between this value and the array, @a data.
 * @param[in] num_data The number of elements in the arrays, @a data,
 * @a edit_mask, and @a result.
 * @param[in] data An input array of size, @a num_data. The bit operation
 * is invoked between this array and @a bit_mask.@n
 * must-be-aligned
 * @param[in] edit_mask A boolean mask array of size, @a num_data. The bit operation
 * is skipped for the elements with the value, false.@n
 * must-be-aligned
 * @param[out] result The output array of size, @a num_data. It stores the result
 * of the bit operation between @a bit_mask and @a data. The bit operation is skipped
 * and the value in array, @a data, is adopted for the elements where corresponding
 * elements in @a edit_mask is false. The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.@n
 * must-be-aligned
 * @return status code
 * @~japanese
 * @brief ビットマスクと一次元配列のビット積を取る。
 * @details 配列の@a i- 番目の要素に対して次の演算を行い、出力@a result を返す:
 * @code
 * result [i] = edit_mask[i] ? (data[i] & bit_maks) : data[i]
 * @endcode
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。
 *
 * @note
 * この関数は、@a data と@a bit_mask の間の非含意のビット演算にも使用できる。
 * 非含意のビット演算を実行するときは、@a bit_mask をビット反転させたもの
 * ( ~@a bit_mask )を関数の入力として与える。非含意のビット演算についての詳細は、例えば次のページを参照:@n
 * http://ja.wikipedia.org/wiki/%E8%AB%96%E7%90%86%E6%BC%94%E7%AE%97
 *
 * @param[in] bit_mask ビットマスク
 * @param[in] num_data @a data, @a edit_mask 及び@a result の要素の数。
 * @param[in] data 入力配列。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] edit_mask データのマスク。要素数は@a num_data でなければならない。
 * この値が true だと、対応する入力配列@a data とビットマスク@a bit_maks のビット積を計算する。
 * この値が false だと、その要素のビット演算は行われず、対応する入力配列@a data の要素がそのまま出力となる。
 * @n must-be-aligned
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。インプレースな変換を許す(@a result == @a data)。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 *
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8And)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/]);
/**
 * @copybrief sakura_OperateBitsUint8And
 * @copydetails sakura_OperateBitsUint8And
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32And)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/]);

/**
 * @~english
 * @brief Invoke bit operation, Converse Nonimplication, between a bit mask and an array.
 * @details Invokes the following bit operation to the @a i- th element of @a result :
 * @code
 * result [i] = edit_mask[i] ? (~data[i] & bit_maks) : data[i]
 * @endcode
 *
 * @note
 * No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @note
 * The function can also be used to invoke bitwise NOR operation of @a data and @a bit_mask .
 * Input the complement of @a bit_mask (~@a bit_mask ) to invoke bitwise NOR operation.
 * For details of bitwise NOR operation, see, e.g.,@n
 * http://en.wikipedia.org/wiki/Truth_table
 *
 * @param[in] bit_mask A bit mask. The bit operation is invoked
 * between this value and the array, @a data.
 * @param[in] num_data The number of elements in the arrays, @a data,
 * @a edit_mask, and @a result.
 * @param[in] data An input array of size, @a num_data. The bit operation
 * is invoked between this array and @a bit_mask.@n
 * must-be-aligned
 * @param[in] edit_mask A boolean mask array of size, @a num_data. The bit operation
 * is skipped for the elements with the value, false.@n
 * must-be-aligned
 * @param[out] result The output array of size, @a num_data. It stores the result
 * of the bit operation between @a bit_mask and @a data. The bit operation is skipped
 * and the value in array, @a data, is adopted for the elements where corresponding
 * elements in @a edit_mask is false. The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.@n
 * must-be-aligned
 * @return status code
 * @~japanese
 * @brief ビットマスクと一次元配列の非逆含意ビット演算を実行する。
 * @details 配列の@a i- 番目の要素に対して次の演算を行い、出力@a result を返す:
 * @code
 * result [i] = edit_mask[i] ? (~data[i] & bit_maks) : data[i]
 * @endcode
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。
 *
 * @note
 * この関数は、@a data と@a bit_mask の間の否定論理和ビット演算(NOR)にも使用できる。
 * 否定論理和のビット演算を実行するときは、@a bit_mask をビット反転させたもの
 * ( ~@a bit_mask )を関数の入力として与える。
 * 否定論理和ビット演算についての詳細は、例えば次のページを参照:@n
 * http://ja.wikipedia.org/wiki/%E8%AB%96%E7%90%86%E6%BC%94%E7%AE%97
 *
 *
 * @param[in] bit_mask ビットマスク
 * @param[in] num_data @a data, @a edit_mask 及び@a result の要素の数。
 * @param[in] data 入力配列。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] edit_mask データのマスク。要素数は@a num_data でなければならない。
 * この値が true だと、対応する入力配列@a data とビットマスク@a bit_maks のビット演算を実行する。
 * この値が false だと、その要素のビット演算は行われず、対応する入力配列@a data の要素がそのまま出力となる。
 * @n must-be-aligned
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。インプレースな変換を許す(@a result == @a data)。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 *
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8ConverseNonImplication)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/]);
/**
 * @copybrief sakura_OperateBitsUint8ConverseNonImplication
 * @copydetails sakura_OperateBitsUint8ConverseNonImplication
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32ConverseNonImplication)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/]);

/**
 * @~english
 * @brief Invoke bit operation, Material Implication, between a bit mask and an array.
 * @details Invokes the following bit operation to the @a i- th element of @a result :
 * @code
 * result [i] = edit_mask[i] ? (~data[i] | bit_maks) : data[i]
 * @endcode
 *
 * @note
 * No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @note
 * The function can also be used to invoke bitwise NAND operation of @a data and @a bit_mask .
 * Input the complement of @a bit_mask (~@a bit_mask ) to invoke bitwise NAND operation.
 * For details of bitwise NAND operation, see, e.g.,@n
 * http://en.wikipedia.org/wiki/Truth_table
 *
 * @param[in] bit_mask A bit mask. The bit operation is invoked
 * between this value and the array, @a data.
 * @param[in] num_data The number of elements in the arrays, @a data,
 * @a edit_mask, and @a result.
 * @param[in] data An input array of size, @a num_data. The bit operation
 * is invoked between this array and @a bit_mask.@n
 * must-be-aligned
 * @param[in] edit_mask A boolean mask array of size, @a num_data. The bit operation
 * is skipped for the elements with the value, false.@n
 * must-be-aligned
 * @param[out] result The output array of size, @a num_data. It stores the result
 * of the bit operation between @a bit_mask and @a data. The bit operation is skipped
 * and the value in array, @a data, is adopted for the elements where corresponding
 * elements in @a edit_mask is false. The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.@n
 * must-be-aligned
 * @return status code
 * @~japanese
 * @brief ビットマスクと一次元配列の含意ビット演算を実行する。
 * @details 配列の@a i- 番目の要素に対して次の演算を行い、出力@a result を返す:
 * @code
 * result [i] = edit_mask[i] ? (~data[i] | bit_maks) : data[i]
 * @endcode
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。
 *
 * @note
 * この関数は、@a data と@a bit_mask の間の否定論理積ビット演算(NAND)にも使用できる。
 * 否定論理積のビット演算を実行するときは、@a bit_mask をビット反転させたもの
 * ( ~@a bit_mask )を関数の入力として与える。
 * 否定論理積ビット演算についての詳細は、例えば次のページを参照:@n
 * http://ja.wikipedia.org/wiki/%E8%AB%96%E7%90%86%E6%BC%94%E7%AE%97
 *
 *
 * @param[in] bit_mask ビットマスク
 * @param[in] num_data @a data, @a edit_mask 及び@a result の要素の数。
 * @param[in] data 入力配列。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] edit_mask データのマスク。要素数は@a num_data でなければならない。
 * この値が true だと、対応する入力配列@a data とビットマスク@a bit_maks のビット演算を実行する。
 * この値が false だと、その要素のビット演算は行われず、対応する入力配列@a data の要素がそのまま出力となる。
 * @n must-be-aligned
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。インプレースな変換を許す(@a result == @a data)。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 *
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8Implication)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/]);
/**
 * @copybrief sakura_OperateBitsUint8Implication
 * @copydetails sakura_OperateBitsUint8Implication
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32Implication)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/]);

/**
 * @~english
 * @brief Invoke bitwise NOT operation of an array.
 * @details Invokes the following bit operation to @a i- th element of @a result :
 * @code
 * result [i] = edit_mask[i] ? ~data[i] : data[i]
 * @endcode
 *
 * @note
 * No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @param[in] num_data The number of elements in the arrays, @a data,
 * @a edit_mask, and @a result.
 * @param[in] data An input array of size, @a num_data. The bit operation
 * is invoked to this array.@n
 * must-be-aligned
 * @param[in] edit_mask A boolean mask array of size, @a num_data. The bit operation
 * is skipped for the elements with the value, false.@n
 * must-be-aligned
 * @param[out] result The output array of size, @a num_data. It stores the result
 * of the bit operation to @a data. The bit operation is skipped
 * and the value in array, @a data, is adopted for the elements where corresponding
 * elements in @a edit_mask is false. The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.@n
 * must-be-aligned
 * @return status code
 * @~japanese
 * @brief 一次元配列をビット反転する。
 * @details 配列の@a i- 番目の要素に対して次の演算を行い、出力@a result を返す:
 * @code
 * result [i] = edit_mask[i] ? ~data[i] : data[i]
 * @endcode
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。
 *
 * @param[in] num_data @a data, @a edit_mask 及び@a result の要素の数。
 * @param[in] data 入力配列。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] edit_mask データのマスク。要素数は@a num_data でなければならない。
 * この値が true だと、対応する入力配列@a data に対するビット反転が実行される。
 * この値が false だと、その要素のビット演算は行われず、対応する入力配列@a data の要素がそのまま出力となる。
 * @n must-be-aligned
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。インプレースな変換を許す(@a result == @a data)。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 *
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8Not)(size_t num_data,
		uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/]);
/**
 * @copybrief sakura_OperateBitsUint8Not
 * @copydetails sakura_OperateBitsUint8Not
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32Not)(
		size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/]);

/**
 * @~english
 * @brief Invoke bit operation OR between a a bit mask and an array.
 * @details Invokes the following bit operation to @a i- th element of @a result :
 * @code
 * result [i] = edit_mask[i] ? (data[i] | bit_maks) : data[i]
 * @endcode
 *
 * @note
 * No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @note
 * The function can also be used to invoke converse implication of @a data and @a bit_mask .
 * Input the complement of @a bit_mask (~@a bit_mask ) to invoke converse implication.
 * For details of converse implication, see, e.g.,@n
 * http://en.wikipedia.org/wiki/Truth_table
 *
 * @param[in] bit_mask A bit mask. The bit operation is invoked
 * between this value and the array, @a data.
 * @param[in] num_data The number of elements in the arrays, @a data,
 * @a edit_mask, and @a result.
 * @param[in] data An input array of size, @a num_data. The bit operation
 * is invoked between this array and @a bit_mask.@n
 * must-be-aligned
 * @param[in] edit_mask A boolean mask array of size, @a num_data. The bit operation
 * is skipped for the elements with the value, false.@n
 * must-be-aligned
 * @param[out] result The output array of size, @a num_data. It stores the result
 * of the bit operation between @a bit_mask and @a data. The bit operation is skipped
 * and the value in array, @a data, is adopted for the elements where corresponding
 * elements in @a edit_mask is false. The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.@n
 * must-be-aligned
 * @return status code
 * @~japanese
 * @brief ビットマスクと一次元配列のビット和を取る。
 * @details 配列の@a i- 番目の要素に対して次の算を行い、出力@a result を返す:
 * @code
 * result [i] = edit_mask[i] ? (data[i] | bit_maks) : data[i]
 * @endcode
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。
 *
 * @note
 * この関数は、@a data と@a bit_mask の間の逆含意のビット演算にも使用できる。
 * 逆含意のビット演算を実行するときは、@a bit_mask をビット反転させたもの
 * ( ~@a bit_mask )を関数の入力として与える。逆含意のビット演算についての詳細は、例えば次のページを参照:@n
 * http://ja.wikipedia.org/wiki/%E8%AB%96%E7%90%86%E6%BC%94%E7%AE%97
 *
 * @param[in] bit_mask ビットマスク
 * @param[in] num_data @a data, @a edit_mask 及び@a result の要素の数。
 * @param[in] data 入力配列。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] edit_mask データのマスク。要素数は@a num_data でなければならない。
 * この値が true だと、対応する入力配列@a data とビットマスク@a bit_maks のビット和を計算する。
 * この値が false だと、その要素のビット演算は行われず、対応する入力配列@a data の要素がそのまま出力となる。
 * @n must-be-aligned
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。インプレースな変換を許す(@a result == @a data)。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 *
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8Or)(uint8_t bit_mask,
		size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/]);
/**
 * @copybrief sakura_OperateBitsUint8Or
 * @copydetails sakura_OperateBitsUint8Or
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32Or)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/]);

/**
 * @~english
 * @brief Invoke bit operation XOR between a bit mask and an array.
 * @details Invokes the following bit operation to @a i- th element of @a result :
 * @code
 * result [i] = edit_mask[i] ? (data[i] ^ bit_maks) : data[i]
 * @endcode
 *
 * @note
 * No operation is done when the data array is zero length, i.e., @a num_data = 0.
 *
 * @note
 * The function can also be used to invoke bitwise XNOR operation of @a data and @a bit_mask .
 * Input the complement of @a bit_mask (~@a bit_mask ) to invoke bitwise XNOR operation.
 * For details of bitwise XNOR operation, see, e.g.,@n
 * http://en.wikipedia.org/wiki/Truth_table
 *
 * @param[in] bit_mask A bit mask. The bit operation is invoked
 * between this value and the array, @a data.
 * @param[in] num_data The number of elements in the arrays, @a data,
 * @a edit_mask, and @a result.
 * @param[in] data An input array of size, @a num_data. The bit operation
 * is invoked between this array and @a bit_mask.@n
 * must-be-aligned
 * @param[in] edit_mask A boolean mask array of size, @a num_data. The bit operation
 * is skipped for the elements with the value, false.@n
 * must-be-aligned
 * @param[out] result The output array of size, @a num_data. It stores the result
 * of the bit operation between @a bit_mask and @a data. The bit operation is skipped
 * and the value in array, @a data, is adopted for the elements where corresponding
 * elements in @a edit_mask is false. The pointer of @a out is allowed to be equal to
 * that of @a in (@a result == @a data), indicating in-place operation.@n
 * must-be-aligned
 * @return status code
 * @~japanese
 * @brief ビットマスクと一次元配列の排他論理和ビット演算(XOR)を実行する。
 * @details 配列の@a i- 番目の要素に対して次の算を行い、出力@a result を返す:
 * @code
 * result [i] = edit_mask[i] ? (data[i] ^ bit_maks) : data[i]
 * @endcode
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。
 *
 * @note
 * この関数は、@a data と@a bit_mask の間の同値ビット演算(XNOR)にも使用できる。
 * 同値ビット演算を実行するときは、@a bit_mask をビット反転させたもの
 * ( ~@a bit_mask )を関数の入力として与える。同値ビット演算についての詳細は、例えば次のページを参照:@n
 * http://ja.wikipedia.org/wiki/%E8%AB%96%E7%90%86%E6%BC%94%E7%AE%97
 *
 * @param[in] bit_mask ビットマスク
 * @param[in] num_data @a data, @a edit_mask 及び@a result の要素の数。
 * @param[in] data 入力配列。要素数は@a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] edit_mask データのマスク。要素数は@a num_data でなければならない。
 * この値が true だと、対応する入力配列@a data とビットマスク@a bit_maks の排他論理和ビット演算を実行する。
 * この値が false だと、その要素のビット演算は行われず、対応する入力配列@a data の要素がそのまま出力となる。
 * @n must-be-aligned
 * @param[out] result 結果の格納先。要素数は@a num_data でなければならない。インプレースな変換を許す(@a result == @a data)。
 * @n must-be-aligned
 * @return 終了ステータス
 *@~
 * MT-safe
 *
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8Xor)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/]);
/**
 * @copybrief sakura_OperateBitsUint8Xor
 * @copydetails sakura_OperateBitsUint8Xor
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32Xor)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/]);

/**
 * @~japanese
 * @brief 補間方法を定義するための列挙型。
 * @~english
 * @brief Enumerations to define interpolation types.
 */
typedef enum {
	/**
	 * @~japanese
	 * @brief 最近接補間法
	 * @~english
	 * @brief Nearest interpolation
	 */LIBSAKURA_SYMBOL(InterpolationMethod_kNearest), /**
	 * @~japanese
	 * @brief 線形補間法
	 * @~english
	 * @brief Linear interpolation
	 */LIBSAKURA_SYMBOL(InterpolationMethod_kLinear), /**
	 * @~japanese
	 * @brief 多項式補間法
	 * @~english
	 * @brief Polynomial interpolation
	 */LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial), /**
	 * @~japanese
	 * @brief スプライン補間法（三次自然スプライン補間）
	 * @~english
	 * @brief Spline interpolation (Natural cubic spline)
	 */LIBSAKURA_SYMBOL(InterpolationMethod_kSpline), /**
	 * @~japanese
	 * @brief 実装されている補間法の個数
	 * @~english
	 * @brief Number of interpolation methods implemented
	 */LIBSAKURA_SYMBOL(InterpolationMethod_kNumMethod)
}LIBSAKURA_SYMBOL(InterpolationMethod);

/**
 * @~japanese
 * @brief 1次元の補間を行う。
 * @details 長さ @a num_base の1次元配列 @a base_position と 長さ @a num_base x @a num_array の1次元
 * 配列 @a base_data で定義される数値データ列をもとにして1次元の補間を行う。
 * @a base_data に @a num_array 個のデータを1次元配列として連結し、一括で補間することができる。
 * @a base_data の各要素が正常値か不正な値かをブール値のマスク @a base_mask で指定することができる。
 * マスク値がtrueならば正常値、falseなら不正な値として扱われ、不正な値は補間処理からは除外される。
 *
 * 補間によって値を得たい点の位置のリストを長さ　@a num_interpolated の配列 @a interpolate_position に
 * 渡すと、補間結果が長さ @a num_interpolated x @a num_array の配列 @a interpolated_data に格納される。
 * 外挿は行わない（データ点が片側にしかない場合にはそのデータ点の値が出力配列 @a interpolated_data にセットされる）。
 * @a interpolated_data に対応する出力マスク配列が @a interpolated_mask である。マスク値がfalseであるような
 * 要素に対応する @a interpolated_data の値は不正であるとみなされる（信用してはいけない）。
 * 出力マスク値がfalseにセットされるケースとしては、必要な入力データが全て不正な値で補間が実行できない場合がある。
 *
 * 戻り値は終了ステータスである。正常終了の場合、
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink
 * を返す。
 * 引数に不正がある場合には
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * を返す。内部で利用するメモリの確保に失敗した場合は、
 * @link sakura_Status::sakura_Status_kNoMemory sakura_Status_kNoMemory @endlink を返す。
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * が返された場合、
 * 考えられる原因は以下の二つである。
 *     - @a interpolation_method が正しくない
 *     - 引数に渡した配列がアラインされていない
 *
 * また、原因不明のエラーでは
 * @link sakura_Status::sakura_Status_kUnknownError sakura_Status_kUnknownError @endlink
 * を返す。
 *
 * @pre @a base_position および @a interpolate_position は昇順または降順にソートされていなければ
 * ならない。また、@a base_position の要素には重複があってはならない。
 *
 * @par sakura_InterpolateXAxisFloat() と sakura_InterpolateYAxisFloat() の違いについて
 * sakura_InterpolateXAxisFloat() と sakura_InterpolateYAxisFloat() の違いは、@a base_data
 * と@a interpolated_data のメモリレイアウトである。たとえば、@a base_data に対して前者が想定する
 * メモリレイアウトは[num_array][num_base]であり、@a base_data に複数のデータを連結する場合には、
 * @verbatim    (データ0の全要素), (データ1の全要素), ... @endverbatim
 * という順序でデータを格納する。
 * 一方後者では、メモリレイアウトとして[num_base][num_array]を
 * 仮定する。すなわち、@a base_data に複数のデータを連結する場合には、
 * @verbatim    (全データの第0要素), (全データの第1要素), ... @endverbatim
 * という順序でデータを格納する。@a interpolated_data のメモリレイアウトもこれに倣う。
 *
 * @par 昇順の場合と降順の場合の速度の違いについて:
 * @a base_position または @a interpolate_position が降順にソートされている場合、
 * 内部では配列要素をコピーして昇順に並べ替えた上で補間を行う。そのため、降順の場合は
 * 昇順よりも処理が遅くなる。
 *
 * @par 多項式補間の動作について:
 * @a polynomial_order はあくまで最大次数を規定するものであり、その次数で必ず
 * 補間が行われるとは限らない。たとえば、@a polynomial_order が2（2次多項式による補間）
 * で@a num_base が2の場合、実際には2点を通る1次多項式が一意に決まるため、2次多項式に
 * よる補間ではなく1次多項式による補間（線形補間）が行われる。
 * @par
 * @a polynomial_order に0を指定した場合、最近接補間が行われる。
 *
 * @par
 * @param[in] interpolation_method 補間方法
 * @param[in] polynomial_order 多項式補間法の場合の最大次数。
 * 実際に適用される次数は、@a num_base との兼ね合いで決まる。
 * @param[in] num_base 補間のためのデータ点の数。
 * @param[in] base_position 補間のための各データの位置。
 * 要素数は@a num_base でなければならない。
 * @a base_position は昇順または降順にソートされていなければならない。
 * must-be-aligned
 * @param[in] num_array 同時に渡すデータ列の数。
 * @param[in] base_data 補間のためのデータ列。
 * 要素数は@a num_base × @a num_array でなければならない。
 * must-be-aligned
 * @param[in] base_mask @a base_data に対応するマスク。
 * 要素数は@a num_base × @a num_array でなければならない。
 * マスク値がfalseである要素に対応する @a base_data の値は補間に使われない。
 * must-be-aligned
 * @param[in] num_interpolated 補間したいデータ点の数。
 * @param[in] interpolate_position 補間したいデータ点のx座標。
 * 要素数は@a num_interpolated でなければならない。
 * @a interpolate_position は昇順または降順にソートされていなければならない。
 * must-be-aligned
 * @param[out] interpolated_data 補間結果。
 * 要素数は@a num_interpolated × @a num_array でなければならない。
 * must-be-aligned
 * @param[out] interpolated_mask 補間結果のマスク。falseの要素に対応する
 * @a interpolated_data の値は不正である（信用してはいけない）。
 * 要素数は@a num_interpolated × @a num_array でなければならない。
 * must-be-aligned
 * @return 終了ステータス。
 *
 * @~english
 * @brief Perform one-dimensional interpolation.
 * @details
 * @param[in] interpolation_method interpolation method.
 * @param[in] polynomial_order maximum polynomial order for polynomial interpolation.
 * Actual order will be determined by a balance
 * between @a polynomial_order and @a num_base.
 * @param[in] num_base number of elements for data points.
 * @param[in] base_position position of data points. Its length must be @a num_base.
 * It must be sorted either ascending or descending.
 * @param[in] num_array number of arrays given in @a base_data.
 * @param[in] base_data value of data points. Its length must be @a num_base times @a num_array.
 * @param[in] num_interpolated number of elements for points that wants to get
 * interpolated value.
 * @param[in] interpolate_position x-coordinate of points that wants to get interpolated
 * value. Its length must be @a num_interpolated.
 * @param[out] interpolated_data storage for interpolation result. Its length must be
 * @a num_interpolated times @a num_array.
 * @return status code.
 *
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InterpolateXAxisFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		double const base_position[/*num_base*/], size_t num_array,
		float const base_data[/*num_base*num_array*/],
		bool const base_mask[/*num_base*num_array*/], size_t num_interpolated,
		double const interpolate_position[/*num_interpolated*/],
		float interpolated_data[/*num_interpolated*num_array*/],
		bool interpolated_mask[/*num_interpolated*num_array*/]);

/**
 * @copybrief sakura_IntepolateXAxisFloat
 * @copydetails sakura_InterpolateXAxisFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InterpolateYAxisFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		double const base_position[/*num_base*/], size_t num_array,
		float const base_data[/*num_base*num_array*/],
		bool const base_mask[/*num_base*num_array*/], size_t num_interpolated,
		double const interpolate_position[/*num_interpolated*/],
		float interpolated_data[/*num_interpolated*num_array*/],
		bool interpolated_mask[/*num_interpolated*num_array*/]);

/**
 * @~japanese
 * @brief position switch calibrationを実行する。
 * @details
 * position switch calibrationを実行する。
 * 具体的には、
 * @verbatim result = scaling_factor * (target - reference) / reference @endverbatim
 * を実行する。これは、position switch観測の温度較正
 * @verbatim calibrated = Tsys * (ON - OFF) / OFF @endverbatim
 * に相当する処理である。ただし、Tsys はシステム雑音温度、ONおよびOFFはそれぞれ
 * on sourceおよびoff sourceの生データである。
 * @n
 * @n
 * 戻り値は終了ステータスである。正常終了の場合、
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink
 * を返す。
 * 引数に不正がある場合には
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * を返す。
 * @n
 * @n
 * インプレースな処理がを許す。すなわち、@a result は@a target もしくは@a reference と同じ配列を渡してよい。
 * その場合、@a target もしくは@a reference は上書きされる。
 * @pre
 * @a num_scaling_factor は1もしくは@a num_data と等しくなければならない。1は周波数に依存しない、または周波数方向に
 * 平均されたシステム雑音温度のみが与えられた場合、@a num_data はシステム雑音温度が周波数に依存する場合に相当する。
 *
 * @param[in] num_scaling_factor @a scaling_factor の要素数。
 * 1または@a num_data のいずれかでなければならない。
 * @param[in] scaling_factor スケーリング因子。システム雑音温度に相当する。
 * 要素数は @a num_scaling_factor でなければならない。
 * must-be-aligned
 * @param[in] num_data データの要素数。
 * @param[in] target ターゲットデータ。on sourceデータに相当する。
 * 要素数は@a num_data でなければならない。
 * must-be-aligned
 * @param[in] reference 参照データ。off sourceデータに相当する。
 * 要素数は@a num_data でなければならない。
 * must-be-aligned
 * @param[out] result 計算結果。較正済みデータに相当する。
 * @a target または@a reference と同じ配列を渡してもよい。
 * 要素数は@a num_data でなければならない。
 * must-be-aligned
 *
 * @return 終了ステータス。
 *
 * @~english
 * @brief Apply position switch calibration.
 * @param[in] num_scaling_factor
 * @param[in] scaling_factor
 * @param[in] num_data
 * @param[in] target
 * @param[in] reference
 * @param[out] result
 *
 * @return status code.
 *
 * @~
 * MT-safe
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ApplyPositionSwitchCalibration)(
		size_t num_scaling_factor,
		float const scaling_factor[/*num_scaling_factor*/], size_t num_data,
		float const target[/*num_data*/], float const reference[/*num_data*/],
		float result[/*num_data*/]);

/**
 * @~japanese
 * @brief コンボリューションに使うカーネルタイプを列挙
 * @~english
 * @brief Enumerations to define kernel types for convolution.
 */
typedef enum {
	/**
	 * @brief Gaussian
	 */LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian), /**
	 * @brief Boxcar
	 */LIBSAKURA_SYMBOL(Convolve1DKernelType_kBoxcar), /**
	 * @brief Hanning
	 */LIBSAKURA_SYMBOL(Convolve1DKernelType_kHanning), /**
	 * @brief Hamming
	 */LIBSAKURA_SYMBOL(Convolve1DKernelType_kHamming), /**
	 * @brief Number of kernel type
	 */LIBSAKURA_SYMBOL(Convolve1DKernelType_kNumType)
}LIBSAKURA_SYMBOL(Convolve1DKernelType);
/**
 * @brief Context struct for convolution
 */
struct LIBSAKURA_SYMBOL(Convolve1DContext);
/**
 * @~japanese
 * @brief コンボリューションに必要なコンテキストを作成する。
 * @details
 * 戻り値は終了ステータスである。正常終了の場合、
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink
 * を返す。
 * 失敗した場合には
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * @link sakura_Status::sakura_Status_kNoMemory sakura_Status_kNoMemory @endlink
 * を返す。
 * @par
 * @param[in] num_data データの要素数。0 < num_data <= INT_MAX
 * @param[in] kernel_type カーネルタイプ
 * Gaussian,BoxCar,Hanning,Hammingを選択可能。各カーネルごとにコンボリューションの結果は異なる。
 * @param[in] kernel_width カーネルの幅. Gaussianカーネルの場合、kernal_widthは半値全幅（FWHM）と解釈される。0 < kernel_width
 * @param[in] use_fft コンボリューションの演算のためにFFTを行うか否かのフラグ。カーネルタイプには無間係。true=行う。false=行わない。
 * FFTを行う場合：
 * 畳み込み定理に従いFFTを利用した演算を行う。
 * 具体的には実数の入力データに対しFFTを行ってできた複素数配列と、事前に作った実数のカーネル配列に対しFFTを行って
 * できた複素数配列とを掛け合せ一つの複素数配列を得る。それを逆FFTし、実数配列である出力データを得る。
 * FFTを行わない場合：
 * 実数の入力データに対して実数のカーネルを用いて演算を行う。
 * @param[out] context コンボリューションのための情報を格納しているコンテキスト. Convolution1Dでの使用後にsakura_DestroyConvolve1DContext
 * により解放されなければならない。終了ステータスが
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink の場合はcontextには
 * コンボリューションに必要な情報が格納されている。終了ステータスが@link sakura_Status::sakura_Status_kNoMemory sakura_Status_kNoMemory @endlink
 * または@link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlinkの場合、contextの値はnullptrである。
 *
 * @return 終了ステータス。
 * @~english
 * @brief create context for convolution
 * @details
 * @param[in] num_data number of data. @num_data must
 * be positive.  0 < num_data < INT32_MAX
 * @param[in] kernel_type type of kernel(Gaussian,BoxCar,Hanning,Hamming).Each kernel can yield different convolution results.
 * @kernel_type is defined as enum.
 * @param[in] kernel_width kernel width. In case of Gaussian kernel, kernel_width will be interpreted as FWHM. 0 < kernel_width
 * @param[in] use_fft true means using FFT, false means not using FFT when convolution is performed. And Independent of the type of kernel.
 * If using FFT, FFT applied kernel which is already included context
 * by CreateConvolve1DContext is multiplied with input data
 * by complex-complex multiplication and then the multiplied complex
 * array is created. Finally inverse FFT is applied against it
 * and then real output data will be created.
 * If not using FFT, it is performed against real input data by real kernel
 * @param[out] context context. It has to be destroyed by sakura_DestroyConvolve1DContext after use by Convolution1D.
 * @return status code.
 * @~
 * MT-unsafe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
		size_t num_data, LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		size_t kernel_width, bool use_fft,
		struct LIBSAKURA_SYMBOL(Convolve1DContext) **context)
				LIBSAKURA_WARN_UNUSED_RESULT;
/**
 * @~japanese
 * @brief コンボリューションを行う。
 * @details sakura_CreateConvolve1DContextで設定した条件に従い、入力データに対してカーネルによるコンボリューションを行う。
 * @n
 * 戻り値は終了ステータスである。正常終了の場合、
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink
 * を返す。
 * 失敗した場合には
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * @link sakura_Status::sakura_Status_kNoMemory sakura_Status_kNoMemory @endlink
 * を返す。
 * @param[in] context コンテキスト
 * @param[in] num_data データの要素数。0 < num_data <= INT_MAX
 * @param[in] input_data 入力データ
 * 配列の長さは @a num_data と同じ。
 * @n must-be-aligned
 * @param[out] output_data 出力データ
 * 配列の長さは @a num_data と同じ。
 * @n must-be-aligned
 * @return 終了ステータス。
 * @~english
 * @brief convolution is performed
 * @details it is performed according to setting of sakura_CreateConvolve1DContext
 * @param[in] context context
 * and @a num_data, @a input_real_array
 * @param[in] num_data number of data
 * @param[in] input_data input data
 * Its length equals to channel number
 * @n must-be-aligned
 * @param[out] output_data
 * Its length equals to channel number
 * @n must-be-aligned
 * @return status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Convolve1DFloat)(
		struct LIBSAKURA_SYMBOL(Convolve1DContext) const *context,
		size_t num_data, float const input_data[/*num_data*/],
		float output_data[/*num_data*/]);
/**
 * @~japanese
 * @brief コンボリューションのために生成したコンテキストを破棄する。
 * @details
 * @param[in] context コンテキスト.
 * @return 終了ステータス。
 * @~english
 * @brief Destroy context
 * @details
 * @param[in] context context.
 * @return status code.
 * @~
 * MT-unsafe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(
		struct LIBSAKURA_SYMBOL(Convolve1DContext) *context);

/**
 * @~japanese
 * @brief 最小二乗フィットを解くための連立方程式の係数値を計算する。
 * @details
 * ( @a num_data ) 個の離散的な点で与えられたデータ yi ( 1 <= i <= @a num_data ) を、各々が同じく ( @a num_data ) 個の離散的な点で与えられる ( @a num_model_bases ) 個の基底関数 ai, bi, ..., ni の線型結合 (A * ai + B * bi + ... + N * ni) で最小二乗フィットし、基底関数の係数値 A, B, C, ... を求めることを考える。この時、これらの数の間には以下のような連立方程式(正規方程式)が成り立つ。
 *
 * @image html GetCoefficientsForLeastSquareFitting.png
 *
 * ここで、総和の記号は、マスクされていない全てのデータについて和を取ることを表す。この関数は、上の連立方程式の係数値、即ち、左辺の行列成分と右辺のベクトル成分を計算する。
 * @par
 * @param[in] num_data 配列 @a data 、 @a mask 、及び、モデルを構成する各基底関数の離散的データ点の要素数。正の数でなければならない。
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] mask 入力データに対するマスク情報。要素数は @a num_mask でなければならない。値がfalseの要素に対応する入力データはフィッティングに用いられない。
 * @n must-be-aligned
 * @param[in] num_model_bases モデルを構成する基底関数の数。正で、且つ、 @a num_data 以下の数でなければならない。
 * @param[in] basis_data モデルを構成する全ての基底関数の離散的な値を格納する１次元配列。関数に対するループはデータに対するループより内側になる。即ち、 @a m 番目のモデル関数の @a n 番目のデータ点の値は、 @a basis_data [ @a num_data * ( @a n -1) + ( @a m -1)]に格納されなければならない。配列の長さは( @a num_model_bases * @a num_data )でなければならない。
 * @n must-be-aligned
 * @param[out] lsq_matrix 求める連立方程式の左辺側の行列成分を格納する１次元配列。この行列は対称行列である。列に対するループは行のループより内側になる。即ち、 @a m 行 @a n 列目の成分値は、 @a lsq_matrix [ @a num_model_bases * ( @a m -1) + ( @a n -1)]に格納される。配列の長さは( @a num_model_bases * @a num_model_bases )となる。
 * @n must-be-aligned
 * @param[out] lsq_vector 求める連立方程式の右辺側のベクトル成分を格納する配列。配列の長さは @a num_model_bases となる。
 * @n must-be-aligned
 * @return 終了ステータス。正常終了時は Status_kOK、パラメータが上記の条件を満さない場合は Status_kInvalidArgument、マスクされていないデータの数が求める連立方程式の本数に満たない場合は Status_kNG、それ以外の例外が内部で発生した場合には Status_kUnknownError となる。
 * @~english
 * @brief Compute coefficients of simultaneous equations used for Least-Square fitting.
 * @details
 * Suppose fitting ( @a num_data ) discrete data points yi with a linear
 * combination of ( @a num_model_bases ) bases (ai, bi, ..., ni), which are
 * also given as ( @a num_data ) discrete points. Assuming the best-fit model
 * is given as (A * ai + B * bi + ... + N * ni), where (A, B, C, ...) are
 * the coefficients to be solved, these values are connected via the following
 * simultaneous equations known as normal equation:
 *
 * @image html GetCoefficientsForLeastSquareFitting.png
 *
 * Note that the summation means all the data points except masked ones
 * are to be added.
 * This function computes the coefficients of the above simultaneous equations, i.e.,
 * the elements of the matrix at the left side and of the vector at the right side.
 * @par
 * @param[in] num_data the number of elements in the arrays @a data
 * and @a mask, and also the number of elements in each model data
 * (i.e., discrete values of basis function) consisting the entire model.
 * it must be a positive number.
 * @param[in] data input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask input mask data with length of @a num_mask .
 * @n must-be-aligned
 * @param[in] num_model_bases number of model basis functions. it must be a
 * positive number, also it must be equal to or less than @a num_data .
 * @param[in] basis_data a 1D array containing values of all basis functions
 * concatenated. loop for basis index must be inside of that for data index,
 * i.e., the @a n -th data of the @a m -th model should be stored at
 * @a basis_data [ @a num_data * @a (n-1) + @a (m-1) ]. its length must be
 * equal to ( @a num_model_bases * @a num_data ).
 * @n must-be-aligned
 * @param[out] lsq_matrix a 1D array containing the values of a matrix
 * at the left side of simultaneous equations for least-square fitting.
 * its length should therefore be equal to ( @a num_model_bases * @a num_model_bases ).
 * loop for columns comes inside that for rows, i.e., the value at the
 * @a m -th row and @a n -th column is stored at @a out [ @a
 * num_model_bases * ( @a m -1) + ( @a n -1)], though @a out is actually
 * symmetric.
 * @n must-be-aligned
 * @param[out] lsq_vector the values of a vector at the right side of
 * simultaneous equations for least-square fitting. its length should be
 * equal to @a num_model_bases.
 * @n must-be-aligned
 * @return status code. Status_kOK if finished successfully,
 * Status_kInvalidArgument in case parameters does not meet the above
 * criteria, Status_kNG in case the number of unmasked data (for which
 * mask is @a false ) is less than the number of simultaneous equations,
 * and Status_kUnknownError in case other exceptions
 * emitted internally.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetLeastSquareFittingCoefficientsDouble)(
		size_t const num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], size_t const num_model_bases,
		double const basis_data[/*num_model_bases*num_data*/],
		double lsq_matrix[/*num_model_bases*num_model_bases*/],
		double lsq_vector[/*num_model_bases*/]) LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 最小二乗フィットを解くための連立方程式の係数値を更新する。
 * @details
 * ( @a num_data ) 個の離散的な点で与えられたデータ yi ( 1 <= i <= @a num_data ) を、各々が同じく ( @a num_data ) 個の離散的な点で与えられる ( @a num_model_bases ) 個の基底関数 ai, bi, ..., ni の線型結合 (A * ai + B * bi + ... + N * ni) で最小二乗フィットし、基底関数の係数値 A, B, C, ... を求めることを考える。この時、これらの数の間には以下のような連立方程式(正規方程式)が成り立つ。
 *
 * @image html GetCoefficientsForLeastSquareFitting.png
 *
 * ここで、総和の記号は、マスクされていない全てのデータについて和を取ることを表す。この関数は、先にGetLeastSquareFittingCoefficients()によって求められた、上の連立方程式の係数値を更新するために用いる。即ち、データ点のうち幾つかを除外して連立方程式を計算し直す際に、先に計算した各成分から除外するデータ点に対応する値を差し引く。除外するデータ数が少ない（前回の計算に用いられたデータ数の半分未満）場合は、一から計算し直すよりも高速である。
 * @par
 * @param[in] num_data 配列 @a data 、及び、モデルを構成する各基底関数の離散的データ点の要素数。正の数でなければならない。
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] num_exclude_indices 先に係数値を計算した状態と比較して、除外したいデータ点の個数をセットする。0以上、且つ、 @a num_data 以下の数でなければならない。
 * @param[in] exclude_indices 除外したいデータ点のインデックス ( @a basis_data の行の添字(0始まり)) を列挙した配列。要素数は @a num_exclude_indices 。
 * @n must-be-aligned
 * @param[in] num_model_bases モデルを構成する基底関数の数。正で、且つ、 @a num_data 以下の数でなければならない。
 * @param[in] basis_data モデルを構成する全ての基底関数の離散的な値を格納する１次元配列。関数に対するループはデータに対するループより内側になる。即ち、 @a m 番目のモデル関数の @a n 番目のデータ点の値は、 @a basis_data [ @a num_data * ( @a n -1) + ( @a m -1)]に格納されなければならない。配列の長さは( @a num_model_bases * @a num_data )でなければならない。
 * @n must-be-aligned
 * @param[in,out] lsq_matrix 更新するべき連立方程式の左辺側の行列成分を格納する１次元配列。この行列は対称行列である。列に対するループは行のループより内側になる。即ち、 @a m 行 @a n 列目の成分値は、 @a lsq_matrix [ @a num_model_bases * ( @a m -1) + ( @a n -1)]に格納される。要素数は必ず ( @a num_model_bases * @a num_model_bases ) でなければならない。
 * @n must-be-aligned
 * @param[in,out] lsq_vector 更新するべき連立方程式の右辺側のベクトル成分を格納する配列。要素数は必ず @a num_model_bases でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス。正常終了時は Status_kOK、パラメータが上記の条件を満さない場合は Status_kInvalidArgument、それ以外の例外が内部で発生した場合には Status_kUnknownError となる。
 * @par 注意:
 * この関数を使うにあたっては、残りのデータ数が連立方程式を解くための必要最低限 ( @a num_model_bases ) を割り込んだり、同じデータ点を重複して差し引いたりすることのないよう、ユーザー自身が気を付ける必要がある。
 * @~english
 * @brief Update coefficients of simultaneous equations used for Least-Square fitting.
 * @details
 * Suppose fitting ( @a num_mask ) discrete data points yi with a linear
 * combination of ( @a num_model_bases ) bases (ai, bi, ..., ni), which are
 * also given as ( @a num_mask ) discrete points. Assuming the best-fit model
 * is given as (A * ai + B * bi + ... + N * ni), where (A, B, C, ...) are
 * the coefficients to be solved, these values are connected via the following
 * simultaneous equations known as normal equation:
 *
 * @image html GetCoefficientsForLeastSquareFitting.png
 *
 * Note that the summation means all the data points except masked ones
 * are to be added.
 * This function updates the coefficients of the above simultaneous equations
 * by subtracting values corresponding to data points which have been used in
 * the previous calculation but not this time. this is faster than newly
 * calculating coefficients if the number of points to be excluded this time
 * is less than half of those previously used.
 * @par
 * @param[in] num_data the number of elements in the arrays @a data and the
 * number of elements in each model data (i.e., discrete values of basis
 * function) consisting the entire model.
 * it must be a positive number.
 * @param[in] data input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] num_exclude_indices the number of data points to be excluded
 * this time. the range of allowed value is between 0 and @a num_data .
 * @param[in] exclude_indices an array containing indices of data points
 * (the row index of @a basis_data ) to be excluded this time. the indices
 * must be stored in the first @a num_exclude_indices elements. its length
 * should be @a num_exclude_indices .
 * @n must-be-aligned
 * @param[in] num_model_bases number of model basis functions. it must be a
 * positive number, also it must be equal to or less than @a num_data .
 * @param[in] in set a 1D array containing the vector components previously
 * calculated.
 * @n must-be-aligned
 * @param[in] basis_data a 1D array containing values of all basis functions
 * concatenated. loop for basis index must be inside of that for data index,
 * i.e., the @a n -th data of the @a m -th model should be stored at
 * @a basis_data [ @a num_data * @a (n-1) + @a (m-1) ]. its length must be
 * equal to ( @a num_model_bases * @a num_data ).
 * @n must-be-aligned
 * @param[in,out] lsq_matrix a 1D array containing the values of a matrix
 * at the left side of simultaneous equations for least-square fitting.
 * its length should therefore be equal to ( @a num_model_bases * @a num_model_bases ).
 * loop for columns comes inside that for rows, i.e., the value at the
 * @a m -th row and @a n -th column is stored at @a out [ @a
 * num_model_bases * ( @a m -1) + ( @a n -1)], though @a out is actually
 * symmetric.
 * @n must-be-aligned
 * @param[in,out] lsq_vector the values of a vector at the right side of
 * simultaneous equations for least-square fitting. its length should be
 * equal to @a num_model_bases.
 * @n must-be-aligned
 * @return status code. Status_kOK if finished successfully,
 * Status_kInvalidArgument in case parameters does not meet the above
 * criteria, and Status_kUnknownError in case other exceptions emitted
 * internally.
 * @par Caution:
 * users must be careful in using this function about which and how many
 * data are to be excluded not to fall into destructive cases that the
 * number of used data becomes less than @a num_model_bases or not to
 * exclude the same data duplicatedly.
 * @~
 * MT-safe
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UpdateLeastSquareFittingCoefficientsDouble)(
		size_t const num_data, float const data[/*num_data*/],
		size_t const num_exclude_indices,
		size_t const exclude_indices[/*num_data*/],
		size_t const num_model_bases,
		double const basis_data[/*num_model_bases*num_data*/],
		double lsq_matrix[/*num_model_bases*num_model_bases*/],
		double lsq_vector[/*num_model_bases*/]) LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 連立方程式をLU分解によって解く。
 * @details
 * 連立方程式 A x = y を解き、ベクトル x の成分を求めることを考える。ここで、A は @a num_equations 行 @a num_equations 列の正方行列であり、x, y は長さ @a num_equations のベクトルである。本関数は与えられた A 及び y に対し、A の LU 分解を行うことによって x の値を求める。
 * @par
 * @param[in] num_equations 連立方程式の本数。
 * @param[in] in_matrix 連立方程式の左辺の行列 (上の方程式の A に対応する) の成分を格納する１次元配列。列に対するループは行のループより内側でなければならない。即ち、 @a m 行 @a n 列目の成分値は、 @a in_matrix [ @a num_equations * ( @a m -1) + ( @a n -1)]に格納されなければならない。配列の長さは( @a num_equations * @a num_equations )となる。
 * @n must-be-aligned
 * @param[in] in_vector 連立方程式の右辺のベクトル (上の方程式の y に対応する) の成分を格納する配列。配列の長さは @a num_equations でなければならない。
 * @n must-be-aligned
 * @param[out] out 連立方程式の解 (上の方程式の x に対応する) を格納する配列。配列の長さは @a num_equations でなければならない。 @a out を指すポインタは @a in_vector と同じでもよい。
 * @n must-be-aligned
 * @return 終了ステータス。
 * @~english
 * @brief Solve simultaneous equations via LU decomposition.
 * @details
 * Suppose solving simultaneous equations A x = y to derive x, where A
 * is a square matrix of @a num_equations rows and columns, and x and y
 * are vectors with length of @a num_equations . Given A and y values,
 * this function computes x values using LU decomposition of A.
 * @par
 * @param[in] num_equations number of equations.
 * @param[in] in_matrix a 1D array containing values of the matrix A in
 * the left side of the above simultaneous equations. loop for columns
 * comes inside that for rows, i.e., the value at the @a m -th row and
 * @a n -th column is stored at
 * @a in_matrix [ @a num_equations * ( @a m -1) + ( @a n -1)].
 * its length must be (@a num_equations * @a num_equations).
 * @n must-be-aligned
 * @param[in] in_vector a 1D array containing values of the vector y in
 * the right side of the above simultaneous equations. its length must be
 * @a num_equations .
 * @n must-be-aligned
 * @param[out] out the solution (x in the above equations). its length
 * must be @a num_equations . the pointer of @a out can be identical with
 * that of @a in_vector .
 * @n must-be-aligned
 * @return status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLUDouble)(
		size_t num_equations,
		double const in_matrix[/*num_equations*num_equations*/],
		double const in_vector[/*num_equations*/],
		double out[/*num_equations*/]) LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief ベースラインフィッティング固有のエラーコードを格納する列挙型。
 * @~english
 * @brief Enumerations to define baseline-specific error code.
 */
typedef enum {
	/**
	 * @~japanese
	 * @brief 成功または正常
	 * @~english
	 * @brief OK
	 */LIBSAKURA_SYMBOL(BaselineStatus_kOK) = 0, /**
	 * @~japanese
	 * @brief 失敗または異常
	 * @~english
	 * @brief NG
	 */LIBSAKURA_SYMBOL(BaselineStatus_kNG) = 1, /**
	 * @~japanese
	 * @brief データ数が不足のため、ベースラインフィッティングを実行できない
	 * @~english
	 * @brief not enough data for baseline fitting
	 */LIBSAKURA_SYMBOL(BaselineStatus_kNotEnoughData) = 2, /**
	 * @~japanese
	 * @brief 実装されているエラーコードの個数
	 * @~english
	 * @brief Number of error codes implemented
	 */LIBSAKURA_SYMBOL(BaselineStatus_kNumStatus)
}LIBSAKURA_SYMBOL(BaselineStatus);

/**
 * @~japanese
 * @brief ベースラインの関数形を格納する列挙型。
 * @~english
 * @brief Enumerations to define baseline type.
 */
typedef enum {
	/**
	 * @~japanese
	 * @brief 羃多項式
	 * @~english
	 * @brief Polynomial
	 */LIBSAKURA_SYMBOL(BaselineType_kPolynomial), /**
	 * @~japanese
	 * @brief チェビシェフ多項式
	 * @~english
	 * @brief Chebyshev Polynomial
	 */LIBSAKURA_SYMBOL(BaselineType_kChebyshev), /**
	 * @~japanese
	 * @brief 三次自然スプライン
	 * @~english
	 * @brief Cubic Spline
	 */LIBSAKURA_SYMBOL(BaselineType_kCubicSpline), /**
	 * @~japanese
	 * @brief 三角多項式
	 * @~english
	 * @brief Sinusoids
	 */LIBSAKURA_SYMBOL(BaselineType_kSinusoid), /**
	 * @~japanese
	 * @brief 実装されている関数形の個数
	 * @~english
	 * @brief Number of baseline functions implemented
	 */LIBSAKURA_SYMBOL(BaselineType_kNumType)
}LIBSAKURA_SYMBOL(BaselineType);

/**
 * @~japanese
 * @brief ベースラインフィッティングに用いるモデルデータと関連する情報を格納する構造体。
 * @~english
 * @brief Context struct for baseline fitting
 */
struct LIBSAKURA_SYMBOL(BaselineContext);

/**
 * @~japanese
 * @brief ベースラインモデル情報を格納するオブジェクトを生成する。
 * @details
 * @param[in] baseline_type ベースラインを表現する関数形。
 * @param[in] order モデルのパラメータ。多項式(poly, chebyshev)では次数、スプラインでは分割数、三角関数では最大の波数。スプラインの場合は正値でなければならない。それ以外のモデルではゼロも許される。このパラメータに基いて構成されるベースラインモデルの基底関数の数はそれぞれ、 @a order+1 （多項式）、 @a order+3 （三次自然スプライン）、 @a order*2+1 （三角関数）となるが、これがデータ点の数 @a num_data より大きな数になってはならない。
 * @param[in] num_data フィットするデータ点の数。
 * @param[out] context ベースラインモデルに関する情報を格納する構造体。
 * @n must-be-aligned
 * @return 終了ステータス。
 * @~english
 * @brief Create an object containing baseline model data.
 * @details
 * @param[in] baseline_type type of basis function.
 * @param[in] order parameter for the specified function.
 * actually it is the order (for polynomial and chebyshev),
 * or number of subsections (for cubic spline), or maximum
 * wave number (for sinusoid). must be positive for cubic
 * spline, while other models accept zero value. the number
 * of model bases, which is @a order+1 for polynomials or
 * @a order+3 for cubic spline or @a order*2+1 for sinusoids,
 * must not exceed @a num_data.
 * @param[in] num_data number of data to fit baseline.
 * @param[out] context an object containing baseline model data.
 * @return status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateBaselineContext)(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order,
		size_t const num_data, LIBSAKURA_SYMBOL(BaselineContext) **context)
				LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief ベースラインモデル情報を格納するオブジェクトを破棄する。
 * @details
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @return 終了ステータス。
 * @~english
 * @brief Destroy an object containing baseline model data.
 * @details
 * @param[in] context an object containing baseline model data.
 * @return status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyBaselineContext)(
LIBSAKURA_SYMBOL(BaselineContext) *context);

/**
 * @~japanese
 * @brief 与えられたデータに対して、同じく与えられたモデル基底関数の線型結合で表されるもののうち最も良く合うものを最小二乗フィットにより求める。
 * @details
 * @param[in] num_data 配列 @a data 、 @a mask 、 @a out 、及び、モデルを構成する各基底関数の離散的データ点の要素数。
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] mask 入力データに対するマスク情報。要素数は @a num_data でなければならない。値がfalseの要素に対応する入力データはフィッティングに用いられない。
 * @n must-be-aligned
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @param[out] out 出力される配列。要素数は @a num_data でなければならない。 @a out を指すポインタは @a data と同じでもよい。
 * @n must-be-aligned
 * @param[out] baseline_status ベースライン固有のエラーコード。
 * @return 終了ステータス。
 * @~english
 * @brief Compute the best-fit model by least-square fitting.
 * @details
 * @param[in] num_data the number of elements in the arrays @a data,
 * @a mask, @a out , and also the number of elements in each model data
 * (i.e., discrete values of basis function) consisting the total model.
 * @param[in] data the input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask the input mask data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] context an object containing baseline model data.
 * @param[out] out the best-fit model with length of @a num_data . the
 * pointer of @a out can be identical with that of @a data .
 * @n must-be-aligned
 * @param[out] baseline_status baseline-specific error code.
 * @return status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitBaselineFloat)(
		size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		float out[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status)
				LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 入力データに対して、与えられたモデル基底関数の線型結合で表されるもののうち最も良く合うものを最小二乗フィットにより求め、差し引く。
 * @details
 * @param[in] num_data 配列 @a data 、 @a mask 、 @a final_mask 、 @a out の要素数。
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] mask 入力データに対するマスク情報。要素数は @a num_data でなければならない。
 * 値がfalseの要素に対応する入力データはフィッティングに用いられない。
 * @n must-be-aligned
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @param[in] clip_threshold_sigma クリッピングの閾値。単位はσ。正値でなければならない。
 * @param[in] num_fitting_max フィッティングを(再帰的に)行う最大回数。
 * 値nが与えられた場合、最初のフィッティング＆差し引きを行った後、
 * 残差データのσを計算し、残差がその値の± @a clip_threshold_sigma
 * 倍を越えるものを除外して再度フィッティング＆差し引きを行うという操作を最大(n-1)回繰り返す。
 * デフォルト値は1、即ち、フィッティング＆差し引きは１回のみ行われ、クリッピングは行わない。
 * もし 0 が渡された場合は、自動的に 1 に変更される。
 * @param[in] get_residual trueの場合、入力データからフィットの結果を差し引いたものを出力として返す。falseの場合は、フィットの結果を出力として返す。
 * @param[out] final_mask 再帰的クリッピングを経た後の最終的なマスク情報。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[out] out 出力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[out] baseline_status ベースライン固有のエラーコード。
 * @return 終了ステータス。
 * @~english
 * @brief Recursively fit a baseline and subtract it from input spectrum.
 * @details
 * @param[in] num_data the number of elements in the arrays @a data,
 * @a mask, @a final_mask, and @a out.
 * @param[in] data the input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask the input mask data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] context an object containing baseline model data.
 * @param[in] clip_threshold_sigma the threshold of clipping in unit of
 * sigma. must be positive.
 * @param[in] num_fitting_max the maximum of total number of times
 * baseline fitting is performed recursively. In case n is given, after
 * the first baseline fitting, subsequent clipping and baseline fitting
 * based on the updated mask are executed (n-1) times at maximum.
 * The default is 1 (i.e., baseline fitting done just once and no
 * clipping applied). In case zero is given, @a num_fitting_max will be
 * automatically changed to 1.
 * @param[in] get_residual set the output to be (input - best-fit) if true,
 * or the best-fit value if false.
 * @param[out] final_mask the final mask data after recursive clipping
 * procedure. its length must be @a num_data .
 * @n must-be-aligned
 * @param[out] out the output data. its length must be @a num_data .
 * @n must-be-aligned
 * @param[out] baseline_status baseline-specific error code.
 * @return status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineFloat)(
		size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		float clip_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual,
		bool final_mask[/*num_data*/], float out[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status)
				LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 入力データに対して、与えられたモデル基底関数の線型結合で表されるもののうち最も良く合うものを最小二乗フィットにより求め、係数を返す。
 * @details
 * @param[in] num_data 配列 @a data 、 @a mask 、and @a final_mask
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] mask 入力データに対するマスク情報。要素数は @a num_data でなければならない。
 * 値がfalseの要素に対応する入力データはフィッティングに用いられない。
 * @n must-be-aligned
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @param[in] clip_threshold_sigma クリッピングの閾値。単位はσ。正値でなければならない。
 * @param[in] num_fitting_max フィッティングを(再帰的に)行う最大回数。
 * 値nが与えられた場合、最初のフィッティング＆差し引きを行った後、
 * 残差データのσを計算し、残差がその値の± @a clip_threshold_sigma
 * 倍を越えるものを除外して再度フィッティング＆差し引きを行うという操作を最大(n-1)回繰り返す。
 * デフォルト値は1、即ち、フィッティング＆差し引きは１回のみ行われ、クリッピングは行わない。
 * もし 0 が渡された場合は、自動的に 1 に変更される。
 * @param[out] final_mask 再帰的クリッピングを経た後の最終的なマスク情報。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[out] coeff 出力データ。要素数は @a num_coeff でなければならない。
 * @n must-be-aligned
 * @param[out] baseline_status ベースライン固有のエラーコード。
 * @return 終了ステータス。
 * @~english
 * @brief Extraction of the coefficients of the polynomial fit.
 * @details
 * @param[in] num_data the number of elements in the arrays @a data,
 * @a mask, and @a final_mask.
 * @param[in] data the input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask the input mask data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] clip_threshold_sigma the threshold of clipping in unit
 * of sigma. must be positive.
 * @param[in] num_fitting_max the maximum of total number of times
 * baseline fitting is performed recursively. In case n is given, after
 * the first baseline fitting, subsequent clipping and baseline fitting
 * based on the updated mask are executed (n-1) times at maximum.
 * The default is 1 (i.e., baseline fitting done just once and no
 * clipping applied). In case zero is given, @a num_fitting_max will be
 * automatically changed to 1.
 * @param[in] num_coeff the number of elements in the arrays @a coeff.
 * @param[out] coeff the coefficients of the polynomial fit. its length must be @a num_coeff.
 * @n must-be-aligned
 * @param[out] final_mask the final mask data after recursive clipping
 * procedure. its length must be @a num_data .
 * @n must-be-aligned
 * @param[out] baseline_status baseline-specific error code.
 * @return status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficentsFloat)(
		size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineContext) const *context,
		float clip_threshold_sigma, uint16_t num_fitting_max, size_t num_coeff,
		double coeff[/*num_coeff*/],
		bool final_mask[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status)
				LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 入力データに多項式ベースラインをフィットし差し引く。
 * @details
 * @param[in] num_data 配列 @a data 、 @a mask 、 @a final_mask 、 @a out の要素数。
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] mask 入力データに対するマスク情報。要素数は @a num_data でなければならない。値がfalseの要素に対応する入力データはフィッティングに用いられない。
 * @n must-be-aligned
 * @param[in] order 多項式モデルの次数。 @a num_data-1 以下の値でなければならない。
 * @param[in] clip_threshold_sigma クリッピングの閾値。単位はσ。正値でなければならない。
 * @param[in] num_fitting_max フィッティングを(再帰的に)行う最大回数。
 * 値nが与えられた場合、最初のフィッティング＆差し引きを行った後、残差データのσを計算し、
 * 残差がその値の± @a clip_threshold_sigma
 * 倍を越えるものを除外して再度フィッティング＆差し引きを行うという操作を最大(n-1)回繰り返す。
 * デフォルト値は1、即ち、フィッティング＆差し引きは１回のみ行われ、クリッピングは行わない。
 * もし 0 が渡された場合は、自動的に 1 に変更される。
 * @param[in] get_residual trueの場合、入力データからフィットの結果を差し引いたものを出力として返す。falseの場合は、フィットの結果を出力として返す。
 * @param[out] final_mask 再帰的クリッピングを経た後の最終的なマスク情報。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[out] out 出力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[out] baseline_status ベースライン固有のエラーコード。
 * @~english
 * @brief Fit a baseline and subtract it from input data.
 * @details
 * @param[in] num_data the number of elements in the arrays @a data,
 * @a mask, @a final_mask, and @a out.
 * @param[in] data the input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask the input mask data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] order order of polynomial model. must be equal or smaller
 * than @a num_data-1 .
 * @param[in] clip_threshold_sigma the threshold of clipping in unit
 * of sigma. must be positive.
 * @param[in] num_fitting_max the maximum of total number of times
 * baseline fitting is performed recursively. In case n is given, after
 * the first baseline fitting, subsequent clipping and baseline fitting
 * based on the updated mask are executed (n-1) times at maximum.
 * The default is 1 (i.e., baseline fitting done just once and no
 * clipping applied). In case zero is given, @a num_fitting_max will be
 * automatically changed to 1.
 * @param[in] get_residual set the output to be (input - best-fit) if true,
 * or the best-fit value if false.
 * @param[out] final_mask the final mask data after recursive clipping
 * procedure. its length must be @a num_data .
 * @n must-be-aligned
 * @param[out] out the output data. its length must be @a num_data .
 * @n must-be-aligned
 * @param[out] baseline_status baseline-specific error code.
 * @return status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselinePolynomialFloat)(
		size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], uint16_t order,
		float clip_threshold_sigma, uint16_t num_fitting_max,
		bool get_residual,
		bool final_mask[/*num_data*/], float out[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status)
				LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 入力データに対して、与えられたベースラインモデルと係数からベースラインを求め、差し引く。
 * @details
 * @param[in] num_data 配列 @a data 、 @a out の要素数。
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] num_coeff 配列 @a coeff の要素数。
 * @param[in] coeff 最小二乗フィットにより得られたベストフィット係数。要素数は @a num_coeff でなければならない。
 * @n must-be-aligned
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @n must-be-aligned
 * @param[out] out 出力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス
 * @~english
 * @brief subtract baseline from input spectrum. baseline is calculated by baseline model and coefficients.
 * @details
 * @param[in] num_data the number of elements in the arrays @a data and @a out.
 * @param[in] data the input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] context an object containing baseline model data.
 * @param[in] num_coeff the number of elements in the arrays @a coeff
 * @param[in] coeff best fit coefficients obtained by least-square fitting. The input data with length of @a num_coeff.
 * @n must-be-aligned
 * @param[out] out the output data. its length must be @a num_data .
 * @n must-be-aligned
 * @param[out] baseline_status baseline-specific error code.
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineUsingCoefficientsFloat)(
		size_t num_data, float const data[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t num_coeff,
		double const coeff[/*num_data*/], float out[/*num_data*/])
				LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 最小二乗フィットにより得られたベストフィット係数を返す
 * @details
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @return 終了ステータス
 * @~english
 * @brief return best fit coefficients obtained by least-square fitting.
 * @details
 * @param[in] context an object containing baseline model data.
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetNumberOfCoefficients)(
LIBSAKURA_SYMBOL(BaselineContext) const *context, size_t * num_coeff)
		LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~english
 * @brief Copy elements in the @a src matrix into the @a dst matrix with flipping elements to reorder
 * as some FFT library expects.
 *
 * @details
 * When you provide @a innerMostUntouched = false, @a elements = {3, 4} and @a src = {
 * @code
     1,   2,   3,
 *   4,   5,   6,
 *   7,   8,   9,
 *  10,  11,  12,
 * @endcode
 * }, then you will get @a dst as below.
 * @code
 *   9,   7,   8,
 *  12,  10,  11,
 *   3,   1,   2,
 *   6,   4,   5,
 * @endcode
 *
 * @param[in] innerMostUntouched If true, the order of the inner most dimension is untouched.
 * @param[in] dims Dimensions of the matrix @a src and @a dst. In other words, a number of elements in @a elements.
 * @param[in] elements Numbers of elements of each dimension of @a src and @a dst with the inner-to-outer order.
 * @param[in] src	Source matrix.
 * @n must-be-aligned
 * @param[in] dst	Destination matrix.
 * @n must-be-aligned
 * @return status code.
 * @~
 * MT-safe
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipMatrixFloat)(
bool innerMostUntouched, size_t dims, size_t const elements[],
		float const src[], float dst[]);

/**
 * @~english
 * @brief Copy elements in the @a src matrix into the @a dst matrix with unflipping elements to the original order.
 *
 * @details
 * When you provide @a innerMostUntouched = false, @a elements = {3, 4} and @a src = {
 * @code
 *   9,   7,   8,
 *  12,  10,  11,
 *   3,   1,   2,
 *   6,   4,   5,
 * @endcode
 * }, then you will get @a dst as below.
 * @code
 *   1,   2,   3,
 *   4,   5,   6,
 *   7,   8,   9,
 *  10,  11,  12,
 * @endcode
 *
 * @param[in] innerMostUntouched If true, the order of the inner most dimension is untouched.
 * @param[in] dims Dimensions of the matrix @a src and @a dst. In other words, a number of elements in @a elements.
 * @param[in] elements Numbers of elements of each dimension of @a src and @a dst with the inner-to-outer order.
 * @param[in] src	Source matrix.
 * @n must-be-aligned
 * @param[in] dst	Destination matrix.
 * @n must-be-aligned
 * @return status code.
 * @~
 * MT-safe
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipMatrixFloat)(
bool innerMostUntouched, size_t dims, size_t const elements[],
		float const src[], float dst[]);

/**
 * @~english
 * @brief Same as @ref sakura_FlipMatrixFloat except the element type of the matrixes.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipMatrixDouble)(
bool innerMostUntouched, size_t dims, size_t const elements[],
		double const src[], double dst[]);

/**
 * @~english
 * @brief Same as @ref sakura_UnflipMatrixFloat except the element type of the matrixes.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipMatrixDouble)(
bool innerMostUntouched, size_t dims, size_t const elements[],
		double const src[], double dst[]);

/**
 * @~english
 * @brief Same as @ref sakura_FlipMatrixFloat except the element type of the matrixes.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipMatrixDouble2)(
bool innerMostUntouched, size_t dims, size_t const elements[],
		double const src[][2], double dst[][2]);

/**
 * @~english
 * @brief Same as @ref sakura_UnflipMatrixFloat except the element type of the matrixes.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipMatrixDouble2)(
bool innerMostUntouched, size_t dims, size_t const elements[],
		double const src[][2], double dst[][2]);

#ifdef __cplusplus
}
/* extern "C" */
#endif

#endif /* LIBSAKURA_LIBSAKURA_SAKURA_H_ */
