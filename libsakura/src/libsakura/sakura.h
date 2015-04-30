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

#define LIBSAKURA_NOEXCEPT /* noexcept */

#ifdef __cplusplus
extern "C" {

#if __cplusplus >= 201103L
# undef LIBSAKURA_NOEXCEPT
# define LIBSAKURA_NOEXCEPT noexcept
#endif

#endif

/**
 * @~english
 * @brief A result of function call.
 *
 * @~japanese
 * @brief 関数の呼び出し結果を表す
 *
 */
typedef enum {
	/**
	 * @~english
	 * @brief Success or normal end
	 *
	 * @~japanese
	 * @brief 成功または正常
	 */LIBSAKURA_SYMBOL
	(Status_kOK) = 0,
	/**
	 * @~english
	 * @brief Failure or abnormal end
	 *
	 * @~japanese
	 * @brief 失敗または異常
	 */LIBSAKURA_SYMBOL(Status_kNG) = 1,
	/**
	 * @~english
	 * @brief Illegal argument(s)
	 *
	 * This includes a violation of the must-be-aligned constraint.
	 *
	 * @~japanese
	 * @brief 引数が不正だった
	 *
	 * @a must-be-aligned 制約に違反している場合も含む。
	 */LIBSAKURA_SYMBOL(Status_kInvalidArgument) = 2,
	/**
	 * @~english
	 * @brief No memory
	 *
	 * @~japanese
	 * @brief メモリーが足りない
	 */LIBSAKURA_SYMBOL(Status_kNoMemory) = 3,
	/**
	 * @~english
	 * @brief Unknown error
	 *
	 * @~japanese
	 * @brief 原因不明のエラー
	 */LIBSAKURA_SYMBOL(Status_kUnknownError) = 99
}LIBSAKURA_SYMBOL (Status);

/**
 * @~english
 * @brief A type of the allocator function called by Sakura Library.
 *
 * Implementation of the function of this type must be reentrant.
 *
 * @note
 * It must return a valid pointer to a memory region of size 0 if 0 is passed as @a size parameter.
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
 * Implementation of the function of this type must be reentrant.
 *
 * @note
 * It must do nothing if NULL is passed as @a pointer parameter.
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
LIBSAKURA_SYMBOL(UserDeallocator) deallocator)
		LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

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
void LIBSAKURA_SYMBOL(CleanUp)() LIBSAKURA_NOEXCEPT;

/**
 * @~english
 * @brief Returns a current time.
 * Precision of the time depends on std::chrono::system_clock.
 * @return Current time in seconds since the Epoch.
 *
 * @~japanese
 * @brief 現在時刻(単位は秒)を返す
 *
 * 精度は std::chrono::system_clock 依存。
 * @return エポックからの経過時間(単位は秒)
 * @~
 * MT-safe
 */
double LIBSAKURA_SYMBOL(GetCurrentTime)() LIBSAKURA_NOEXCEPT;

/*
 * memory alignment(for SIMD)
 */
/**
 * @~english
 * @brief Checks if @a ptr points the aligned address Sakura Library requires.
 *
 * @param[in] ptr An address to be checked. NULL is allowed.
 * @return true if the address is aligned, otherwise false
 * @~japanese
 * @brief Sakuraライブラリが想定するアライメントに、@a ptr が合っているか調べる
 *
 * @param[in] ptr アラインされているか調べたいアドレス。NULL も受け付ける。
 * @return アラインされているなら true , そうでないなら false
 * @~
 * MT-safe
 */

bool LIBSAKURA_SYMBOL(IsAligned)(void const *ptr) LIBSAKURA_NOEXCEPT;

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
size_t LIBSAKURA_SYMBOL (GetAlignment)() LIBSAKURA_NOEXCEPT;

/**
 * @~english
 * @brief Returns an aligned address close to @a arena by adding 0 or minimum offset.
 *
 * It returns @a arena if @a arena is already aligned.
 *
 * @param[in] arena start address of a memory region
 * @param[in] size_of_arena size of the memory region pointed by @a arena
 * @param[in] size_required required size after alignment
 * @return aligned address if at least @a size_required bytes are available in @a arena after alignment,
 * otherwise NULL.
 *
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
		size_t size_required) LIBSAKURA_NOEXCEPT;
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
 *
 * @~english
 * @brief Returns an aligned address close to @a arena by adding 0 or minimum offset.
 *
 * It returns @a arena if @a arena is already aligned.
 *
 * @param[in] arena start address of an array
 * @param[in] elements_in_arena The number of elements in @a arena , not a size in bytes.
 * @param[in] elements_required required number of elements after alignment
 * @return aligned address if at least @a elements_required are available in @a arena after alignment,
 * otherwise NULL.
 *
 * @~
 * MT-safe
 */
float *LIBSAKURA_SYMBOL(AlignFloat)(size_t elements_in_arena, float *arena,
		size_t elements_required) LIBSAKURA_NOEXCEPT;

/**
 * @~japanese
 * @copydoc sakura_AlignFloat()
 *
 * @~english
 * @copydoc sakura_AlignFloat()
 *
 */
double *LIBSAKURA_SYMBOL(AlignDouble)(size_t elements_in_arena, double *arena,
		size_t elements_required) LIBSAKURA_NOEXCEPT;

/**
 * @~japanese
 * @brief @ref sakura_ComputeStatisticsFloat と @ref sakura_ComputeAccurateStatisticsFloat の結果を格納する構造体
 *
 * @~english
 * @brief A structure to which the result of @ref sakura_ComputeStatisticsFloat and @ref sakura_ComputeAccurateStatisticsFloat is stored.
 *
 * You can also figure out following statistics from the members of this struct if count > 0:
 *  - mean = sum / count
 *  - rms = sqrt(square_sum / count)
 *  - variance = abs(square_sum / count - mean * mean)
 *  - stddev = sqrt(variance)
 */
typedef struct {
	/**
	 * @~
	 * number of valid data
	 */
	size_t count;
	/**
	 * @~
	 * sum of valid data
	 */
	double sum;
	/**
	 * @~
	 * sum of squared valid data
	 */
	double square_sum;

	float min;
	/**
	 * @~
	 * max value of valid data. NaN if no valid data.
	 */
	float max;
	/**
	 * @~
	 * index for one of min value. -1 if there is no valid data.
	 */
	int index_of_min;
	/**
	 * @~
	 * index for one of max value. -1 if there is no valid data.
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
 * @param[in] num_data The number of elements in @a data and @a is_valid . @a num_data <= INT32_MAX
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
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_ComputeStatisticsFloat
 * @copydetails sakura_ComputeStatisticsFloat
 * @~
 * The result of this function is more accurate than that of @ref sakura_ComputeStatisticsFloat if
 * num_data is large. This function is slower than @ref sakura_ComputeStatisticsFloat .
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeAccurateStatisticsFloat)(
		size_t num_data, float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResultFloat) *result) LIBSAKURA_NOEXCEPT;

/**
 * @~japanese
 * @brief validな値のみを先頭に詰めて昇順にソートする
 *
 * @param[in] num_data @a data 及び@a is_valid の要素の数
 * @param[in] is_valid データのマスク。この値が false だと、
 * 対応する@a data の要素が無視される
 * @param[in,out] data ソート対象のデータ。In placeでソートするので、この配列内の順序は変更される。
 * 対応する@a is_valid がtrueの場合、InfやNaNであってはならない。
 * @param[out] new_num_data (validでないデータを除いた)ソートされた要素数( <= @a num_data ) の格納先
 * @return 終了ステータス
 *
 * @~english
 * @brief Sorts only valid data in ascending order.
 *
 * @param[in] num_data The number of elements in @a data and @a is_valid .
 * @param[in] is_valid Masks of @a data. If a value of element is false,
 * the corresponding element in @a data is ignored.
 * @param[in,out] data Data to be sorted. Since data is sorted in place, contents of this array are not preserved.
 * If corresponding element in @a is_valid is true, the element in @a data must not be Inf nor NaN.
 * @param[out] new_num_data The number of sorted elements that don't include invalid data( <= @a num_data ) is stored here.
 * @return status code
 *
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SortValidValuesDenselyFloat)(
		size_t num_data, bool const is_valid[], float data[],
		size_t *new_num_data) LIBSAKURA_NOEXCEPT;

 /**
 * @~english
 * @brief Computes median absolute deviation.
 *
 * This function applies abs(@a data[i] - median) for each element in @a data and
 * stores results to @a new_data. Then it sorts @a new_data in ascending order.
 * @a data must be sorted in advance using e.g. @ref sakura_SortValidValuesDenselyFloat .
 * The median is a center value if @a num_data is odd. Otherwise, the median is
 * an average of two center values.
 * The caller is responsible to take a median value from @a new_data as
 * a median absolute deviation.
 *
 *  mad = (@a new_data[@a num_data/2] + @a new_data[@a num_data/2 - @a num_data%2]) / 2 if @a num_data > 0
  *
 * @param[in] num_data	The number of elements in @a data and @a new_data .
 * @param[in] data	Each value in data must not be Inf nor NaN. @a data must be sorted.
 * <br/>must-be-aligned
 * @param[in] new_data	An array where the results are stored. @a new_data may points @a data.
 * <br/>must-be-aligned
 * @return status code
 * @~
 * MT-safe
  */
 LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeMedianAbsoluteDeviationFloat)(
 		size_t num_data, float const data[], float new_data[]) LIBSAKURA_NOEXCEPT;

 /**
 * @~english
 * @brief Grids data with convolution.
 *
 * This function plot @a value from @a start_spectrum to @a end_spectrum whose location is at (@a x , @a y)
 * onto @a grid as points which have width represented by @a convolution_table.
 * The results of 2 * @a support pixels from the edges of the grids are not reliable due to
 * implementation to gain speed.
 * Polarizations and channels are mapped by @a polarization_map and @a channel_map when gridding.
 * All floating-point values passed to this function must not be NaN nor +-Inf.
 *
 * @param[in] num_spectra	The number of spectra. It should be  0 <= @a start_spectrum <= @a end_spectrum <= @a num_spectra .
 * @param[in] start_spectrum	Starting index of spectrum to be processed.
 * @param[in] end_spectrum	An index next to the last index of spectrum to be processed.
 * @param[in] spectrum_mask	Masks to each spectrum. The number of elements must be @a num_spectra.
 * If a value of element is false, the corresponding spectrum is ignored.
 * <br/>must-be-aligned
 * @param[in] x	X position of the point projected to 2D plane of the grid.
 * The number of elements must be @a num_spectra . Each element must be INT32_MIN < @a x[i] < INT32_MAX .
 * <br/>must-be-aligned
 * @param[in] y	Y position of the point projected to 2D plane of the grid.
 * The number of elements must be @a num_spectra . Each element must be INT32_MIN < @a y[i] < INT32_MAX .
 * <br/>must-be-aligned
 * @param[in] support	A half width(the number of pixels from center to edge) of convolution kernel on the grid plane of size @a width x @a height. It must be 0 < @a support <= 256 .<br/>
 * Note that this function may cause stack overflow because it allocates memory of a size proportional to @a support * @a sampling on stack
 * when @a support * @a sampling is too large.
 * @param[in] sampling	A resolution of convolution kernel(samples/grid aka samples/pixel). It must be 0 < @a sampling <= INT32_MAX .
 * @param[in] num_polarizations	The number of polarizations. It must be 0 < @a num_polarizations <= INT32_MAX .
 * @param[in] polarization_map	The number of elements must be @a num_polarizations. Each element must be in range [0,@a num_polarizations_for_grid).
 * <br/>must-be-aligned
 * @param[in] num_channels	The number of channels. It must be 0 < @a num_channels <= INT32_MAX .
 * @param[in] channel_map	The number of elements must be @a num_channels. Each element must be in range [0,@a num_channels_for_grid).
 * <br/>must-be-aligned
 * @param[in] mask	Masks to each @a value . Its memory layout should be [@a num_spectra][@a num_polarizations][@a num_channels].
 * If a value of element is false, the corresponding combination of the spectrum, polarization and channel is ignored.
 * <br/>must-be-aligned
 * @param[in] value	Values to be gridded. Its memory layout should be [@a num_spectra][@a num_polarizations][@a num_channels].
 * <br/>must-be-aligned
 * @param[in] weight Weights for @a value. Its memory layout should be [@a num_spectra][@a num_channels].
 * <br/>must-be-aligned
 * @param[in] weight_only	True if you want to get a grid of weight itself rather than production of @a value and @a weight. Otherwise false.
 * @param[in] num_convolution_table	The number of elements of @a convolution_table. It should be ceil(sqrt(2.)*(@a support+1)*@a sampling) <= @a num_convolution_table <= INT32_MAX / 32 .
 * @param[in] convolution_table	An array which represents convolution kernel. The number of elements must be @a num_convolution_table. The first element corresponds to center of the point.
 * <br/>must-be-aligned
 * @param[in] num_polarizations_for_grid	The number of polarizations on the grid. It should be 0 < @a num_polarizations_for_grid <= INT32_MAX .
 * @param[in] num_channels_for_grid	The number of channels on the grid. It should be 0 < @a num_channels_for_grid <= INT32_MAX .
 * @param[in] width	Width of the grid . It should be 0 < @a width <= INT32_MAX .
 * @param[in] height Height of the grid . It should be 0 < @a height <= INT32_MAX .
 * @param[out] weight_sum	Sum of weights. Its memory layout should be [@a num_polarizations_for_grid][@a num_channels_for_grid].
 * <br/>must-be-aligned
 * @param[out] weight_of_grid	Weight for each grid. Its memory layout should be [@a height][@a width][@a num_polarizations_for_grid][@a num_channels_for_grid].
 * <br/>must-be-aligned
 * @param[out] grid	The resulting grid. Its memory layout should be [@a height][@a width][@a num_polarizations_for_grid][@a num_channels_for_grid].
 * <br/>must-be-aligned
 *
 * @~japanese
 * @brief 畳み込みしながらグリッドする
 *
 * @a start_spectrum から @a end_spectrum までの、@a x , @a y 座標上の @a value の値を
 * @a convolution_table で表される広がりを持った点として、@a grid 上にプロットする。
 * エッジから 2 * @a support ピクセル分の値は信頼できない(理由:処理の高速化のため)。
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
 * @param[in] sampling 畳み込みカーネルの精度(/pixel)。範囲は、0 < @a sampling <= INT32_MAX
 * @param[in] num_polarizations 範囲は、0 < @a num_polarizations <= INT32_MAX
 * @param[in] polarization_map	要素数は、@a num_polarizations 。各要素の値は、[0,@a num_polarizations_for_grid)でなければならない。<br/>must-be-aligned
 * @param[in] num_channels 範囲は、0 < @a num_channels <= INT32_MAX
 * @param[in] channel_map	要素数は、@a num_channels 。各要素の値は、[0,@a num_channels_for_grid)でなければならない。<br/>must-be-aligned
 * @param[in] mask	要素のレイアウトは、[@a num_spectra][@a num_polarizations][@a num_channels]。falseの場合は、該当するスペクトル、偏波、チ ャネルの組み合わせのデータは無視される。<br/>must-be-aligned
 * @param[in] value	要素のレイアウトは、[@a num_spectra][@a num_polarizations][@a num_channels]。グリッディングする値。<br/>must-be-aligned
 * @param[in] weight 要素のレイアウトは、[@a num_spectra][@a num_channels]。重み。<br/>must-be-aligned
 * @param[in] weight_only @a value に重みを掛けたものではなく、重み自体をグリッディングする場合は、true。
 * @param[in] num_convolution_table @a convolution_table の要素数。 範囲は、ceil(sqrt(2.)*(@a support+1)*@a sampling) <= num_convolution_table <= INT32_MAX / 32
 * @param[in] convolution_table	要素数は、@a num_convolution_table 。畳み込みに使用する重みカーブ。要素0が中心を表す。<br/>must-be-aligned
 * @param[in] num_polarizations_for_grid 範囲は、0 < @a num_polarizations_for_grid <= INT32_MAX
 * @param[in] num_channels_for_grid 範囲は、0 < @a num_channels_for_grid <= INT32_MAX
 * @param[in] width 範囲は、0 < @a width <= INT32_MAX
 * @param[in] height 範囲は、0 < @a height <= INT32_MAX
 * @param[out] weight_sum	要素のレイアウトは、[@a num_polarizations_for_grid][@a num_channels_for_grid]。重みの合計。<br/>must-be-aligned
 * @param[out] weight_of_grid	要素のレイアウトは、[@a height][@a width][@a num_polarizations_for_grid][@a num_channels_for_grid]。グリッドの重み。<br/>must-be-aligned
 * @param[out] grid	要素のレイアウトは、[@a height][@a width][@a num_polarizations_for_grid][@a num_channels_for_grid]。グリッディング結果。<br/>must-be-aligned
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
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

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
 * @param[in] num_data The number of elements in the arrays, @a data and @a result .
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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesInclusiveFloat)(
		size_t num_data, float const data[/*num_data*/], size_t num_condition,
		float const lower_bounds[/*num_condition*/],
		float const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_SetTrueIfInRangesInclusiveFloat
 * @copydetails sakura_SetTrueIfInRangesInclusiveFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesInclusiveInt)(
		size_t num_data, int const data[/*num_data*/], size_t num_condition,
		int const lower_bounds[/*num_condition*/],
		int const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveFloat)(
		size_t num_data, float const data[/*num_data*/], size_t num_condition,
		float const lower_bounds[/*num_condition*/],
		float const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_SetTrueIfInRangesExclusiveFloat
 * @copydetails sakura_SetTrueIfInRangesExclusiveFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfInRangesExclusiveInt)(
		size_t num_data, int const data[/*num_data*/], size_t num_condition,
		int const lower_bounds[/*num_condition*/],
		int const upper_bounds[/*num_condition*/],
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_SetTrueIfGreaterThanFloat
 * @copydetails sakura_SetTrueIfGreaterThanFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_SetTrueIfGreaterThanOrEqualsFloat
 * @copydetails sakura_SetTrueIfGreaterThanOrEqualsFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfGreaterThanOrEqualsInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_SetTrueIfLessThanFloat
 * @copydetails sakura_SetTrueIfLessThanFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsFloat)(
		size_t num_data, float const data[/*num_data*/], float threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_SetTrueIfLessThanOrEqualsFloat
 * @copydetails sakura_SetTrueIfLessThanOrEqualsFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetTrueIfLessThanOrEqualsInt)(
		size_t num_data, int const data[/*num_data*/], int threshold,
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SetFalseIfNanOrInfFloat)(
		size_t num_data, float const data[/*num_data*/],
		bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

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
		uint8_t const data[/*num_data*/], bool result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_Uint8ToBool
 * @copydetails sakura_Uint8ToBool
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Uint32ToBool)(size_t num_data,
		uint32_t const data[/*num_data*/], bool result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

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
bool const data[/*num_data*/], bool result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

/**
 * @~english
 * @brief Invoke bit operation AND between a bit mask and an array.
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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseAndUint8)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_OperateBitwiseAndUint8
 * @copydetails sakura_OperateBitwiseAndUint8
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseAndUint32)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint8)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_OperateBitwiseConverseNonImplicationUint8
 * @copydetails sakura_OperateBitwiseConverseNonImplicationUint8
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseConverseNonImplicationUint32)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint8)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_OperateBitwiseImplicationUint8
 * @copydetails sakura_OperateBitwiseImplicationUint8
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseImplicationUint32)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseNotUint8)(
		size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_OperateBitwiseNotUint8
 * @copydetails sakura_OperateBitwiseNotUint8
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseNotUint32)(
		size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @~english
 * @brief Invoke bit operation OR between a bit mask and an array.
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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseOrUint8)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_OperateBitwiseOrUint8
 * @copydetails sakura_OperateBitwiseOrUint8
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseOrUint32)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseXorUint8)(
		uint8_t bit_mask, size_t num_data, uint8_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint8_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;
/**
 * @copybrief sakura_OperateBitwiseXorUint8
 * @copydetails sakura_OperateBitwiseXorUint8
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitwiseXorUint32)(
		uint32_t bit_mask, size_t num_data, uint32_t const data[/*num_data*/],
		bool const edit_mask[/*num_data*/], uint32_t result[/*num_data*/])
				LIBSAKURA_NOEXCEPT;

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
	 */LIBSAKURA_SYMBOL(InterpolationMethod_kNumElements)
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
 * @a polynomial_order は @a interpolation_method が
 * @link sakura_InterpolationMethod::sakura_InterpolationMethod_kPolynomial sakura_InterpolationMethod_kPolynomial @endlink
 * である場合のみ有効であり、その他の場合は無視される。
 * @param[in] num_base 補間のためのデータ点の数。@a num_base は正値でなければならない。
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
 * @details It performs one-dimensional interpolation based on two input arrays @a base_position
 * and @a base_data. Size of input array is @a num_base for @a base_position and @a num_interpolated
 * times @a num_array for @a base_data, where @a num_array is a number of arrays that are passed to the
 * function so that interpolation on multiple arrays can be performed simultaneously.
 * One can set boolean mask for @a base_data using @a base_mask, which has same array shape as
 * @a base_data. Mask value is true if data is valid while the value is false if the data is invalid.
 * Invalid data will be excluded from interpolation.
 *
 * List of locations where interpolated value is evaluated have to be specified via @a interpolate_position
 * whose length is @a num_interpolated. Interpolation result on each point specified by
 * @a interpolate_position is stored to @a interpolated_data. No extrapolation will be performed. Instead,
 * out of range points are filled by the value of nearest points.
 * Output boolean mask is @a interpolated_mask. Data will be invalid if corresponding mask is false, i.e.,
 * the value is just a nominal one but a result of actual interpolation. The mask will be false when
 * interpolation is skipped due to insufficient number of valid data elements.
 *
 * The function returns result status. In the successful run, returned value is
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink while appropriate
 * error status will be returned for failure:
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * for invalid input arguments,
 * @link sakura_Status::sakura_Status_kNoMemory sakura_Status_kNoMemory @endlink
 * for memory allocation error for internal variables, and
 * @link sakura_Status::sakura_Status_kUnknownError sakura_Status_kUnknownError @endlink
 * for other unknown error.
 * Possible reason for
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * is either of (1) invalid @a interpolation_method or (2) any of input/output array is not aligned.
 *
 * @pre @a base_position and @a interpolate_position must be sorted. Also, these arrays should not
 * have duplicate values.
 *
 * @par Difference between sakura_InterpolateXAxisFloat() and sakura_InterpolateYAxisFloat():
 * Difference between these two similar functions is a memory layout of @a base_data and
 * @a interpolated_data. The former assumes the layout like [num_array][num_base] so that
 * @a base_data should store the data in the following order,
 * @verbatim data0[0], data0[1], ..., data1[0], data1[1], ... @endverbatim
 * On the other hand, the latter requires [num_base][num_array], i.e,
 * @verbatim data0[0], data1[0], ..., data0[1], data1[1], ... @endverbatim
 * Result array, @a interpolated_data, follows the memory layout required for @a base_data.
 *
 * @par Impact of sort order on performance:
 * When input arrays, @a base_position and/or @a base_data, are sorted in descending order,
 * the arrays are internally reversed in ascending order and store them to working array.
 * Therefore, descending inputs may cause degradation of performance compared with ascending inputs.
 *
 * @par Note on polynomial interpolation:
 * Note that @a polynomial_order defines maximum order for polynomial interpolation.
 * In other words, it doesn't assure the interpolation to be specified order.
 * For example, suppose that @a polynomial_order is 2 and @a num_base is also 2.
 * In this case, effective polynomial order is 1 since we can obtain unique polynomial
 * with order 1, not 2, that passes through given two points.
 * Note also that @a polynomial_order 0 is equivalent to nearest interpolation.
 *
 * @par
 * @param[in] interpolation_method interpolation method.
 * @param[in] polynomial_order maximum polynomial order for polynomial interpolation.
 * Actual order will be determined by a balance
 * between @a polynomial_order and @a num_base.
 * This parameter is effective only when @a interpolation_method is
 * @link sakura_InterpolationMethod::sakura_InterpolationMethod_kPolynomial sakura_InterpolationMethod_kPolynomial @endlink.
 * In other interpolation methods, it is ignored.
 * @param[in] num_base number of elements for data points. Its value must be greater than 0.
 * @param[in] base_position position of data points. Its length must be @a num_base.
 * It must be sorted either ascending or descending.
 * must-be-aligned
 * @param[in] num_array number of arrays given in @a base_data.
 * @param[in] base_data value of data points. Its length must be @a num_base times @a num_array.
 * must-be-aligned
 * @param[in] base_mask boolean mask for data. Its length must be @a num_base times @a num_array.
 * False points will be excluded from the interpolation
 * must-be-aligned
 * @param[in] num_interpolated number of elements for points that wants to get
 * interpolated value.
 * @param[in] interpolated_position location of points that wants to get interpolated
 * value. Its length must be @a num_interpolated.
 * must-be-aligned
 * @param[out] interpolated_data storage for interpolation result. Its length must be
 * @a num_interpolated times @a num_array.
 * must-be-aligned
 * @param[out] interpolated_mask boolean mask for interpolation result. Its length must be
 * @a num_interpolated times @a num_array.
 * must-be-aligned
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
		double const interpolated_position[/*num_interpolated*/],
		float interpolated_data[/*num_interpolated*num_array*/],
		bool interpolated_mask[/*num_interpolated*num_array*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @copybrief sakura_IntepolateXAxisFloat
 * @copydetails sakura_InterpolateXAxisFloat
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(InterpolateYAxisFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_base,
		double const base_position[/*num_base*/], size_t num_array,
		float const base_data[/*num_base*num_array*/],
		bool const base_mask[/*num_base*num_array*/], size_t num_interpolated,
		double const interpolated_position[/*num_interpolated*/],
		float interpolated_data[/*num_interpolated*num_array*/],
		bool interpolated_mask[/*num_interpolated*num_array*/])
				LIBSAKURA_NOEXCEPT;

/**
 * @~
 * @brief Normalize data against reference value with scaling factor.
 * @details
 * Normalize the data against reference value with scaling factor. The function normalizes the @a data
 * by the assumption that
 * the value in @a reference is normalized to @a scaling_factor. Specifically, it will calculate,
 * @verbatim result = scaling_factor * (data - reference) / reference @endverbatim
 *
 * @n
 * The function returns result status. For successful run, return value will be
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink.
 * On the other hand, it will return
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * if any invalid values are passed to arguments.
 * @n
 * @n
 * The function allows in-place calculation, i.e., @a result array can be either @a data or
 * @a reference. In this case, @a data or @a reference will be overwritten.
 * @pre
 * @a num_scaling_factor should be 1 or equal to @a num_data. These values corresponds to the case of
 * frequency-independent (channel averaged) or frequency-dependent system temperature,
 * respectively.
 *
 * @param[in] num_scaling_factor number of elements of @a scaling_factor. Its value should be
 * 1 or equal to @a num_data. If @a num_scaling_factor is 1, @a scaling_factor[0] will be
 * applied to all array elements.
 * @param[in] scaling_factor scaling factor corresponding to system temperature. number of
 * elements must be @a num_scaling_factor.
 * must-be-aligned
 * @param[in] num_data number of data
 * @param[in] data data to be normalized. number of elements must be @a num_data.
 * must-be-aligned
 * @param[in] reference reference data. number of elements must be @a num_data.
 * must-be-aligned
 * @param[out] result resulting normalized data. one can give same array with either
 * @a data or @a reference for in-place calculation. number of elements must be @a num_data.
 * must-be-aligned
 *
 * @return status code.
 *
 * MT-safe
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(NormalizeDataAgainstReferenceFloat)(
		size_t num_scaling_factor,
		float const scaling_factor[/*num_scaling_factor*/], size_t num_data,
		float const data[/*num_data*/], float const reference[/*num_data*/],
		float result[/*num_data*/]) LIBSAKURA_NOEXCEPT;

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
	 */LIBSAKURA_SYMBOL(Convolve1DKernelType_kNumElements)
}LIBSAKURA_SYMBOL(Convolve1DKernelType);
/**
 * @brief Context struct for convolution
 */
struct LIBSAKURA_SYMBOL(Convolve1DContextFloat);
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
 * @param[in] num_data The number of data. @a num_data must
 * be positive.  0 < num_data < INT32_MAX
 * @param[in] kernel_type type of kernel(Gaussian,BoxCar,Hanning,Hamming).Each kernel can yield different convolution results.
 * @a kernel_type is defined as enum.
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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateConvolve1DContextFloat)(
		size_t num_data, LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		size_t kernel_width, bool use_fft,
		struct LIBSAKURA_SYMBOL(Convolve1DContextFloat) **context)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

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
 * @brief Convolution is performed.
 * @details Convolution operations are performed by shifting a kernel over the input data.
 * The kernel is created by the settings of sakura_CreateConvolve1DContext.
 * @param[in] context
 * Context.
 * @param[in] num_data
 * The number of elements in @a input_data and @a output_data. (0 < @a num_data <= INT_MAX)
 * @param[in] input_data
 * Input data. Its length equals to the number of elements in @a input_data.
 * @n must-be-aligned
 * @param[out] output_data
 * Output data. Its length equals to the number of elements in @a output_data.
 * @n must-be-aligned
 * @return status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Convolve1DFloat)(
		struct LIBSAKURA_SYMBOL(Convolve1DContextFloat) const *context,
		size_t num_data, float const input_data[/*num_data*/],
		float output_data[/*num_data*/]) LIBSAKURA_NOEXCEPT;
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
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyConvolve1DContextFloat)(
		struct LIBSAKURA_SYMBOL(Convolve1DContextFloat) *context)
				LIBSAKURA_NOEXCEPT;

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
 * @param[in] num_lsq_bases 実際に係数の計算に用いる基底関数の数。正で、且つ、@a num_model_basis以下の数でなければならない。
 * @param[out] lsq_matrix 求める連立方程式の左辺側の行列成分を格納する１次元配列。この行列は対称行列である。列に対するループは行のループより内側になる。即ち、 @a m 行 @a n 列目の成分値は、 @a lsq_matrix [ @a num_lsq_bases * ( @a m -1) + ( @a n -1)]に格納される。配列の長さは( @a num_lsq_bases * @a num_lsq_bases )となる。
 * @n must-be-aligned
 * @param[out] lsq_vector 求める連立方程式の右辺側のベクトル成分を格納する配列。配列の長さは @a num_lsq_bases となる。
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
 * @param[in] num_data The number of elements in the arrays @a data
 * and @a mask, and also the number of elements in each model data
 * (i.e., discrete values of basis function) consisting the entire model.
 * It must be a positive number.
 * @param[in] data Input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask Input mask data with length of @a num_mask .
 * @n must-be-aligned
 * @param[in] num_model_bases Number of model basis functions. It must be a
 * positive number, also it must be equal to or less than @a num_data .
 * @param[in] basis_data A 1D array containing values of all basis functions
 * concatenated. Loop for basis index must be inside of that for data index,
 * i.e., the @a n -th data of the @a m -th model should be stored at
 * @a basis_data [ @a num_data * @a (n-1) + @a (m-1) ]. Its length must be
 * equal to ( @a num_model_bases * @a num_data ).
 * @n must-be-aligned
 * @param[in] num_lsq_bases The number of model basis functions to be used
 * in actual fitting. It must be a positive number and must not exceed
 * @a num_model_bases .
 * @param[out] lsq_matrix A 1D array containing the values of a matrix
 * at the left side of simultaneous equations for least-square fitting.
 * Its length should therefore be equal to ( @a num_lsq_bases * @a num_lsq_bases ).
 * Loop for columns comes inside that for rows, i.e., the value at the
 * @a m -th row and @a n -th column is stored at @a out [ @a
 * num_lsq_bases * ( @a m -1) + ( @a n -1)], though @a out is actually
 * symmetric.
 * @n must-be-aligned
 * @param[out] lsq_vector The values of a vector at the right side of
 * simultaneous equations for least-square fitting. Its length should be
 * equal to @a num_model_bases .
 * @n must-be-aligned
 * @return Status code. sakura_Status_kOK if finished successfully,
 * sakura_Status_kInvalidArgument in case parameters does not meet the
 * above criteria, sakura_Status_kNG in case the number of unmasked data
 * (for which mask is @a false ) is less than the number of simultaneous
 * equations, and sakura_Status_kUnknownError in case other exceptions
 * emitted internally.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetLSQCoefficientsDouble)(
		size_t const num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], size_t const num_model_bases,
		double const basis_data[/*num_model_bases*num_data*/],
		size_t const num_lsq_bases,
		double lsq_matrix[/*num_lsq_bases*num_lsq_bases*/],
		double lsq_vector[/*num_lsq_bases*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 最小二乗フィットを解くための連立方程式の係数値を更新する。
 * @details
 * ( @a num_data ) 個の離散的な点で与えられたデータ yi ( 1 <= i <= @a num_data ) を、各々が同じく ( @a num_data ) 個の離散的な点で与えられる ( @a num_model_bases ) 個の基底関数 ai, bi, ..., ni の線型結合 (A * ai + B * bi + ... + N * ni) で最小二乗フィットし、基底関数の係数値 A, B, C, ... を求めることを考える。この時、これらの数の間には以下のような連立方程式(正規方程式)が成り立つ。
 *
 * @image html GetCoefficientsForLeastSquareFitting.png
 *
 * ここで、総和の記号は、マスクされていない全てのデータについて和を取ることを表す。この関数は、先にGetLSQCoefficients()によって求められた、上の連立方程式の係数値を更新するために用いる。即ち、データ点のうち幾つかを除外して連立方程式を計算し直す際に、先に計算した各成分から除外するデータ点に対応する値を差し引く。除外するデータ数が少ない（前回の計算に用いられたデータ数の半分未満）場合は、一から計算し直すよりも高速である。
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
 * @param[in] num_lsq_bases 実際に係数の計算に用いる基底関数の数。正で、且つ、@a num_model_basis以下の数でなければならない。
 * @param[in,out] lsq_matrix 更新するべき連立方程式の左辺側の行列成分を格納する１次元配列。この行列は対称行列である。列に対するループは行のループより内側になる。即ち、 @a m 行 @a n 列目の成分値は、 @a lsq_matrix [ @a num_lsq_bases * ( @a m -1) + ( @a n -1)]に格納される。要素数は必ず ( @a num_lsq_bases * @a num_lsq_bases ) でなければならない。
 * @n must-be-aligned
 * @param[in,out] lsq_vector 更新するべき連立方程式の右辺側のベクトル成分を格納する配列。要素数は必ず @a num_lsq_bases でなければならない。
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
 * @param[in] num_data The number of elements in the arrays @a data and the
 * number of elements in each model data (i.e., discrete values of basis
 * function) consisting the entire model.
 * It must be a positive number.
 * @param[in] data Input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] num_exclude_indices The number of data points to be excluded
 * this time. The range of allowed value is between 0 and @a num_data .
 * @param[in] exclude_indices An array containing indices of data points
 * (the row index of @a basis_data ) to be excluded this time. The indices
 * must be stored in the first @a num_exclude_indices elements. Its length
 * should be @a num_exclude_indices .
 * @n must-be-aligned
 * @param[in] num_model_bases Number of model basis functions. It must be a
 * positive number, also it must be equal to or less than @a num_data .
 * @param[in] basis_data A 1D array containing values of all basis functions
 * concatenated. Loop for basis index must be inside of that for data index,
 * i.e., the @a n -th data of the @a m -th model should be stored at
 * @a basis_data [ @a num_data * @a (n-1) + @a (m-1) ]. Its length must be
 * equal to ( @a num_model_bases * @a num_data ).
 * @n must-be-aligned
 * @param[in] num_lsq_bases The number of model basis functions to be used
 * in actual fitting. It must be a positive number and must not exceed
 * @a num_model_bases .
 * @param[in,out] lsq_matrix A 1D array containing the values of a matrix
 * at the left side of simultaneous equations for least-square fitting.
 * Its length should therefore be equal to ( @a num_lsq_bases * @a num_lsq_bases ).
 * Loop for columns comes inside that for rows, i.e., the value at the
 * @a m -th row and @a n -th column is stored at @a out [ @a
 * num_lsq_bases * ( @a m -1) + ( @a n -1)], though @a out is actually
 * symmetric.
 * @param[in,out] lsq_vector The values of a vector at the right side of
 * simultaneous equations for least-square fitting. Its length should be
 * equal to @a num_lsq_bases .
 * @n must-be-aligned
 * @return Status code. sakura_Status_kOK if finished successfully,
 * sakura_Status_kInvalidArgument in case parameters does not meet the above
 * criteria, and sakura_Status_kUnknownError in case other exceptions emitted
 * internally.
 * @par Caution:
 * Users must be careful in using this function about which and how many
 * data are to be excluded not to fall into destructive cases that the
 * number of used data becomes less than @a num_model_bases or not to
 * exclude the same data duplicatedly.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UpdateLSQCoefficientsDouble)(
		size_t const num_data, float const data[/*num_data*/],
		size_t const num_exclude_indices,
		size_t const exclude_indices[/*num_data*/],
		size_t const num_model_bases,
		double const basis_data[/*num_model_bases*num_data*/],
		size_t const num_lsq_bases,
		double lsq_matrix[/*num_lsq_bases*num_lsq_bases*/],
		double lsq_vector[/*num_lsq_bases*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

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
 * @param[in] num_equations Number of equations.
 * @param[in] in_matrix A 1D array containing values of the matrix A in
 * the left side of the above simultaneous equations. Loop for columns
 * comes inside that for rows, i.e., the value at the @a m -th row and
 * @a n -th column is stored at
 * @a in_matrix [ @a num_equations * ( @a m -1) + ( @a n -1)].
 * Its length must be (@a num_equations * @a num_equations).
 * @n must-be-aligned
 * @param[in] in_vector A 1D array containing values of the vector y in
 * the right side of the above simultaneous equations. Its length must be
 * @a num_equations .
 * @n must-be-aligned
 * @param[out] out The solution (x in the above equations). Its length
 * must be @a num_equations . The pointer of @a out can be identical with
 * that of @a in_vector .
 * @n must-be-aligned
 * @return Status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLUDouble)(
		size_t num_equations,
		double const in_matrix[/*num_equations*num_equations*/],
		double const in_vector[/*num_equations*/],
		double out[/*num_equations*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

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
	 * @brief Not enough data for baseline fitting
	 */LIBSAKURA_SYMBOL(BaselineStatus_kNotEnoughData) = 2, /**
	 * @~japanese
	 * @brief 実装されているエラーコードの個数
	 * @~english
	 * @brief Number of error codes implemented
	 */LIBSAKURA_SYMBOL(BaselineStatus_kNumElements)
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
	 */LIBSAKURA_SYMBOL(BaselineType_kNumElements)
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
 * @param[in] order モデルのパラメータ。多項式(poly, chebyshev)では最大の次数、スプラインでは最大の分割数、三角関数では最大の波数。スプラインの場合は正値でなければならない。それ以外のモデルではゼロも許される。このパラメータに基いて構成されるベースラインモデルの基底関数の数はそれぞれ、 @a order+1 （多項式）、 @a order+3 （三次自然スプライン）、 @a order*2+1 （三角関数）となるが、これがデータ点の数 @a num_data より大きな数になってはならない。
 * @param[in] num_data フィットするデータ点の数。
 * @param[out] context ベースラインモデルに関する情報を格納する構造体。
 * @n must-be-aligned
 * @return 終了ステータス。
 * @~english
 * @brief Create an object containing baseline model data.
 * @details
 * @param[in] baseline_type Type of basis function.
 * @param[in] order Parameter for the specified function.
 * It is the maximum order (for sakura_BaselineType_kPolynomial
 * and sakura_BaselineType_kChebyshev), or the maximum number of
 * spline pieces (for sakura_BaselineType_kCubicSpline), or the
 * maximum wave number (for sakura_BaselineType_kSinusoid). It
 * must be positive for sakura_BaselineType_kCubicSpline, while
 * other models accept zero value. The number of model bases,
 * which is @a order+1 for sakura_BaselineType_kPolynomial and
 * sakura_BaselineType_kChebyshev, or @a order+3 for
 * sakura_BaselineType_kCubicSpline or @a order*2+1 for
 * sakura_BaselineType_kSinusoid, must not exceed @a num_data.
 * @param[in] num_data Number of data to fit baseline.
 * @param[out] context An object containing baseline model data.
 * @return Status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateBaselineContext)(
LIBSAKURA_SYMBOL(BaselineType) const baseline_type, uint16_t const order,
		size_t const num_data,
		struct LIBSAKURA_SYMBOL(BaselineContext) **context)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief ベースラインモデル情報を格納するオブジェクトを破棄する。
 * @details
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @return 終了ステータス。
 * @~english
 * @brief Destroy an object containing baseline model data.
 * @details
 * @param[in] context An object containing baseline model data.
 * @return Status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyBaselineContext)(
		struct LIBSAKURA_SYMBOL(BaselineContext) *context) LIBSAKURA_NOEXCEPT;

/**
 * @~japanese
 * @brief 入力データに対して、与えられたモデル基底関数の線型結合で表されるもののうち最も良く合うものを最小二乗フィットにより求め、差し引く。
 * @details
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @param[in] order モデルのパラメータ。多項式(polynomial, chebyshev)では次数、三角関数では最大の波数。値は@a context を生成したときに指定した@a order より大きな数になってはならない。このパラメータに基いて構成されるベースラインモデルの基底関数の数はそれぞれ、 @a order+1 （多項式）、 @a order*2+1 （三角関数）となるが、これがデータ点の数 @a num_data より大きな数になってはならない。
 * @param[in] num_data 配列 @a data 、 @a mask 、 @a final_mask 、 @a out の要素数。
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] mask 入力データに対するマスク情報。要素数は @a num_data でなければならない。
 * 値がfalseの要素に対応する入力データはフィッティングに用いられない。
 * @n must-be-aligned
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
 * @param[in] context An object containing baseline model data.
 * @param[in] order Parameter for the specified function.
 * It is the order for sakura_BaselineType_kPolynomial and
 * sakura_BaselineType_kChebyshev, or the maximum wave number
 * for sakura_BaselineType_kSinusoid. The value should not exceed
 * the @a order specified in creation of @a context .
 * The number of model bases, which is @a order+1 for
 * sakura_BaselineType_kPolynomial and
 * sakura_BaselineType_kChebyshev, or @a order*2+1 for
 * sakura_BaselineType_kSinusoid, must not exceed @a num_data.
 * @param[in] num_data The number of elements in the arrays @a data,
 * @a mask, @a final_mask, and @a out.
 * @param[in] data The input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask The input mask data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] clip_threshold_sigma The threshold of clipping in unit of
 * sigma. must be positive.
 * @param[in] num_fitting_max The maximum of total number of times
 * baseline fitting is performed recursively. In case n is given, after
 * the first baseline fitting, subsequent clipping and baseline fitting
 * based on the updated mask are executed (n-1) times at maximum.
 * The default is 1 (i.e., baseline fitting done just once and no
 * clipping applied). In case zero is given, @a num_fitting_max will be
 * automatically changed to 1.
 * @param[in] get_residual Set the output to be (input - best-fit) if true,
 * or the best-fit value if false.
 * @param[out] final_mask The final mask data after recursive clipping
 * procedure. its length must be @a num_data .
 * @n must-be-aligned
 * @param[out] out The output data. Its length must be @a num_data .
 * @n must-be-aligned
 * @param[out] baseline_status Baseline-specific error code.
 * @return Status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineFloat)(
		struct LIBSAKURA_SYMBOL(BaselineContext) const *context,
		uint16_t const order, size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], float clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual,
		bool final_mask[/*num_data*/], float out[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 入力データに対して、与えられた個数の区分からなる3次スプライン曲線を最小二乗フィットし、差し引く。
 * @details
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @param[in] num_pieces スプライン曲線の区分の数。正値、かつ、( @a num_data / 4) 以下の数でなければならない。
 * @param[in] num_data 配列 @a data 、 @a mask 、 @a final_mask 、 @a out の要素数。
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] mask 入力データに対するマスク情報。要素数は @a num_data でなければならない。
 * 値がfalseの要素に対応する入力データはフィッティングに用いられない。
 * @n must-be-aligned
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
 * @brief Recursively fit a cubic spline baseline and subtract it from input spectrum.
 * @details
 * @param[in] context An object containing baseline model data.
 * @param[in] num_pieces Number of spline pieces. It must be positive
 * and also must not exceed ( @a num_data * 4 ).
 * @param[in] num_data The number of elements in the arrays @a data,
 * @a mask, @a final_mask, and @a out.
 * @param[in] data The input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask The input mask data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] clip_threshold_sigma The threshold of clipping in unit of
 * sigma. must be positive.
 * @param[in] num_fitting_max The maximum of total number of times
 * baseline fitting is performed recursively. In case n is given, after
 * the first baseline fitting, subsequent clipping and baseline fitting
 * based on the updated mask are executed (n-1) times at maximum.
 * The default is 1 (i.e., baseline fitting done just once and no
 * clipping applied). In case zero is given, @a num_fitting_max will be
 * automatically changed to 1.
 * @param[in] get_residual Set the output to be (input - best-fit) if true,
 * or the best-fit value if false.
 * @param[out] final_mask The final mask data after recursive clipping
 * procedure. its length must be @a num_data .
 * @n must-be-aligned
 * @param[out] out The output data. Its length must be @a num_data .
 * @n must-be-aligned
 * @param[out] baseline_status Baseline-specific error code.
 * @return Status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineFloat)(
		struct LIBSAKURA_SYMBOL(BaselineContext) const *context,
		size_t num_pieces, size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], float clip_threshold_sigma,
		uint16_t num_fitting_max, bool get_residual,
		bool final_mask[/*num_data*/], float out[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 入力データに対して、与えられたモデル基底関数の線型結合で表されるもののうち最も良く合うものを最小二乗フィットにより求め、係数を返す。
 * @details
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @param[in] num_data 配列 @a data 、 @a mask 、 @a final_mask の要素数。
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] mask 入力データに対するマスク情報。要素数は @a num_data でなければならない。
 * 値がfalseの要素に対応する入力データはフィッティングに用いられない。
 * @n must-be-aligned
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
 * @param[in] context An object containing baseline model data.
 * @param[in] num_data The number of elements in the arrays @a data,
 * @a mask, and @a final_mask.
 * @param[in] data The input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask The input mask data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] clip_threshold_sigma The threshold of clipping in unit
 * of sigma. must be positive.
 * @param[in] num_fitting_max The maximum of total number of times
 * baseline fitting is performed recursively. In case n is given, after
 * the first baseline fitting, subsequent clipping and baseline fitting
 * based on the updated mask are executed (n-1) times at maximum.
 * The default is 1 (i.e., baseline fitting done just once and no
 * clipping applied). In case zero is given, @a num_fitting_max will be
 * automatically changed to 1.
 * @param[in] num_coeff The number of elements in the arrays @a coeff.
 * @param[out] coeff The coefficients of the polynomial fit. Its length
 * must be @a num_coeff .
 * @n must-be-aligned
 * @param[out] final_mask The final mask data after recursive clipping
 * procedure. Its length must be @a num_data .
 * @n must-be-aligned
 * @param[out] baseline_status Baseline-specific error code.
 * @return Status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsFloat)(
		struct LIBSAKURA_SYMBOL(BaselineContext) const *context,
		size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], float clip_threshold_sigma,
		uint16_t num_fitting_max, size_t num_coeff, double coeff[/*num_coeff*/],
		bool final_mask[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 入力データに対して、与えられた個数の区分からなる3次スプライン曲線を最小二乗フィットし、係数を返す。
 * @details
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @param[in] num_data 配列 @a data 、 @a mask 、 @a final_mask の要素数。
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] mask 入力データに対するマスク情報。要素数は @a num_data でなければならない。
 * 値がfalseの要素に対応する入力データはフィッティングに用いられない。
 * @n must-be-aligned
 * @param[in] clip_threshold_sigma クリッピングの閾値。単位はσ。正値でなければならない。
 * @param[in] num_fitting_max フィッティングを(再帰的に)行う最大回数。
 * 値nが与えられた場合、最初のフィッティング＆差し引きを行った後、
 * 残差データのσを計算し、残差がその値の± @a clip_threshold_sigma
 * 倍を越えるものを除外して再度フィッティング＆差し引きを行うという操作を最大(n-1)回繰り返す。
 * デフォルト値は1、即ち、フィッティング＆差し引きは１回のみ行われ、クリッピングは行わない。
 * もし 0 が渡された場合は、自動的に 1 に変更される。
 * @param[in] num_pieces スプライン曲線の区分の数。正値、かつ、( @a num_data / 4 ) 以下の数でなければならない。
 * @param[out] coeff 出力データ。左端の区分から、0次、1次、2次、3次の順に係数値が入る。要素数は ( @a num_pieces * 4 ) でなければならない。
 * @n must-be-aligned
 * @param[out] final_mask 再帰的クリッピングを経た後の最終的なマスク情報。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[out] baseline_status ベースライン固有のエラーコード。
 * @return 終了ステータス。
 * @~english
 * @brief Extraction of the coefficients of cubic spline fit.
 * @details
 * @param[in] context An object containing baseline model data.
 * @param[in] num_data The number of elements in the arrays @a data,
 * @a mask, and @a final_mask.
 * @param[in] data The input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] mask The input mask data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] clip_threshold_sigma The threshold of clipping in unit
 * of sigma. must be positive.
 * @param[in] num_fitting_max The maximum of total number of times
 * baseline fitting is performed recursively. In case n is given, after
 * the first baseline fitting, subsequent clipping and baseline fitting
 * based on the updated mask are executed (n-1) times at maximum.
 * The default is 1 (i.e., baseline fitting done just once and no
 * clipping applied). In case zero is given, @a num_fitting_max will be
 * automatically changed to 1.
 * @param[in] num_pieces The number of spline pieces.
 * @param[out] coeff The coefficients of the best-fit cubic spline curve.
 * Its length must be @a num_pieces * 4 .
 * @n must-be-aligned
 * @param[out] final_mask The final mask data after recursive clipping
 * procedure. Its length must be @a num_data .
 * @n must-be-aligned
 * @param[out] baseline_status Baseline-specific error code.
 * @return Status code.
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitBaselineCoefficientsCubicSplineFloat)(
		struct LIBSAKURA_SYMBOL(BaselineContext) const *context,
		size_t num_data, float const data[/*num_data*/],
		bool const mask[/*num_data*/], float clip_threshold_sigma,
		uint16_t num_fitting_max, size_t num_pieces,
		double coeff[/*4*num_piece*/], bool final_mask[/*num_data*/],
		LIBSAKURA_SYMBOL(BaselineStatus) *baseline_status)
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 入力データに対して、与えられたベースラインモデルと係数からベースラインを求め、差し引く。
 * @details
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @param[in] num_data 配列 @a data 、 @a out の要素数。
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] num_coeff 配列 @a coeff の要素数。値は正値でなくてはならない。
 * また、値は@a context が持つベースラインモデルの基底関数の数より大きな数になってはならない。
 * @param[in] coeff ベースラインの係数。要素数は @a num_coeff でなければならない。
 * @n must-be-aligned
 * @param[out] out 出力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス
 * @~english
 * @brief Subtract baseline from input data. Baseline is calculated by baseline model and given coefficients.
 * @details
 * @param[in] context An object containing baseline model data.
 * @param[in] num_data The number of elements in @a data and @a out.
 * @param[in] data The input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] num_coeff The number of elements in @a coeff .
 * The value must be positive and must not exceed the number of
 * model bases defined in @a context .
 * @param[in] coeff Coefficients of model data. Its length must
 * be @a num_coeff.
 * @n must-be-aligned
 * @param[out] out The output data. Its length must be @a num_data .
 * @n must-be-aligned
 * @return Status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineUsingCoefficientsFloat)(
		struct LIBSAKURA_SYMBOL(BaselineContext) const *context,
		size_t num_data, float const data[/*num_data*/], size_t num_coeff,
		double const coeff[/*num_data*/], float out[/*num_data*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief 入力データに対して、与えられた係数から3次スプライン曲線によりベースラインを求め、差し引く。
 * @details
 * @param[in] context 3次曲線モデルに関する情報を格納する構造体。
 * @param[in] num_data 配列 @a data 、 @a out の要素数。
 * @param[in] data 入力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @param[in] num_pieces スプライン曲線の区間の数。0の場合は差し引きを行わない。
 * @param[in] coeff 3次スプライン曲線の係数。左端の区分から、0次、1次、2次、3次の順に係数値が入る。要素数は @a num_pieces * 4 でなければならない。
 * @n must-be-aligned
 * @param[in] boundary num_pieces個の区間の0番目、1番目、…の先頭の位置。要素数は @a num_pieces でなければならない。
 * @n must-be-aligned
 * @param[out] out 出力データ。要素数は @a num_data でなければならない。
 * @n must-be-aligned
 * @return 終了ステータス
 * @~english
 * @brief Subtract cubic spline baseline from input data. Baseline is calculated by cubic curve model and given coefficients.
 * @details
 * @param[in] context An object containing model data of cubic curve.
 * @param[in] num_data The number of elements in @a data and @a out .
 * @param[in] data The input data with length of @a num_data .
 * @n must-be-aligned
 * @param[in] num_pieces The number of spline pieces. If zero is
 * given, no subtraction executed.
 * @param[in] coeff Coefficients of cubic spline curve. Its length
 * must be @a num_pieces * 4.
 * @n must-be-aligned
 * @param[in] boundary A 1D array containing the left edge positions of
 * pieces of the spline curve. The values should be stored in
 * left-to-right order. Its length must be @a num_pieces .
 * @n must-be-aligned
 * @param[out] out The output data. Its length must be @a num_data .
 * @n must-be-aligned
 * @return Status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselineCubicSplineUsingCoefficientsFloat)(
		struct LIBSAKURA_SYMBOL(BaselineContext) const *context,
		size_t num_data, float const data[/*num_data*/], size_t num_pieces,
		double const coeff[/*4*num_pieces*/],
		double const boundary[/*num_pieces*/], float out[/*num_data*/])
				LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~japanese
 * @brief ベースラインのモデル基底関数の個数を返す
 * @details
 * @param[in] context ベースラインモデルに関する情報を格納する構造体。
 * @param[in] order モデルのパラメータ。多項式(polynomial, chebyshev)では次数、スプラインでは分割数、三角関数では最大の波数。スプラインの場合は正値でなければならない。それ以外のモデルではゼロも許される。値は@a context を生成したときに指定した@a order より大きな数になってはならない。
 * @param[in] num_coeff ベースラインのモデル基底関数の個数。
 * @return 終了ステータス
 * @~english
 * @brief Return the number of basis functions composing baseline model.
 * @details
 * @param[in] context An object containing baseline model data.
 * @param[in] order Parameter for the specified function.
 * It is the order (for sakura_BaselineType_kPolynomial and
 * sakura_BaselineType_kChebyshev), or the number of pieces (for
 * sakura_BaselineType_kCubicSpline), or the maximum wave number
 * (for sakura_BaselineType_kSinusoid). It must be positive for
 * sakura_BaselineType_kCubicSpline, while other models accept
 * zero value. The value should not exceed the @a order specified
 * in creation of @a context .
 * @param[in] num_coeff Number of basis functions.
 * @return Status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetNumberOfCoefficients)(
		struct LIBSAKURA_SYMBOL(BaselineContext) const *context, uint16_t order,
		size_t *num_coeff) LIBSAKURA_NOEXCEPT LIBSAKURA_WARN_UNUSED_RESULT;

/**
 * @~english
 * @brief Copy elements in the @a src matrix into the @a dst matrix with flipping elements to reorder
 * as some FFT library expects.
 *
 * @details
 * When you provide @a inner_most_untouched = false, @a elements = {3, 4} and @a src = {
 * @code
 *   1,   2,   3,
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
 * If @a inner_most_untouched = true, then you will get @a dst as below.
 * @code
 *  7,   8,   9,
 * 10,  11,  12,
 *  1,   2,   3,
 *  4,   5,   6,
 * @endcode
 *
 * @param[in] inner_most_untouched If true, the order of the inner most dimension is untouched.
 * @param[in] dims Dimensions of the matrix @a src and @a dst. In other words, the number of elements in @a elements.
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
bool inner_most_untouched, size_t dims, size_t const elements[],
		float const src[], float dst[]) LIBSAKURA_NOEXCEPT;

/**
 * @~english
 * @brief Copy elements in the @a src matrix into the @a dst matrix with unflipping elements to the original order.
 *
 * @details
 * When you provide @a inner_most_untouched = false, @a elements = {3, 4} and @a src = {
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
 * @param[in] inner_most_untouched If true, the order of the inner most dimension is untouched.
 * @param[in] dims Dimensions of the matrix @a src and @a dst. In other words, the number of elements in @a elements.
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
bool inner_most_untouched, size_t dims, size_t const elements[],
		float const src[], float dst[]) LIBSAKURA_NOEXCEPT;

/**
 * @~english
 * @brief Same as @ref sakura_FlipMatrixFloat except the element type of the matrixes.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipMatrixDouble)(
bool inner_most_untouched, size_t dims, size_t const elements[],
		double const src[], double dst[]) LIBSAKURA_NOEXCEPT;

/**
 * @~english
 * @brief Same as @ref sakura_UnflipMatrixFloat except the element type of the matrixes.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipMatrixDouble)(
bool inner_most_untouched, size_t dims, size_t const elements[],
		double const src[], double dst[]) LIBSAKURA_NOEXCEPT;

/**
 * @~english
 * @brief Same as @ref sakura_FlipMatrixFloat except the element type of the matrixes.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(FlipMatrixDouble2)(
bool inner_most_untouched, size_t dims, size_t const elements[],
		double const src[][2], double dst[][2]) LIBSAKURA_NOEXCEPT;

/**
 * @~english
 * @brief Same as @ref sakura_UnflipMatrixFloat except the element type of the matrixes.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(UnflipMatrixDouble2)(
bool inner_most_untouched, size_t dims, size_t const elements[],
		double const src[][2], double dst[][2]) LIBSAKURA_NOEXCEPT;

#ifdef __cplusplus
}
/* extern "C" */
#endif

#endif /* LIBSAKURA_LIBSAKURA_SAKURA_H_ */
