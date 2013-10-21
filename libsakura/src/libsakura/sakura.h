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
	 */LIBSAKURA_SYMBOL(Status_kInvalidArgument) = 2, /**
	 * @~japanese
	 * @brief メモリーが足りない
	 */LIBSAKURA_SYMBOL(Status_kNoMemory) = 3, /**
	 * @~japanese
	 * @brief 原因不明のエラー
	 */LIBSAKURA_SYMBOL(Status_kUnknownError) = 99
}LIBSAKURA_SYMBOL(Status);

/**
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
 * @brief Initializes Sakura Library
 * @return Only when sakura_Status_kOK is returned, you can use Sakura Library.
 * @~japanese
 * @brief Sakuraライブラリを初期化する
 *
 * 他の全てのSakuraライブラリAPIの呼び出しに先立って、呼び出すこと。
 * マルチスレッドセーフではないので、単一のスレッドから呼び出すこと。
 * @ref sakura_CleanUp() の呼び出しを挟まず、複数回この関数を呼び出してはならない。
 * @param[in]	allocator	Sakuraライブラリ内で、メモリーを確保するときに呼び出されるアロケーター。NULLの場合はmalloc(3)が使用される。 @ref sakura_UserAllocator
 * @param[in]	deallocator	Sakuraライブラリ内で、メモリーを開放するときに呼び出されるデアロケーター。NULLの場合はfree(3)が使用される。 @ref sakura_UserDeallocator
 * @~
 * MT-unsafe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Initialize)(
LIBSAKURA_SYMBOL(UserAllocator) allocator,
LIBSAKURA_SYMBOL(UserDeallocator) deallocator);

/**
 * @~english
 * @brief Cleans up Sakura Library
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
 * @brief Returns the current time.
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
 * @~japanese
 * @brief Sakuraライブラリがベクトル演算を行う配列に期待するアライメントを返す
 *
 * @return Sakuraライブラリは、戻り値の倍数アドレスにアラインされることを期待する
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
 * @brief @ref sakura_ComputeStatistics の結果を格納する構造体
 */
typedef struct {
	size_t count; /**< 個数 */
	float sum; /**< 合計 */
	float mean; /**< 平均 */
	float rms; /**< 二乗平均平方根 */
	float stddev; /**< 分散 */
	float min; /**< 最小値 */
	float max; /**< 最大値 */
	int index_of_min; /**< 最小値のインデックス(有効な値がなければ-1) */
	int index_of_max; /**< 最大値のインデックス(有効な値がなければ-1) */
}LIBSAKURA_SYMBOL(StatisticsResult);
/**
 * @~japanese
 * @brief 統計値を計算する。どのような統計値を算出するかは
 * @ref sakura_StatisticsResult を参照
 *
 * @param[in] num_data @a data 及び@a is_valid の要素の数
 * @param[in] data 対象となるデータ。対応する@a is_valid がtrueの場合、Inf/NaNであってはならない。
 * <br/>must-be-aligned
 * @param[in] is_valid データのマスク。この値が false だと、
 * 対応する@a data の要素が無視される。
 * <br/>must-be-aligned
 * @param[out] result 結果の格納先。計算不可能な場合は、構造体内のメンバーの値にNaNが設定される。
 * 同じ値が複数あった場合、どの値のインデックスが@a index_of_min, @a index_of_maxに格納されるかは不定である。
 * @return 終了ステータス
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(ComputeStatistics)(size_t num_data,
		float const data[], bool const is_valid[],
		LIBSAKURA_SYMBOL(StatisticsResult) *result);
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
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SortValidValuesDensely)(
		size_t num_data, bool const is_valid[], float data[],
		size_t *new_num_data);
/**
 * @~japanese
 * @brief 畳み込みしながらグリッドする
 *
 * @a start_spectrum から @a end_spectrum までの、@a x , @a y 座標上の @a value の値を
 * @a convolution_table で表される広がりを持った点として、@a grid 上にプロットする。
 * @a grid にプロットする際は、@a polarization_map , @a channel_map によって偏波とチャネルのマッピングが行われる。
 *
 * 各浮動小数点の数値はNaN/+-Infであってはならない。
 * @param[in] num_spectra 次の関係でなければならない。 0 <= start_spectrum <= end_spectrum <= num_spectra
 * @param[in] start_spectrum 開始spectrumの添字
 * @param[in] end_spectrum 終了spectrumの添字+1
 * @param[in] spectrum_mask	要素数はnum_spectra。falseのスペクトルは無視される。<br/>must-be-aligned
 * @param[in] x 要素数は@a num_spectra 。2次元平面に投射済みのx座標。<br/>must-be-aligned
 * @param[in] y 要素数は@a num_spectra 。2次元平面に投射済みのy座標。<br/>must-be-aligned
 * @param[in] support @a width x @a height 平面における畳み込みカーネルの広がり(中心か らのpixel数)。範囲は、0 < support <= (INT32_MAX - 1) / 2<br/>
 * ただし、@a support * @a sampling に比例するサイズの領域をスタック上に確保するので、@a support * @a sampling が大きな値となる場合、スタックオーバーフローを起こす。
 * @param[in] sampling 畳み込みカーネルの精度(/pixel)。範囲は、0 < sampling <= INT32_MAX
 * @param[in] num_polarizations 範囲は、0 < num_polarization <= INT32_MAX
 * @param[in] polarization_map	要素数は、num_polarization。各要素の値は、[0,num_polarization_for_grid)でなければならない。<br/>must-be-aligned
 * @param[in] num_channels 範囲は、0 < num_channels <= INT32_MAX
 * @param[in] channel_map	要素数は、num_channels。各要素の値は、[0,num_channels_for_grid)でなければならない。<br/>must-be-aligned
 * @param[in] mask	要素のレイアウトは、[num_spectra][num_polarization][num_channels]。falseの場合は、該当するスペクトル、偏波、チ ャネルの組み合わせのデータは無視される。<br/>must-be-aligned
 * @param[in] value	要素のレイアウトは、[num_spectra][num_polarization][num_channels]。グリッディングする値。<br/>must-be-aligned
 * @param[in] weight 要素のレイアウトは、[num_spectra][num_channels]。重み。<br/>must-be-aligned
 * @param[in] weight_only @a value に重みを掛けたものではなく、重み自体をグリッディングする場合は、true。
 * @param[in] num_convolution_table @a convolution_table の要素数。 範囲は、ceil(sqrt(2.)*(support+1)*sampling) <= convolution_table <= INT32_MAX / 32
 * @param[in] convolution_table	要素数は、@a num_convolution_table 。畳み込みに使用する重みカーブ。各要素の値は、NaN,Infであってはならない。要素0が中心を表す。<br/>must-be-aligned
 * @param[in] num_polarization_for_grid 範囲は、0 < num_polarization_for_grid <= INT32_MAX
 * @param[in] num_channels_for_grid 範囲は、0 < num_channels_for_grid <= INT32_MAX
 * @param[in] width 範囲は、0 < width <= INT32_MAX
 * @param[in] height 範囲は、0 < height <= INT32_MAX
 * @param[out] weight_sum	要素のレイアウトは、[num_polarization_for_grid][num_channels_for_grid]。重みの合計。<br/>must-be-aligned
 * @param[out] weight_of_grid	要素のレイアウトは、[height][width][num_polarization_for_grid][num_channels_for_grid]。グリッドの重み。<br/>must-be-aligned
 * @param[out] grid	要素のレイアウトは、[height][width][num_polarization_for_grid][num_channels_for_grid]。グリッディング結果。<br/>must-be-aligned
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GridConvolving)(size_t num_spectra,
		size_t start_spectrum, size_t end_spectrum,
		bool const spectrum_mask[/*num_spectra*/],
		double const x[/*num_spectra*/], double const y[/*num_spectra*/],
		size_t support, size_t sampling, size_t num_polarization,
		uint32_t const polarization_map[/*num_polarization*/],
		size_t num_channels, uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const value/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const weight/*[num_spectra]*/[/*num_channels*/], bool weight_only,
		size_t num_convolution_table/*= ceil(sqrt(2.)*(support+1)*sampling)*/,
		float const convolution_table[/*num_convolution_table*/],
		size_t num_polarization_for_grid, size_t num_channels_for_grid,
		size_t width, size_t height,
		double weight_sum/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]);

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
 * Array elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] num_condition The number of elements in the arrays, @a lower_bounds
 * and @a upper_bounds.
 * @param[in] lower_bounds The input array of size, @a num_condition.
 * Array elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] upper_bounds The input array of size, @a num_condition.
 * Array elements should not contain Inf nor NaN.
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
 * 要素はInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[in] num_condition 一次元配列@a lower_bounds 及び@a upper_bounds の要素の数。
 * 下限値と上限値の組の数を表す。
 * @param[in] lower_bounds 入力一次元配列。検定条件の下限値を格納する。要素数は@a num_condition でなければならない。
 * 要素はInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[in] upper_bounds 入力一次元配列。検定条件の上限値を格納する。要素数は@a num_condition でなければならない。
 * 要素はInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
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
 * @brief Returns if the values in input array are greater than a threshold.
 * @details Returns true if the corresponding element in the input array
 * is greater than a @a threshold,
 * @par
 * @a data[i] >  @a threshold .
 *
 * @note
 * No evaluation is done when the data array is zero length, i.e., @a num_data = 0.@n
 *
 * @param[in] num_data The number of elements in the arrays, @a data
 * and @a result
 * @param[in] data An input array of size, @a num_data.
 * Array elements should not contain Inf nor NaN.
 * @n must-be-aligned
 * @param[in] threshold the threshold of evaluation. The value should not be Inf nor NaN.
 * @param[out] result The output array of size, @a num_data.
 * @n must-be-aligned
 * @return status code
 * @~japanese
 * @brief 入力配列の値が、与えられたしきい値より大きいかどうかを検定する。
 * @details 入力配列の要素の値がしきい値( @a threshold )より大きければ真を返す。@n
 * すなわち、
 * @par
 * @a data[i] >  @a threshold .
 *
 * を検定する。@n
 *
 * @note
 * 入力配列の要素数が0 (@a num_data = 0)の時は、演算は実行されない。@n
 *
 * @param[in] num_data 一次元配列@a data 及び@a result の要素の数。
 * @param[in] data 入力一次元配列。検定の対象となる値を格納する。要素数は@a num_data でなければならない。
 * 要素はInfやNaNを含んではならない。InfやNaNの場合の動作は不定。
 * @n must-be-aligned
 * @param[in] threshold しきい値。値はInfやNaNであってはならない。InfやNaNの場合の動作は不定。
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
 * @brief 一次元の補間を行う。
 * @details 長さ@a num_base の配列@a x_base と@a y_base で定義されるデータ点を
 * もとにして1次元の補間を行う。補間したい点のx座標のリストを長さ
 * @a num_interpolated の配列@a x_interpolated に渡すと、補間結果が
 * 長さ@a num_interpolated の配列@a y_interpolated に格納される。
 * 外挿は行わない（データ点が片側にしかない場合にはそのデータ点の値が出力配列
 * @a y_interpolated にセットされる）。
 *
 * 戻り値は終了ステータスである。正常終了の場合、
 * @link sakura_Status::sakura_Status_kOK sakura_Status_kOK @endlink
 * を返す。
 * 引数に不正がある場合には
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * を返す。
 * @link sakura_Status::sakura_Status_kInvalidArgument sakura_Status_kInvalidArgument @endlink
 * が返された場合、
 * 考えられる原因は以下の三つである。
 *     - @a interpolatin_method が正しくない
 *     - 引数に渡した配列がアラインされていない
 *
 * また、原因不明のエラーでは
 * @link sakura_Status::sakura_Status_kUnknownError sakura_Status_kUnknownError @endlink
 * を返す。
 *
 * @par
 * @pre @a x_base および@a x_interpolated は昇順または降順にソートされていなければ
 * ならない。また、@a x_base の要素には重複があってはならない。
 *
 * @par 昇順の場合と降順の場合の速度の違いについて:
 * @a x_base または @a x_interpolated が降順にソートされている場合、
 * 内部では配列要素をコピーして昇順に並べ替えた上で補間を行う。そのため、降順の場合は
 * 昇順よりも処理が遅くなる。
 *
 * @par 多項式補間の動作について:
 * @a polynomial_order はあくまで最大次数を規定するものであり、その次数で必ず
 * 補間が行われるとは限らない。たとえば、@a polynomial_order が2（二次多項式による補間）
 * で@a num_base が2の場合、実際には2点を通る一次多項式が一意に決まるため、二次多項式に
 * よる補間ではなく一次多項式による補間（線形補間）が行われる。
 * @par
 * @a polynomial_order に0を指定した場合、最近接補間が行われる。
 *
 * @par
 * @param[in] interpolation_method 補間方法
 * @param[in] polynomial_order 多項式補間法の場合の最大次数。
 * 実際の次数は、@a num_base との兼ね合いで決まる。
 * @param[in] num_base 補間のためのデータ点の数。
 * @param[in] x_base 補間のためのデータ点のx座標。
 * 要素数は@a num_base でなければならない。
 * @a x_base は昇順または降順にソートされていなければならない。
 * must-be-aligned
 * @param[in] num_array
 * @param[in] y_base 補間のためのデータ点のy座標。
 * 要素数は@a num_base * @a num_array でなければならない。
 * must-be-aligned
 * @param[in] num_interpolated 補間したいデータ点の数。
 * @param[in] x_interpolated 補間したいデータ点のx座標。
 * 要素数は@a num_interpolated でなければならない。
 * @a x_interpolated は昇順または降順にソートされていなければならない。
 * must-be-aligned
 * @param[out] y_interpolated 補間結果。
 * 要素数は@a num_interpolated * @a num_array でなければならない。
 * must-be-aligned
 * @return 終了ステータス。
 *
 * @~english
 * @brief Perform one-dimensional interpolation
 * @details
 * @param[in] interpolation_method interpolation method.
 * @param[in] polynomial_order maximum polynomial order for polynomial interpolation.
 * Actual order will be determined by a balance
 * between @a polynomial_order and @a num_base.
 * @param[in] num_base number of elements for data points.
 * @param[in] x_base x-coordinate of data points. Its length must be @a num_base.
 * It must be sorted either ascending or descending.
 * @param[in] num_array
 * @param[in] y_base y-coordinate of data points. Its length must be @a num_base.
 * @param[in] num_interpolated number of elements for points that wants to get
 * interpolated value.
 * @param[in] x_interpolated x-coordinate of points that wants to get interpolated
 * value. Its length must be @a num_interpolated.
 * @param[out] y_interpolated storage for interpolation result. Its length must be
 * @a num_interpolated.
 * @return status code.
 *
 * @~
 * MT-safe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Interpolate1DXFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_x_base,
		double const x_base[/*num_x_base*/], size_t num_y,
		float const data_base[/*num_x_base*num_y*/], size_t num_x_interpolated,
		double const x_interpolated[/*num_x_interpolated*/],
		float data_interpolated[/*num_x_interpolated*num_y*/]);

/**
 * @brief Perform one-dimensional interpolation along interpolation axis
 * of two-dimensional array.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Interpolate1DYFloat)(
LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
		uint8_t polynomial_order, size_t num_y_base,
		double const y_base[/*num_y_base*/], size_t num_x,
		float const data_base[/*num_x*num_y_base*/], size_t num_y_interpolated,
		double const y_interpolated[/*num_y_interpolated*/],
		float data_interpolated[/*num_x*num_y_interpolated*/]);

/**
 * @japanese
 * @brief position switch calibrationを実行する。
 * @details
 * position switch calibrationを実行する。
 * 具体的には、
 * @f[
 *  {\rm result} = {\rm scaling\_factor} \times \frac{{\rm target} - {\rm reference}}
 *               {\rm reference}
 * @f]
 * を実行する。これは、position switch観測の温度較正
 * @f[
 *  \rm{calibrated} = T_{\rm sys} \times \frac{{\rm ON} - {\rm OFF}}{\rm OFF}
 * @f]
 * に相当する処理である。ただし、@f$ T_{\rm{sys}}@f$ はシステム雑音温度、ONおよびOFFはそれぞれ
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
 * @english
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
 * @brief スムージングに使うカーネルタイプを列挙
 * @~english
 * @brief Enumerations to define kernel types.
 */
typedef enum {
	/**
	 * @brief Gaussian
	 */LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian), /**
	 * @brief BoxCar
	 */LIBSAKURA_SYMBOL(Convolve1DKernelType_kBoxcar), /**
	 * @brief Hanning
	 */LIBSAKURA_SYMBOL(Convolve1DKernelType_kHanning), /**
	 * @brief Hamming
	 */LIBSAKURA_SYMBOL(Convolve1DKernelType_kHamming)
}LIBSAKURA_SYMBOL(Convolve1DKernelType);
/**
 * @brief Context struct for Convolution
 */
struct LIBSAKURA_SYMBOL(Convolve1DContext);
/**
 * @~japanese
 * @brief コンテキストを作成する。
 * @details
 * @param[in] num_channel チャネル数
 * @param[in] kernel_type カーネルタイプ
 * Gaussian,BoxCar,Hanning,Hammingを選択可能。
 * @param[in] kernel_width カーネル幅
 * カーネルのシグマ＝カーネル幅／√(８ln2) により計算する。
 * @param[in] use_fft FFTを行うか否かのフラグ。true=行う。false=行わない。
 * @param[in,out] context コンテキスト
 * FFT済みカーネル、作成済み実数複素数FFTプラン、複素数実数FFTプラン、チャネル数、実数配列
 * を持つ。
 * @return 終了ステータス。
 * @~english
 * @brief Create context
 * @details
 * @param[in] num_channel number of channel of input spectrum. @num_channel must
 * be positive.
 * @param[in] kernel_type type of kernel(Gaussian,BoxCar,Hanning,Hamming)
 * @kernel_type is defined as enum.
 * @param[in] kernel_width kernel width which proposion to sigma
 * @param[in] use_fft if use fft then true, if not, faulse.
 * @param[in,out] context context to store number of channel,fftwf_plan,
 * fftwf_complex and flexible real array.
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateConvolve1DContext)(
		size_t num_channel, LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
		size_t kernel_width, bool use_fft,
		struct LIBSAKURA_SYMBOL(Convolve1DContext) **context);
/**
 * @~japanese
 * @brief コンボリューションを行う。
 * @details 入力スペクトルに対するコンボリューションによるスムージング処理を行う。
 * FFTを使用する場合：
 * 入力スペクトルに対しFFTを行った複素数配列と、事前に作った複素数のFFT済みカーネルとを
 * 掛け合せ一つの複素数配列を得る。それをIFFTし、実数配列である出力スペクトルを得る。
 * FFTを使用しない場合：
 * FFTを使用せず、実数配列のまま入力スペクトルとカーネルとでコンボリューションを行う。
 * @param[in,out] context コンテキスト
 * FFT済みカーネル、作成済み実数複素数FFTプラン、複素数実数FFTプラン、チャネル数、実数配列
 * を持つ。
 * @param[in] input_spectrum 入力スペクトル
 * 配列の長さは @a num_channel と同じ。
 * @param[in] input_flag 入力フラグ
 * 配列の長さは @a num_channel と同じ。
 * @param[out] output_spectrum 出力スペクトル
 * 配列の長さは @a num_channel と同じ。
 * @return 終了ステータス。
 * @~english
 * @brief Do Convolution
 * @details It can do smoothing input spectrum by doing convolution
 * with using fft or not. If using fft, fft applied kernl which is
 * already included context by CreateConvolve1DContext will be multiplied
 * with input spectrum by complex-complex multiplication and then
 * the multiplied complex array will be created. Finally IFFT will be
 * applied against it and then real output spectrum will be created.
 * @param[in,out] context context which contain @a fftwf_plan, @a fftw_complex
 * and @a num_channel, @a input_real_array
 * @param[in] input_spectrum input spectrum
 * Its length equals to number of channel
 * @param[in] input_flag
 * Its length equals to number of channel
 * @param[out] output_spectrum
 * Its length equals to number of channel
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Convolve1D)(
		struct LIBSAKURA_SYMBOL(Convolve1DContext) **context,
		float input_spectrum[/*num_in*/], bool const input_flag[/*num_in*/],
		float output_spectrum[/*num_in*/]);
/**
 * @~japanese
 * @brief コンテキストを作成する。
 * @details
 * @param[in] num_channel チャネル数
 * @param[in] kernel_type カーネルタイプ
 * Gaussian,BoxCar,Hanning,Hammingを選択可能。
 * @param[in] kernel_width カーネル幅
 * カーネルのシグマ＝カーネル幅／√(８ln2) により計算する。
 * @param[in] use_fft FFTを行うか否かのフラグ。true=行う。false=行わない。
 * @param[in,out] context コンテキスト
 * FFT済みカーネル、作成済み実数複素数FFTプラン、複素数実数FFTプラン、チャネル数、実数配列
 * を持つ。
 * @return 終了ステータス。
 * @~english
 * @brief Create context
 * @details
 * @param[in] num_channel number of channel of input spectrum. @num_channel must
 * be positive.
 * @param[in] kernel_type type of kernel(Gaussian,BoxCar,Hanning,Hamming)
 * @kernel_type is defined as enum.
 * @param[in] kernel_width kernel width which proposion to sigma
 * @param[in] use_fft if use fft then true, if not, faulse.
 * @param[in,out] context context to store number of channel,fftwf_plan,
 * fftwf_complex and flexible real array.
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DestroyConvolve1DContext)(
		struct LIBSAKURA_SYMBOL(Convolve1DContext) *context);
/**
 * @brief Logical operation AND between two boolean arrays.
 * @details
 * @param[in] num_in channel number.
 * @param[in] in1 the first input data.
 * @param[in] in2 the second input data.
 * @param[out] out the output data.
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateLogicalAnd)(size_t num_in,
bool const in1[/*num_in*/], bool const in2[/*num_in*/],
bool out[/*num_in*/]);
/**
 * @brief Compute subtraction between two float arrays (in1 - in2).
 * @details
 * @param[in] num_in channel number.
 * @param[in] in1 the first input data.
 * @param[in] in2 the second input data.
 * @param[out] out the output data.
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateFloatSubtraction)(
		size_t num_in, float const in1[/*num_in*/], float const in2[/*num_in*/],
		float out[/*num_in*/]);
/**
 * @brief Compute values for Least-Square fitting from input data and a set of model data.
 * @details
 * @param[in] num_in channel number.
 * @param[in] in_data input data.
 * @param[in] in_mask channel mask for the input data.
 * @param[in] num_model number of model functions.
 * @param[in] model model data.
 * @param[out] out the least-square matrix
 * @param[out] out_vector the right side value of the simultaneous equations of least-square fitting.
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetLeastSquareMatrix)(size_t num_in,
		float const in_data[/*num_in*/], bool const in_mask[/*num_in*/],
		size_t num_model, double const model[/*num_model * num_in*/],
		double out[/*num_model * num_model*/],
		double out_vector[/*num_model*/]);
/**
 * @brief Solve simultaneous equations via LU decomposition.
 * @details
 * @param[in] num_eqn number of equations.
 * @param[in] lsq_matrix0 the matrix in the left side of simultaneous equations.
 * @param[in] lsq_vector0 the right side value of simultaneous equations.
 * @param[out] out the solution.
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLU)(
		size_t num_eqn, double const lsq_matrix0[/*num_eqn * num_eqn*/],
		double const lsq_vector0[/*num_eqn*/], double out[/*num_eqn*/]);
/**
 * @brief Compute the best-fit model spectrum using model spectra and coefficients.
 * @details
 * @param[in] num_chan channel number.
 * @param[in] num_eqn number of models.
 * @param[in] model model data.
 * @param[in] coeff coefficients.
 * @param[out] out the best-fit model.
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DoGetBestFitModel)(size_t num_chan,
		size_t num_eqn, double const model[/*num_eqn * num_chan*/],
		double const coeff[/*num_eqn*/], float out[/*num_in*/]);
/**
 * @brief Compute the best-fit model spectrum by least-square fitting.
 * @details
 * @param[in] num_in channel number.
 * @param[in] in_data the input data.
 * @param[in] in_mask channel mask for the input data.
 * @param[in] num_model number of model functions.
 * @param[in] model model data.
 * @param[out] out the best-fit model.
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitModel)(size_t num_in,
		float const in_data[/*num_in*/], bool const in_mask[/*num_in*/],
		size_t num_model, double const model[/*num_model * num_in*/],
		float out[/*num_in*/]);
/**
 * @~japanese
 * @brief 入力スペクトルに多項式ベースラインをフィットし差し引く。
 * @~english
 * @brief Fit a baseline and subtract it from input spectrum.
 * @details
 * @param[in] num_chan channel number.
 * @param[in] in_data the input data.
 * @param[in] in_mask channel mask for the input data.
 * @param[in] order order of polynomial model.
 * @param[in] clipping_threshold_sigma the threshold of clipping.
 * @param[in] num_fitting_max the maximum number of recursive baseline fitting.
 * In case n is given, after the first baseline fitting, clipping and baseline
 * fitting based on the updated channel mask are executed for (n-1) times at maximum.
 * Default is 1 (i.e., baseline fitting done just once and no clipping applied)
 * @param[in] get_residual set the output to be (input - best-fit) if true, or the best-fit value if false.
 * @param[out] out the output data.
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SubtractBaselinePolynomial)(
		size_t num_chan, float const in_data[/*num_chan*/],
		bool const in_mask[/*num_chan*/], unsigned int order,
		float clipping_threshold_sigma, unsigned int num_fitting_max,
		bool get_residual, float out[/*num_chan*/]);
/**
 * @brief Compute a set of model spectra.
 * @details
 * @param[in] num_chan channel number.
 * @param[in] order order of polynomial model.
 * @param[out] out the set of baseline model data.
 * @return status code.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBaselineModel)(size_t num_chan,
		unsigned int order, double out[/*(order+1)*num_chan*/]);
/**
 * @brief Actually fit a baseline and subtract it from input spectrum.
 * @details
 * @param[in] num_chan channel number.
 * @param[in] in_data the input data.
 * @param[in] in_mask channel mask for the input data.
 * @param[in] num_model number of model functions.
 * @param[in] model_data the model data.
 * @param[in] clipping_threshold_sigma the threshold of clipping.
 * @param[in] num_fitting_max the maximum number of iterative clipping.
 * @param[in] get_residual set the output to be (input - best-fit) if true, or the best-fit value if false.
 * @param[out] out the output data.
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DoSubtractBaseline)(size_t num_chan,
		float const in_data[/*num_chan*/], bool const in_mask[/*num_chan*/],
		unsigned int num_model, double const model_data[/*num_model * num_chan*/],
		float clipping_threshold_sigma, unsigned int num_fitting_max,
		bool get_residual, float out[/*num_chan*/]);

#ifdef __cplusplus
}
/* extern "C" */
#endif

#endif /* LIBSAKURA_LIBSAKURA_SAKURA_H_ */
