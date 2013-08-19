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
 * @brief 関数の呼び出し結果を表す。
 *
 */
typedef enum {
	/**
	 * @~japanese
	 * @brief 成功または正常
	 */LIBSAKURA_SYMBOL(Status_kOK) = 0,
	/**
	 * @~japanese
	 * @brief 失敗または異常
	 */LIBSAKURA_SYMBOL(Status_kNG) = 1,
	/**
	 * @~japanese
	 * @brief 引数が不正だった。
	 */LIBSAKURA_SYMBOL(Status_kInvalidArgument) = 2,
	/**
	 * @~japanese
	 * @brief 原因不明のエラー。
	 */LIBSAKURA_SYMBOL(Status_kUnknownError) = 99
}LIBSAKURA_SYMBOL(Status);

/**
 * @~english
 * @brief Initializes Sakura Library
 * @return Only when sakura_Status_kOK is returned, you can use Sakura Library.
 * @~japanese
 * @brief Sakuraライブラリを初期化する。
 *
 * 他の全てのSakuraライブラリAPIの呼び出しに先立って、呼び出すこと。
 * マルチスレッドセーフではないので、単一のスレッドから呼び出すこと。
 * @ref sakura_CleanUp() の呼び出しを挟まず、複数回この関数を呼び出してはならない。
 * @~
 * MT-unsafe
 */LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Initialize)();

/**
 * @~english
 * @brief Cleans up Sakura Library
 * @~japanese
 * @brief Sakuraライブラリをクリーンアップする。
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
 * @brief 現在時刻(単位は秒)を返す。
 *
 * 精度はgettimeofday(2)依存。
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
 * @param ptr アラインされているか調べたいアドレス。NULL も受け付ける。
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
 * アドレスを必要最小限増加させ、アラインされたアドレスを返す。
 * @a arena がアラインされていれば@a arena を返す。
 *
 * @param arena 確保されている領域の先頭アドレス
 * @param size_of_arena 確保されている領域のサイズ
 * @param size_required アライン後も利用可能でなければならないサイズ
 * @return アラインされたアドレス。もし、 @a size_required を格納するのに
 * 十分な大きさの@a size_of_arena が無いならば、 nullptr を返す。
 * @~
 * MT-safe
 */
void const *LIBSAKURA_SYMBOL(AlignAny)(size_t size_of_arena, void const *arena,
		size_t size_required);
float const *LIBSAKURA_SYMBOL(AlignFloat)(size_t elements_in_arena,
		float const *arena, size_t elements_required);
double const *LIBSAKURA_SYMBOL(AlignDouble)(size_t elements_in_arena,
		double const *arena, size_t elements_required);


/**
 * @~japanese
 * Sakuraライブラリが動的にメモリーを確保するときに呼び出す関数の型。
 * 関数はリエントラントな実装でなければならない。
 * @param size
 */
typedef void *(*LIBSAKURA_SYMBOL(UserAllocator))(size_t size);

/**
 * @~japanese
 * Sakuraライブラリが動的に確保したメモリーを開放するときに呼び出す関数の型。
 * 関数はリエントラントな実装でなければならない。
 * @param pointer
 */
typedef void (*LIBSAKURA_SYMBOL(UserDeallocator))(void *pointer);

/**
 * @~japanese
 * @ref sakura_ComputeStatistics の結果を格納する構造体。
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
 * @ref sakura_StatisticsResult を参照。The bit operation is invoked
 *
 * @param num_data @a data 及び@a is_valid の要素の数
 * @param data 対象となるデータ。対応する@a is_valid がtrueの場合、NaNであってはならない。
 * @param is_valid データのマスク。この値が false だと、
 * 対応する@a data の要素が無視される
 * @param result 結果の格納先。計算不可能な場合は、構造体内のメンバーの値にNaNが設定される。
 * 同じ値が複数あった場合、どの値のインデックスが@a index_of_min, @a index_of_maxに格納されるかは不定である。
 *
 * @~
 * MT-safe
 */
void LIBSAKURA_SYMBOL(ComputeStatistics)(size_t num_data, float const data[],
		bool const is_valid[], LIBSAKURA_SYMBOL(StatisticsResult) *result);

/**
 * @~japanese
 * @brief validな値のみを先頭に詰めて昇順にソートする。
 *
 * @param is_valid データのマスク。この値が false だと、
 * 対応する@a data の要素が無視される
 * @param num_data @a data 及び@a is_valid の要素の数
 * @param data ソート対象のデータ。In placeでソートするので、この配列内の順序は変更される。
 * 対応する@a is_valid がtrueの場合、NaNであってはならない。
 * @return (validでないデータを除いた)ソートされた要素数( <= @a num_data)
 * @~
 * MT-safe
 */
size_t LIBSAKURA_SYMBOL(SortValidValuesDensely)(size_t num_data,
		bool const is_valid[], float data[]);

/**
 * @~japanese
 * @brief 畳み込みしながらグリッドする。
 * 各浮動小数点の数値はNaN/+-Infであってはならない。
 * @param num_spectra 次の関係でなければならない。 0 <= start_spectrum <= end_spectrum <= num_spectra
 * @param start_spectrum
 * @param end_spectrum
 * @param spectrum_mask	要素数はnum_spectra。falseのスペクトルは無視される。
 * @param x
 * @param y
 * @param support
 * @param sampling
 * @param num_polarizations
 * @param polarization_map	各要素の値は、[0,num_polarization_for_grid)でなければならない。要素数は、num_polarizationでなければならない。
 * @param num_channels
 * @param channel_map	各要素の値は、[0,num_channels_for_grid)でなければならない。要素数は、num_channelsでなければならない。
 * @param mask
 * @param value
 * @param weight
 * @param do_weight
 * @param num_convolution_table >= ceil(sqrt(2.)*(support+1)*sampling)
 * @param convolution_table
 * @param num_polarization_for_grid
 * @param num_channels_for_grid
 * @param width
 * @param height
 * @param weight_sum
 * @param weight_of_grid
 * @param grid
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GridConvolving)(size_t num_spectra,
		size_t start_spectrum, size_t end_spectrum,
		bool const spectrum_mask[/*num_spectra*/],
		double const x[/*num_spectra*/],
		double const y[/*num_spectra*/],
		size_t support, size_t sampling,
		size_t num_polarization,
		uint32_t const polarization_map[/*num_polarization*/],
		size_t num_channels,
		uint32_t const channel_map[/*num_channels*/],
		bool const mask/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const value/*[num_spectra][num_polarization]*/[/*num_channels*/],
		float const weight/*[num_spectra]*/[/*num_channels*/],
		bool do_weight,
		size_t num_convolution_table/*= ceil(sqrt(2.)*(support+1)*sampling)*/,
		float const convolution_table[/*num_convolution_table*/],
		size_t num_polarization_for_grid, size_t num_channels_for_grid,
		size_t width, size_t height,
		double weight_sum/*[num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float weight_of_grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/],
		float grid/*[height][width][num_polarization_for_grid]*/[/*num_channels_for_grid*/]
		);

/**
 * @~english
 * @brief Invoke bit operation AND between an utint8_t value and array.
 * @details Invokes the following bit operation to @a i- th element of @a out : @n
 * @a out [i] = ( @a edit_mask [i] ? (@a bit_maks & @a in [i]) : @a in [i] )
 *
 * @param bit_mask An uint8_t value. The bit operation is invoked
 * between this value and the array, @a in.
 * @param num_in The number of elements in the arrays, @a in,
 * @a edit_mask, and @a out.
 * @param in An input array (uint8_t) of size, @a num_in. The bit operation
 * is invoked between this array and @a bit_mask.
 * @param edit_mask A boolean mask array of size, @a num_in. The bit operation
 * is skipped for the elements with the value, false.
 * @param out The output array (uint8_t) of size, @a num_in. It stores the result
 * of the bit operation between @a bit_mask and @a in. The bit operation is skipped
 * and the value in array, @a in, is adopted for the elements where corresponding
 * elements in @a edit_mask is false.
 * @return @a sakura_Status
 * @~japanese
 * @brief TBD
 * @details Invokes the following bit operation to @a i- th element of @a out : @n
 * @a out [i] = ( @a edit_mask[i] ? (@a bit_maks & @a in [i]) : @a in [i] )
 *
 * @param bit_mask TBD
 * @param num_in @a in, @a edit_mask 及び@a out の要素の数。
 * @param in TBD
 * @param edit_mask データのマスク。この値が false だと、
 * 対応する要素のTBDが無視される
 * @param out 結果の格納先。
 * @return @a sakura_Status
 *
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint8And)(uint8_t bit_mask, size_t num_in,
		uint8_t const in[/*num_in*/], bool const edit_mask[/*num_in*/], uint8_t out[/*num_in*/]);

/**
 * @~
 * @brief Invoke bit operation AND between an utint32_t value and array.
 * @details Invokes the following bit operation to @a i- th element of @a out : @n
 * @a out [i] = ( @a edit_mask [i] ? (@a bit_maks & @a in [i]) : @a in [i] )
 *
 * @param bit_mask An uint32_t value. The bit operation is invoked
 * between this value and the array, @a in.
 * @param num_in The number of elements in the arrays, @a in,
 * @a edit_mask, and @a out.
 * @param in An input array (uint32_t) of size, @a num_in. The bit operation
 * is invoked between this array and @a bit_mask.
 * @param edit_mask A boolean mask array of size, @a num_in. The bit operation
 * is skipped for the elements with the value, false.
 * @param out The output array (uint32_t) of size, @a num_in. It stores the result
 * of the bit operation between @a bit_mask and @a in. The bit operation is skipped
 * and the value in array, @a in, is adopted for the elements where corresponding
 * elements in @a edit_mask is false.
 * @return @a sakura_Status
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateBitsUint32And)(uint32_t bit_mask, size_t num_in,
		uint32_t const in[/*num_in*/], bool const edit_mask[/*num_in*/], uint32_t out[/*num_in*/]);


/**
 * @~japanese
 * @brief 補間方法を定義するための列挙型。
 * @~english
 * @brief Enumerations to define interpolation types.
 **/
typedef enum {
	/**
	 * @brief Nearest interpolation
	 **/
	LIBSAKURA_SYMBOL(InterpolationMethod_kNearest),
	/**
	 * @brief Linear interpolation
	 **/
	LIBSAKURA_SYMBOL(InterpolationMethod_kLinear),
	/**
	 * @brief Polynomial interpolation
	 **/
	LIBSAKURA_SYMBOL(InterpolationMethod_kPolynomial),
	/**
	 * @brief Spline interpolation
	 **/
	LIBSAKURA_SYMBOL(InterpolationMethod_kSpline),
	/**
	 * @brief Number of interpolation methods implemented
	 **/
	LIBSAKURA_SYMBOL(InterpolationMethod_kNumMethod)
} LIBSAKURA_SYMBOL(InterpolationMethod);


/**
 * @brief Perform one-dimensional interpolation
 **/
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(Interpolate1dFloat)(LIBSAKURA_SYMBOL(
		InterpolationMethod) interoplation_method, int polynomial_order,
		size_t num_base, double const x_base[], float const y_base[],
		size_t num_interpolated, double x_interpolated[], float y_interpolated[]);


/**
 * @brief Perform pseudo two-dimensional interpolation
 **/
//LIBSAKURA_SYMBOL(Status) InterpolatePseudo2dFloat(LIBSAKURA_SYMBOL(InterpolationMethod) interpolation_method,
//		int polynomial_order, double x_interpolated, size_t num_base, double x_base[],
//		size_t num_array, float *y_base[], float y_interpolated[]);

/**
 * @~japanese
 * @brief スムージングに使うカーネルタイプを列挙
 * @~english
 * @brief Enumerations to define kernel types.
 **/
typedef enum {
	/**
	 * @brief Gaussian
	 **/
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kGaussian),
	/**
	 * @brief BoxCar
	 **/
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kBoxcar),
	/**
	 * @brief Hanning
	 **/
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kHanning),
	/**
	 * @brief Hamming
	 **/
	LIBSAKURA_SYMBOL(Convolve1DKernelType_kHamming)
} LIBSAKURA_SYMBOL(Convolve1DKernelType);

/**
 * @brief Context struct for Convolution
 **/
typedef struct {
  float fft_applied_kernel[0];
}LIBSAKURA_SYMBOL(Convole1DContext);

/**
 * @brief Creating 1D Kernel with FFT or without FFT
 **/

LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(CreateConvole1DContext)(
             size_t num_channel,LIBSAKURA_SYMBOL(Convolve1DKernelType) kernel_type,
             size_t kernel_width,bool use_fft,LIBSAKURA_SYMBOL(Convole1DContext) **context);

/**
 * @brief Logical operation AND between two boolean arrays.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateLogicalAnd)(size_t num_in,
		bool const in1[/*num_in*/], bool const in2[/*num_in*/], bool out[/*num_in*/]);

/**
 * @brief Compute subtraction between two float arrays (in1 - in2).
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(OperateFloatSubtraction)(size_t num_in,
		float const in1[/*num_in*/], float const in2[/*num_in*/], float out[/*num_in*/]);

/**
 * @brief Compute values for Least-Square fitting from input data and a set of model data.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetLeastSquareMatrix)(size_t num_in,
		float const in_data[/*num_in*/], bool const in_mask[/*num_in*/],
		size_t num_model, double const model[/*num_model * num_in*/],
		double out[/*num_model * num_model*/], double out_vector[/*num_model*/]);

/**
 * @brief Solve simultaneous equations via LU decomposition.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(SolveSimultaneousEquationsByLU)(size_t num_eqn,
		double const lsq_matrix0[/*num_eqn * num_eqn*/],
		double const lsq_vector0[/*num_eqn*/], double out[/*num_eqn*/]);

/**
 * @brief Compute the best-fit model spectrum using model spectra and coefficients.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(DoGetBestFitModel)(size_t num_chan,
		size_t num_eqn, double const model[/*num_eqn * num_chan*/],
		double const coeff[/*num_eqn*/], float out[/*num_in*/]);

/**
 * @brief Compute the best-fit model spectrum by least-square fitting.
 */
LIBSAKURA_SYMBOL(Status) LIBSAKURA_SYMBOL(GetBestFitModel)(size_t num_in,
		float const in_data[/*num_in*/], bool const in_mask[/*num_in*/],
		size_t num_model, double const model[/*num_model * num_in*/],
		float out[/*num_in*/]);

#ifdef __cplusplus
}
/* extern "C" */
#endif

#endif /* LIBSAKURA_LIBSAKURA_SAKURA_H_ */
