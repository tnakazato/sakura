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
#include <sys/time.h>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <memory>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <complex>
#include <functional>

#include <libsakura/sakura.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"

#ifdef DISABLE_ALIGNAS
#undef SIMD_ALIGN
#define SIMD_ALIGN
#endif

//#define AVX __attribute__((aligned (32)))
#define AVX

#define ELEMENTSOF(x) (sizeof(x) / sizeof((x)[0]))
#define STATIC_ASSERT(x) static_assert((x), # x)

using namespace std;

namespace {
inline double CurrentTime() {
	return LIBSAKURA_SYMBOL(GetCurrentTime)();
}

template<typename T>
class Finalizer {
	T f;
public:
	Finalizer(T func) :
			f(func) {
	}
	~Finalizer() {
		f();
	}
};

template<typename T>
T square(T n) {
	return n * n;
}

typedef int32_t integer;

constexpr double Sqrt_(double a, double x, unsigned n) {
	return ((n > 0 ? Sqrt_(a, x, n - 1) : a)
			+ a / (n > 0 ? Sqrt_(a, x, n - 1) : a)) / 2.;
}

constexpr double Sqrt(double v) {
	return Sqrt_(v, v, 5);
}

constexpr int64_t Ceil(double v) {
	return static_cast<int64_t>(v) >= v ?
			static_cast<int64_t>(v) : static_cast<int64_t>(v) + 1;
}

template<typename T>
inline bool NearlyEqual(T aa, T bb) {
	T threshold = 0.00007;
	//T threshold = 0.0000003;
	return aa == bb || (abs(aa - bb) / (abs(aa) + abs(bb))) < threshold;
}

template<size_t CYCLE>
struct InitFuncs {
	static bool Mask1D(size_t i) {
		return i % CYCLE == 0;
	}
	static bool Mask2D(size_t i, size_t j) {
		return i % CYCLE == 0 || j % CYCLE == 0;
	}
	static size_t Map(size_t i, size_t max_size) {
		return i % max_size;
	}
};

template<size_t CYCLE>
struct InitFuncs2 {
	static bool Mask1D(size_t i) {
		return i % CYCLE == 0;
	}
	static bool Mask2D(size_t i, size_t j) {
		return i % CYCLE == 0 || j % CYCLE == 0;
	}
	static size_t Map(size_t i, size_t max) {
		return 0;
	}
};

template<size_t NROW, size_t ROW_FACTOR, size_t NVISCHAN, size_t NVISPOL>
struct SIMD_ALIGN RowBase {
	static constexpr size_t kNROW = NROW;
	static constexpr size_t kROW_FACTOR = ROW_FACTOR;
	static constexpr size_t kNVISCHAN = NVISCHAN;
	static constexpr size_t kNVISPOL = NVISPOL;

	SIMD_ALIGN
	float values[NROW][NVISPOL][NVISCHAN]; //
	SIMD_ALIGN
	double x[NROW]; //
	SIMD_ALIGN
	double y[NROW]; //
	SIMD_ALIGN bool sp_mask[NROW]; //
	SIMD_ALIGN bool mask[NROW][NVISPOL][NVISCHAN]; //
	SIMD_ALIGN
	float weight[NROW][NVISCHAN];
	void Init() {
		{
			auto a = &values;
			auto value = 0.0f;
			for (size_t i = 0; i < ELEMENTSOF(*a); ++i) {
				for (size_t j = 0; j < ELEMENTSOF((*a)[0]); ++j) {
					for (size_t k = 0; k < ELEMENTSOF((*a)[0][0]); ++k) {
						(*a)[i][j][k] = value;
					}
				}
			}
		}
		{
			auto a = &x;
			auto b = &y;
			auto value = 0.0;
			for (size_t i = 0; i < ELEMENTSOF(*a); ++i) {
				(*a)[i] = value;
				(*b)[i] = value;
			}
		}
		{
			auto a = &sp_mask;
			auto value = true;
			for (size_t i = 0; i < ELEMENTSOF(*a); ++i) {
				(*a)[i] = value;
			}
		}
		{
			auto a = &mask;
			auto value = true;
			for (size_t i = 0; i < ELEMENTSOF(*a); ++i) {
				for (size_t j = 0; j < ELEMENTSOF((*a)[0]); ++j) {
					for (size_t k = 0; k < ELEMENTSOF((*a)[0][0]); ++k) {
						(*a)[i][j][k] = value;
					}
				}
			}
		}
		{
			auto a = &weight;
			auto value = 1.0f;
			for (size_t i = 0; i < ELEMENTSOF(*a); ++i) {
				for (size_t j = 0; j < ELEMENTSOF((*a)[0]); ++j) {
					(*a)[i][j] = value;
				}
			}
		}
	}

	void SetRandomXY(double xmin, double xmax, double ymin, double ymax) {
		for (size_t i = 0; i < ELEMENTSOF(x); ++i) {
			x[i] = xmin + (xmax- xmin) * drand48();
			y[i] = ymin + (ymax- ymin) * drand48();
		}
	}

	void SetRandomXYInt(double xmin, double xmax, double ymin, double ymax) {
		for (size_t i = 0; i < ELEMENTSOF(x); ++i) {
			x[i] = int(xmin + (xmax- xmin) * drand48());
			y[i] = int(ymin + (ymax- ymin) * drand48());
		}
	}

	void SetRandomValues(float min, float max) {
		for (size_t i = 0; i < ELEMENTSOF(values); ++i) {
			for (size_t j = 0; j < ELEMENTSOF(values[0]); ++j) {
				for (size_t k = 0; k < ELEMENTSOF(values[0][0]); ++k) {
					values[i][j][k] = min + (max-min) * drand48();
				}
			}
		}
	}
};

template<size_t NVISCHAN, size_t NVISPOL, size_t NCHAN, size_t NPOL,
		size_t SAMPLING, size_t SUPPORT, size_t NX, size_t NY,
		typename InitFuncs>
struct SIMD_ALIGN GridBase {
	static constexpr size_t kNVISCHAN = NVISCHAN;
	static constexpr size_t kNVISPOL = NVISPOL;
	static constexpr size_t kNCHAN = NCHAN;
	static constexpr size_t kNPOL = NPOL;
	static constexpr size_t kNX = NX;
	static constexpr size_t kNY = NY;
	static constexpr size_t kSAMPLING = SAMPLING;
	static constexpr size_t kSUPPORT = SUPPORT;
	static constexpr size_t kCONV_TABLE_SIZE = (size_t(
			Ceil(Sqrt(2.) * ((SUPPORT + 1) * SAMPLING))));

	// inputs
	SIMD_ALIGN
	float convTab[kCONV_TABLE_SIZE]; //
	SIMD_ALIGN
	uint32_t chanmap[NVISCHAN]; //
	SIMD_ALIGN
	uint32_t polmap[NVISPOL];

// inputs/outputs
	typedef float GridType[NY][NX][NPOL][NCHAN];
	typedef float GridTypeRef[NPOL][NCHAN][NY][NX];
	typedef double SumType[NPOL][NCHAN];

	SIMD_ALIGN
	GridType grid; //
	SIMD_ALIGN
	GridType wgrid; //
	SIMD_ALIGN
	SumType sumwt;

protected:
	void ResetMap() {
		for (size_t i = 0; i < ELEMENTSOF(chanmap); ++i) {
			chanmap[i] = InitFuncs::Map(i, NCHAN);
		}
		for (size_t i = 0; i < ELEMENTSOF(polmap); ++i) {
			polmap[i] = InitFuncs::Map(i, NPOL);
		}
	}

	bool CmpGrid(char const msg[], GridType *a,
			GridType *b) {
		cout << "comparing " << msg << "\n";
		size_t count = 0;
		size_t count2 = 0;
		bool differ = false;
		for (size_t i = 0; i < ELEMENTSOF(*a); ++i) {
			for (size_t j = 0; j < ELEMENTSOF((*a)[0]); ++j) {
				for (size_t k = 0; k < ELEMENTSOF((*a)[0][0]); ++k) {
					for (size_t l = 0; l < ELEMENTSOF((*a)[0][0][0]); ++l) {
						float aa = (*a)[i][j][k][l];
						float bb = (*b)[i][j][k][l];
						if (NearlyEqual<float>(aa, bb)) {
							//if (bb != 0.) cout << aa << ", " << bb << endl;
							if (differ) {
								cout << "... just before [" << i << "][" << j
										<< "][" << k << "][" << l << "]\n";
							}
							count = 0;
							differ = false;
						} else {
							if (count < 20) {
								cout << "[" << i << "][" << j << "][" << k
										<< "][" << l << "]: "
										<< setprecision(10) << aa << " <> "
										<< bb << endl;
								differ = true;
							}
							count++;
							count2++;
							if (count2 > 200) {
								goto exit;
							}
						}
					}
				}
			}
		}
		exit: if (differ) {
			cout << "... END\n";
		}
		return count2 > 0 ? false : true;
	}

	bool CmpSumwt(SumType *a, SumType *b) {
		cout << "comparing sumwt\n";
		size_t count = 0;
		size_t count2 = 0;
		bool differ = false;
		for (size_t i = 0; i < ELEMENTSOF(*a); ++i) {
			for (size_t j = 0; j < ELEMENTSOF((*a)[0]); ++j) {
				double aa = (*a)[i][j];
				double bb = (*b)[i][j];
				if (NearlyEqual<double>(aa, bb)) {
					if (differ) {
						cout << "... just before [" << i << "][" << j << "]\n";
					}
					count = 0;
					differ = false;
				} else {
					if (count < 20) {
						cout << "[" << i << "][" << j << "]: "
								<< setprecision(10) << aa << " <> " << bb
								<< endl;
						differ = true;
					}
					count++;
					count2++;
				}
			}
		}
		if (differ) {
			cout << "... END\n";
		}
		return count2 > 0 ? false : true;
	}

	void InitTabs() {
#if 0
		static float ct[] = { 2.0f, 1.0f, 0.f };
		SetConvTab(ELEMENTSOF(ct), ct);
#else
		for (size_t i = 0; i < ELEMENTSOF(convTab); ++i) {
			convTab[i] = exp(-square(2.0 * i / ELEMENTSOF(convTab)));
		}
#endif
		ResetMap();
		ClearResults();
	}

public:
	void SetUp() {
		LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(
				Initialize)(nullptr, nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
		srand48(0);
		assert(LIBSAKURA_SYMBOL(IsAligned(this)));
	}

	void TearDown() {
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	void SetConvTab(size_t num, float src[]) {
		auto a = &convTab;
		assert(num <= ELEMENTSOF(*a));
		size_t i;
		for (i = 0; i < num; ++i) {
			(*a)[i] = src[i];
		}
		for (; i < ELEMENTSOF(*a); ++i) {
			(*a)[i] = 0;
		}
	}

	template<typename RT>
	void TestInvalidArg() {
		STATIC_ASSERT(RT::kNVISCHAN == kNVISCHAN);
		STATIC_ASSERT(RT::kNVISPOL == kNVISPOL);

		SIMD_ALIGN RT rows;
		// OK
		LIBSAKURA_SYMBOL(
				Status) result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW,
				rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL,
				polmap, RT::kNVISCHAN, chanmap, rows.mask[0][0],
				rows.values[0][0], rows.weight[0], false, ELEMENTSOF(convTab),
				convTab, NPOL, NCHAN, NX, NY, sumwt[0], wgrid[0][0][0],
				grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

		result = sakura_GridConvolving(0, 0, 0, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], true, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

		{ // nullptr
		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, nullptr,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				nullptr, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, nullptr, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, nullptr,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, nullptr, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, nullptr, rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], nullptr,
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], nullptr,
				rows.weight[0], true, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				nullptr, false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), nullptr, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, nullptr, wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], nullptr, grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
	}

		{ // alignment
		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask+1,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x+1, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y+1, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap+1,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap+1, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, &rows.mask[0][0][1], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], &rows.values[0][0][1],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], &rows.values[0][0][1],
				rows.weight[0], true, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				&rows.weight[0][1], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab+1, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, &sumwt[0][1], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], &wgrid[0][0][0][1], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], &grid[0][0][0][1]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
	}

		{ // out of range
			constexpr uint32_t int32max = INT32_MAX;
			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW + 1, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, 0, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, (int32max - 1) / 2 + 1, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, 0, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, int32max+1, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, 0, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, int32max+1, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					0, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					int32max+1, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, 0,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, int32max+1,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					0, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					int32max+1, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, 0, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, int32max+1, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, NX, 0, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, NX, int32max+1, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, 0, convTab, NPOL,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab)+1, convTab, NPOL,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab)-1, convTab, NPOL,
					NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
					rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
					RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
					NCHAN, int32max / 2, 3, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		}
#if 0
		result = sakura_GridConvolving(RT::kNROW, 0, RT::kNROW, rows.sp_mask,
				rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVISPOL, polmap,
				RT::kNVISCHAN, chanmap, rows.mask[0][0], rows.values[0][0],
				rows.weight[0], false, ELEMENTSOF(convTab), convTab, NPOL,
				NCHAN, NX, NY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
#endif

	}

	template<typename T, typename F>
	static void Dump(T *a, F func, int width, int prec) {
		auto rowSep = "";
		for (size_t i = 0; i < ELEMENTSOF(*a); ++i) {
			cout << rowSep;
			auto colSep = "{";
			for (size_t j = 0; j < ELEMENTSOF((*a)[0]); ++j) {
				cout << colSep << setw(width) << fixed << setprecision(prec)
						<< func(a, i, j);
				colSep = ", ";
			}
			cout << "}";
			rowSep = ",\n";
		}
		cout << endl;
	}
	static void Dump(SumType *a, int width, int prec = 1) {
		Dump(a, [](SumType *a, size_t i, size_t j)->double {
			return (*a)[i][j];
		}, width, prec);
	}
	static void Dump(GridType *a, int width, int prec = 1) {
		for (size_t pol = 0; pol < ELEMENTSOF((*a)[0][0]); ++pol) {
			for (size_t ch = 0; ch < ELEMENTSOF((*a)[0][0][0]); ++ch) {
				cout << "(pol, ch) = (" << pol << ", " << ch << ")\n";
				Dump(a, [=](GridType *a, size_t i, size_t j)->float {
					return (*a)[i][j][pol][ch];
				}, width, prec);
			}
		}
	}

	static void Transform(GridTypeRef *src, GridType *dst) {
		for (size_t i = 0; i < ELEMENTSOF(*src); ++i) {
			for (size_t j = 0; j < ELEMENTSOF((*src)[0]); ++j) {
				for (size_t k = 0; k < ELEMENTSOF((*src)[0][0]); ++k) {
					for (size_t l = 0; l < ELEMENTSOF((*src)[0][0][0]); ++l) {
						(*dst)[k][l][i][j] = (*src)[i][j][k][l];
					}
				}
			}
		}
	}

	static void ClearResults(SumType *sumwt, GridType *wgrid, GridType *grid) {
		memset(grid, 0, sizeof(*grid));
		memset(wgrid, 0, sizeof(*wgrid));
		memset(sumwt, 0, sizeof(*sumwt));
	}

	void ClearResults() {
		ClearResults(&sumwt, &wgrid, &grid);
	}

	template<typename RT>
	void TrySpeed(bool do_verify, bool weightOnly, RT *rows,
			double (*sumwt2)[NPOL][NCHAN],
			float (*wgrid2)[NY][NX][NPOL][NCHAN],
			float (*grid2)[NY][NX][NPOL][NCHAN]) {
		STATIC_ASSERT(RT::kNVISCHAN == kNVISCHAN);
		STATIC_ASSERT(RT::kNVISPOL == kNVISPOL);
		{
			cout << "Gridding by C++ ... " << flush;
			double start = CurrentTime();
			for (size_t i = 0; i < RT::kROW_FACTOR; ++i) {
				LIBSAKURA_SYMBOL(
						Status) result = sakura_GridConvolving(RT::kNROW, 0,
						RT::kNROW, rows->sp_mask, rows->x, rows->y, SUPPORT,
						SAMPLING, RT::kNVISPOL, polmap,
						RT::kNVISCHAN, chanmap, rows->mask[0][0],
						rows->values[0][0], rows->weight[0], weightOnly,
						ELEMENTSOF(convTab), convTab, NPOL, NCHAN, NX, NY,
						sumwt[0], wgrid[0][0][0], grid[0][0][0]);
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
			}
			double end = CurrentTime();
			double new_time = end - start;
			cout << "done. " << new_time << " sec\n";
		}
		if (do_verify) {
			EXPECT_TRUE(CmpSumwt(&sumwt, sumwt2));
			EXPECT_TRUE(CmpGrid("wgrid", &wgrid, wgrid2));
			EXPECT_TRUE(CmpGrid("grid", &grid, grid2));
		}
	}

	template<typename RT>
	void TestGrid(
			std::function<
					void(GridBase *, RT *, SumType *sumwt, GridType *wgrid,
							GridType *grid)> func) {
		STATIC_ASSERT(RT::kNVISCHAN == kNVISCHAN);
		STATIC_ASSERT(RT::kNVISPOL == kNVISPOL);
		InitTabs();

		RT *rows = nullptr;
		unique_ptr<void, DefaultAlignedMemory> rowsMem(
				DefaultAlignedMemory::AlignedAllocateOrException(sizeof(RT),
						&rows));
		rows->Init();

		GridType *wgrid2 = nullptr;
		unique_ptr<void, DefaultAlignedMemory> wgrid2Mem(
				DefaultAlignedMemory::AlignedAllocateOrException(
						sizeof(*wgrid2), &wgrid2));
		// init wgrid2;

		GridType *grid2 = nullptr;
		unique_ptr<void, DefaultAlignedMemory> grid2Mem(
				DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*grid2),
						&grid2));

		SumType *sumwt2 = nullptr;
		unique_ptr<void, DefaultAlignedMemory> sumwt2Mem(
				DefaultAlignedMemory::AlignedAllocateOrException(
						sizeof(*sumwt2), &sumwt2));

		func(this, rows, sumwt2, wgrid2, grid2);
	}
};

typedef GridBase<1, 1, 1, 1, 2, 2, 7, 7, InitFuncs<4> > TestMinimal;
typedef GridBase<1024, 4, 1024, 4, 100, 10, 200, 180, InitFuncs<1024> > TestTypical;

#if 0
typedef TestBase<512, 2, 1024, 4, 1024, 4, 100, 10, 200, 180, InitFuncs<2> > TestCase1;
typedef TestBase<512 - 1, 10, 1024 - 1, 4 - 1, 1024 - 1, 4, 100, 10, 100, 90,
InitFuncs<7> > TestCase2;
typedef TestBase<512, 1, 512, 4, 1024, 4, 100, 10, 200, 180, InitFuncs2<1> > TestNto1;
typedef TestBase<512, 10000, 3, 1, 3, 1, 100, 10, 3, 3, InitFuncs<1> > TestSmall;
#endif
}
			// namespace

TEST(Gridding, Basic) {
typedef TestMinimal TestCase;
typedef RowBase<1, 1, TestCase::kNVISCHAN, TestCase::kNVISPOL> RowType;

TestCase *test_case;
unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
test_case->SetUp();
auto finFunc = [=]() {test_case->TearDown();};
auto fin = Finalizer<decltype(finFunc)>(finFunc);

test_case->TestInvalidArg<RowType>();
test_case->TestGrid<RowType>(
		[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 8, 7, 6, 5, 4, 3, 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				bool weightOnly = true;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, wgrid2, grid2);
			}
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 8, 7, 6, 5, 4, 3, 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->x[0] = 3;
				rows->y[0] = 3;
				rows->values[0][0][0] = 2;
				{
					auto p = *sumwt2;
					p[0][0] = 116.f;
				}
				TestCase::GridType wgrid = {
					{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
					{	0.0, 3.0, 4.0, 4.0, 4.0, 3.0, 0.0},
					{	0.0, 4.0, 6.0, 6.0, 6.0, 4.0, 0.0},
					{	0.0, 4.0, 6.0, 8.0, 6.0, 4.0, 0.0},
					{	0.0, 4.0, 6.0, 6.0, 6.0, 4.0, 0.0},
					{	0.0, 3.0, 4.0, 4.0, 4.0, 3.0, 0.0},
					{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
				};
				TestCase::GridType grid = {
					{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
					{	0.0, 6.0, 8.0, 8.0, 8.0, 6.0, 0.0},
					{	0.0, 8.0, 12.0, 12.0, 12.0, 8.0, 0.0},
					{	0.0, 8.0, 12.0, 16.0, 12.0, 8.0, 0.0},
					{	0.0, 8.0, 12.0, 12.0, 12.0, 8.0, 0.0},
					{	0.0, 6.0, 8.0, 8.0, 8.0, 6.0, 0.0},
					{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
				};
				bool weightOnly = false;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, &wgrid, &grid);
				TestCase::Dump(&tc->sumwt, 5);
				TestCase::Dump(&tc->wgrid, 5);
				TestCase::Dump(&tc->grid, 5);
			}
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 8, 7, 6, 5, 4, 3, 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->x[0] = 2.5;
				rows->y[0] = 2.5;
				rows->values[0][0][0] = 2;
				{
					auto p = *sumwt2;
					p[0][0] = 109.f;
				}
				TestCase::GridType wgrid = {
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0},
						{  0.0,   4.0,   5.0,   5.0,   4.0,   3.0,   0.0},
						{  0.0,   5.0,   7.0,   7.0,   5.0,   3.0,   0.0},
						{  0.0,   5.0,   7.0,   7.0,   5.0,   3.0,   0.0},
						{  0.0,   4.0,   5.0,   5.0,   4.0,   3.0,   0.0},
						{  0.0,   3.0,   3.0,   3.0,   3.0,   1.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}
				};
				TestCase::GridType grid = {
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0},
						{  0.0,   8.0,  10.0,  10.0,   8.0,   6.0,   0.0},
						{  0.0,  10.0,  14.0,  14.0,  10.0,   6.0,   0.0},
						{  0.0,  10.0,  14.0,  14.0,  10.0,   6.0,   0.0},
						{  0.0,   8.0,  10.0,  10.0,   8.0,   6.0,   0.0},
						{  0.0,   6.0,   6.0,   6.0,   6.0,   2.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}
				};
				bool weightOnly = false;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, &wgrid, &grid);
				TestCase::Dump(&tc->sumwt, 5);
				TestCase::Dump(&tc->wgrid, 5);
				TestCase::Dump(&tc->grid, 5);
			}
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 8, 7, 6, 5, 4, 3, 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->x[0] = 2.49;
				rows->y[0] = 2.49;
				rows->values[0][0][0] = 2;
				{
					auto p = *sumwt2;
					p[0][0] = 109.f;
				}
				TestCase::GridType wgrid = {
						{  1.0,   3.0,   3.0,   3.0,   3.0,   0.0,   0.0},
						{  3.0,   4.0,   5.0,   5.0,   4.0,   0.0,   0.0},
						{  3.0,   5.0,   7.0,   7.0,   5.0,   0.0,   0.0},
						{  3.0,   5.0,   7.0,   7.0,   5.0,   0.0,   0.0},
						{  3.0,   4.0,   5.0,   5.0,   4.0,   0.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}
				};
				TestCase::GridType grid = {
						{  2.0,   6.0,   6.0,   6.0,   6.0,   0.0,   0.0},
						{  6.0,   8.0,  10.0,  10.0,   8.0,   0.0,   0.0},
						{  6.0,  10.0,  14.0,  14.0,  10.0,   0.0,   0.0},
						{  6.0,  10.0,  14.0,  14.0,  10.0,   0.0,   0.0},
						{  6.0,   8.0,  10.0,  10.0,   8.0,   0.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}
				};
				bool weightOnly = false;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, &wgrid, &grid);
				TestCase::Dump(&tc->sumwt, 5);
				TestCase::Dump(&tc->wgrid, 5);
				TestCase::Dump(&tc->grid, 5);
			}
			{ // incremental grid
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 8, 7, 6, 5, 4, 3, 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->x[0] = 3;
				rows->y[0] = 3;
				rows->values[0][0][0] = 2;
				{
					auto p = *sumwt2;
					p[0][0] = 225.f;
				}
				TestCase::GridType wgrid = {
						{  1.0,   3.0,   3.0,   3.0,   3.0,   0.0,   0.0},
						{  3.0,   7.0,   9.0,   9.0,   8.0,   3.0,   0.0},
						{  3.0,   9.0,  13.0,  13.0,  11.0,   4.0,   0.0},
						{  3.0,   9.0,  13.0,  15.0,  11.0,   4.0,   0.0},
						{  3.0,   8.0,  11.0,  11.0,  10.0,   4.0,   0.0},
						{  0.0,   3.0,   4.0,   4.0,   4.0,   3.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}
				};
				TestCase::GridType grid = {
						{  2.0,   6.0,   6.0,   6.0,   6.0,   0.0,   0.0},
						{  6.0,  14.0,  18.0,  18.0,  16.0,   6.0,   0.0},
						{  6.0,  18.0,  26.0,  26.0,  22.0,   8.0,   0.0},
						{  6.0,  18.0,  26.0,  30.0,  22.0,   8.0,   0.0},
						{  6.0,  16.0,  22.0,  22.0,  20.0,   8.0,   0.0},
						{  0.0,   6.0,   8.0,   8.0,   8.0,   6.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}
				};
				bool weightOnly = false;
				tc->TrySpeed(false, weightOnly, rows, sumwt2, &wgrid, &grid);
				rows->x[0] = 2.49;
				rows->y[0] = 2.49;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, &wgrid, &grid);
				TestCase::Dump(&tc->sumwt, 5);
				TestCase::Dump(&tc->wgrid, 5);
				TestCase::Dump(&tc->grid, 5);
			}
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 8, 7, 6, 5, 4, 3, 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->x[0] = 3;
				rows->y[0] = 3;
				rows->values[0][0][0] = 100;
				{
					auto p = *sumwt2;
					p[0][0] = 348.f;
				}
				TestCase::GridType wgrid = {
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0},
						{  0.0,   9.0,  12.0,  12.0,  12.0,   9.0,   0.0},
						{  0.0,  12.0,  18.0,  18.0,  18.0,  12.0,   0.0},
						{  0.0,  12.0,  18.0,  24.0,  18.0,  12.0,   0.0},
						{  0.0,  12.0,  18.0,  18.0,  18.0,  12.0,   0.0},
						{  0.0,   9.0,  12.0,  12.0,  12.0,   9.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}
				};
				bool weightOnly = true;
				tc->TrySpeed(false, weightOnly, rows, sumwt2, wgrid2, grid2);

				rows->values[0][0][0] = 0;
				weightOnly = false;
				tc->TrySpeed(false, weightOnly, rows, sumwt2, wgrid2, grid2);

				rows->values[0][0][0] = 2;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, &wgrid, &wgrid);
			}
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 8, 7, 6, 5, 4, 3, 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->x[0] = 3;
				rows->y[0] = 3;
				rows->values[0][0][0] = 100;
				rows->weight[0][0] = 0.5;
				{
					auto p = *sumwt2;
					p[0][0] = 58.f;
				}
				TestCase::GridType wgrid = {
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0},
						{  0.0,   1.5,   2.0,   2.0,   2.0,   1.5,   0.0},
						{  0.0,   2.0,   3.0,   3.0,   3.0,   2.0,   0.0},
						{  0.0,   2.0,   3.0,   4.0,   3.0,   2.0,   0.0},
						{  0.0,   2.0,   3.0,   3.0,   3.0,   2.0,   0.0},
						{  0.0,   1.5,   2.0,   2.0,   2.0,   1.5,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}
				};
				bool weightOnly = true;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, &wgrid, &wgrid);
			}
		});
}

TEST(Gridding, Clipping) {
typedef TestMinimal TestCase;
typedef RowBase<1, 1, TestCase::kNVISCHAN, TestCase::kNVISPOL> RowType;

TestCase *test_case;
unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
test_case->SetUp();
auto finFunc = [=]() {test_case->TearDown();};
auto fin = Finalizer<decltype(finFunc)>(finFunc);

test_case->TestGrid<RowType>(
		[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 8, 7, 6, 5, 4, 3, 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->x[0] = 1.5;
				rows->y[0] = 1.5;
				rows->values[0][0][0] = 2;
				{
					auto p = *sumwt2;
					p[0][0] = 109.f;
				}
				TestCase::GridType wgrid = {
						{  4.0,   5.0,   5.0,   4.0,   3.0,   0.0,   0.0},
						{  5.0,   7.0,   7.0,   5.0,   3.0,   0.0,   0.0},
						{  5.0,   7.0,   7.0,   5.0,   3.0,   0.0,   0.0},
						{  4.0,   5.0,   5.0,   4.0,   3.0,   0.0,   0.0},
						{  3.0,   3.0,   3.0,   3.0,   1.0,   0.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}
				};
				TestCase::GridType grid = {
						{  8.0,  10.0,  10.0,   8.0,   6.0,   0.0,   0.0},
						{ 10.0,  14.0,  14.0,  10.0,   6.0,   0.0,   0.0},
						{ 10.0,  14.0,  14.0,  10.0,   6.0,   0.0,   0.0},
						{  8.0,  10.0,  10.0,   8.0,   6.0,   0.0,   0.0},
						{  6.0,   6.0,   6.0,   6.0,   2.0,   0.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}
				};
				bool weightOnly = false;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, &wgrid, &grid);
			}
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 8, 7, 6, 5, 4, 3, 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->x[0] = 4.5;
				rows->y[0] = 4.5;
				rows->values[0][0][0] = 2;
				bool weightOnly = false;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, wgrid2, grid2);
			}
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 8, 7, 6, 5, 4, 3, 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->x[0] = 4.49;
				rows->y[0] = 4.49;
				rows->values[0][0][0] = 2;
				{
					auto p = *sumwt2;
					p[0][0] = 109.f;
				}
				TestCase::GridType wgrid = {
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0},
						{  0.0,   0.0,   1.0,   3.0,   3.0,   3.0,   3.0},
						{  0.0,   0.0,   3.0,   4.0,   5.0,   5.0,   4.0},
						{  0.0,   0.0,   3.0,   5.0,   7.0,   7.0,   5.0},
						{  0.0,   0.0,   3.0,   5.0,   7.0,   7.0,   5.0},
						{  0.0,   0.0,   3.0,   4.0,   5.0,   5.0,   4.0}
				};
				TestCase::GridType grid = {
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0},
						{  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0},
						{  0.0,   0.0,   2.0,   6.0,   6.0,   6.0,   6.0},
						{  0.0,   0.0,   6.0,   8.0,  10.0,  10.0,   8.0},
						{  0.0,   0.0,   6.0,  10.0,  14.0,  14.0,  10.0},
						{  0.0,   0.0,   6.0,  10.0,  14.0,  14.0,  10.0},
						{  0.0,   0.0,   6.0,   8.0,  10.0,  10.0,   8.0}
				};
				bool weightOnly = false;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, &wgrid, &grid);
				TestCase::Dump(&tc->sumwt, 5);
				TestCase::Dump(&tc->wgrid, 5);
				TestCase::Dump(&tc->grid, 5);
			}
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 8, 7, 6, 5, 4, 3, 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->values[0][0][0] = 2;
				rows->x[0] = 1.49;
				rows->y[0] = 1.5;
				bool weightOnly = false;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, wgrid2, grid2);

				rows->x[0] = 1.5;
				rows->y[0] = 1.49;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, wgrid2, grid2);

				rows->x[0] = 4.49;
				rows->y[0] = 4.5;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, wgrid2, grid2);

				rows->x[0] = 4.5;
				rows->y[0] = 4.49;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, wgrid2, grid2);
			}
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 8, 7, 6, 5, 4, 3, 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->x[0] = -1;
				rows->y[0] = 3;
				rows->values[0][0][0] = 2;
				bool weightOnly = false;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, wgrid2, grid2);

				rows->x[0] = 8;
				rows->y[0] = 3;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, wgrid2, grid2);

				rows->x[0] = 3;
				rows->y[0] = -1;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, wgrid2, grid2);

				rows->x[0] = 3;
				rows->y[0] = 8;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, wgrid2, grid2);
			}
		});
}

TEST(Gridding, ChMapping) {
typedef GridBase<2, 1, 2, 1, 1, 1, 5, 5, InitFuncs<4> > TestCase;
typedef RowBase<1, 1, TestCase::kNVISCHAN, TestCase::kNVISPOL> RowType;

TestCase *test_case;
unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
test_case->SetUp();
auto finFunc = [=]() {test_case->TearDown();};
auto fin = Finalizer<decltype(finFunc)>(finFunc);

test_case->TestInvalidArg<RowType>();
test_case->TestGrid<RowType>(
		[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->x[0] = 2;
				rows->y[0] = 2;
				rows->values[0][0][0] = 2;
				rows->values[0][0][1] = 3;
				rows->weight[0][1] = 0.5;
				{
					auto p = *sumwt2;
					p[0][0] = 10;
					p[0][1] = 5;
				}
				TestCase::GridTypeRef rwgrid = { {
						{
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 1.0, 1.0, 1.0, 0.0},
							{	0.0, 1.0, 2.0, 1.0, 0.0},
							{	0.0, 1.0, 1.0, 1.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						},
						{	{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.5, 0.5, 0.5, 0.0},
							{	0.0, 0.5, 1.0, 0.5, 0.0},
							{	0.0, 0.5, 0.5, 0.5, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						}
					}
				};
				TestCase::GridTypeRef rgrid = { {
						{
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 2.0, 2.0, 2.0, 0.0},
							{	0.0, 2.0, 4.0, 2.0, 0.0},
							{	0.0, 2.0, 2.0, 2.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						},
						{
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 1.5, 1.5, 1.5, 0.0},
							{	0.0, 1.5, 3.0, 1.5, 0.0},
							{	0.0, 1.5, 1.5, 1.5, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						}
					}
				};
				TestCase::GridType wgrid;
				TestCase::GridType grid;
				TestCase::Transform(&rwgrid, &wgrid);
				TestCase::Transform(&rgrid, &grid);
				bool weightOnly = false;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, &wgrid, &grid);
			}
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				tc->chanmap[0] = 1;
				tc->chanmap[1] = 1;
				rows->x[0] = 2;
				rows->y[0] = 2;
				rows->values[0][0][0] = 2;
				rows->values[0][0][1] = 3;
				rows->weight[0][1] = 0.5;
				{
					auto p = *sumwt2;
					p[0][0] = 0;
					p[0][1] = 15;
				}
				TestCase::GridTypeRef rwgrid = { { {
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						},
						{
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 1.5, 1.5, 1.5, 0.0},
							{	0.0, 1.5, 3.0, 1.5, 0.0},
							{	0.0, 1.5, 1.5, 1.5, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						}
					}
				};
				TestCase::GridTypeRef rgrid = { { {
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						},
						{
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 3.5, 3.5, 3.5, 0.0},
							{	0.0, 3.5, 7.0, 3.5, 0.0},
							{	0.0, 3.5, 3.5, 3.5, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						}
					}
				};
				TestCase::GridType wgrid;
				TestCase::GridType grid;
				TestCase::Transform(&rwgrid, &wgrid);
				TestCase::Transform(&rgrid, &grid);
				bool weightOnly = false;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, &wgrid, &grid);
			}
});
}

TEST(Gridding, PolMapping) {
typedef GridBase<1, 2, 1, 2, 1, 1, 5, 5, InitFuncs<4> > TestCase;
typedef RowBase<1, 1, TestCase::kNVISCHAN, TestCase::kNVISPOL> RowType;

TestCase *test_case;
unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
test_case->SetUp();
auto finFunc = [=]() {test_case->TearDown();};
auto fin = Finalizer<decltype(finFunc)>(finFunc);

test_case->TestInvalidArg<RowType>();
test_case->TestGrid<RowType>(
		[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				rows->x[0] = 2;
				rows->y[0] = 2;
				rows->values[0][0][0] = 2;
				rows->values[0][1][0] = 3;
				{
					auto p = *sumwt2;
					p[0][0] = 10;
					p[1][0] = 10;
				}
				TestCase::GridTypeRef rwgrid = { { {
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 1.0, 1.0, 1.0, 0.0},
							{	0.0, 1.0, 2.0, 1.0, 0.0},
							{	0.0, 1.0, 1.0, 1.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						}}, { {
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 1.0, 1.0, 1.0, 0.0},
							{	0.0, 1.0, 2.0, 1.0, 0.0},
							{	0.0, 1.0, 1.0, 1.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						}}
				};
				TestCase::GridTypeRef rgrid = { { {
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 2.0, 2.0, 2.0, 0.0},
							{	0.0, 2.0, 4.0, 2.0, 0.0},
							{	0.0, 2.0, 2.0, 2.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						}}, { {
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 3.0, 3.0, 3.0, 0.0},
							{	0.0, 3.0, 6.0, 3.0, 0.0},
							{	0.0, 3.0, 3.0, 3.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						}}
				};
				TestCase::GridType wgrid;
				TestCase::GridType grid;
				TestCase::Transform(&rwgrid, &wgrid);
				TestCase::Transform(&rgrid, &grid);
				bool weightOnly = false;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, &wgrid, &grid);
			}
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				static float ct[] = { 2, 1, -.9f };
				tc->SetConvTab(ELEMENTSOF(ct), ct);
				tc->polmap[0] = 1;
				tc->polmap[1] = 1;
				rows->x[0] = 2;
				rows->y[0] = 2;
				rows->values[0][0][0] = 2;
				rows->values[0][1][0] = 3;
				rows->weight[0][0] = 0.5;
				{
					auto p = *sumwt2;
					p[0][0] = 0;
					p[0][1] = 10;
				}
				TestCase::GridTypeRef rwgrid = { { {
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						}}, { {
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 1.0, 1.0, 1.0, 0.0},
							{	0.0, 1.0, 2.0, 1.0, 0.0},
							{	0.0, 1.0, 1.0, 1.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						}}
				};
				TestCase::GridTypeRef rgrid = { { {
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						}}, { {
							{	0.0, 0.0, 0.0, 0.0, 0.0},
							{	0.0, 2.5, 2.5, 2.5, 0.0},
							{	0.0, 2.5, 5.0, 2.5, 0.0},
							{	0.0, 2.5, 2.5, 2.5, 0.0},
							{	0.0, 0.0, 0.0, 0.0, 0.0}
						}}
				};
				TestCase::GridType wgrid;
				TestCase::GridType grid;
				TestCase::Transform(&rwgrid, &wgrid);
				TestCase::Transform(&rgrid, &grid);
				bool weightOnly = false;
				tc->TrySpeed(true, weightOnly, rows, sumwt2, &wgrid, &grid);
			}
});
}

TEST(Gridding, Typical) {
typedef TestTypical TestCase;
typedef RowBase<512, 1, TestCase::kNVISCHAN, TestCase::kNVISPOL> RowType;

TestCase *test_case;
unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
test_case->SetUp();
auto finFunc = [=]() {test_case->TearDown();};
auto fin = Finalizer<decltype(finFunc)>(finFunc);

test_case->TestInvalidArg<RowType>();
test_case->TestGrid<RowType>(
		[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
			{
				TestCase::ClearResults(sumwt2, wgrid2, grid2);
				tc->ClearResults();
				rows->SetRandomXYInt(2 * TestCase::kSUPPORT, TestCase::kNX - 2 * TestCase::kSUPPORT,
						2 * TestCase::kSUPPORT, TestCase::kNY - 2 * TestCase::kSUPPORT);
				bool weightOnly = false;
				tc->TrySpeed(false, weightOnly, rows, sumwt2, wgrid2, grid2);
				//TestCase::Dump(&tc->sumwt, 5);
			}
		});
}

#if 0
TEST(Gridding, Generic) {
typedef TestCase1 TestCase;
TestCase *test_case;
unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
test_case->SetUp();
test_case->TearDown();
}

TEST(Gridding, Odd) {
typedef TestCase2 TestCase;
TestCase *test_case;
unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
test_case->SetUp();
test_case->TearDown();
}

TEST(Gridding, Small) {
typedef TestSmall TestCase;
TestCase *test_case;
unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
test_case->SetUp();
test_case->TearDown();
}
#endif
