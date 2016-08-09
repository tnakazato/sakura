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
#include <sys/time.h>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <memory>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <complex>
#include <functional>

#include <libsakura/sakura.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"
#include "testutil.h"

#define ELEMENTSOF(x) (sizeof(x) / sizeof((x)[0]))
#define STATIC_ASSERT(x) static_assert((x), # x)

#if defined(DISABLE_ALIGNAS)
# define SIMD_ALIGN /* nothing */
#endif

using namespace std;

namespace {
inline double CurrentTime() {
	return GetCurrentTime();
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
T Square(T n) {
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

template<size_t kNRowArg, size_t kRowFactorArg, size_t kNVisChanArg, size_t kNVisPolArg>
struct SIMD_ALIGN RowBase {
	static constexpr size_t kNRow = kNRowArg;
	static constexpr size_t kRowFactor = kRowFactorArg;
	static constexpr size_t kNVisChan = kNVisChanArg;
	static constexpr size_t kNVisPol = kNVisPolArg;

	SIMD_ALIGN
	float values[kNRow][kNVisPol][kNVisChan]; //
	SIMD_ALIGN
	double x[kNRow]; //
	SIMD_ALIGN
	double y[kNRow]; //
	SIMD_ALIGN bool sp_mask[kNRow]; //
	SIMD_ALIGN bool mask[kNRow][kNVisPol][kNVisChan]; //
	SIMD_ALIGN
	float weight[kNRow][kNVisChan];
	void Init() {
		SetValues(0.f);
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
		SetWeight(1.f);
	}

	void SetRandomXY(double xmin, double xmax, double ymin, double ymax) {
		for (size_t i = 0; i < ELEMENTSOF(x); ++i) {
			x[i] = xmin + (xmax - xmin) * drand48();
			y[i] = ymin + (ymax - ymin) * drand48();
		}
	}

	void SetRandomXYInt(double xmin, double xmax, double ymin, double ymax) {
		for (size_t i = 0; i < ELEMENTSOF(x); ++i) {
			x[i] = int(xmin + (xmax - xmin) * drand48());
			y[i] = int(ymin + (ymax - ymin) * drand48());
		}
	}

	void SetXY(double xx, double yy) {
		for (size_t i = 0; i < ELEMENTSOF(x); ++i) {
			x[i] = xx;
			y[i] = yy;
		}
	}

	void SetRandomValues(float min, float max) {
		for (size_t i = 0; i < ELEMENTSOF(values); ++i) {
			for (size_t j = 0; j < ELEMENTSOF(values[0]); ++j) {
				for (size_t k = 0; k < ELEMENTSOF(values[0][0]); ++k) {
					values[i][j][k] = min + (max - min) * drand48();
				}
			}
		}
	}
	void SetValues(float value) {
		for (size_t i = 0; i < ELEMENTSOF(values); ++i) {
			for (size_t j = 0; j < ELEMENTSOF(values[0]); ++j) {
				for (size_t k = 0; k < ELEMENTSOF(values[0][0]); ++k) {
					values[i][j][k] = value;
				}
			}
		}
	}

	void SetWeight(float value) {
		for (size_t i = 0; i < ELEMENTSOF(weight); ++i) {
			for (size_t j = 0; j < ELEMENTSOF(weight[0]); ++j) {
				weight[i][j] = value;
			}
		}
	}
};

template<typename MessageType>
void ReportBenchmark(MessageType const &key, double sec) {
	std::cout << std::setprecision(5) << "#x# benchmark " << key << " " << sec
			<< std::endl;
}

template<size_t kNVisChanArg, size_t kNVisPolArg, size_t kNChanArg, size_t kNPolArg,
		size_t SAMPLING, size_t SUPPORT, size_t kNXArg, size_t kNYArg,
		typename InitFuncs>
struct SIMD_ALIGN GridBase {
	static constexpr size_t kNVisChan = kNVisChanArg;
	static constexpr size_t kNVisPol = kNVisPolArg;
	static constexpr size_t kNChan = kNChanArg;
	static constexpr size_t kNPol = kNPolArg;
	static constexpr size_t kNX = kNXArg;
	static constexpr size_t kNY = kNYArg;
	static constexpr size_t kSampling = SAMPLING;
	static constexpr size_t kSupport = SUPPORT;
	static constexpr size_t kConvTableSize = (size_t(
			Ceil(Sqrt(2.) * ((SUPPORT + 1) * SAMPLING))));

	// inputs
	SIMD_ALIGN
	float conv_tab[kConvTableSize]; //
	SIMD_ALIGN
	uint32_t chanmap[kNVisChan]; //
	SIMD_ALIGN
	uint32_t polmap[kNVisPol];

	// inputs/outputs
	typedef float GridType[kNY][kNX][kNPol][kNChan];
	typedef float GridTypeRef[kNPol][kNChan][kNY][kNX];
	typedef double SumType[kNPol][kNChan];

	SIMD_ALIGN
	GridType grid; //
	SIMD_ALIGN
	GridType wgrid; //
	SIMD_ALIGN
	SumType sumwt;

protected:
	void ResetMap() {
		for (size_t i = 0; i < ELEMENTSOF(chanmap); ++i) {
			chanmap[i] = InitFuncs::Map(i, kNChan);
		}
		for (size_t i = 0; i < ELEMENTSOF(polmap); ++i) {
			polmap[i] = InitFuncs::Map(i, kNPol);
		}
	}

	template<typename GridType>
	bool CmpGrid(char const msg[], GridType *a, GridType *b) {
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
							++count;
							++count2;
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

	template<typename SumType>
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
					++count;
					++count2;
				}
			}
		}
		if (differ) {
			cout << "... END\n";
		}
		return count2 > 0 ? false : true;
	}

	void InitTabs() {
		InitConvTab();
		ResetMap();
		ClearResults();
	}

public:
	void SetUp() {
		LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
				nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
		srand48(0);
		assert(LIBSAKURA_SYMBOL(IsAligned(this)));
		assert(
				kConvTableSize
						== static_cast<size_t>(ceil(
								sqrt(2.) * (kSupport + 1) * kSampling)));
	}

	void TearDown() {
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	void SetConvTab(size_t num, float src[]) {
		auto a = &conv_tab;
		assert(num <= ELEMENTSOF(*a));
		size_t i;
		for (i = 0; i < num; ++i) {
			(*a)[i] = src[i];
		}
		for (; i < ELEMENTSOF(*a); ++i) {
			(*a)[i] = 0;
		}
	}

	void InitConvTab() {
		for (size_t i = 0; i < ELEMENTSOF(conv_tab); ++i) {
			conv_tab[i] = exp(-Square(2.0 * i / ELEMENTSOF(conv_tab)));
		}
	}

	template<typename RT>
	void TestInvalidArg() {
		STATIC_ASSERT(RT::kNVisChan == kNVisChan);
		STATIC_ASSERT(RT::kNVisPol == kNVisPol);

		InitTabs();
		RT *rows_ptr;
		unique_ptr<void, DefaultAlignedMemory> rows_storage(
				DefaultAlignedMemory::AlignedAllocateOrException(
						sizeof(*rows_ptr), &rows_ptr));
		rows_ptr->Init();
		SIMD_ALIGN RT &rows = *rows_ptr;
		// OK
		LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow,
				0, RT::kNRow, rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
				RT::kNVisPol, polmap, RT::kNVisChan, chanmap, rows.mask[0][0],
				rows.values[0][0], rows.weight[0], false, ELEMENTSOF(conv_tab),
				conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0], wgrid[0][0][0],
				grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

		result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(0, 0, 0, rows.sp_mask, rows.x,
				rows.y, SUPPORT, SAMPLING, RT::kNVisPol, polmap, RT::kNVisChan,
				chanmap, rows.mask[0][0], rows.values[0][0], rows.weight[0],
				true, ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY,
				sumwt[0], wgrid[0][0][0], grid[0][0][0]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

		{ // nullptr
			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					nullptr, rows.x, rows.y, SUPPORT, SAMPLING, RT::kNVisPol,
					polmap, RT::kNVisChan, chanmap, rows.mask[0][0],
					rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, nullptr, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, nullptr, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, nullptr, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, nullptr,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap, nullptr,
					rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], nullptr, rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], nullptr, rows.weight[0], true,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], nullptr, false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), nullptr, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, nullptr,
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					nullptr, grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], nullptr);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
		}

		{ // alignment
			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask + 1, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x + 1, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y + 1, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap + 1, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap + 1,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					&rows.mask[0][0][1], rows.values[0][0], rows.weight[0],
					false, ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY,
					sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], &rows.values[0][0][1], rows.weight[0],
					false, ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY,
					sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], &rows.values[0][0][1], rows.weight[0],
					true, ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY,
					sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], &rows.weight[0][1],
					false, ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY,
					sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab + 1, kNPol, kNChan, kNX, kNY,
					sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY,
					&sumwt[0][1], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					&wgrid[0][0][0][1], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], &grid[0][0][0][1]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
		}

		{ // out of range
			constexpr uint32_t int32max = INT32_MAX;
			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow + 1,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, 0, SAMPLING, RT::kNVisPol,
					polmap, RT::kNVisChan, chanmap, rows.mask[0][0],
					rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, (int32max - 1) / 2 + 1,
					SAMPLING, RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, 0, RT::kNVisPol,
					polmap, RT::kNVisChan, chanmap, rows.mask[0][0],
					rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, int32max + 1,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING, 0, polmap,
					RT::kNVisChan, chanmap, rows.mask[0][0], rows.values[0][0],
					rows.weight[0], false, ELEMENTSOF(conv_tab), conv_tab, kNPol,
					kNChan, kNX, kNY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					int32max + 1, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, 0, chanmap, rows.mask[0][0],
					rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, int32max + 1, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, 0, kNChan, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, int32max + 1, kNChan, kNX, kNY,
					sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, 0, kNX, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, int32max + 1, kNX, kNY,
					sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, 0, kNY, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, int32max + 1, kNY,
					sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, 0, sumwt[0],
					wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, kNX, int32max + 1,
					sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					0, conv_tab, kNPol, kNChan, kNX, kNY, sumwt[0], wgrid[0][0][0],
					grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab) + 1, conv_tab, kNPol, kNChan, kNX, kNY,
					sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab) - 1, conv_tab, kNPol, kNChan, kNX, kNY,
					sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

			result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow, 0, RT::kNRow,
					rows.sp_mask, rows.x, rows.y, SUPPORT, SAMPLING,
					RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
					rows.mask[0][0], rows.values[0][0], rows.weight[0], false,
					ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan, int32max / 2, 3,
					sumwt[0], wgrid[0][0][0], grid[0][0][0]);
			EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);

		}
	}

	template<typename T, typename F>
	static void Dump(T *a, F func, int width, int prec) {
		auto row_sep = "";
		for (size_t i = 0; i < ELEMENTSOF(*a); ++i) {
			cout << row_sep;
			auto col_sep = "{";
			for (size_t j = 0; j < ELEMENTSOF((*a)[0]); ++j) {
				cout << col_sep << setw(width) << fixed << setprecision(prec)
						<< func(a, i, j);
				col_sep = ", ";
			}
			cout << "}";
			row_sep = ",\n";
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
	void TrySpeed(bool do_verify, bool weight_only, RT *rows,
			double (*sumwt2)[kNPol][kNChan], float (*wgrid2)[kNY][kNX][kNPol][kNChan],
			float (*grid2)[kNY][kNX][kNPol][kNChan], char const *key = nullptr) {
		STATIC_ASSERT(RT::kNVisChan == kNVisChan);
		STATIC_ASSERT(RT::kNVisPol == kNVisPol);
		{
			cout << "Gridding ... " << flush;
			double start = CurrentTime();
			for (size_t i = 0; i < RT::kRowFactor; ++i) {
				LIBSAKURA_SYMBOL(
						Status) result = LIBSAKURA_SYMBOL(GridConvolvingFloat)(RT::kNRow,
						0, RT::kNRow, rows->sp_mask, rows->x, rows->y, SUPPORT,
						SAMPLING, RT::kNVisPol, polmap, RT::kNVisChan, chanmap,
						rows->mask[0][0], rows->values[0][0], rows->weight[0],
						weight_only, ELEMENTSOF(conv_tab), conv_tab, kNPol, kNChan,
						kNX, kNY, sumwt[0], wgrid[0][0][0], grid[0][0][0]);
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
			}
			double end = CurrentTime();
			double time_diff = end - start;
			cout << "done. \n";
			if (key) {
				ReportBenchmark(key, time_diff);
			}
		}
		if (do_verify) {
			EXPECT_TRUE(CmpSumwt(&sumwt, sumwt2));
			EXPECT_TRUE(CmpGrid("wgrid", &wgrid, wgrid2));
			EXPECT_TRUE(CmpGrid("grid", &grid, grid2));
		}
	}

	template<typename RT>
	void TestGrid(
			function<
					void(GridBase *, RT *, SumType *sumwt, GridType *wgrid,
							GridType *grid)> func) {
		STATIC_ASSERT(RT::kNVisChan == kNVisChan);
		STATIC_ASSERT(RT::kNVisPol == kNVisPol);
		InitTabs();

		RT *rows = nullptr;
		unique_ptr<void, DefaultAlignedMemory> rows_mem(
				DefaultAlignedMemory::AlignedAllocateOrException(sizeof(RT),
						&rows));
		rows->Init();

		GridType *wgrid2 = nullptr;
		unique_ptr<void, DefaultAlignedMemory> wgrid2_mem(
				DefaultAlignedMemory::AlignedAllocateOrException(
						sizeof(*wgrid2), &wgrid2));
		// init wgrid2;

		GridType *grid2 = nullptr;
		unique_ptr<void, DefaultAlignedMemory> grid2_mem(
				DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*grid2),
						&grid2));

		SumType *sumwt2 = nullptr;
		unique_ptr<void, DefaultAlignedMemory> sumwt2_mem(
				DefaultAlignedMemory::AlignedAllocateOrException(
						sizeof(*sumwt2), &sumwt2));

		func(this, rows, sumwt2, wgrid2, grid2);
	}
};

typedef GridBase<1, 1, 1, 1, 2, 2, 7, 7, InitFuncs<4> > TestMinimal;
typedef GridBase<1024, 4, 1024, 4, 100, 10, 200, 180, InitFuncs<1024> > TestTypical;
typedef GridBase<1023, 4, 17, 2, 100, 10, 200, 180, InitFuncs<1023> > TestOdd;

}
// namespace
#pragma GCC diagnostic ignored "-Wmissing-braces"

TEST(Gridding, Basic) {
	typedef TestMinimal TestCase;
	typedef RowBase<2, 1, TestCase::kNVisChan, TestCase::kNVisPol> RowType;

	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
					&test_case));
	test_case->SetUp();
	auto fin_func = [=]() {test_case->TearDown();};
	auto fin = Finalizer<decltype(fin_func)>(fin_func);

	test_case->TestInvalidArg<RowType>();
	test_case->TestGrid<RowType>(
			[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
				rows->x[1] = INT64_MIN;
				rows->y[1] = INT64_MAX;
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					bool weight_only = true;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
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
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &grid);
					TestCase::Dump(&tc->sumwt, 5);
					TestCase::Dump(&tc->wgrid, 5);
					TestCase::Dump(&tc->grid, 5);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->x[0] = 2.5;
					rows->y[0] = 2.5;
					rows->values[0][0][0] = 2;
					{
						auto p = *sumwt2;
						p[0][0] = 109.f;
					}
					TestCase::GridType wgrid = {
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 4.0, 5.0, 5.0, 4.0, 3.0, 0.0},
						{	0.0, 5.0, 7.0, 7.0, 5.0, 3.0, 0.0},
						{	0.0, 5.0, 7.0, 7.0, 5.0, 3.0, 0.0},
						{	0.0, 4.0, 5.0, 5.0, 4.0, 3.0, 0.0},
						{	0.0, 3.0, 3.0, 3.0, 3.0, 1.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
					};
					TestCase::GridType grid = {
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 8.0, 10.0, 10.0, 8.0, 6.0, 0.0},
						{	0.0, 10.0, 14.0, 14.0, 10.0, 6.0, 0.0},
						{	0.0, 10.0, 14.0, 14.0, 10.0, 6.0, 0.0},
						{	0.0, 8.0, 10.0, 10.0, 8.0, 6.0, 0.0},
						{	0.0, 6.0, 6.0, 6.0, 6.0, 2.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
					};
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &grid);
					TestCase::Dump(&tc->sumwt, 5);
					TestCase::Dump(&tc->wgrid, 5);
					TestCase::Dump(&tc->grid, 5);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->x[0] = 2.49;
					rows->y[0] = 2.49;
					rows->values[0][0][0] = 2;
					{
						auto p = *sumwt2;
						p[0][0] = 109.f;
					}
					TestCase::GridType wgrid = {
						{	1.0, 3.0, 3.0, 3.0, 3.0, 0.0, 0.0},
						{	3.0, 4.0, 5.0, 5.0, 4.0, 0.0, 0.0},
						{	3.0, 5.0, 7.0, 7.0, 5.0, 0.0, 0.0},
						{	3.0, 5.0, 7.0, 7.0, 5.0, 0.0, 0.0},
						{	3.0, 4.0, 5.0, 5.0, 4.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
					};
					TestCase::GridType grid = {
						{	2.0, 6.0, 6.0, 6.0, 6.0, 0.0, 0.0},
						{	6.0, 8.0, 10.0, 10.0, 8.0, 0.0, 0.0},
						{	6.0, 10.0, 14.0, 14.0, 10.0, 0.0, 0.0},
						{	6.0, 10.0, 14.0, 14.0, 10.0, 0.0, 0.0},
						{	6.0, 8.0, 10.0, 10.0, 8.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
					};
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &grid);
					TestCase::Dump(&tc->sumwt, 5);
					TestCase::Dump(&tc->wgrid, 5);
					TestCase::Dump(&tc->grid, 5);
				}
				{ // incremental grid
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->x[0] = 3;
					rows->y[0] = 3;
					rows->values[0][0][0] = 2;
					{
						auto p = *sumwt2;
						p[0][0] = 225.f;
					}
					TestCase::GridType wgrid = {
						{	1.0, 3.0, 3.0, 3.0, 3.0, 0.0, 0.0},
						{	3.0, 7.0, 9.0, 9.0, 8.0, 3.0, 0.0},
						{	3.0, 9.0, 13.0, 13.0, 11.0, 4.0, 0.0},
						{	3.0, 9.0, 13.0, 15.0, 11.0, 4.0, 0.0},
						{	3.0, 8.0, 11.0, 11.0, 10.0, 4.0, 0.0},
						{	0.0, 3.0, 4.0, 4.0, 4.0, 3.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
					};
					TestCase::GridType grid = {
						{	2.0, 6.0, 6.0, 6.0, 6.0, 0.0, 0.0},
						{	6.0, 14.0, 18.0, 18.0, 16.0, 6.0, 0.0},
						{	6.0, 18.0, 26.0, 26.0, 22.0, 8.0, 0.0},
						{	6.0, 18.0, 26.0, 30.0, 22.0, 8.0, 0.0},
						{	6.0, 16.0, 22.0, 22.0, 20.0, 8.0, 0.0},
						{	0.0, 6.0, 8.0, 8.0, 8.0, 6.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
					};
					bool weight_only = false;
					tc->TrySpeed(false, weight_only, rows, sumwt2, &wgrid, &grid);
					rows->x[0] = 2.49;
					rows->y[0] = 2.49;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &grid);
					TestCase::Dump(&tc->sumwt, 5);
					TestCase::Dump(&tc->wgrid, 5);
					TestCase::Dump(&tc->grid, 5);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->x[0] = 3;
					rows->y[0] = 3;
					rows->values[0][0][0] = 100;
					{
						auto p = *sumwt2;
						p[0][0] = 348.f;
					}
					TestCase::GridType wgrid = {
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 9.0, 12.0, 12.0, 12.0, 9.0, 0.0},
						{	0.0, 12.0, 18.0, 18.0, 18.0, 12.0, 0.0},
						{	0.0, 12.0, 18.0, 24.0, 18.0, 12.0, 0.0},
						{	0.0, 12.0, 18.0, 18.0, 18.0, 12.0, 0.0},
						{	0.0, 9.0, 12.0, 12.0, 12.0, 9.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
					};
					bool weight_only = true;
					tc->TrySpeed(false, weight_only, rows, sumwt2, wgrid2, grid2);

					rows->values[0][0][0] = 0;
					weight_only = false;
					tc->TrySpeed(false, weight_only, rows, sumwt2, wgrid2, grid2);

					rows->values[0][0][0] = 2;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &wgrid);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
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
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 1.5, 2.0, 2.0, 2.0, 1.5, 0.0},
						{	0.0, 2.0, 3.0, 3.0, 3.0, 2.0, 0.0},
						{	0.0, 2.0, 3.0, 4.0, 3.0, 2.0, 0.0},
						{	0.0, 2.0, 3.0, 3.0, 3.0, 2.0, 0.0},
						{	0.0, 1.5, 2.0, 2.0, 2.0, 1.5, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
					};
					bool weight_only = true;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &wgrid);
				}
			});
}

TEST(Gridding, Masking) {
	typedef TestMinimal TestCase;
	typedef RowBase<1, 1, TestCase::kNVisChan, TestCase::kNVisPol> RowType;

	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
					&test_case));
	test_case->SetUp();
	auto fin_func = [=]() {test_case->TearDown();};
	auto fin = Finalizer<decltype(fin_func)>(fin_func);

	test_case->TestGrid<RowType>(
			[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->x[0] = 3;
					rows->y[0] = 3;
					rows->values[0][0][0] = 2;
					rows->sp_mask[0] = false;
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->x[0] = 3;
					rows->y[0] = 3;
					rows->values[0][0][0] = 2;
					rows->sp_mask[0] = true;
					rows->mask[0][0][0] = false;
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);
				}
			});
}

TEST(Gridding, Clipping) {
	typedef TestMinimal TestCase;
	typedef RowBase<1, 1, TestCase::kNVisChan, TestCase::kNVisPol> RowType;

	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
					&test_case));
	test_case->SetUp();
	auto fin_func = [=]() {test_case->TearDown();};
	auto fin = Finalizer<decltype(fin_func)>(fin_func);

	test_case->TestGrid<RowType>(
			[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->x[0] = 1.5;
					rows->y[0] = 1.5;
					rows->values[0][0][0] = 2;
					{
						auto p = *sumwt2;
						p[0][0] = 109.f;
					}
					TestCase::GridType wgrid = {
						{	4.0, 5.0, 5.0, 4.0, 3.0, 0.0, 0.0},
						{	5.0, 7.0, 7.0, 5.0, 3.0, 0.0, 0.0},
						{	5.0, 7.0, 7.0, 5.0, 3.0, 0.0, 0.0},
						{	4.0, 5.0, 5.0, 4.0, 3.0, 0.0, 0.0},
						{	3.0, 3.0, 3.0, 3.0, 1.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
					};
					TestCase::GridType grid = {
						{	8.0, 10.0, 10.0, 8.0, 6.0, 0.0, 0.0},
						{	10.0, 14.0, 14.0, 10.0, 6.0, 0.0, 0.0},
						{	10.0, 14.0, 14.0, 10.0, 6.0, 0.0, 0.0},
						{	8.0, 10.0, 10.0, 8.0, 6.0, 0.0, 0.0},
						{	6.0, 6.0, 6.0, 6.0, 2.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
					};
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &grid);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->x[0] = 4.5;
					rows->y[0] = 4.5;
					rows->values[0][0][0] = 2;
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->x[0] = 4.49;
					rows->y[0] = 4.49;
					rows->values[0][0][0] = 2;
					{
						auto p = *sumwt2;
						p[0][0] = 109.f;
					}
					TestCase::GridType wgrid = {
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 0.0, 1.0, 3.0, 3.0, 3.0, 3.0},
						{	0.0, 0.0, 3.0, 4.0, 5.0, 5.0, 4.0},
						{	0.0, 0.0, 3.0, 5.0, 7.0, 7.0, 5.0},
						{	0.0, 0.0, 3.0, 5.0, 7.0, 7.0, 5.0},
						{	0.0, 0.0, 3.0, 4.0, 5.0, 5.0, 4.0}
					};
					TestCase::GridType grid = {
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 0.0, 2.0, 6.0, 6.0, 6.0, 6.0},
						{	0.0, 0.0, 6.0, 8.0, 10.0, 10.0, 8.0},
						{	0.0, 0.0, 6.0, 10.0, 14.0, 14.0, 10.0},
						{	0.0, 0.0, 6.0, 10.0, 14.0, 14.0, 10.0},
						{	0.0, 0.0, 6.0, 8.0, 10.0, 10.0, 8.0}
					};
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &grid);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->values[0][0][0] = 2;
					rows->x[0] = 1.49;
					rows->y[0] = 1.5;
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);

					rows->x[0] = 1.5;
					rows->y[0] = 1.49;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);

					rows->x[0] = 4.49;
					rows->y[0] = 4.5;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);

					rows->x[0] = 4.5;
					rows->y[0] = 4.49;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {8, 7, 6, 5, 4, 3, 2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->x[0] = -1;
					rows->y[0] = 3;
					rows->values[0][0][0] = 2;
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);

					rows->x[0] = 8;
					rows->y[0] = 3;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);

					rows->x[0] = 3;
					rows->y[0] = -1;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);

					rows->x[0] = 3;
					rows->y[0] = 8;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);
				}
			});
}

TEST(Gridding, ChMapping) {
	typedef GridBase<2, 1, 2, 1, 1, 1, 5, 5, InitFuncs<4> > TestCase;
	typedef RowBase<1, 1, TestCase::kNVisChan, TestCase::kNVisPol> RowType;

	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
					&test_case));
	test_case->SetUp();
	auto fin_func = [=]() {test_case->TearDown();};
	auto fin = Finalizer<decltype(fin_func)>(fin_func);

	test_case->TestInvalidArg<RowType>();
	test_case->TestGrid<RowType>(
			[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {2, 1, -.9f};
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
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &grid);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {2, 1, -.9f};
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
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &grid);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->x[0] = 2;
					rows->y[0] = 2;
					rows->values[0][0][0] = 2;
					rows->values[0][0][1] = 3;
					rows->weight[0][1] = 0.5;
					rows->mask[0][0][0] = false;
					{
						auto p = *sumwt2;
						p[0][0] = 0;
						p[0][1] = 5;
					}
					TestCase::GridTypeRef rwgrid = { {
							{
								{	0.0, 0.0, 0.0, 0.0, 0.0},
								{	0.0, 0.0, 0.0, 0.0, 0.0},
								{	0.0, 0.0, 0.0, 0.0, 0.0},
								{	0.0, 0.0, 0.0, 0.0, 0.0},
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
					TestCase::GridType wgrid;
					TestCase::GridType grid;
					TestCase::Transform(&rwgrid, &wgrid);
					TestCase::Transform(&rgrid, &grid);
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &grid);
				}
			});
}

TEST(Gridding, PolMapping) {
	typedef GridBase<1, 2, 1, 2, 1, 1, 5, 5, InitFuncs<4> > TestCase;
	typedef RowBase<1, 1, TestCase::kNVisChan, TestCase::kNVisPol> RowType;

	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
					&test_case));
	test_case->SetUp();
	auto fin_func = [=]() {test_case->TearDown();};
	auto fin = Finalizer<decltype(fin_func)>(fin_func);

	test_case->TestInvalidArg<RowType>();
	test_case->TestGrid<RowType>(
			[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {2, 1, -.9f};
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
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &grid);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {2, 1, -.9f};
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
						p[1][0] = 10;
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
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &grid);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					static float ct[] = {2, 1, -.9f};
					tc->SetConvTab(ELEMENTSOF(ct), ct);
					rows->x[0] = 2;
					rows->y[0] = 2;
					rows->values[0][0][0] = 2;
					rows->values[0][1][0] = 3;
					rows->mask[0][0][0] = false;
					{
						auto p = *sumwt2;
						p[0][0] = 0;
						p[1][0] = 5;
					}
					TestCase::GridTypeRef rwgrid = { { {
								{	0.0, 0.0, 0.0, 0.0, 0.0},
								{	0.0, 0.0, 0.0, 0.0, 0.0},
								{	0.0, 0.0, 0.0, 0.0, 0.0},
								{	0.0, 0.0, 0.0, 0.0, 0.0},
								{	0.0, 0.0, 0.0, 0.0, 0.0}
							}}, { {
								{	0.0, 0.0, 0.0, 0.0, 0.0},
								{	0.0, 0.5, 0.5, 0.5, 0.0},
								{	0.0, 0.5, 1.0, 0.5, 0.0},
								{	0.0, 0.5, 0.5, 0.5, 0.0},
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
								{	0.0, 1.5, 1.5, 1.5, 0.0},
								{	0.0, 1.5, 3.0, 1.5, 0.0},
								{	0.0, 1.5, 1.5, 1.5, 0.0},
								{	0.0, 0.0, 0.0, 0.0, 0.0}
							}}
					};
					TestCase::GridType wgrid;
					TestCase::GridType grid;
					TestCase::Transform(&rwgrid, &wgrid);
					TestCase::Transform(&rgrid, &grid);
					bool weight_only = false;
					tc->TrySpeed(true, weight_only, rows, sumwt2, &wgrid, &grid);
				}
			});
}

#define ENDOFARRAY(x) (&(x)[ELEMENTSOF(x)])

TEST(Gridding, Typical) {
	typedef TestTypical TestCase;
	typedef RowBase<512, 1, TestCase::kNVisChan, TestCase::kNVisPol> RowType;

	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
					&test_case));
	test_case->SetUp();
	auto fin_func = [=]() {test_case->TearDown();};
	auto fin = Finalizer<decltype(fin_func)>(fin_func);

	test_case->TestInvalidArg<RowType>();
	test_case->TestGrid<RowType>(
			[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
				{
					constexpr int pos = 50;
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					rows->SetWeight(7);
					rows->SetValues(5);
					rows->SetXY(pos, pos);
					fill_n(tc->conv_tab, ELEMENTSOF(tc->conv_tab), 11);
					bool weight_only = false;
					tc->TrySpeed(false, weight_only, rows, sumwt2, wgrid2, grid2, "Gridding_Typical_1");

					double wsum = rows->kNRow;
					wsum *= rows->weight[0][0] * tc->conv_tab[0];
					wsum *= Square(2 * tc->kSupport + 1);
					bool result = true;
					for (size_t i = 0; i < ELEMENTSOF(tc->sumwt); ++i) {
						if (!all_of(tc->sumwt[i], ENDOFARRAY(tc->sumwt[i]), [wsum](decltype(tc->sumwt[i][0]) v) -> bool {
											return v == wsum;
										})) {
							result = false;
						}
					}
					ASSERT_TRUE(result);

					float const grid_value = wsum * rows->values[0][0][0] / Square(2 * tc->kSupport + 1);
					float const wgrid_value = wsum / Square(2 * tc->kSupport + 1);

					auto is_inside = [=](size_t x, size_t y) -> bool {
						return pos - tc->kSupport <= x && x <= pos + tc->kSupport &&
						pos - tc->kSupport <= y && y <= pos + tc->kSupport;
					};
					bool fail = false;
					auto a = &tc->wgrid;
					for (size_t y = 0; y < ELEMENTSOF(*a); ++y) {
						for (size_t x = 0; x < ELEMENTSOF((*a)[0]); ++x) {
							for (size_t pol = 0; pol < ELEMENTSOF((*a)[0][0]); ++pol) {
								for (size_t ch = 0; ch < ELEMENTSOF((*a)[0][0][0]); ++ch) {
									auto rw = is_inside(x, y) ? wgrid_value : 0;
									auto w = tc->wgrid[y][x][pol][ch];
									auto rv = is_inside(x, y) ? grid_value : 0;
									auto v = tc->grid[y][x][pol][ch];
									if (!(rw == w && rv == v)) {
										cout << pol << ", " << ch << ": " << rv << ", " << v << endl;
										fail = true;
									}
								}
							}
						}
					}
					EXPECT_TRUE(fail == false);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					tc->InitConvTab();
					rows->Init();
					rows->sp_mask[rows->kNRow -1] = false;
					for (size_t i = 0; i < rows->kNRow; ++i) {
						rows->mask[i][rows->kNVisPol-2][rows->kNVisChan-2] = false;
					}
					rows->weight[0][0] = 2.f;
					rows->SetValues(1);
					rows->SetRandomXYInt(2 * TestCase::kSupport, TestCase::kNX - 2 * TestCase::kSupport,
							2 * TestCase::kSupport, TestCase::kNY - 2 * TestCase::kSupport);
					bool weight_only = false;
					tc->TrySpeed(false, weight_only, rows, sumwt2, wgrid2, grid2, "Gridding_Typical_2");
					ASSERT_DOUBLE_EQ(86829.7941513061523438, tc->sumwt[0][0]);
					ASSERT_DOUBLE_EQ(86660.2047096043825150, tc->sumwt[0][1]);
					EXPECT_TRUE(all_of(&tc->sumwt[0][2], &tc->sumwt[0][ELEMENTSOF(tc->sumwt[0])],
									[tc](decltype(tc->sumwt[0][0]++) x) -> bool {return x == tc->sumwt[0][1];}));
					for (size_t i = 1; i < rows->kNVisPol-2; ++i) {
						EXPECT_EQ(0, memcmp(tc->sumwt[0], tc->sumwt[i], sizeof(tc->sumwt[0])));
					}
					EXPECT_EQ(0, memcmp(tc->sumwt[0], tc->sumwt[rows->kNVisPol-1], sizeof(tc->sumwt[0])));
					bool fail = false;
					auto a = &tc->wgrid;
					for (size_t y = 0; y < ELEMENTSOF(*a); ++y) {
						for (size_t x = 0; x < ELEMENTSOF((*a)[0]); ++x) {
							for (size_t pol = 0; pol < ELEMENTSOF((*a)[0][0]); ++pol) {
								for (size_t ch = 0; ch < ELEMENTSOF((*a)[0][0][0]); ++ch) {
									auto rw = tc->wgrid[y][x][0][1];
									auto w = tc->wgrid[y][x][pol][ch];
									auto rv = tc->grid[y][x][0][1];
									auto v = tc->grid[y][x][pol][ch];
									if (pol == rows->kNVisPol-2 && ch ==rows->kNVisChan-2) {
										EXPECT_EQ(0., w);
										EXPECT_EQ(0., v);
									} else if (ch == 0) {
										EXPECT_LE(rw, w);
										EXPECT_LE(rv, v);
									} else if (rw == w && rv == v) {
									} else {
										cout << pol << ", " << ch << ": " << rv << ", " << v << endl;
										fail = true;
									}
								}
							}
						}
					}
					EXPECT_TRUE(fail == false);
				}
			});
}

TEST(Gridding, Odd) {
	typedef TestTypical TestCase;
	typedef RowBase<256, 1, TestCase::kNVisChan, TestCase::kNVisPol> RowType;

	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
					&test_case));
	test_case->SetUp();
	auto fin_func = [=]() {test_case->TearDown();};
	auto fin = Finalizer<decltype(fin_func)>(fin_func);

	test_case->TestInvalidArg<RowType>();
	test_case->TestGrid<RowType>(
			[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
				{
					constexpr int pos = 50;
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					rows->SetWeight(7);
					rows->SetValues(5);
					rows->SetXY(pos, pos);
					tc->chanmap[0] = 1;
					fill_n(tc->conv_tab, ELEMENTSOF(tc->conv_tab), 11);
					bool weight_only = false;
					tc->TrySpeed(false, weight_only, rows, sumwt2, wgrid2, grid2, "Gridding_Odd_1");

					double wsum = rows->kNRow;
					wsum *= rows->weight[0][0] * tc->conv_tab[0];
					wsum *= Square(2 * tc->kSupport + 1);
					bool result = true;
					for (size_t i = 0; i < ELEMENTSOF(tc->sumwt); ++i) {
						if (!(tc->sumwt[i][0] == 0 && tc->sumwt[i][1]) == wsum * 2) {
							result = false;
						}
						if (!all_of(&tc->sumwt[i][2], ENDOFARRAY(tc->sumwt[i]), [wsum](decltype(tc->sumwt[i][0]) v) -> bool {
											return v == wsum;
										})) {
							cout << i << endl;
							result = false;
						}
					}
					ASSERT_TRUE(result);

					float const grid_value = wsum * rows->values[0][0][0] / Square(2 * tc->kSupport + 1);
					float const wgrid_value = wsum / Square(2 * tc->kSupport + 1);

					auto is_inside = [=](size_t x, size_t y) -> bool {
						return pos - tc->kSupport <= x && x <= pos + tc->kSupport &&
						pos - tc->kSupport <= y && y <= pos + tc->kSupport;
					};
					bool fail = false;
					auto a = &tc->wgrid;
					for (size_t y = 0; y < ELEMENTSOF(*a); ++y) {
						for (size_t x = 0; x < ELEMENTSOF((*a)[0]); ++x) {
							for (size_t pol = 0; pol < ELEMENTSOF((*a)[0][0]); ++pol) {
								for (size_t ch = 0; ch < ELEMENTSOF((*a)[0][0][0]); ++ch) {
									float k = ch == 0 ? 0 : (ch == 1 ? 2 : 1);
									auto rw = is_inside(x, y) ? wgrid_value * k : 0;
									auto w = tc->wgrid[y][x][pol][ch];
									auto rv = is_inside(x, y) ? grid_value * k : 0;
									auto v = tc->grid[y][x][pol][ch];
									if (!(rw == w && rv == v)) {
										cout << pol << ", " << ch << ": " << rv << ", " << v << endl;
										fail = true;
									}
								}
							}
						}
					}
					EXPECT_TRUE(fail == false);
				}
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					tc->InitConvTab();
					rows->Init();
					tc->chanmap[0] = 1;
					rows->sp_mask[rows->kNRow -1] = false;
					for (size_t i = 0; i < rows->kNRow; ++i) {
						rows->mask[i][rows->kNVisPol-2][rows->kNVisChan-2] = false;
					}
					rows->weight[0][0] = 2.f;
					rows->SetValues(1);
					rows->SetRandomXYInt(2 * TestCase::kSupport, TestCase::kNX - 2 * TestCase::kSupport,
							2 * TestCase::kSupport, TestCase::kNY - 2 * TestCase::kSupport);
					bool weight_only = false;
					tc->TrySpeed(false, weight_only, rows, sumwt2, wgrid2, grid2, "Gridding_Odd_2");
					ASSERT_DOUBLE_EQ(0, tc->sumwt[0][0]);
					ASSERT_DOUBLE_EQ(86660.2047096043825150, tc->sumwt[0][1]);
					ASSERT_DOUBLE_EQ(43245.307633951306, tc->sumwt[0][2]);
					EXPECT_TRUE(all_of(&tc->sumwt[0][3], &tc->sumwt[0][ELEMENTSOF(tc->sumwt[0])],
									[tc](decltype(tc->sumwt[0][0]++) x) -> bool {return x == tc->sumwt[0][2];}));
					for (size_t i = 1; i < rows->kNVisPol-2; ++i) {
						EXPECT_EQ(0, memcmp(tc->sumwt[0], tc->sumwt[i], sizeof(tc->sumwt[0])));
					}
					EXPECT_EQ(0, memcmp(tc->sumwt[0], tc->sumwt[rows->kNVisPol-1], sizeof(tc->sumwt[0])));
					bool fail = false;
					auto a = &tc->wgrid;
					for (size_t y = 0; y < ELEMENTSOF(*a); ++y) {
						for (size_t x = 0; x < ELEMENTSOF((*a)[0]); ++x) {
							for (size_t pol = 0; pol < ELEMENTSOF((*a)[0][0]); ++pol) {
								for (size_t ch = 0; ch < ELEMENTSOF((*a)[0][0][0]); ++ch) {
									auto rw = tc->wgrid[y][x][0][2];
									auto w = tc->wgrid[y][x][pol][ch];
									auto rv = tc->grid[y][x][0][2];
									auto v = tc->grid[y][x][pol][ch];
									if (pol == rows->kNVisPol-2 && ch ==rows->kNVisChan-2) {
										EXPECT_EQ(0., w);
										EXPECT_EQ(0., v);
									} else if (ch == 0) {
										EXPECT_EQ(0, w);
										EXPECT_EQ(0, v);
									} else if (ch == 1) {
										EXPECT_LE(rw, w);
										EXPECT_LE(rv, v);
									} else if (rw == w && rv == v) {
									} else {
										cout << pol << ", " << ch << ": " << rv << ", " << v << endl;
										fail = true;
									}
								}
							}
						}
					}
					EXPECT_TRUE(fail == false);
				}
			});
}

TEST(Gridding, VectorizedWeightOnly) {
	typedef GridBase<64, 1, 64, 1, 2, 2, 10, 10, InitFuncs<4> > TestCase;
	typedef RowBase<1, 1, TestCase::kNVisChan, TestCase::kNPol> RowType;

	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
					&test_case));
	test_case->SetUp();
	auto fin_func = [=]() {test_case->TearDown();};
	auto fin = Finalizer<decltype(fin_func)>(fin_func);

	test_case->TestGrid<RowType>(
			[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					rows->SetValues(1000);
					rows->SetXY(4.5, 4.4);
					rows->SetWeight(2);
					fill_n(rows->mask[0][0], TestCase::kNChan, false);
					fill_n(rows->mask[0][0], TestCase::kNChan / 2, true);
					static float ct[] = {2, 1, 0.5f, 0.1f, -.9f, -10, -15, -30 };
					tc->SetConvTab(ELEMENTSOF(ct), ct);

					double ref[TestCase::kNY][TestCase::kNX] = {
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, -20.0, -20.0, -20.0, -20.0, -60.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, -1.8, 0.2, 0.2, -1.8, -20.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.2, 2.0, 2.0, 0.2, -20.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.2, 2.0, 2.0, 0.2, -20.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, -1.8, 0.2, 0.2, -1.8, -20.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
						{	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
					double sumwt = 0;
					for (size_t y = 0; y < TestCase::kNY; ++y) {
						for (size_t x = 0; x < TestCase::kNX; ++x) {
							sumwt += ref[y][x];
							for (size_t ch = 0; ch < TestCase::kNChan/2; ++ch) {
								(*grid2)[y][x][0][ch] = ref[y][x];
								(*wgrid2)[y][x][0][ch] = ref[y][x];
								//sumwt2[y][x][0][ch] = ref[y][x];
							}
						}
					}
					fill_n((*sumwt2)[0], TestCase::kNChan/2, sumwt);

					bool weight_only = true;
					tc->TrySpeed(true, weight_only, rows, sumwt2, wgrid2, grid2);
				}
			});
}

TEST(Gridding, PerformanceVectorized) {
	typedef TestTypical TestCase;
	typedef RowBase<512, 2, TestCase::kNVisChan, TestCase::kNVisPol> RowType;

	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
					&test_case));
	test_case->SetUp();
	auto fin_func = [=]() {test_case->TearDown();};
	auto fin = Finalizer<decltype(fin_func)>(fin_func);

	test_case->TestGrid<RowType>(
			[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					rows->SetValues(1);
					rows->SetRandomXY(-2. * TestCase::kSupport,
							TestCase::kNX + 2 * TestCase::kSupport,
							-2. * TestCase::kSupport,
							TestCase::kNY + 2 * TestCase::kSupport);
					bool weight_only = false;
					tc->TrySpeed(false, weight_only, rows, sumwt2, wgrid2, grid2, "Gridding_VectorizedWeighted");
				}
			});
}

TEST(Gridding, PerformanceScalar) {
	typedef TestTypical TestCase;
	typedef RowBase<256, 1, TestCase::kNVisChan, TestCase::kNVisPol> RowType;

	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
					&test_case));
	test_case->SetUp();
	auto fin_func = [=]() {test_case->TearDown();};
	auto fin = Finalizer<decltype(fin_func)>(fin_func);

	test_case->TestGrid<RowType>(
			[](TestCase *tc, RowType *rows, TestCase::SumType *sumwt2, TestCase::GridType *wgrid2, TestCase::GridType *grid2) {
				{
					TestCase::ClearResults(sumwt2, wgrid2, grid2);
					tc->ClearResults();
					rows->SetValues(1);
					rows->SetRandomXY(-2. * TestCase::kSupport, TestCase::kNX + 2 * TestCase::kSupport,
							-2. * TestCase::kSupport, TestCase::kNY + 2 * TestCase::kSupport);
					tc->chanmap[0] = 1;
					bool weight_only = false;
					tc->TrySpeed(false, weight_only, rows, sumwt2, wgrid2, grid2, "Gridding_ScalarWeighted");
				}
			});
}

