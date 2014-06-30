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

#include <libsakura/sakura.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"

#define DO_VERIFY 1
#define DO_FORTRAN (1 || DO_VERIFY)

//#define AVX __attribute__((aligned (32)))
#define AVX
#define elementsof(x) (sizeof(x) / sizeof((x)[0]))

using namespace std;

namespace {
inline double CurrentTime() {
	return LIBSAKURA_SYMBOL(GetCurrentTime)();
}

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

extern "C" {
void ggridsd_(double const xy[][2], const complex<float> *values,
		integer *nvispol, integer *nvischan, integer *dowt, integer const *flag,
		integer const *rflag, float const *weight, integer *nrow, integer *irow,
		complex<float> *grid, float *wgrid, integer *nx, integer *ny,
		integer *npol, integer *nchan, integer *support, integer *sampling,
		float *convFunc, integer *chanmap, integer *polmap, double *sumwt);
}

template<typename T>
inline bool NearlyEqual(T aa, T bb) {
	T threshold = 0.00007;
	//T threshold = 0.0000003;
	return aa == bb || (abs(aa - bb) / (abs(aa) + abs(bb))) < threshold;
}

template<size_t NROW, size_t ROW_FACTOR, size_t NVISCHAN, size_t NVISPOL,
		size_t NCHAN, size_t NPOL, size_t SAMPLING, size_t SUPPORT, size_t NX,
		size_t NY, typename InitFuncs>
class SIMD_ALIGN TestBase {
	static constexpr size_t CONV_TABLE_SIZE = (size_t(
			Ceil(Sqrt(2.) * ((SUPPORT + 1) * SAMPLING))));

	SIMD_ALIGN float convTab[CONV_TABLE_SIZE];SIMD_ALIGN uint32_t chanmap[NVISCHAN];SIMD_ALIGN uint32_t polmap[NVISPOL];

	SIMD_ALIGN integer flag[NROW][NVISCHAN][NVISPOL];SIMD_ALIGN bool flag_[NROW][NVISPOL][NVISCHAN];

	SIMD_ALIGN integer rflag[NROW];SIMD_ALIGN bool rflag_[NROW];

	SIMD_ALIGN float weight[NROW][NVISCHAN];SIMD_ALIGN double sumwt[NCHAN][NPOL];SIMD_ALIGN double sumwt2[NPOL][NCHAN];

	double fortran_time;

protected:
	void clearTabs(complex<float> (*grid)[NCHAN][NPOL][NY][NX],
			float (*wgrid)[NCHAN][NPOL][NY][NX], double (*sumwt)[NCHAN][NPOL]) {
		memset(*grid, 0, sizeof(*grid));
		memset(*wgrid, 0, sizeof(*wgrid));
		memset(*sumwt, 0, sizeof(*sumwt));
	}

	template<typename A, typename B>
	void clearTabs(A *grid, A *wgrid, B *sumwt) {
		memset(*grid, 0, sizeof(*grid));
		memset(*wgrid, 0, sizeof(*wgrid));
		memset(*sumwt, 0, sizeof(*sumwt));
	}

	bool cmpwgrid(float (*a)[NCHAN][NPOL][NY][NX],
			float (*b)[NY][NX][NPOL][NCHAN]) {
		cout << "comparing wgrid\n";
		size_t count = 0;
		size_t count2 = 0;
		bool differ = false;
		for (size_t i = 0; i < elementsof(*a); i++) {
			for (size_t j = 0; j < elementsof((*a)[0]); j++) {
				for (size_t k = 0; k < elementsof((*a)[0][0]); k++) {
					for (size_t l = 0; l < elementsof((*a)[0][0][0]); l++) {
						float aa = (*a)[i][j][k][l];
						float bb = (*b)[k][l][j][i];
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

	bool cmpsumwt(double (*a)[NCHAN][NPOL], double (*b)[NPOL][NCHAN]) {
		cout << "comparing sumwt\n";
		size_t count = 0;
		size_t count2 = 0;
		bool differ = false;
		for (size_t i = 0; i < elementsof(*a); i++) {
			for (size_t j = 0; j < elementsof((*a)[0]); j++) {
				double aa = (*a)[i][j];
				double bb = (*b)[j][i];
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

	bool cmpgrid(complex<float> (*a)[NCHAN][NPOL][NY][NX],
			float (*b)[NY][NX][NPOL][NCHAN]) {
		cout << "comparing grid\n";
		size_t count = 0;
		size_t count2 = 0;
		bool differ = false;
		for (size_t i = 0; i < elementsof(*a); i++) {
			for (size_t j = 0; j < elementsof((*a)[0]); j++) {
				for (size_t k = 0; k < elementsof((*a)[0][0]); k++) {
					for (size_t l = 0; l < elementsof((*a)[0][0][0]); l++) {
						float aa = (*a)[i][j][k][l].real();
						float bb = (*b)[k][l][j][i];
						if (NearlyEqual<float>(aa, bb)) {
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

	void trySpeed(bool dowt, float (*values_)[NROW][NVISPOL][NVISCHAN],
			double (*x)[NROW], double (*y)[NROW],
			complex<float> (*grid)[NCHAN][NPOL][NY][NX],
			float (*wgrid)[NCHAN][NPOL][NY][NX],
			float (*grid2)[NY][NX][NPOL][NCHAN],
			float (*wgrid2)[NY][NX][NPOL][NCHAN]) {
		{
			//cout << "Zeroing values... " << flush;
			//double start = CurrentTime();
			clearTabs(grid2, wgrid2, &sumwt2);
			//double end = CurrentTime();
			//cout << "done. " << end - start << " sec\n";
		}

		{
			cout << "Gridding by C++(Speed) ... " << flush;
			double start = CurrentTime();
			for (size_t i = 0; i < ROW_FACTOR; ++i) {
				LIBSAKURA_SYMBOL(Status) result = sakura_GridConvolving(NROW, 0,
						NROW, rflag_, *x, *y, SUPPORT, SAMPLING, NVISPOL,
						(uint32_t*) polmap, NVISCHAN, (uint32_t*) chanmap,
						flag_[0][0], (*values_)[0][0], weight[0], dowt,
						elementsof(convTab), convTab, NPOL, NCHAN, NX, NY,
						sumwt2[0], (*wgrid2)[0][0][0], (*grid2)[0][0][0]);
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
			}
			double end = CurrentTime();
			double new_time = end - start;
			cout << "done. " << new_time << " sec\n";
			cout << "ratio: " << (new_time / fortran_time) << ", "
					<< (fortran_time / new_time) << " times faster" << endl;
		}
#if DO_VERIFY
		EXPECT_TRUE(cmpgrid(grid, grid2));
		EXPECT_TRUE(cmpwgrid(wgrid, wgrid2));
#endif
		EXPECT_TRUE(cmpsumwt(&sumwt, &sumwt2));
	}

	void initTabs() {
		for (size_t i = 0; i < elementsof(convTab); i++) {
			convTab[i] = exp(-square(2.0 * i / elementsof(convTab)));
		}
		for (size_t i = 0; i < elementsof(chanmap); i++) {
			chanmap[i] = InitFuncs::Map(i, NCHAN);
		}
		for (size_t i = 0; i < elementsof(polmap); i++) {
			polmap[i] = InitFuncs::Map(i, NPOL);
		}
		for (size_t i = 0; i < elementsof(flag); i++) {
			for (size_t j = 0; j < elementsof(flag[0]); j++) {
				for (size_t k = 0; k < elementsof(flag[0][0]); k++) {
					flag_[i][k][j] = InitFuncs::Mask2D(j, k);
					flag[i][j][k] = flag_[i][k][j] == false;
				}
			}
		}
		for (size_t i = 0; i < elementsof(rflag); i++) {
			rflag_[i] = InitFuncs::Mask1D(i);
			rflag[i] = rflag_[i] == false;
		}
		for (size_t i = 0; i < elementsof(weight); i++) {
			for (size_t j = 0; j < elementsof(weight[0]); j++) {
				weight[i][j] = drand48(); // - 0.1;
			}
		}
	}

public:
	void SetUp() {
		LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
				nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
		srand48(0);
		assert(LIBSAKURA_SYMBOL(IsAligned(this)));
	}

	void TearDown() {
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	void TestInvalidArg() {
		LIBSAKURA_SYMBOL(Status) result = sakura_GridConvolving(NROW, 0, NROW,
				0, 0, 0, SUPPORT, SAMPLING, NVISPOL, (uint32_t*) polmap,
				NVISCHAN, (uint32_t*) chanmap, flag_[0][0], 0, weight[0], false,
				elementsof(convTab), convTab, NPOL, NCHAN, NX, NY, 0, 0, 0);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
	}

	void TestGrid() {
		initTabs();
		double (*xy_)[NROW][2] = nullptr;
		unique_ptr<void, DefaultAlignedMemory> xy(
				DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*xy_),
						&xy_));

		double (*x_)[NROW] = nullptr;
		unique_ptr<void, DefaultAlignedMemory> x(
				DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*x_),
						&x_));
		double (*y_)[NROW] = nullptr;
		unique_ptr<void, DefaultAlignedMemory> y(
				DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*y_),
						&y_));
		//cout << "Initializing xy... " << flush;
		{
			//double start = CurrentTime();
			for (size_t i = 0; i < elementsof(*xy_); i++) {
				(*xy_)[i][0] = (*x_)[i] = NX * drand48();
				(*xy_)[i][1] = (*y_)[i] = NY * drand48();
			}
			//double end = CurrentTime();
			//cout << "done. " << end - start << " sec\n";
		}

		complex<float> (*values)[NROW][NVISCHAN][NVISPOL] = nullptr;
		unique_ptr<void, DefaultAlignedMemory> valueMem(
				DefaultAlignedMemory::AlignedAllocateOrException(
						sizeof(*values), &values));

		float (*values2)[NROW][NVISPOL][NVISCHAN] = nullptr;
		unique_ptr<void, DefaultAlignedMemory> valueMem2(
				DefaultAlignedMemory::AlignedAllocateOrException(
						sizeof(*values2), &values2));
		cout << sizeof(*values) << endl;
		//cout << "Initializing values... " << flush;
		{
			//double start = CurrentTime();
			for (size_t i = 0; i < elementsof(*values); i++) {
				for (size_t j = 0; j < elementsof((*values)[0]); j++) {
					for (size_t k = 0; k < elementsof((*values)[0][0]); k++) {
						(*values)[i][j][k] = (*values2)[i][k][j] = drand48();
					}
				}
			}
			//double end = CurrentTime();
			//cout << "done. " << end - start << " sec\n";
		}

		complex<float> (*grid)[NCHAN][NPOL][NY][NX] = nullptr;
		unique_ptr<void, DefaultAlignedMemory> gridMem(
				DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*grid),
						&grid));

		float (*wgrid)[NCHAN][NPOL][NY][NX] = nullptr;
		unique_ptr<void, DefaultAlignedMemory> wgridMem(
				DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*wgrid),
						&wgrid));
		cout << "grid, wgrid size: " << sizeof(*grid) << ", " << sizeof(*wgrid)
				<< endl;
		{
			//cout << "Zeroing values... " << flush;
			//double start = CurrentTime();
			clearTabs(grid, wgrid, &sumwt);
			//double end = CurrentTime();
			//cout << "done. " << end - start << " sec\n";
		}

		assert(sizeof(int) == sizeof(uint32_t));
		bool dowt = true;
		{
			cout << "Gridding by Fortran ... " << flush;
			double start = CurrentTime();
			integer dowt_ = integer(dowt);
			integer nvispol = NVISPOL;
			integer nvischan = NVISCHAN;
			integer nrow = NROW;
			integer nx = NX, ny = NY;
			integer npol = NPOL, nchan = NCHAN;
			integer support = SUPPORT, sampling = SAMPLING;
			for (size_t i = 0; i < ROW_FACTOR * DO_FORTRAN; ++i) {
				integer irow = -1;
#if 1
				ggridsd_((*xy_), (*values)[0][0], &nvispol, &nvischan, &dowt_,
						flag[0][0], rflag, weight[0], &nrow, &irow,
						(*grid)[0][0][0], (*wgrid)[0][0][0], &nx, &ny, &npol,
						&nchan, &support, &sampling, convTab, (int*) chanmap,
						(int*) polmap, sumwt[0]);
#endif
			}
			double end = CurrentTime();
			fortran_time = end - start;
			cout << "done. " << fortran_time << " sec\n";
		}
#if !DO_VERIFY
		gridMem.reset();
		wgridMem.reset();
		grid = nullptr;
		wgrid = nullptr;
#endif

		valueMem.reset();

		float (*grid2)[NY][NX][NPOL][NCHAN] = nullptr;
		unique_ptr<void, DefaultAlignedMemory> grid2Mem(
				DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*grid2),
						&grid2));

		float (*wgrid2)[NY][NX][NPOL][NCHAN] = nullptr;
		unique_ptr<void, DefaultAlignedMemory> wgrid2Mem(
				DefaultAlignedMemory::AlignedAllocateOrException(
						sizeof(*wgrid2), &wgrid2));

		trySpeed(dowt, values2, x_, y_, grid, wgrid, grid2, wgrid2);
	}
};

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

typedef TestBase<512, 1, 1024, 4, 1024, 4, 100, 10, 200, 180, InitFuncs<1> > TestTypical;
typedef TestBase<512, 2, 1024, 4, 1024, 4, 100, 10, 200, 180, InitFuncs<2> > TestCase1;
typedef TestBase<512 - 1, 10, 1024 - 1, 4 - 1, 1024 - 1, 4, 100, 10, 100, 90,
		InitFuncs<7> > TestCase2;
typedef TestBase<512, 1, 512, 4, 1024, 4, 100, 10, 200, 180, InitFuncs2<1> > TestNto1;
typedef TestBase<512, 10000, 3, 1, 3, 1, 100, 10, 3, 3, InitFuncs<1> > TestSmall;

}

TEST(Gridding, Typical) {
	typedef TestTypical TestCase;
	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
	test_case->SetUp();
	test_case->TestInvalidArg();
	test_case->TestGrid();
	test_case->TearDown();
}

TEST(Gridding, Generic) {
	typedef TestCase1 TestCase;
	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
	test_case->SetUp();
	test_case->TestInvalidArg();
	test_case->TestGrid();
	test_case->TearDown();
}

TEST(Gridding, Odd) {
	typedef TestCase2 TestCase;
	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
	test_case->SetUp();
	test_case->TestInvalidArg();
	test_case->TestGrid();
	test_case->TearDown();
}

TEST(Gridding, Nto1) {
	typedef TestNto1 TestCase;
	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
	test_case->SetUp();
	test_case->TestInvalidArg();
	test_case->TestGrid();
	test_case->TearDown();
}

TEST(Gridding, Small) {
	typedef TestSmall TestCase;
	TestCase *test_case;
	unique_ptr<void, DefaultAlignedMemory> test_case_storage(
		DefaultAlignedMemory::AlignedAllocateOrException(sizeof(*test_case),
				&test_case));
	test_case->SetUp();
	test_case->TestInvalidArg();
	test_case->TestGrid();
	test_case->TearDown();
}
