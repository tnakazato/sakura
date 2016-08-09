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
#include <iostream>
#include <iomanip>
#include <memory>
#include <complex>
#include <functional>
#include <numeric>

#include <libsakura/sakura.h>

#include <libsakura/localdef.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"
#include "testutil.h"

using namespace std;

namespace {

constexpr bool verbose = false;

template<typename T>
void InitArray(size_t size, T p[]) {
	for (size_t i = 0; i < size; ++i) {
		p[i] = T(i + 1);
	}
}

template<typename T>
bool IsEqual(size_t size, T const p[], T const q[]) {
	for (size_t i = 0; i < size; ++i) {
		if (p[i] != q[i]) {
			cout << i << ": " << p[i] << " <> " << q[i] << endl;
			return false;
		}
	}
	return true;
}

template<typename T, size_t COL>
void Print(size_t size, T const p[]) {
	if (!verbose)
		return;
	for (size_t i = 0; i < size; ++i) {
		cout << setw(3) << p[i] << ", ";
		if ((i + 1) % COL == 0) {
			cout << endl;
		}
	}
}

template<typename T>
T Product(size_t n, T const data[]) {
	if (n == 0)
		return 0;
	return accumulate(&data[0], &data[n], T(1), multiplies<T>());
}

template<typename T>
struct TestTarget {
	static LIBSAKURA_SYMBOL(Status) Flip(
	bool inner_most_untouched, size_t dims, size_t const elements[],
			T const src[], T dst[]) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kNG);
	}
	static LIBSAKURA_SYMBOL(Status) Unflip(
	bool inner_most_untouched, size_t dims, size_t const elements[],
			T const src[], T dst[]) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kNG);
	}
};

template<>
struct TestTarget<float> {
	static LIBSAKURA_SYMBOL(Status) Flip(
	bool inner_most_untouched, size_t dims, size_t const elements[],
			float const src[], float dst[]) {
		return LIBSAKURA_SYMBOL(FlipArrayFloat)(inner_most_untouched, dims,
				elements, src, dst);
	}
	static LIBSAKURA_SYMBOL(Status) Unflip(
	bool inner_most_untouched, size_t dims, size_t const elements[],
			float const src[], float dst[]) {
		return LIBSAKURA_SYMBOL(UnflipArrayFloat)(inner_most_untouched, dims,
				elements, src, dst);
	}
};

template<>
struct TestTarget<double> {
	static LIBSAKURA_SYMBOL(Status) Flip(
	bool inner_most_untouched, size_t dims, size_t const elements[],
			double const src[], double dst[]) {
		return LIBSAKURA_SYMBOL(FlipArrayDouble)(inner_most_untouched, dims,
				elements, src, dst);
	}
	static LIBSAKURA_SYMBOL(Status) Unflip(
	bool inner_most_untouched, size_t dims, size_t const elements[],
			double const src[], double dst[]) {
		return LIBSAKURA_SYMBOL(UnflipArrayDouble)(inner_most_untouched, dims,
				elements, src, dst);
	}
};

template<>
struct TestTarget<complex<double> > {
	static LIBSAKURA_SYMBOL(Status) Flip(
	bool inner_most_untouched, size_t dims, size_t const elements[],
			complex<double> const src[], complex<double> dst[]) {
		return LIBSAKURA_SYMBOL(FlipArrayDouble2)(inner_most_untouched, dims,
				elements, reinterpret_cast<double const (*)[2]>(src),
				reinterpret_cast<double (*)[2]>(dst));
	}
	static LIBSAKURA_SYMBOL(Status) Unflip(
	bool inner_most_untouched, size_t dims, size_t const elements[],
			complex<double> const src[], complex<double> dst[]) {
		return LIBSAKURA_SYMBOL(UnflipArrayDouble2)(inner_most_untouched, dims,
				elements, reinterpret_cast<double const (*)[2]>(src),
				reinterpret_cast<double (*)[2]>(dst));
	}
};

void ReportNone(double sec, char const action[]) {
}

template<typename T, size_t COL, bool kTiming = false, int kRepeat = 1,
		typename Func = decltype(ReportNone)>
void TestGeneric(bool inner_most_untouched, size_t dims,
		size_t const elements[], T const ref[],
		Func timingReporter = ReportNone) {
	const size_t prod = Product(dims, elements);
	T *data = nullptr;
	unique_ptr<void, DefaultAlignedMemory> data_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(data[0]) * prod, &data));
	InitArray(prod, data);
	if (!kTiming) {
		if (verbose) {
			cout << "src:\n";
		}
		Print<T, COL>(prod, data);
	}

	T *flipped_data = nullptr;
	unique_ptr<void, DefaultAlignedMemory> flipped_data_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(flipped_data[0]) * prod, &flipped_data));
	LIBSAKURA_SYMBOL(Status) result;
	{
		double start = GetCurrentTime();
		for (int i = 0; i < kRepeat; ++i) {
			result = TestTarget<T>::Flip(inner_most_untouched, dims, elements,
					data, flipped_data);
		}
		double end = GetCurrentTime();
		if (kTiming) {
			timingReporter(end - start, "Flip");
		}
	}
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
	if (!kTiming) {
		if (verbose) {
			cout << "flipped(inner_most_untouched=" << inner_most_untouched
					<< "):\n";
		}
		Print<T, COL>(prod, flipped_data);
	}
	if (ref) {
		EXPECT_TRUE(IsEqual(prod, ref, flipped_data));
	}

	T *rev_data = nullptr;
	unique_ptr<void, DefaultAlignedMemory> rev_data_storage(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(rev_data[0]) * prod, &rev_data));

	{
		double start = GetCurrentTime();
		for (int i = 0; i < kRepeat; ++i) {
			result = TestTarget<T>::Unflip(inner_most_untouched, dims, elements,
					flipped_data, rev_data);
		}
		double end = GetCurrentTime();
		if (kTiming) {
			timingReporter(end - start, "Unflip");
		}
	}
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
	if (!kTiming) {
		if (verbose) {
			cout << "unflipped:\n";
		}
		Print<T, COL>(prod, rev_data);
	}
	EXPECT_TRUE(IsEqual(prod, data, rev_data));

}

template<typename T>
void Tests() {
	{
		static size_t const elements[] = { 4 };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { T(3), T(4), T(1), T(2) };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 4>(false, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { 4 };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { T(1), T(2), T(3), T(4) };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 4>(true, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { 5 };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { T(4), T(5), T(1), T(2), T(3) };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 5>(false, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { 5 };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { T(1), T(2), T(3), T(4), T(5) };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 5>(true, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { 5, 3 };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { T(14), T(15), T(11), T(12), T(13), T(4), T(5),
				T(1), T(2), T(3), T(9), T(10), T(6), T(7), T(8) };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 5>(false, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { 5, 3 };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { T(11), T(12), T(13), T(14), T(15), T(1), T(2),
				T(3), T(4), T(5), T(6), T(7), T(8), T(9), T(10) };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 5>(true, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { 3, 4 };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { T(9), T(7), T(8), T(12), T(10), T(11), T(3), T(
				1), T(2), T(6), T(4), T(5) };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 3>(false, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { 3, 4 };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { T(7), T(8), T(9), T(10), T(11), T(12), T(1), T(
				2), T(3), T(4), T(5), T(6) };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 3>(true, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { 2, 3, 2 };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { T(12), T(11), T(8), T(7), T(10), T(9),

		T(6), T(5), T(2), T(1), T(4), T(3) };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 2>(false, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { 2, 3, 2 };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { T(11), T(12), T(7), T(8), T(9), T(10),

		T(5), T(6), T(1), T(2), T(3), T(4) };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 2>(true, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { 1 };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { T(1) };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 1>(false, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { 1 };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { T(1) };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 1>(true, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 1>(false, dims, elements, ref);
	}
	if (verbose)
		cout << "---\n";
	{
		static size_t const elements[] = { };
		const size_t dims = ELEMENTSOF(elements);
		static T const ref[] = { };
		assert(ELEMENTSOF(ref) == Product(dims, elements));
		TestGeneric<T, 1>(true, dims, elements, ref);
	}
}

template<typename T, bool Flag>
void TestsError() {
	{
		const size_t dims = 0;
		SIMD_ALIGN
		T const src[] = { };
		SIMD_ALIGN
		T dst[] = { };
		auto result = TestTarget<T>::Flip(Flag, dims, nullptr, src, dst);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
	}
	{
		static size_t const elements[] = { };
		const size_t dims = ELEMENTSOF(elements);
		SIMD_ALIGN
		T dst[] = { };
		auto result = TestTarget<T>::Flip(Flag, dims, elements, nullptr, dst);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
	}
	{
		static size_t const elements[] = { };
		const size_t dims = ELEMENTSOF(elements);
		SIMD_ALIGN
		T const src[] = { };
		auto result = TestTarget<T>::Flip(Flag, dims, elements, src, nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
	}
	{
		static size_t const elements[] = { 1 };
		const size_t dims = ELEMENTSOF(elements);
		SIMD_ALIGN
		struct {
			char dummy[1];
			T array[1];
		} src;
		SIMD_ALIGN
		T dst[1];
		auto result = TestTarget<T>::Flip(Flag, dims, elements,
				reinterpret_cast<T const *>(&src.dummy[1]), dst);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
	}
	{
		static size_t const elements[] = { 1 };
		const size_t dims = ELEMENTSOF(elements);
		SIMD_ALIGN
		T const src[] = { 1 };
		SIMD_ALIGN
		struct {
			char dummy[1];
			T array[1];
		} dst;
		auto result = TestTarget<T>::Flip(Flag, dims, elements, src,
				reinterpret_cast<T *>(&dst.dummy[1]));
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
	}
}

template<typename MessageType>
void ReportBenchmark(MessageType const &key, double sec) {
	std::cout << std::setprecision(5) << "#x# benchmark FFT_" << key << " " << sec
			<< std::endl;
}

template<typename T>
void TestsSpeed(char const key_base[]) {
	auto inner_most_untouched = "_touched";
	auto reporter =
			[key_base, &inner_most_untouched](double sec, char const action[]) {
				string key(action);
				ReportBenchmark(key + key_base + inner_most_untouched, sec);
			};
	{
		static size_t const elements[] = { 128, 1024, 512 };
		const size_t dims = ELEMENTSOF(elements);
		TestGeneric<T, 512, true, 8>(false, dims, elements, nullptr, reporter);
	}
	inner_most_untouched = "_untouched";
	{
		static size_t const elements[] = { 128, 1024, 512 };
		const size_t dims = ELEMENTSOF(elements);
		TestGeneric<T, 512, true, 8>(true, dims, elements, nullptr, reporter);
	}
}

} // namespace

TEST(FFT, Basic) {
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
			nullptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

	cout << "float:\n";
	Tests<float>();
	cout << "double:\n";
	Tests<double>();
	cout << "complex<double>:\n";
	Tests<complex<double> >();

	LIBSAKURA_SYMBOL(CleanUp)();
}

TEST(FFT, ErrorCase) {
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
			nullptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

	cout << "float:\n";
	TestsError<float, false>();
	TestsError<float, true>();
	cout << "double:\n";
	TestsError<double, false>();
	TestsError<double, true>();
	cout << "complex<double>:\n";
	TestsError<complex<double>, false>();
	TestsError<complex<double>, true>();

	LIBSAKURA_SYMBOL(CleanUp)();
}

TEST(FFT, Performance) {
	LIBSAKURA_SYMBOL(Status) result = LIBSAKURA_SYMBOL(Initialize)(nullptr,
			nullptr);
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);

	cout << "float:\n";
	TestsSpeed<float>("MatrixFloat");
	cout << "double:\n";
	TestsSpeed<double>("MatrixDouble");
	cout << "complex<double>:\n";
	TestsSpeed<complex<double> >("MatrixDouble2");

	LIBSAKURA_SYMBOL(CleanUp)();
}
