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
#include <iostream>
#include <iomanip>
#include <memory>
#include <complex>

#include <libsakura/sakura.h>

#include <libsakura/localdef.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"

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
	//cout << endl;
}

template<typename T>
T Product(size_t n, T const data[]) {
	if (n == 0)
		return 0;
	T result = 1;
	for (size_t i = 0; i < n; ++i) {
		result *= data[i];
	}
	return result;
}

template<typename T>
struct TestTarget {
	static LIBSAKURA_SYMBOL(Status) Flip(
	bool innerMostUntouched, size_t dims, size_t const elements[],
			T const src[], T dst[]) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kNG);
	}
	static LIBSAKURA_SYMBOL(Status) Unflip(
	bool innerMostUntouched, size_t dims, size_t const elements[],
			T const src[], T dst[]) {
		assert(false);
		return LIBSAKURA_SYMBOL(Status_kNG);
	}
};

template<>
struct TestTarget<float> {
	static LIBSAKURA_SYMBOL(Status) Flip(
	bool innerMostUntouched, size_t dims, size_t const elements[],
			float const src[], float dst[]) {
		return LIBSAKURA_SYMBOL(FlipMatrixFloat)(innerMostUntouched, dims,
				elements, src, dst);
	}
	static LIBSAKURA_SYMBOL(Status) Unflip(
	bool innerMostUntouched, size_t dims, size_t const elements[],
			float const src[], float dst[]) {
		return LIBSAKURA_SYMBOL(UnflipMatrixFloat)(innerMostUntouched, dims,
				elements, src, dst);
	}
};

template<>
struct TestTarget<double> {
	static LIBSAKURA_SYMBOL(Status) Flip(
	bool innerMostUntouched, size_t dims, size_t const elements[],
			double const src[], double dst[]) {
		return LIBSAKURA_SYMBOL(FlipMatrixDouble)(innerMostUntouched, dims,
				elements, src, dst);
	}
	static LIBSAKURA_SYMBOL(Status) Unflip(
	bool innerMostUntouched, size_t dims, size_t const elements[],
			double const src[], double dst[]) {
		return LIBSAKURA_SYMBOL(UnflipMatrixDouble)(innerMostUntouched, dims,
				elements, src, dst);
	}
};

template<>
struct TestTarget<complex<double> > {
	static LIBSAKURA_SYMBOL(Status) Flip(
	bool innerMostUntouched, size_t dims, size_t const elements[],
			complex<double> const src[], complex<double> dst[]) {
		return LIBSAKURA_SYMBOL(FlipMatrixDouble2)(innerMostUntouched, dims,
				elements, reinterpret_cast<double const (*)[2]>(src),
				reinterpret_cast<double (*)[2]>(dst));
	}
	static LIBSAKURA_SYMBOL(Status) Unflip(
	bool innerMostUntouched, size_t dims, size_t const elements[],
			complex<double> const src[], complex<double> dst[]) {
		return LIBSAKURA_SYMBOL(UnflipMatrixDouble2)(innerMostUntouched, dims,
				elements, reinterpret_cast<double const (*)[2]>(src),
				reinterpret_cast<double (*)[2]>(dst));
	}
};

template<typename T, size_t COL, bool timing = false>
void TestGeneric(bool innerMostUntouched, size_t dims, size_t const elements[],
		T const ref[]) {
	const size_t prod = Product(dims, elements);
	T *data = nullptr;
	unique_ptr<void, DefaultAlignedMemory> dataStorage(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(data[0]) * prod, &data));
	InitArray(prod, data);
	if (!timing) {
		Print<T, COL>(prod, data);
	}

	T *flippedData = nullptr;
	unique_ptr<void, DefaultAlignedMemory> flippedDataStorage(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(flippedData[0]) * prod, &flippedData));
	int repeat = 1;
	LIBSAKURA_SYMBOL(Status) result;
	{
		double start = LIBSAKURA_SYMBOL(GetCurrentTime)();
		for (int i = 0; i < repeat; ++i) {
			result = TestTarget<T>::Flip(innerMostUntouched, dims, elements,
					data, flippedData);
		}
		double end = LIBSAKURA_SYMBOL(GetCurrentTime)();
		if (timing) {
			cout << "Flip time: " << end - start << " sec" << endl;
		}
	}
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
	if (!timing) {
		Print<T, COL>(prod, flippedData);
	}
	if (ref) {
		EXPECT_TRUE(IsEqual(prod, ref, flippedData));
	}

	T *revData = nullptr;
	unique_ptr<void, DefaultAlignedMemory> revDataStorage(
			DefaultAlignedMemory::AlignedAllocateOrException(
					sizeof(revData[0]) * prod, &revData));

	{
		double start = LIBSAKURA_SYMBOL(GetCurrentTime)();
		for (int i = 0; i < repeat; ++i) {
			result = TestTarget<T>::Unflip(innerMostUntouched, dims, elements,
					flippedData, revData);
		}
		double end = LIBSAKURA_SYMBOL(GetCurrentTime)();
		if (timing) {
			cout << "Unflip time: " << end - start << " sec" << endl;
		}
	}
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), result);
	if (!timing) {
		Print<T, COL>(prod, revData);
	}
	EXPECT_TRUE(IsEqual(prod, data, revData));

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
		T const src[] = { 1, 1 };
		SIMD_ALIGN
		T dst[] = { 1, 1 };
		auto result = TestTarget<T>::Flip(Flag, dims, elements, &src[1], dst);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
	}
	{
		static size_t const elements[] = { 1 };
		const size_t dims = ELEMENTSOF(elements);
		SIMD_ALIGN
		T const src[] = { 1, 1 };
		SIMD_ALIGN
		T dst[] = { 1, 1 };
		auto result = TestTarget<T>::Flip(Flag, dims, elements, src, &dst[1]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kInvalidArgument), result);
	}
}

template<typename T>
void TestsSpeed() {
	{
		static size_t const elements[] = { 1024, 1024, 512 };
		const size_t dims = ELEMENTSOF(elements);
		TestGeneric<T, 1024, true>(false, dims, elements, nullptr);
	}
	{
		static size_t const elements[] = { 1024, 1024, 512 };
		const size_t dims = ELEMENTSOF(elements);
		TestGeneric<T, 1024, true>(true, dims, elements, nullptr);
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
	TestsSpeed<float>();
	cout << "double:\n";
	TestsSpeed<double>();
	cout << "complex<double>:\n";
	TestsSpeed<complex<double> >();

	LIBSAKURA_SYMBOL(CleanUp)();
}
