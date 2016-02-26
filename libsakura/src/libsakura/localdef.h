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
/**
 * @file
 * Contains utility stuffs for internal use.
 *
 * This file is an internal header file for libsakura.
 * This is not a part of libsakura API.
 *
 */

#ifndef LIBSAKURA_LIBSAKURA_LOCALDEF_H_
#define LIBSAKURA_LIBSAKURA_LOCALDEF_H_

/** Don't include this header file from other header files.
 * This file includes the definitions only for internal use of libsakura.
 */

#define MARK_AS_USED(x) ((void)(0 && (x)))

#if defined(__AVX2__)
# define ARCH_AFTER_HASWELL 1
#endif

#if defined(__AVX__)
# define ARCH_AFTER_SANDY_BRIDGE 1
#endif

#if !defined(ARCH_AFTER_SANDY_BRIDGE) && !defined(ARCH_SCALAR)
# define ARCH_DEFAULT 1
#endif

#ifndef ARCH_SUFFIX
# define ARCH_SUFFIX Default
#endif
#define CONCAT_SYM(A, B) A ## B
#define ADDSUFFIX(A, B) CONCAT_SYM(A, B)

#define ELEMENTSOF(x) (sizeof(x) / sizeof((x)[0]))
#define STATIC_ASSERT(x) static_assert((x), # x)

#undef likely
#undef unlikely
#if defined(__GNUC__) && (__GNUC__ > 2) && defined(__OPTIMIZE__)
# define likely(expr) (__builtin_expect(!!(expr), 1))
# define unlikely(expr) (__builtin_expect(!!(expr), 0))
#else
# define likely(expr) (expr)
# define unlikely(expr) (expr)
#endif

/*
 * ALIGNMENT must be a power of 2.
 */
#define LIBSAKURA_ALIGNMENT (64u/* double(or mmx) 64bits */ / 8u)
#if defined(__SSE__)
# undef LIBSAKURA_ALIGNMENT
# define LIBSAKURA_ALIGNMENT (128u/* sse 128bits */ / 8u)
#endif
#if defined(__AVX__)
# undef LIBSAKURA_ALIGNMENT
# define LIBSAKURA_ALIGNMENT (256u/* avx 256bits */ / 8u)
#endif

#if defined(__cplusplus)

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>

#if defined(DISABLE_ALIGNAS) || defined(__INTEL_COMPILER) && __INTEL_COMPILER < 1600
#define SIMD_ALIGN /* nothing */
#else
#define SIMD_ALIGN alignas(LIBSAKURA_ALIGNMENT)
#endif

namespace {

inline bool IsAligned(void const*ptr, size_t alignemnt) {
	uintptr_t addr = (uintptr_t) ptr;
	return addr % alignemnt == 0;
}

class ScopeGuard {
	typedef std::function<void(void) noexcept> Func;
public:
	ScopeGuard() = delete;
	explicit ScopeGuard(Func clean_up, bool enabled = true) :
	clean_up_(clean_up), engaged_(enabled), called_(false) {
	}
	ScopeGuard(ScopeGuard const &other) = delete;
	ScopeGuard &operator =(ScopeGuard const &other) = delete;
	ScopeGuard const &operator =(ScopeGuard const &other) const = delete;
	void *operator new(std::size_t) = delete;

	/**
	 * @~
	 * @brief Calls @ clean_up parameter provided to the constructor
	 * if the ScopeGuard is enabled.
	 */
	~ScopeGuard() {
		if (engaged_) {
			assert(! called_);
			clean_up_();
			// called_ = true;
		}
	}
	/**
	 * @~
	 * @brief Disables the ScopeGuard
	 *
	 * Destructor won't call @a clean_up parameter provided to the constructor
	 * if the ScopeGuard is disabled.
	 */
	void Disable() {
		engaged_ = false;
	}
	/**
	 * @~
	 * @brief Enables the ScopeGuard if it has not cleaned up yet
	 *
	 * Don't call @a Enable() after explicit call of @ref CleanUpNow().
	 */
	void Enable() {
		assert(! called_);
		engaged_ = true;
	}

	/**
	 * @~
	 * @brief Calls @a clean_up parameter provided to the constructor
	 * if the ScopeGuard is enabled
	 *
	 * Calling this method more than once is not allowed.
	 */
	void CleanUpNow() {
		if (engaged_) {
			assert(! called_);
			clean_up_();
			called_ = true;
			engaged_ = false;
		}
	}
private:
	Func clean_up_;
	bool engaged_;
	bool called_;
};

#ifdef __has_builtin
  /* clang has this macro */
#else
  #define __has_builtin(x) 0  // Compatibility with non-clang compilers.
#endif

#if defined(__GNUG__) && (!defined(__clang__) || __has_builtin(__builtin_assume_aligned))
/**
 * @~japanese
 * @brief @a ptr がsakuraのアライメント要件を満たしていると見なし、そのアドレスを返す。
 *
 * コンパイラがサポートしていれば、
 * コンパイラは、戻り値がsakuraのアライメント要件を満たしているものとして最適化を行う。
 *
 * @param ptr sakuraのアライメント要件を満たしているアドレス
 * @param alignment アライメントサイズ
 * @return @a ptr (sakuraのアライメント要件を満たしているというコンパイラ依存の属性付き)
 */
template<typename T, size_t ALIGNMENT = LIBSAKURA_ALIGNMENT>
inline T AssumeAligned(T ptr) {
	return reinterpret_cast<T>(__builtin_assume_aligned(ptr, ALIGNMENT));
}
#else /* defined(__GNUG__) && (!defined(__clang__) || __has_builtin(__builtin_assume_aligned)) */

//#warning __builtin_assume_aligned is not supported.

/**
 * @~japanese
 * @brief @a ptr がsakuraのアライメント要件を満たしていると見なし、そのアドレスを返す。
 *
 * コンパイラがサポートしていれば、
 * コンパイラは、戻り値がsakuraのアライメント要件を満たしているものとして最適化を行う。
 *
 * @param ptr sakuraのアライメント要件を満たしているアドレス
 * @param alignment アライメントサイズ
 * @return @a ptr (sakuraのアライメント要件を満たしているというコンパイラ依存の属性付き)
 */
template<typename T, size_t ALIGNMENT = LIBSAKURA_ALIGNMENT>
inline /*alignas(LIBSAKURA_ALIGNMENT)*/T AssumeAligned(T ptr) {
	return ptr;
}
#endif /* defined(__GNUG__) && (!defined(__clang__) || __has_builtin(__builtin_assume_aligned)) */
}

#else /* defined(__cplusplus) */

/*#warning __builtin_assume_aligned is not supported.*/
#ifdef defined(__GNUG__) && (!defined(__clang__) || __has_builtin(__builtin_assume_aligned))
# define AssumeAligned(ptr) ((__typeof__(ptr))__builtin_assume_aligned((ptr), LIBSAKURA_ALIGNMENT))
#else /* defined(__GNUG__) && (!defined(__clang__) || __has_builtin(__builtin_assume_aligned)) */
# define AssumeAligned(ptr) (ptr)
#endif /* defined(__GNUG__) && (!defined(__clang__) || __has_builtin(__builtin_assume_aligned)) */

#endif /* defined(__cplusplus) */

#endif /* LIBSAKURA_LIBSAKURA_LOCALDEF_H_ */
