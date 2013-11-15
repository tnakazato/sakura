/*
 * localdef.h
 *
 *  Created on: 2013/02/22
 *      Author: kohji
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

#ifndef ARCH_SUFFIX
# define ARCH_SUFFIX Default
#endif
#define CONCAT_SYM(A, B) A ## B
#define ADDSUFFIX(A, B) CONCAT_SYM(A, B)

#define ELEMENTSOF(x) (sizeof(x) / sizeof((x)[0]))
#define STATIC_ASSERT(x) static_assert((x), # x)

/*
 * ALIGNMENT must be a power of 2.
 */
#define LIBSAKURA_ALIGNMENT (256u/* avx 256bits */ / 8u)

#if defined(__cplusplus)

#include <cassert>
#include <cstddef>
#include <functional>

namespace {

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

#ifdef __GNUG__
template<typename T>
inline T AssumeAligned(T ptr) {
	return reinterpret_cast<T>(__builtin_assume_aligned(ptr,
	LIBSAKURA_ALIGNMENT));
}
#else /* __GNUG__ */
template<typename T>
inline /*alignas(LIBSAKURA_ALIGNMENT)*/T *AssumeAligned(T *ptr) {
	return ptr;
}
#endif /* __GNUG__ */
}

#else /* defined(__cplusplus) */

#ifdef __GNUC__
# define AssumeAligned(ptr) ((__typeof__(ptr))__builtin_assume_aligned((ptr), LIBSAKURA_ALIGNMENT))
#else /* __GNUC__ */
# define AssumeAligned(ptr) (ptr)
#endif /* __GNUC__ */

#endif /* defined(__cplusplus) */

#endif /* LIBSAKURA_LIBSAKURA_LOCALDEF_H_ */
