/*
 * localdef.h
 *
 *  Created on: 2013/02/22
 *      Author: kohji
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

#if defined(__cplusplus)

#include <cassert>
#include <cstddef>
#include <functional>

namespace {

class ScopeGuard {
	typedef std::function<void(void) noexcept> Func;
public:
	ScopeGuard() = delete;
	explicit ScopeGuard(Func clean_up, bool engaged = true) :
			clean_up_(clean_up), engaged_(engaged), called_(false) {
	}
	ScopeGuard(ScopeGuard const &other) = delete;
	ScopeGuard &operator =(ScopeGuard const &other) = delete;
	ScopeGuard const &operator =(ScopeGuard const &other) const = delete;
	void *operator new(std::size_t) = delete;

	~ScopeGuard() {
		if (engaged_) {
			assert(! called_);
			clean_up_();
			// called_ = true;
		}
	}
	void Disable() {
		engaged_ = false;
	}
	void Enable() {
		assert(! called_);
		engaged_ = true;
	}

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

}

#endif /* defined(__cplusplus) */

#endif /* LIBSAKURA_LIBSAKURA_LOCALDEF_H_ */
