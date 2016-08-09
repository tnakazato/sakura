/*
 * testutil.h
 *
 *  Created on: Aug 8, 2016
 *      Author: nakazato
 */

#ifndef _LIBSAKURA_TEST_TESTUTIL_H_
#define _LIBSAKURA_TEST_TESTUTIL_H_

#ifdef __cplusplus
extern "C" {

#if __cplusplus >= 201103L
# undef LIBSAKURA_NOEXCEPT
# define LIBSAKURA_NOEXCEPT noexcept
#endif

#endif

double GetCurrentTime() LIBSAKURA_NOEXCEPT;

#ifdef __cplusplus
}
/* extern "C" */
#endif

#endif /* _LIBSAKURA_TEST_TESTUTIL_H_ */
