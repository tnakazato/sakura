/*
 * localdef.h
 *
 *  Created on: 2013/02/22
 *      Author: kohji
 */

#ifndef LOCALDEF_H_
#define LOCALDEF_H_

/* Don't include this header file from other header files.
 * This file includes the definitions only for internal use of libsakura.
 */

#ifndef ARCH_SUFFIX
# define ARCH_SUFFIX Default
#endif
#define CONCAT_SYM(A, B) A ## B
#define ADDSUFFIX(A, B) CONCAT_SYM(A, B)

#define elementsof(x) (sizeof(x) / sizeof((x)[0]))

#endif /* LOCALDEF_H_ */
