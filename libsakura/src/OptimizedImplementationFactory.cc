// Author: Kohji Nakamura <k.nakamura@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution

#include <cstdio>
#include <stdint.h>
#include "libsakura/OptimizedImplementationFactory.h"
#include "libsakura/OptimizedImplementationFactoryImpl.h"

#define elementsof(x) (sizeof(x) / sizeof((x)[0]))

using namespace std;
using namespace sakura;

namespace {
  struct reg_t {
    uint32_t eax, ebx, ecx, edx;
  };

  void get_cpuid(reg_t &reg) {
    __asm__ volatile("cpuid"
		     : "=a" (reg.eax),
		       "=b" (reg.ebx),
		       "=c" (reg.ecx),
		       "=d" (reg.edx)
		     : "a" (reg.eax)
		     );
  }

  struct simd_feature_t {
    bool mmx, sse, sse2, sse3, ssse3, sse4_1, sse4_2, avx;
  };

  struct simd_flags_t {
    uint32_t ecx, edx;
    bool simd_feature_t::*flag;
  } simd_flags[] = {
    {0, 0x0800000, &simd_feature_t::mmx},
    {0, 0x2000000, &simd_feature_t::sse},
    {0, 0x4000000, &simd_feature_t::sse2},
    {0x00000001, 0, &simd_feature_t::sse3},
    {0x00000200, 0, &simd_feature_t::ssse3},
    {0x00080000, 0, &simd_feature_t::sse4_1},
    {0x00100000, 0, &simd_feature_t::sse4_2},
    {0x10000000, 0, &simd_feature_t::avx},
  };

  void get_cpu_feature(simd_feature_t &simd_feature) {
    reg_t reg;
    reg.eax = 1;
    get_cpuid(reg);
#if 0
    printf("eax: %08x\n", reg.eax);
    printf("ebx: %08x\n", reg.ebx);
    printf("ecx: %08x\n", reg.ecx);
    printf("edx: %08x\n", reg.edx);
#endif

    static simd_feature_t const simd_all_false = {};
    simd_feature = simd_all_false;
    for (int i = 0; i < elementsof(simd_flags); ++i) {
      if ((reg.ecx & simd_flags[i].ecx) == simd_flags[i].ecx &&
	  (reg.edx & simd_flags[i].edx) == simd_flags[i].edx) {
	simd_feature.*simd_flags[i].flag = true;
      }
    }
  }

  GriddingDefault griddingDefault;

  class OptimizedImplementationFactoryDefault:
    public OptimizedImplementationFactory {
  public:
    virtual Gridding *getGriddingImpl() {
      return &griddingDefault;
    }
  } defaultFactory;

  GriddingAfterSandyBridge griddingAfterSandyBridge;

  class OptimizedImplementationFactoryAfterSandyBridge:
    public OptimizedImplementationFactory {
  public:
    virtual Gridding *getGriddingImpl() {
      return &griddingAfterSandyBridge;
    }
  } afterSandyBridge;
}

namespace sakura {
  OptimizedImplementationFactory *
  OptimizedImplementationFactory::getFactory() {
    simd_feature_t simd_feature;
    get_cpu_feature(simd_feature);

    if (simd_feature.avx) {
      return &afterSandyBridge;
    }
    return &defaultFactory;
  }
}
