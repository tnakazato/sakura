#include <iostream>

#include <libsakura/sakura.h>

#include "gtest/gtest.h"

TEST(BasicInterpolationTest, VerifyCalled) {
	double *x_base;
	float *y_base;
	double *x_interpolated;
	float *y_interpolated;
	sakura_Status result = sakura_Interpolate1dFloat(sakura_InterpolationMethod_kNearest,
			0, 0, x_base, y_base, 0, x_interpolated, y_interpolated);
	EXPECT_EQ(result, sakura_Status_kOK);

}
