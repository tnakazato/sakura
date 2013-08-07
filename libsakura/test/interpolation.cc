#include <iostream>

#include <libsakura/sakura.h>

#include "gtest/gtest.h"

TEST(Interpolate1dFloatTest, InvalidType) {
	const size_t num_base = 2;
	const size_t num_interpolated = 5;
	double x_base[num_base] = {0.0, 1.0};
	float y_base[num_base] = {1.0, -1.0};
	double x_interpolated[num_interpolated] = {-1.0, 0.0, 0.1, 0.7, 1.5};
	float y_interpolated[num_interpolated];
	sakura_Status result = sakura_Interpolate1dFloat(sakura_InterpolationMethod_kNumMethod,
			0, num_base, x_base, y_base, num_interpolated, x_interpolated, y_interpolated);

	// Should return InvalidArgument status
	EXPECT_EQ(result, sakura_Status_kInvalidArgument) << "Interpolate1dFloat should fail!";
}

TEST(Interpolate1dFloatTest, Nearest) {
	const size_t num_base = 2;
	const size_t num_interpolated = 6;
	double x_base[num_base] = {0.0, 1.0};
	float y_base[num_base] = {1.0, -1.0};
	double x_interpolated[num_interpolated] = {-1.0, 0.0, 0.1, 0.5, 0.7, 1.5};
	float y_interpolated[num_interpolated];
	sakura_Status result = sakura_Interpolate1dFloat(sakura_InterpolationMethod_kNearest,
			0, num_base, x_base, y_base, num_interpolated, x_interpolated, y_interpolated);


	// Basic check whether function is completed or not
	EXPECT_EQ(result, sakura_Status_kOK) << "Interpolate1dFloat had any problems during execution.";

	// Value check
	for (size_t index = 0; index < num_interpolated; ++index) {
		float reference;
		if (x_interpolated[index] <= 0.5 * (x_base[0] + x_base[1])) {
			reference = y_base[0];
		}
		else {
			reference = y_base[1];
		}
		std::cout << "Expected value at index " << index << ": " << reference << std::endl;
		EXPECT_EQ(y_interpolated[index], reference)
		<< "interpolated value differs from expected value at " << index << ": "
		<< reference << ", " << y_interpolated[index];
	}
}
