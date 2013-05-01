/*
 * init.cc
 *
 *  Created on: 2013/05/01
 *      Author: kohji
 */

#include <libsakura/sakura.h>
#include "gtest/gtest.h"

TEST(Global, Init) {
	sakura_Status result = sakura_Initialize();
	EXPECT_EQ(result, sakura_Status_kOK);

	sakura_CleanUp();

	result = sakura_Initialize();
	EXPECT_EQ(result, sakura_Status_kOK);

	sakura_CleanUp();
}
