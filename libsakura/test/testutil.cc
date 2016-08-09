/*
 * testutil.cc
 *
 *  Created on: Aug 9, 2016
 *      Author: nakazato
 */

#include <chrono>

extern "C" double GetCurrentTime() noexcept {
	using namespace std::chrono;
	auto now = system_clock::now().time_since_epoch();
	return duration_cast<duration<double, seconds::period>>(now).count();
}

