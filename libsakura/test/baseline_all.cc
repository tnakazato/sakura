
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <sys/time.h>
#include <random>
#include <stdio.h>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"
#include "baseline.h"

/* the number of elements in input/output array to test */
#define NUM_DATA 5
#define NUM_DATA2 15
#define NUM_DATA3 4096
#define NUM_MODEL 3
#define NUM_REPEAT 20000
#define NUM_MODEL2 4
#define NUM_REPEAT2 2000000

using namespace std;


/*
 * A super class to test baseline functions
 */
class Baseline: public ::testing::Test {
protected:

	Baseline() :
			verbose(false) {
	}

	virtual void SetUp() {
		// Initialize sakura
		LIBSAKURA_SYMBOL (Status) status = LIBSAKURA_SYMBOL(Initialize)(nullptr,
				nullptr);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}

	virtual void TearDown() {
		LIBSAKURA_SYMBOL(CleanUp)();
	}





	// Set (1+x+x*x) float values into an array
	void SetFloatPolynomial(size_t const num_data, float *data) {
		for (size_t i = 0; i < num_data; ++i) {
			double x = (double) i;
			data[i] = (float) (1.0 + x + x * x);
		}
	}

	//Set (A[0]+A[1]*x+A[2]*x*x+A[3]*x*x*x) float values into an array
	void SetFloatPolynomial(size_t num_data, float *data,
			double *coeff_answer) {
		for (size_t i = 0; i < num_data; ++i) {
			double x = (double) i;
			data[i] = (float) (coeff_answer[0] + coeff_answer[1] * x
					+ coeff_answer[2] * x * x + coeff_answer[3] * x * x * x);
		}
	}


	//Set constat any type values into an array
	template<typename T>
	void Set_Value_Constant(T value, size_t const num_data, T *data) {
		for (size_t i = 0; i < num_data; ++i) {
			data[i] = value;
		}
	}


	// Check if the expected and actual values are enough close to each other
	void CheckAlmostEqual(double expected, double actual, double tolerance) {
		double deviation = fabs(actual - expected);
		double val = max(fabs(actual), fabs(expected)) * tolerance + tolerance;
		ASSERT_LE(deviation, val);
	}

	template<typename T>
	void Print1DArray(char const *name, size_t print_length, T const *data,
			size_t start_idx = 0, bool print_name = true, bool newline = true) {
		if (print_name)
			cout << name << " = [";
		size_t index = start_idx + print_length - 1;
		if (typeid(float const*) == typeid(data)) {
			for (size_t i = start_idx; i < index; ++i)
				cout << data[i] << ", ";
		} else if (typeid(double const*) == typeid(data)) {
			for (size_t i = start_idx; i < index; ++i)
				cout << data[i] << ", ";
		} else if (typeid(bool const*) == typeid(data)) {
			for (size_t i = start_idx; i < index; ++i)
				cout << (data[i] ? "T" : "F") << ", ";
		}
		cout << data[index] << " ]";
		if (newline)
			cout << endl;
	}

	//given as 1D float array but actually stores (num_row * num_column) 2D data
	//for which column loop comes inside row loop.
	template<typename T>
	void PrintArray(char const *name, size_t num_row, size_t num_column,
			T const *data) {
		cout << name << " = [";
		for (size_t i = 0; i < num_row; ++i) {
			Print1DArray(name, num_column, data, num_column * i, false, false);
			if (i < num_row - 1)
				cout << ", " << endl;
		}
		cout << " ]" << endl;
	}

	bool verbose;
};

//Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
void Destroy(LIBSAKURA_SYMBOL(BaselineContextFloat) *context, size_t status) {
	LIBSAKURA_SYMBOL (Status) destroy_status =
			sakura_DestroyBaselineContextFloat(context);
	EXPECT_EQ(status, destroy_status);
	//std::cout << "Destroy Status : " << destroy_status << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
// TEST
////////////////////////////////////////////////////////////////////////////////

/*
 * Test sakura_DestroyBaselineContextFloat
 * successful case
 */
TEST_F(Baseline, DestroyBaselineContextFloat) {
	uint16_t const order(20);
	size_t const num_chan(4096);

	double start, end;
	double elapsed_time = 0.0;
	size_t const num_repeat(NUM_REPEAT);
	for (size_t i = 0; i < num_repeat; ++i) {
		LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;
		LIBSAKURA_SYMBOL (Status) create_status =
				sakura_CreateBaselineContextFloat(
						LIBSAKURA_SYMBOL(BaselineType_kPolynomial), order,
						num_chan, &context);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);

		start = LIBSAKURA_SYMBOL(GetCurrentTime)();

		//Destroy
		Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));

		end = LIBSAKURA_SYMBOL(GetCurrentTime)();
		elapsed_time += (end - start);
	}
	cout << "Elapsed Time: " << elapsed_time << " sec." << endl;
}

/*
 * Test sakura_DestroyBaselineContextFloat
 * failure case : context is a null pointer
 * returned value must be Status_kInvalidArgument
 */
TEST_F(Baseline, DestroyBaselineContextFloatWithContextNullPointer) {
	LIBSAKURA_SYMBOL(BaselineContextFloat) * context = nullptr;

	Destroy(context, LIBSAKURA_SYMBOL(Status_kInvalidArgument));
}


