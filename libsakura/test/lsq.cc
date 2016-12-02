/*
 * @SAKURA_LICENSE_HEADER_START@
 * Copyright (C) 2013-2016
 * National Astronomical Observatory of Japan
 * 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
 * 
 * This file is part of Sakura.
 * 
 * Sakura is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * Sakura is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
 * @SAKURA_LICENSE_HEADER_END@
 */
/*
 * lsq.cc
 *
 *  Created on: 2016/09/26
 *      Author: Wataru Kawasaki
 */

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits.h>
#include <map>
#include <memory>
#include <string>
#include <sys/time.h>
#include <random>
#include <stdio.h>

#include <libsakura/sakura.h>
#include <libsakura/localdef.h>
#include <libsakura/memory_manager.h>
#include "loginit.h"
#include "aligned_memory.h"
#include "gtest/gtest.h"
#include "baseline.h"
#include "testutil.h"

/* the number of elements in input/output array to test */
#define NUM_DATA 5
#define NUM_MODEL 3
#define NUM_REPEAT 20000

using namespace std;

typedef enum {
	ApiName_kLSQFitPolynomialFloat,
	ApiName_kLSQFitCubicSplineFloat,
	ApiName_kLSQFitSinusoidFloat,
	ApiName_kNumElements
} ApiName;

map<ApiName, string> ApiNameStr = { { ApiName_kLSQFitPolynomialFloat,
		"LSQFitPolynomialFloat" }, { ApiName_kLSQFitCubicSplineFloat,
		"LSQFitCubicSplineFloat" }, { ApiName_kLSQFitSinusoidFloat,
		"LSQFitSinusoidFloat" }, { ApiName_kNumElements, "" } };

typedef enum {
	TestCategory_kSuccessful,
	TestCategory_kPerformance,
	TestCategory_kFailure,
	TestCategory_kNumElements
} TestCategory;

typedef enum {
	ParamDataType_kInt,
	ParamDataType_kSizeT,
	ParamDataType_kUInt16T,
	ParamDataType_kFloat,
	ParamDataType_kDouble,
	ParamDataType_kBool,
	ParamDataType_kFloatPointer,
	ParamDataType_kDoublePointer,
	ParamDataType_kBoolPointer,
	ParamDataType_kContextPointer,
	ParamDataType_kFitStatusPointer,
	ParamDataType_kNumElements
} ParamDataType;

bool IsPointer(ParamDataType const &data_type) {
	bool res = false;
	if ((data_type == ParamDataType_kFloatPointer)
			|| (data_type == ParamDataType_kDoublePointer)
			|| (data_type == ParamDataType_kBoolPointer)
			|| (data_type == ParamDataType_kContextPointer)
			|| (data_type == ParamDataType_kFitStatusPointer))
		res = true;
	return res;
}

typedef enum {
	ParamCategory_kValid,
	ParamCategory_kInvalid,
	ParamCategory_kUndefined,
	ParamCategory_kNumElements
} ParamCategory;

typedef enum {
	ParamValueType_kLowerBound,
	ParamValueType_kUpperBound,
	ParamValueType_kInRange,
	ParamValueType_kAligned,
	ParamValueType_kValidPointer,
	ParamValueType_kTooLess,
	ParamValueType_kTooGreat,
	ParamValueType_kOutofRange,
	ParamValueType_kNaN,
	ParamValueType_kInf,
	ParamValueType_kNegativeInf,
	ParamValueType_kNull,
	ParamValueType_kNotAligned,
	ParamValueType_kUndefined,
	ParamValueType_kNumElements
} ParamValueType;

struct ParamAttr {
	//size_t param_id;
	string name;
	ParamDataType data_type;
	ParamCategory category;
	ParamValueType value_type;
};

struct TestCase {
	size_t test_id;
	string title;
	string desc;
	ApiName api_name;
	TestCategory category;
	size_t num_repeat;
	vector<ParamAttr> param;
	LIBSAKURA_SYMBOL(Status) expect_status;
	//Func compare_func;
	void AddParamAttr(string const name, ParamDataType const data_type,
			ParamCategory const category, ParamValueType const value_type) {
		param.push_back( { name, data_type, category, value_type });
	}
	size_t GetParamAttrIndex(string const param_name) {
		bool found = false;
		size_t idx;
		for (size_t i = 0; i < param.size(); ++i) {
			if (param[i].name == param_name) {
				idx = i;
				found = true;
				break;
			}
		}
		if (!found) {
			throw std::runtime_error("parameter not found.");
		}
		return idx;
	}
	ParamAttr& GetParamAttrByName(string const param_name) {
		return param[GetParamAttrIndex(param_name)];
	}
};

TestCase CreateDefaultTestCase(ApiName const api_name) {
	TestCase test_case;
	test_case.test_id = 0; //UINT_MAX;
	test_case.title = ApiNameStr[api_name];
	test_case.desc = "default case";
	test_case.api_name = api_name;
	test_case.category = TestCategory_kSuccessful;
	test_case.num_repeat = 1;
	test_case.expect_status = LIBSAKURA_SYMBOL(Status_kOK);

	test_case.param.clear();
	size_t num_params = 0;
	vector<string> param_name(0);
	vector<ParamDataType> param_dtype(0);
	vector<ParamCategory> param_category(0);
	switch (api_name) {
	case ApiName_kLSQFitPolynomialFloat:
		num_params = 14;
		param_name = {"context",
			"order", "num_data",
			"data", "mask",
			"clip_threshold_sigma", "num_fitting_max",
			"num_coeff", "coeff",
			"best_fit", "residual",
			"final_mask", "rms",
			"lsqfit_status"};
		param_dtype = {ParamDataType_kContextPointer,
			ParamDataType_kUInt16T, ParamDataType_kSizeT,
			ParamDataType_kFloatPointer, ParamDataType_kBoolPointer,
			ParamDataType_kFloat, ParamDataType_kUInt16T,
			ParamDataType_kSizeT, ParamDataType_kDoublePointer,
			ParamDataType_kFloatPointer, ParamDataType_kFloatPointer,
			ParamDataType_kBoolPointer, ParamDataType_kFloatPointer,
			ParamDataType_kFitStatusPointer};
		break;
		case ApiName_kLSQFitCubicSplineFloat:
		break;
		case ApiName_kLSQFitSinusoidFloat:
		break;
		default:
		assert(false);
		break;
	}
	assert(param_name.size() == num_params);
	assert(param_dtype.size() == num_params);
	for (size_t i = 0; i < num_params; ++i) {
		test_case.AddParamAttr(param_name[i], param_dtype[i],
				ParamCategory_kValid,
				(IsPointer(param_dtype[i]) ?
						ParamValueType_kValidPointer : ParamValueType_kInRange));
	}
	return test_case;
}

vector<TestCase> CreateTestCases(ApiName const api_name) {
	vector<TestCase> test_cases(0);
	TestCase default_test_case = CreateDefaultTestCase(api_name);
	//first one is for successful case with default values for all parameters
	test_cases.push_back(default_test_case);
	TestCase* tc;
	auto get_working_test_case = [&]() {
		test_cases.push_back(default_test_case);
		return &test_cases[test_cases.size()-1];};
	auto add_test_case = [&](string const title_suffix,
			string const test_description,
			TestCategory const test_category,
			string const param_name,
			ParamValueType const param_vtype,
			LIBSAKURA_SYMBOL(Status) const expect_status) {
		tc = get_working_test_case();
		tc->test_id = test_cases.size() - 1;
		tc->title = ApiNameStr[api_name] + "_" + title_suffix;
		tc->desc = test_description;
		tc->category = test_category;
		tc->expect_status = expect_status;
		tc->GetParamAttrByName(param_name).value_type = param_vtype;

	};
	switch (api_name) {
	case ApiName_kLSQFitPolynomialFloat:
		//#000b - contextVP (context for chebyshev)
		add_test_case("contextVP", "context for chebyshev",
				TestCategory_kSuccessful, "context",
				ParamValueType_kValidPointer, LIBSAKURA_SYMBOL(Status_kOK));
		//#001 - contextNULL
		add_test_case("contextNULL", "", TestCategory_kFailure, "context",
				ParamValueType_kNull,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		//#002 - contextOR (context for cspline)
		add_test_case("contextOR", "context for cspline", TestCategory_kFailure,
				"context", ParamValueType_kOutofRange,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		//#003 - contextOR (context for sinusoids)
		add_test_case("contextOR", "context for sinusoid",
				TestCategory_kFailure, "context", ParamValueType_kOutofRange,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		//#004 - orderLB
		add_test_case("orderLB", "order=0", TestCategory_kSuccessful, "order",
				ParamValueType_kLowerBound, LIBSAKURA_SYMBOL(Status_kOK));
		//#005 - orderTG (order is larger than the one used to create context)
		add_test_case("orderTG", "order=4", TestCategory_kFailure, "order",
				ParamValueType_kTooGreat,
				LIBSAKURA_SYMBOL(Status_kInvalidArgument));
		break;
	default:
		assert(false);
		break;
	}
	return test_cases;
}

template<typename T>
void AllocateAligned(size_t const length, T **data, void **data_storage) {
	*data = nullptr;
	unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(T) * length, data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(*data));
	*data_storage = storage_for_data.release();
}

template<typename T>
void Allocate(size_t const length, T **data) {
	*data = nullptr;
	unique_ptr<T, LIBSAKURA_PREFIX::Memory> work_data(
			static_cast<T *>(LIBSAKURA_PREFIX::Memory::Allocate(
					sizeof(T) * length)), LIBSAKURA_PREFIX::Memory());
	if (work_data == nullptr) {
		throw bad_alloc();
	}
	*data = work_data.release();
}

struct ParamSet {
	map<string, size_t> size;
	map<string, uint16_t> ui16;
	map<string, float> fval;
	map<string, float *> fptr;
	map<string, double *> dptr;
	map<string, bool *> bptr;
	map<string, LIBSAKURA_SYMBOL(LSQFitContextFloat) *> ctxt;
	map<string, LIBSAKURA_SYMBOL(LSQFitStatus) *> fsta;
	map<string, void *> sto;

	void AllocateAlignedFloat(string const &name, string const &length_name) {
		AllocateAligned(size[length_name], &fptr[name], &sto[name]);
	}
	void AllocateAlignedDouble(string const &name, string const &length_name) {
		AllocateAligned(size[length_name], &dptr[name], &sto[name]);
	}
	void AllocateAlignedBool(string const &name, string const &length_name) {
		AllocateAligned(size[length_name], &bptr[name], &sto[name]);
	}
	void AllocateFloat(string const &name) {
		Allocate(1, &fptr[name]);
	}
	void AllocateStatus(string const &name) {
		Allocate(1, &fsta[name]);
	}
	void AllocateContexts(string const &name, string const &defaultvalue_name) {
		LIBSAKURA_SYMBOL(Status) create_status;
		uint16_t order = ui16["order"];
		size_t num_data = size["num_data"];

		//polynomial
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_poly = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextPolynomialFloat)(
				LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, num_data,
				&context_poly);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		ctxt["poly"] = context_poly;

		//chebyshev
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_chebyshev = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextPolynomialFloat)(
				LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
				&context_chebyshev);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		ctxt["chebyshev"] = context_chebyshev;

		//cspline
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_cspline = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextCubicSplineFloat)(order, num_data,
				&context_cspline);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		ctxt["cspline"] = context_cspline;

		//sinusoid
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_sinusoid = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextSinusoidFloat)(order, num_data,
				&context_sinusoid);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
		ctxt["sinusoid"] = context_sinusoid;

		ctxt[name] = ctxt[defaultvalue_name];
	}
	void Free(string const &name) {
		if (sto.count(name) == 1) {
			if (sto[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(sto[name]);
			}
		} else if (fptr.count(name) == 1) {
			if (fptr[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(fptr[name]);
			}
		} else if (dptr.count(name) == 1) {
			if (dptr[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(dptr[name]);
			}
		} else if (bptr.count(name) == 1) {
			if (bptr[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(bptr[name]);
			}
		} else if (fsta.count(name) == 1) {
			if (fsta[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(fsta[name]);
			}
		} else if (ctxt.count(name) == 1) {
			if (ctxt[name] != nullptr) {
				LIBSAKURA_SYMBOL(Status) status;
				status = LIBSAKURA_SYMBOL(DestroyLSQFitContextFloat)(
						ctxt[name]);
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
			}
		}
	}
};

ParamSet CreateDefaultParameterSet(ApiName const api_name) {
	ParamSet ps;
	switch (api_name) {
	case ApiName_kLSQFitPolynomialFloat:
		ps.ui16["order"] = 3;
		ps.size["num_data"] = 15;
		ps.AllocateContexts("context", "poly");
		ps.AllocateAlignedFloat("data", "num_data");
		ps.AllocateAlignedBool("mask", "num_data");
		for (size_t i = 0; i < ps.size["num_data"]; ++i) {
			ps.fptr["data"][i] = 10.0;
			ps.bptr["mask"][i] = true;
		}
		ps.fval["clip_threshold_sigma"] = 3.0;
		ps.ui16["num_fitting_max"] = 2;
		ps.size["num_coeff"] = ps.ui16["order"] + 1;
		ps.AllocateAlignedDouble("coeff", "num_coeff");
		ps.AllocateAlignedFloat("best_fit", "num_data");
		ps.AllocateAlignedFloat("residual", "num_data");
		ps.AllocateAlignedBool("final_mask", "num_data");
		ps.AllocateFloat("rms");
		ps.AllocateStatus("lsqfit_status");
		break;
	default:
		assert(false);
		break;
	}
	return ps;
}

//interpret test cases, setup all parameter values to be used for testing
void Prologue(TestCase const &tc, TestCase const &default_tc, ParamSet &ps) {
	cout << "    {" << tc.title << ": " << tc.desc << "}   ";
	//cout << endl;
	for (size_t j = 0; j < tc.param.size(); ++j) {
		//cout << "(" << tc.param[j].name << ")" << "[" << tc.param[j].value_type << "]   ";
	}
	//cout << endl;

	ps = CreateDefaultParameterSet(tc.api_name);
	ParamAttr param;
	ParamAttr default_param;
	string name;
	ParamDataType dtype;
	ParamValueType vtype;
	for (size_t i = 0; i < tc.param.size(); ++i) {
		param = tc.param[i];
		name = param.name;
		dtype = param.data_type;
		default_param = default_tc.param[i];
		assert(name == default_param.name);
		assert(dtype == default_param.data_type);
		vtype = param.value_type;
		if (name == "context") {
			assert(dtype == ParamDataType_kContextPointer);
			if (vtype == ParamValueType_kNull) {
				ps.ctxt[name] = nullptr;
			} else {
				if (tc.desc.find("poly") != string::npos) {
					ps.ctxt[name] = ps.ctxt["poly"];
				} else if (tc.desc.find("chebyshev") != string::npos) {
					ps.ctxt[name] = ps.ctxt["chebyshev"];
				} else if (tc.desc.find("cspline") != string::npos) {
					ps.ctxt[name] = ps.ctxt["cspline"];
				} else if (tc.desc.find("sinusoid") != string::npos) {
					ps.ctxt[name] = ps.ctxt["sinusoid"];
				}
			}
		} else if (name == "order") {
			if (vtype == ParamValueType_kLowerBound) {
				ps.ui16[name] = 0;
			} else if (vtype == ParamValueType_kTooGreat) {
				ps.ui16[name] = 4;
			}
			ps.size["num_coeff"] = ps.ui16["order"] + 1;
			ps.Free("coeff");
			ps.AllocateAlignedDouble("coeff", "num_coeff");
		}
	}
}

void RunApi(TestCase const &tc, ParamSet &ps,
LIBSAKURA_SYMBOL(Status) &run_status) {
	switch (tc.api_name) {
	case ApiName_kLSQFitPolynomialFloat:
		run_status = LIBSAKURA_SYMBOL(LSQFitPolynomialFloat)(
				ps.ctxt[tc.param[0].name], ps.ui16[tc.param[1].name],
				ps.size[tc.param[2].name], ps.fptr[tc.param[3].name],
				ps.bptr[tc.param[4].name], ps.fval[tc.param[5].name],
				ps.ui16[tc.param[6].name], ps.size[tc.param[7].name],
				ps.dptr[tc.param[8].name], ps.fptr[tc.param[9].name],
				ps.fptr[tc.param[10].name], ps.bptr[tc.param[11].name],
				ps.fptr[tc.param[12].name], ps.fsta[tc.param[13].name]);
		break;
	default:
		break;
	}
}

void CheckStatus(TestCase const &tc, LIBSAKURA_SYMBOL(Status) const &status) {
	EXPECT_EQ(tc.expect_status, status);
}

void CheckValues(TestCase const &tc, ParamSet &ps) {
	cout << "          ****[";
	for (size_t j = 0; j < ps.size["num_coeff"]; ++j) {
		cout << ps.dptr["coeff"][j] << ", ";
	}
	cout << "]****" << endl;
}

void Execute(TestCase const &tc, ParamSet &ps) {
	double time_elapsed = 0.0;
	LIBSAKURA_SYMBOL (Status) run_status;
	for (size_t i = 0; i < tc.num_repeat; ++i) {
		cout << "####[" << tc.num_repeat << "]" << "{" << ps.ui16["order"]
				<< "}" << "<" << ps.ctxt["context"] << ">" << "####" << endl;
		double time_start = GetCurrentTime();
		RunApi(tc, ps, run_status);
		double time_end = GetCurrentTime();
		time_elapsed += (time_end - time_start);

		CheckStatus(tc, run_status);
		if ((tc.category == TestCategory_kSuccessful)
				|| (tc.category == TestCategory_kPerformance)) {
			CheckValues(tc, ps);
		}
	}

	if (tc.category == TestCategory_kPerformance) {
		cout << setprecision(5) << "#x# benchmark Lsq_"
				<< ApiNameStr[tc.api_name] << " " << time_elapsed << endl;
	}
}

void FreeParamSet(ParamSet &ps) {
	for (auto it = ps.sto.begin(); it != ps.sto.end(); ++it) {
		auto key = it->first;
		if (ps.sto[key] != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(ps.sto[key]);
		}
	}
	for (auto it = ps.fptr.begin(); it != ps.fptr.end(); ++it) {
		auto key = it->first;
		if ((ps.fptr[key] != nullptr) && (ps.sto.count(key) == 0)) {
			LIBSAKURA_PREFIX::Memory::Free(ps.fptr[key]);
		}
	}
	for (auto it = ps.dptr.begin(); it != ps.dptr.end(); ++it) {
		auto key = it->first;
		if ((ps.dptr[key] != nullptr) && (ps.sto.count(key) == 0)) {
			LIBSAKURA_PREFIX::Memory::Free(ps.dptr[key]);
		}
	}
	for (auto it = ps.bptr.begin(); it != ps.bptr.end(); ++it) {
		auto key = it->first;
		if ((ps.bptr[key] != nullptr) && (ps.sto.count(key) == 0)) {
			LIBSAKURA_PREFIX::Memory::Free(ps.bptr[key]);
		}
	}
	for (auto it = ps.fsta.begin(); it != ps.fsta.end(); ++it) {
		auto key = it->first;
		if (ps.fsta[key] != nullptr) {
			LIBSAKURA_PREFIX::Memory::Free(ps.fsta[key]);
		}
	}
	if (ps.ctxt.size() > 1) {
		for (auto it = ps.ctxt.begin(); it != ps.ctxt.end(); ++it) {
			//ctxt["context"] need not to be free since it is a copy of one of the others.
			auto key = it->first;
			auto p = ps.ctxt[key];
			if ((p != nullptr) && (key != "context")) {
				LIBSAKURA_SYMBOL(Status) status =
				LIBSAKURA_SYMBOL(DestroyLSQFitContextFloat)(p);
				EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
				p = nullptr;
			}
		}
	} else if (ps.ctxt.size() == 1) {
		LIBSAKURA_SYMBOL(Status) status =
		LIBSAKURA_SYMBOL(DestroyLSQFitContextFloat)(ps.ctxt["context"]);
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}
}

void Epilogue(ParamSet &ps) {
	FreeParamSet(ps);
}

void RunTest(ApiName const api_name) {
	vector<TestCase> test_cases = CreateTestCases(api_name);
	for (size_t i = 0; i < test_cases.size(); ++i) {
		ParamSet ps;
		Prologue(test_cases[i], test_cases[0], ps);
		Execute(test_cases[i], ps);
		Epilogue(ps);
	}
}

/*
 * A super class to test lsq functions
 */
class Lsq: public ::testing::Test {
protected:

	Lsq() :
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

	//Set (A[0]+A[1]*x+A[2]*x*x+A[3]*x*x*x) float values into an array
	void SetFloatPolynomial(size_t num_data, float *data,
			double *coeff_answer) {
		for (size_t i = 0; i < num_data; ++i) {
			double x = (double) i;
			data[i] = (float) (coeff_answer[0] + coeff_answer[1] * x
					+ coeff_answer[2] * x * x + coeff_answer[3] * x * x * x);
		}
	}

	template<typename T>
	void SetConstant(T value, size_t const num_data, T *data) {
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

	bool verbose;
};

//Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
void Destroy(LIBSAKURA_SYMBOL(LSQFitContextFloat) *context, size_t status) {
	LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyLSQFitContextFloat(
			context);
	EXPECT_EQ(status, destroy_status);
	cout << "Destroy Status : " << destroy_status << endl;
}

/*
 * Test sakura_CreateLSQFitContextFloatWithPolynomialPerformanceTest
 * successful case (with normal polynomial model)
 */
/*
 TEST_F(Baseline, CreateLSQFitContextFloatWithPolynomialPerformanceTest) {
 uint16_t const order(1000);
 size_t const num_chan(65535);

 LIBSAKURA_SYMBOL(LSQFitContextFloat) * context = nullptr;

 size_t const num_repeat(1);
 double start, end;
 double elapsed_time = 0.0;
 for (size_t i = 0; i < num_repeat; ++i) {
 start = GetCurrentTime();
 LIBSAKURA_SYMBOL (Status) create_status =
 sakura_CreateLSQFitContextPolynomialFloat(
 LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order,
 num_chan, &context);
 end = GetCurrentTime();
 elapsed_time += (end - start);
 EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), create_status);
 Destroy(context, LIBSAKURA_SYMBOL(Status_kOK));
 }
 cout << setprecision(5)
 << "#x# benchmark Baseline_CreateLSQFitContextPolynomialFloatWithPolynomialPerformanceTest"
 << " " << elapsed_time << endl;
 }
 */

/*
 * Tests for LSQFitPolynomialFloat()
 */
TEST_F(Lsq, LSQFitPolynomialFloat) {
	RunTest(ApiName_kLSQFitPolynomialFloat);
}
