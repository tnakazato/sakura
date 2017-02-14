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
#include <random>
#include <sstream>
#include <stdio.h>
#include <string>
#include <sys/time.h>

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

//Aliases for Returned Status
LIBSAKURA_SYMBOL(Status) Ret_kOK = LIBSAKURA_SYMBOL(Status_kOK);
LIBSAKURA_SYMBOL(Status) Ret_kNG = LIBSAKURA_SYMBOL(Status_kNG);
LIBSAKURA_SYMBOL(Status) Ret_kIA = LIBSAKURA_SYMBOL(Status_kInvalidArgument);
LIBSAKURA_SYMBOL(Status) Ret_kNM = LIBSAKURA_SYMBOL(Status_kNoMemory);
LIBSAKURA_SYMBOL(Status) Ret_kUE = LIBSAKURA_SYMBOL(Status_kUnknownError);

typedef enum {
	ApiName_kCreateLSQFitContextPolynomialFloat,
	ApiName_kCreateLSQFitContextCubicSplineFloat,
	ApiName_kCreateLSQFitContextSinusoidFloat,
	ApiName_kDestroyLSQFitContextFloat,
	ApiName_kGetNumberOfCoefficientsFloat,
	ApiName_kLSQFitPolynomialFloat,
	ApiName_kLSQFitCubicSplineFloat,
	ApiName_kLSQFitSinusoidFloat,
	ApiName_kSubtractPolynomialFloat,
	ApiName_kSubtractCubicSplineFloat,
	ApiName_kSubtractSinusoidFloat,
	ApiName_kNumElements
} ApiName;

map<ApiName, string> ApiNameStr = { {
		ApiName_kCreateLSQFitContextPolynomialFloat,
		"CreateLSQFitContextPolynomialFloat" }, {
		ApiName_kCreateLSQFitContextCubicSplineFloat,
		"CreateLSQFitContextCubicSplineFloat" }, {
		ApiName_kCreateLSQFitContextSinusoidFloat,
		"CreateLSQFitContextSinusoidFloat" }, {
		ApiName_kDestroyLSQFitContextFloat, "DestroyLSQFitContextFloat" }, {
		ApiName_kGetNumberOfCoefficientsFloat, "GetNumberOfCoefficientsFloat" },
		{ ApiName_kLSQFitPolynomialFloat, "LSQFitPolynomialFloat" }, {
				ApiName_kLSQFitCubicSplineFloat, "LSQFitCubicSplineFloat" }, {
				ApiName_kLSQFitSinusoidFloat, "LSQFitSinusoidFloat" }, {
				ApiName_kSubtractPolynomialFloat, "SubtractPolynomialFloat" }, {
				ApiName_kSubtractCubicSplineFloat, "SubtractCubicSplineFloat" },
		{ ApiName_kSubtractSinusoidFloat, "SubtractSinusoidFloat" }, {
				ApiName_kNumElements, "" } };

//Test Category
typedef enum {
	TCat_kOK, TCat_kPF, TCat_kNG, TCat_kNumElements
} TCat;

typedef enum {
	PDType_kInt,
	PDType_kSizeT,
	PDType_kUInt16T,
	PDType_kFloat,
	PDType_kDouble,
	PDType_kBool,
	PDType_kFloatPointer,
	PDType_kDoublePointer,
	PDType_kBoolPointer,
	PDType_kContextPointer,
	PDType_kFitStatusPointer,
	PDType_kNumElements
} PDType;

vector<string> LsqFuncTypeStr = { "poly", "chebyshev", "cspline", "sinusoid" };

bool IsPointer(PDType const &data_type) {
	bool res = false;
	if ((data_type == PDType_kFloatPointer)
			|| (data_type == PDType_kDoublePointer)
			|| (data_type == PDType_kBoolPointer)
			|| (data_type == PDType_kContextPointer)
			|| (data_type == PDType_kFitStatusPointer))
		res = true;
	return res;
}

typedef enum {
	PCat_kValid, PCat_kInvalid, PCat_kUndefined, PCat_kNumElements
} PCat;

typedef enum {
	PVType_kLB, //LowerBound
	PVType_kUB, //UpperBound
	PVType_kIR, //InRange
	PVType_kAL, //Aligned
	PVType_kVP, //ValidPointer
	PVType_kTL, //TooLess
	PVType_kTG, //TooGreat
	PVType_kOR, //OutofRange
	PVType_kNAN, //NaN
	PVType_kINF, //Inf
	PVType_kNINF, //NegativeInf
	PVType_kNULL, //Null
	PVType_kNAL, //NotAligned
	PVType_kUNDEF, //Undefined
	PVType_kNumElements
} PVType;

struct PAttr {
	//size_t param_id;
	string name;
	PDType data_type;
	PCat category;
	PVType value_type;
};

struct TestCase {
	size_t test_id;
	string title;
	string desc;
	ApiName api_name;
	TCat category;
	size_t num_repeat;
	vector<PAttr> param;LIBSAKURA_SYMBOL(Status) expect_status;

	//Func compare_func;
	void AddParamAttr(string const name, PDType const data_type,
			PCat const param_category, PVType const value_type) {
		param.push_back( { name, data_type, param_category, value_type });
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
	PAttr& GetParamAttrByName(string const param_name) {
		return param[GetParamAttrIndex(param_name)];
	}
};

TestCase CreateDefaultTestCase(ApiName const api_name) {
	TestCase test_case;
	test_case.test_id = 0; //UINT_MAX;
	test_case.title = ApiNameStr[api_name];
	test_case.desc = "default case";
	test_case.api_name = api_name;
	test_case.category = TCat_kOK;
	test_case.num_repeat = 1;
	test_case.expect_status = Ret_kOK;
	test_case.param.clear();
	size_t num_params = 0;
	vector<string> param_name(0);
	vector<PDType> param_dtype(0);
	vector<PCat> param_category(0);

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
		param_dtype = {PDType_kContextPointer,
			PDType_kUInt16T, PDType_kSizeT,
			PDType_kFloatPointer, PDType_kBoolPointer,
			PDType_kFloat, PDType_kUInt16T,
			PDType_kSizeT, PDType_kDoublePointer,
			PDType_kFloatPointer, PDType_kFloatPointer,
			PDType_kBoolPointer, PDType_kFloatPointer,
			PDType_kFitStatusPointer};
		break;
		case ApiName_kLSQFitCubicSplineFloat:
		num_params = 13;
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
		test_case.AddParamAttr(param_name[i], param_dtype[i], PCat_kValid,
				(IsPointer(param_dtype[i]) ? PVType_kVP : PVType_kIR));
	}
	return test_case;
}

vector<string> Split(string const &s, char const delim) {
	vector<string> elems;
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		if (!item.empty()) {
			elems.push_back(item);
		}
	}
	return elems;
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
			TCat const test_category,
			string const param_name,
			PVType const param_vtype,
			LIBSAKURA_SYMBOL(Status) const expect_status) {
		tc = get_working_test_case();
		tc->test_id = test_cases.size() - 1;
		tc->title = ApiNameStr[api_name] + "_" + title_suffix;
		tc->desc = test_description;
		tc->category = test_category;
		vector<string> names = Split(param_name, ',');
		for (size_t i = 0; i < names.size(); ++i) {
			tc->GetParamAttrByName(names[i]).value_type = param_vtype;
		}
		tc->expect_status = expect_status;
	};
	switch (api_name) {
	case ApiName_kLSQFitPolynomialFloat:
		//#000b - contextVP (context for chebyshev)
		add_test_case("contextVP", "context for chebyshev", TCat_kOK, "context",
				PVType_kVP, Ret_kOK);
		//#001 - contextNULL
		add_test_case("contextNULL", "", TCat_kNG, "context", PVType_kNULL,
				Ret_kIA);
		//#002 - contextOR (context for cspline)
		add_test_case("contextOR", "context for cspline", TCat_kNG, "context",
				PVType_kOR, Ret_kIA);
		//#003 - contextOR (context for sinusoids)
		add_test_case("contextOR", "context for sinusoid", TCat_kNG, "context",
				PVType_kOR, Ret_kIA);
		//#004 - orderLB
		add_test_case("orderLB", "order=0", TCat_kOK, "order", PVType_kLB,
				Ret_kOK);
		//#005 - orderTG (order is larger than the one used to create context)
		add_test_case("orderTG", "order=4", TCat_kNG, "order", PVType_kTG,
				Ret_kIA);
		//#006 - num_dataTL
		//add_test_case("num_data", "=14", PVType_kTL, TCat_kNG, Ret_kIA);
		add_test_case("num_dataTL", "num_data=14", TCat_kNG, "num_data",
				PVType_kTL, Ret_kIA);
		//#007 - num_dataTG
		add_test_case("num_dataTG", "num_data=16", TCat_kNG, "num_data",
				PVType_kTG, Ret_kIA);
		//#008 - num_dataLB
		add_test_case("num_dataLB", "num_data=2", TCat_kOK, "num_data",
				PVType_kLB, Ret_kOK);
		//#009 - dataNULL
		add_test_case("dataNULL", "data null pointer", TCat_kNG, "data",
				PVType_kNULL, Ret_kIA);
		//#010 - dataNA
		add_test_case("dataNA", "data not aligned", TCat_kNG, "data",
				PVType_kNAL, Ret_kIA);
		//#011 - maskNULL
		add_test_case("maskNULL", "mask null pointer", TCat_kNG, "mask",
				PVType_kNULL, Ret_kIA);
		//#012 - maskNA
		add_test_case("maskNA", "mask not aligned", TCat_kNG, "mask",
				PVType_kNAL, Ret_kIA);
		//#013 - minimum required effective data
		add_test_case("effdataLB", "minimum required effective data", TCat_kOK,
				"mask", PVType_kVP, Ret_kOK);
		//#014 - less than minimum required effective data
		add_test_case("effdataTL", "less than minimum required effective data",
				TCat_kNG, "mask", PVType_kVP, Ret_kNG);
		//#015 - no effective data
		add_test_case("noeffdata", "no effective data", TCat_kNG, "mask",
				PVType_kVP, Ret_kNG);
		//#016 - clip_threshold_sigma is negative
		add_test_case("clip_threshold_sigmaOR", "clip_threshold_sigma<0",
				TCat_kNG, "clip_threshold_sigma", PVType_kOR, Ret_kIA);
		//#017 - clip_threshold_sigma is zero
		add_test_case("clip_threshold_sigmaOR", "clip_threshold_sigma=0",
				TCat_kNG, "clip_threshold_sigma", PVType_kOR, Ret_kIA);
		//#018 - clip_threshold_sigmaUB
		add_test_case("clip_threshold_sigmaUB", "clip_threshold_sigma=FLT_MAX",
				TCat_kOK, "clip_threshold_sigma", PVType_kUB, Ret_kOK);
		//#019 - clip_threshold_sigma is very small positive
		add_test_case("clip_threshold_sigmaTL",
				"clip_threshold_sigma is positive but too small", TCat_kNG,
				"clip_threshold_sigma", PVType_kTL, Ret_kNG);
		//#020 - num_fitting_maxLB
		add_test_case("num_fitting_maxLB", "", TCat_kOK, "num_fitting_max",
				PVType_kLB, Ret_kOK);
		//#021 - num_fitting_maxUB
		add_test_case("num_fitting_maxUB", "", TCat_kOK, "num_fitting_max",
				PVType_kUB, Ret_kOK);
		//#022 - num_coeffTL
		add_test_case("num_coeffTL", "", TCat_kNG, "num_coeff", PVType_kTL,
				Ret_kIA);
		//#023 - num_coeffLB
		add_test_case("num_coeffLB", "", TCat_kOK, "num_coeff", PVType_kLB,
				Ret_kOK);
		//#024 - num_coeffUB
		add_test_case("num_coeffUB", "", TCat_kOK, "num_coeff", PVType_kUB,
				Ret_kOK);
		//#025 - num_coeffTG
		add_test_case("num_coeffTG", "", TCat_kNG, "num_coeff", PVType_kTG,
				Ret_kIA);
		//#026 - coeffNULL
		add_test_case("coeffNULL", "", TCat_kOK, "coeff", PVType_kNULL,
				Ret_kOK);
		//#027 - best_fitNULL
		add_test_case("best_fitNULL", "", TCat_kOK, "best_fit", PVType_kNULL,
				Ret_kOK);
		//#028 - residualNULL
		add_test_case("residualNULL", "", TCat_kOK, "residual", PVType_kNULL,
				Ret_kOK);
		//#029 - coeffNULLbest_fitNULL
		add_test_case("coeffNULLbest_fitNULL", "", TCat_kOK, "coeff,best_fit",
				PVType_kNULL, Ret_kOK);
		//#030 - coeffNULLresidualNULL
		add_test_case("coeffNULLresidualNULL", "", TCat_kOK, "coeff,residual",
				PVType_kNULL, Ret_kOK);
		//#031 - best_fitNULLresidualNULL
		add_test_case("best_fitNULLresidualNULL", "", TCat_kOK,
				"best_fit,residual", PVType_kNULL, Ret_kOK);
		//#032 - coeffNULLbest_fitNULLresidualNULL
		add_test_case("coeffNULLbest_fitNULLresidualNULL", "", TCat_kOK,
				"coeff,best_fit,residual", PVType_kNULL, Ret_kOK);
		//#033 - coeffNA
		add_test_case("coeffNA", "", TCat_kNG, "coeff", PVType_kNAL, Ret_kIA);
		//#034 - best_fitNA
		add_test_case("best_fitNA", "", TCat_kNG, "best_fit", PVType_kNAL,
				Ret_kIA);
		//#035 - best_fit=data
		add_test_case("best_fitVP", "best_fit=data", TCat_kOK, "best_fit",
				PVType_kVP, Ret_kOK);
		//#036 - residualNA
		add_test_case("residualNA", "", TCat_kNG, "residual", PVType_kNAL,
				Ret_kIA);
		//#037 - residual=data
		add_test_case("residualVP", "residual=data", TCat_kOK, "residual",
				PVType_kVP, Ret_kOK);
		//#038 - final_maskNULL
		add_test_case("final_maskNULL", "", TCat_kNG, "final_mask",
				PVType_kNULL, Ret_kIA);
		//#039 - final_maskNA
		add_test_case("final_maskNA", "", TCat_kNG, "final_mask", PVType_kNAL,
				Ret_kIA);
		//#040 - final_mask=mask
		add_test_case("final_maskVP", "final_mask=mask", TCat_kOK, "final_mask",
				PVType_kVP, Ret_kOK);
		//#041 - rmsNULL
		add_test_case("rmsNULL", "", TCat_kNG, "rms", PVType_kNULL, Ret_kIA);
		//#042 - lsqfit_statusNULL
		add_test_case("lsqfit_statusNULL", "", TCat_kNG, "lsqfit_status",
				PVType_kNULL, Ret_kIA);
		//#043 - performance_order
		add_test_case("performance_order", "order=99", TCat_kPF, "order",
				PVType_kIR, Ret_kOK);
		//#044 - performance_num_data
		add_test_case("performance_num_data", "num_data=10000", TCat_kPF,
				"order", PVType_kIR, Ret_kOK);
		//#045 - performance_num_fitting_max
		add_test_case("performance_num_fitting_max", "num_fitting_max=100",
				TCat_kPF, "order", PVType_kIR, Ret_kOK);
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
void AllocateNotAligned(size_t const length, T **data, void **data_storage) {
	*data = nullptr;
	unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(T) * (length + 1), data));
	++(*data);
	assert(!LIBSAKURA_SYMBOL(IsAligned)(*data));
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

	template<typename T> void DoSetNullPointer(string const &name, T &ptr_map) {
		if (ptr_map.count(name) == 1) {
			ptr_map[name] = nullptr;
			assert(ptr_map[name] == nullptr);
		}
	}
	void SetNullPointer(string const &name) {
		Free(name);
		DoSetNullPointer(name, sto);
		DoSetNullPointer(name, fptr);
		DoSetNullPointer(name, dptr);
		DoSetNullPointer(name, bptr);
		DoSetNullPointer(name, ctxt);
		DoSetNullPointer(name, fsta);
	}
	void AllocateAlignedFloat(string const &name, string const &length_name) {
		Free(name);
		AllocateAligned(size[length_name], &fptr[name], &sto[name]);
	}
	void AllocateNotAlignedFloat(string const &name,
			string const &length_name) {
		Free(name);
		AllocateNotAligned(size[length_name], &fptr[name], &sto[name]);
	}
	void AllocateAlignedDouble(string const &name, string const &length_name) {
		Free(name);
		AllocateAligned(size[length_name], &dptr[name], &sto[name]);
	}
	void AllocateNotAlignedDouble(string const &name,
			string const &length_name) {
		Free(name);
		AllocateNotAligned(size[length_name], &dptr[name], &sto[name]);
	}
	void AllocateAlignedBool(string const &name, string const &length_name) {
		Free(name);
		AllocateAligned(size[length_name], &bptr[name], &sto[name]);
	}
	void AllocateNotAlignedBool(string const &name, string const &length_name) {
		Free(name);
		AllocateNotAligned(size[length_name], &bptr[name], &sto[name]);
	}
	void AllocateFloat(string const &name) {
		Free(name);
		Allocate(1, &fptr[name]);
	}
	void AllocateStatus(string const &name) {
		Free(name);
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
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["poly"] = context_poly;

		//polynomial(num_data=2,order=1)
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_polynd2o1 = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextPolynomialFloat)(
				LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), 1, 2,
				&context_polynd2o1);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["polynd2o1"] = context_polynd2o1;

		//polynomial(num_data=100,order=99)
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_polynd100o99 = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextPolynomialFloat)(
				LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), 99, 100,
				&context_polynd100o99);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["polynd100o99"] = context_polynd100o99;

		//polynomial(num_data=10000)
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_polynd10000 = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextPolynomialFloat)(
				LIBSAKURA_SYMBOL(LSQFitType_kPolynomial), order, 10000,
				&context_polynd10000);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["polynd10000"] = context_polynd10000;

		//chebyshev
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_chebyshev = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextPolynomialFloat)(
				LIBSAKURA_SYMBOL(LSQFitType_kChebyshev), order, num_data,
				&context_chebyshev);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["chebyshev"] = context_chebyshev;

		//cspline
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_cspline = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextCubicSplineFloat)(order, num_data,
				&context_cspline);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["cspline"] = context_cspline;

		//sinusoid
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_sinusoid = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextSinusoidFloat)(order, num_data,
				&context_sinusoid);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["sinusoid"] = context_sinusoid;

		ctxt[name] = ctxt[defaultvalue_name];
	}
	void Free(string const &name) {
		if (sto.count(name) == 1) {
			if (sto[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(sto[name]);
				sto[name] = nullptr;
			}
		} else if (fptr.count(name) == 1) {
			if (fptr[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(fptr[name]);
				fptr[name] = nullptr;
			}
		} else if (dptr.count(name) == 1) {
			if (dptr[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(dptr[name]);
				dptr[name] = nullptr;
			}
		} else if (bptr.count(name) == 1) {
			if (bptr[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(bptr[name]);
				bptr[name] = nullptr;
			}
		} else if (fsta.count(name) == 1) {
			if (fsta[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(fsta[name]);
				fsta[name] = nullptr;
			}
		} else if (ctxt.count(name) == 1) {
			if (ctxt[name] != nullptr) {
				if ((name != "context") || (ctxt.size() == 1)) {
					LIBSAKURA_SYMBOL(Status) status;
					status = LIBSAKURA_SYMBOL(DestroyLSQFitContextFloat)(
							ctxt[name]);
					EXPECT_EQ(Ret_kOK, status);
				}
				ctxt[name] = nullptr;
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

bool HasString(string const &str, string const &key) {
	return (str.find(key) != string::npos);
}

//interpret test cases, setup all parameter values to be used for testing
void Prologue(TestCase &tc, TestCase const &default_tc, ParamSet &ps) {
	auto title_has = [&](string const &s) {return HasString(tc.title, s);};
	auto desc_has = [&](string const &s) {return HasString(tc.desc, s);};
	cout << "    {" << tc.title << ": " << tc.desc << "}   " << flush << endl;
	/*
	 cout << endl;
	 for (size_t j = 0; j < tc.param.size(); ++j) {
	 cout << "(" << tc.param[j].name << ")" << "[" << tc.param[j].value_type << "]   ";
	 }
	 cout << endl;
	 */

	ps = CreateDefaultParameterSet(tc.api_name);
	PAttr param;
	PAttr default_param;
	string name;
	PDType dtype;
	PVType vtype;
	if (tc.category == TCat_kPF) {
		if (tc.api_name == ApiName_kLSQFitPolynomialFloat) {
			if (desc_has("order=99")) {
				tc.num_repeat = 40;
				ps.ui16["order"] = 99;
				ps.size["num_coeff"] = ps.ui16["order"] + 1;
				ps.ctxt["context"] = ps.ctxt["polynd100o99"];
				ps.size["num_data"] = ps.size["num_coeff"];
				ps.AllocateAlignedFloat("data", "num_data");
				ps.AllocateAlignedBool("mask", "num_data");
				for (size_t i = 0; i < ps.size["num_data"]; ++i) {
					ps.fptr["data"][i] = 10.0f;
					ps.bptr["mask"][i] = true;
				}
				ps.ui16["num_fitting_max"] = 1;
				ps.AllocateAlignedDouble("coeff", "num_coeff");
				ps.AllocateAlignedFloat("best_fit", "num_data");
				ps.AllocateAlignedFloat("residual", "num_data");
				ps.AllocateAlignedBool("final_mask", "num_data");
			} else if (desc_has("num_data=10000")) {
				tc.num_repeat = 700;
				ps.ctxt["context"] = ps.ctxt["polynd10000"];
				ps.size["num_data"] = 10000;
				ps.AllocateAlignedFloat("data", "num_data");
				ps.AllocateAlignedBool("mask", "num_data");
				for (size_t i = 0; i < ps.size["num_data"]; ++i) {
					ps.fptr["data"][i] = 10.0f;
					ps.bptr["mask"][i] = true;
				}
				ps.AllocateAlignedDouble("coeff", "num_coeff");
				ps.AllocateAlignedFloat("best_fit", "num_data");
				ps.AllocateAlignedFloat("residual", "num_data");
				ps.AllocateAlignedBool("final_mask", "num_data");
			} else if (desc_has("num_fitting_max=100")) {
				tc.num_repeat = 25000;
				ps.fval["clip_threshold_sigma"] = FLT_MAX;
				ps.ui16["num_fitting_max"] = 100;
				ps.fptr["data"][0] = 0.0;
			}
		}
		return;
	}
	for (size_t i = 0; i < tc.param.size(); ++i) {
		param = tc.param[i];
		name = param.name;
		dtype = param.data_type;
		default_param = default_tc.param[i];
		assert(name == default_param.name);
		assert(dtype == default_param.data_type);
		vtype = param.value_type;
		if (vtype == PVType_kNULL) {
			ps.SetNullPointer(name);
		} else if (vtype == PVType_kNAL) {
			if ((name == "data") || (name == "best_fit")
					|| (name == "residual")) {
				ps.AllocateNotAlignedFloat(name, "num_data");
			} else if (name == "coeff") {
				ps.AllocateNotAlignedDouble(name, "num_coeff");
			} else if ((name == "mask") || (name == "final_mask")) {
				ps.AllocateNotAlignedBool(name, "num_data");
			}
		} else {
			if (name == "context") {
				assert(dtype == PDType_kContextPointer);
				for (size_t j = 0; j < LsqFuncTypeStr.size(); ++j) {
					if (desc_has(LsqFuncTypeStr[j])) {
						ps.ctxt[name] = ps.ctxt[LsqFuncTypeStr[j]];
						break;
					}
				}
			} else if (name == "order") {
				if (vtype == PVType_kLB) {
					ps.ui16[name] = 0;
					ps.size["num_coeff"] = ps.ui16[name] + 1;
					ps.AllocateAlignedDouble("coeff", "num_coeff");
				} else if (vtype == PVType_kTG) {
					ps.ui16[name] = 4;
					ps.size["num_coeff"] = ps.ui16[name] + 1;
					ps.AllocateAlignedDouble("coeff", "num_coeff");
				}
			} else if (name == "num_data") {
				if (vtype == PVType_kTL) {
					ps.size[name]--;
				} else if (vtype == PVType_kTG) {
					ps.size[name]++;
				} else if (vtype == PVType_kLB) {
					ps.ui16["order"] = 1;
					ps.size["num_coeff"] = ps.ui16["order"] + 1;
					ps.size[name] = ps.size["num_coeff"];
					ps.ctxt["context"] = ps.ctxt["polynd2o1"];
				}
			} else if (name == "mask") {
				if (title_has("effdataLB")) {
					for (size_t j = ps.ui16["order"] + 1;
							j < ps.size["num_data"]; ++j) {
						ps.bptr[name][j] = false;
					}
				} else if (title_has("effdataTL")) {
					for (size_t j = ps.ui16["order"]; j < ps.size["num_data"];
							++j) {
						ps.bptr[name][j] = false;
					}
				} else if (title_has("noeffdata")) {
					for (size_t j = 0; j < ps.size["num_data"]; ++j) {
						ps.bptr[name][j] = false;
					}
				}
			} else if (name == "clip_threshold_sigma") {
				if (desc_has(name + "<0")) {
					ps.fval[name] = -1.0f;
				} else if (desc_has(name + "=0")) {
					ps.fval[name] = 0.0f;
				} else if (desc_has("FLT_MAX")) {
					ps.fval[name] = FLT_MAX;
				} else if (desc_has("small")) {
					ps.fval[name] = 0.1f;
					ps.ui16["num_fitting_max"] = 5;
					ps.fptr["data"][0] = 1000.0;
				}
			} else if (name == "num_fitting_max") {
				if (vtype == PVType_kLB) {
					ps.ui16[name] = 0;
				} else if (vtype == PVType_kUB) {
					ps.ui16[name] = UINT16_MAX;
				}
			} else if (name == "num_coeff") {
				if (vtype == PVType_kTL) {
					ps.ui16["order"] = 2;
					ps.size[name] = 2;
					ps.AllocateAlignedDouble("coeff", name);
				} else if (vtype == PVType_kLB) {
					ps.ui16["order"] = 2;
					ps.size[name] = 3;
					ps.AllocateAlignedDouble("coeff", name);
				} else if (vtype == PVType_kUB) {
					ps.ui16["order"] = 3;
					ps.size[name] = 4;
					ps.AllocateAlignedDouble("coeff", name);
				} else if (vtype == PVType_kTG) {
					ps.ui16["order"] = 3;
					ps.size[name] = 5;
					ps.AllocateAlignedDouble("coeff", name);
				}
			} else if ((name == "best_fit") || (name == "residual")) {
				if (desc_has(name + "=data")) {
					ps.Free(name);
					ps.fptr[name] = ps.fptr["data"];
				}
			} else if (name == "final_mask") {
				if (desc_has(name + "=mask")) {
					ps.Free(name);
					ps.bptr[name] = ps.bptr["mask"];
				}
			}
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

template<typename T> bool HasKey(T const &m, string const &key) {
	return (m.count(key) > 0);
}

void CheckValues(TestCase const &tc, ParamSet &ps) {
	if (!HasString(tc.title, "coeffNULL") && HasKey(ps.dptr, "num_coeff")
			&& HasKey(ps.size, "num_coeff")) {
		//cout << "          ****coeff=[";
		for (size_t j = 0; j < ps.size["num_coeff"]; ++j) {
			//cout << ps.dptr["coeff"][j] << ", ";
		}
		//cout << "]****" << endl;
	}
	if (!HasString(tc.title, "best_fitNULL") && (HasKey(ps.fptr, "best_fit"))
			&& (HasKey(ps.size, "num_data"))) {
		//cout << "          ****best_fit=[";
		for (size_t j = 0; j < ps.size["num_data"]; ++j) {
			//cout << ps.fptr["best_fit"][j] << ", ";
		}
		//cout << "]****" << endl;
	}
	if (!HasString(tc.title, "residualNULL") && (HasKey(ps.fptr, "residual"))
			&& (HasKey(ps.size, "num_data"))) {
		//cout << "          ****residual=[";
		for (size_t j = 0; j < ps.size["num_data"]; ++j) {
			//cout << ps.fptr["residual"][j] << ", ";
		}
		//cout << "]****" << endl;
	}
}

void Execute(TestCase const &tc, ParamSet &ps) {
	double time_elapsed = 0.0;
	LIBSAKURA_SYMBOL (Status) run_status;
	for (size_t i = 0; i < tc.num_repeat; ++i) {
		/*
		 cout << "####[" << tc.num_repeat << "]" << "{" << ps.ui16["order"]
		 << "}" << "<" << ps.ctxt["context"] << ">" << "####" << endl;
		 */
		double time_start = GetCurrentTime();
		RunApi(tc, ps, run_status);
		double time_end = GetCurrentTime();
		time_elapsed += (time_end - time_start);
	}

	if (tc.category == TCat_kPF) {
		//cout << endl;
		cout << setprecision(5) << "#x# benchmark Lsq_" << tc.title << " "
				<< time_elapsed << endl;
	}

	CheckStatus(tc, run_status);
	if ((tc.category == TCat_kOK) || (tc.category == TCat_kPF)) {
		CheckValues(tc, ps);
	}
}

template<typename T> void FreeParamSet(ParamSet &ps, T &ptr) {
	for (auto it = ptr.begin(); it != ptr.end(); ++it) {
		ps.Free(it->first);
	}
}

void Epilogue(ParamSet &ps) {
	FreeParamSet(ps, ps.sto);
	FreeParamSet(ps, ps.fptr);
	FreeParamSet(ps, ps.dptr);
	FreeParamSet(ps, ps.bptr);
	FreeParamSet(ps, ps.ctxt);
	FreeParamSet(ps, ps.fsta);
	//cout << endl;
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
		EXPECT_EQ(Ret_kOK, status);
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

//Destroy(context, Ret_kOK);
/*
 void Destroy(LIBSAKURA_SYMBOL(LSQFitContextFloat) *context, size_t status) {
 LIBSAKURA_SYMBOL (Status) destroy_status = sakura_DestroyLSQFitContextFloat(
 context);
 EXPECT_EQ(status, destroy_status);
 cout << "Destroy Status : " << destroy_status << endl;
 }
 */

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
 EXPECT_EQ(Ret_kOK, create_status);
 Destroy(context, Ret_kOK);
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
