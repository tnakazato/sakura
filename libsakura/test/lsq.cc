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

//Aliases for returned status
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
	TCat_kOK, //Success
	TCat_kPF, //Performance test
	TCat_kNG, //Fail
	TCat_kNumElements
} TCat;

//Parameter Data Type
typedef enum {
	PDType_kInt,
	PDType_kSizeT,
	PDType_kUInt16T,
	PDType_kFloat,
	PDType_kDouble,
	PDType_kBool,
	PDType_kLSQFitType,
	PDType_kSizeTPointer,
	PDType_kFloatPointer,
	PDType_kDoublePointer,
	PDType_kBoolPointer,
	PDType_kContextPointer,
	PDType_kContextPointerPointer,
	PDType_kCsplineCoeffPointer,
	PDType_kFitStatusPointer,
	PDType_kNumElements
} PDType;

vector<string> LsqFuncTypeStr = { "poly", "chebyshev", "cspline", "sinusoid" };

bool IsPointer(PDType const &data_type) {
	bool res = false;
	if ((data_type == PDType_kSizeTPointer)
			|| (data_type == PDType_kFloatPointer)
			|| (data_type == PDType_kDoublePointer)
			|| (data_type == PDType_kBoolPointer)
			|| (data_type == PDType_kContextPointer)
			|| (data_type == PDType_kContextPointerPointer)
			|| (data_type == PDType_kCsplineCoeffPointer)
			|| (data_type == PDType_kFitStatusPointer))
		res = true;
	return res;
}

//Parameter Category
typedef enum {
	PCat_kValid, //Valid value
	PCat_kInvalid, //Invalid value
	PCat_kUndefined, //Parameter value is not defined
	PCat_kNumElements
} PCat;

//Parameter Value Type
typedef enum {
	PVType_kLB, //Lower bound
	PVType_kUB, //Upper bound
	PVType_kIR, //In range
	PVType_kAL, //Aligned
	PVType_kVP, //Valid pointer
	PVType_kTL, //Too less
	PVType_kTG, //Too great
	PVType_kOR, //Out of range
	PVType_kDUP, //Duplicate
	PVType_kIVO, //Invalid order
	PVType_kNAN, //NaN
	PVType_kINF, //Infinity
	PVType_kNINF, //Negative Infinity
	PVType_kNULL, //Null
	PVType_kNAL, //Not aligned
	PVType_kUNDEF, //Undefined
	PVType_kNumElements
} PVType;

//Parameter Attribute
struct PAttr {
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
		//cout << "@@@@ GetParamAttrIndex"
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
	test_case.test_id = 0;
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
	case ApiName_kCreateLSQFitContextPolynomialFloat:
		num_params = 4;
		param_name = {"lsqfit_type", "order", "num_data", "context"};
		param_dtype = {PDType_kLSQFitType, PDType_kUInt16T, PDType_kSizeT, PDType_kContextPointerPointer};
		break;
		case ApiName_kCreateLSQFitContextCubicSplineFloat:
		num_params = 3;
		param_name = {"npiece", "num_data", "context"};
		param_dtype = {PDType_kUInt16T, PDType_kSizeT, PDType_kContextPointerPointer};
		break;
		case ApiName_kCreateLSQFitContextSinusoidFloat:
		num_params = 3;
		param_name = {"nwave", "num_data", "context"};
		param_dtype = {PDType_kUInt16T, PDType_kSizeT, PDType_kContextPointerPointer};
		break;
		case ApiName_kDestroyLSQFitContextFloat:
		num_params = 1;
		param_name = {"context"};
		param_dtype = {PDType_kContextPointer};
		break;
		case ApiName_kGetNumberOfCoefficientsFloat:
		num_params = 3;
		param_name = {"context", "order", "num_coeff"};
		param_dtype = {PDType_kContextPointer, PDType_kUInt16T, PDType_kSizeTPointer};
		break;
		case ApiName_kLSQFitPolynomialFloat:
		num_params = 14;
		param_name = {"context", "order", "num_data",
			"data", "mask", "clip_threshold_sigma",
			"num_fitting_max", "num_coeff", "coeff",
			"best_fit", "residual", "final_mask",
			"rms", "lsqfit_status"};
		param_dtype = {PDType_kContextPointer, PDType_kUInt16T, PDType_kSizeT,
			PDType_kFloatPointer, PDType_kBoolPointer, PDType_kFloat,
			PDType_kUInt16T, PDType_kSizeT, PDType_kDoublePointer,
			PDType_kFloatPointer, PDType_kFloatPointer, PDType_kBoolPointer,
			PDType_kFloatPointer, PDType_kFitStatusPointer};
		break;
		case ApiName_kLSQFitCubicSplineFloat:
		num_params = 14;
		param_name = {"context", "num_pieces", "num_data",
			"data", "mask", "clip_threshold_sigma",
			"num_fitting_max", "coeff", "best_fit",
			"residual", "final_mask", "rms",
			"boundary", "lsqfit_status"};
		param_dtype = {PDType_kContextPointer, PDType_kSizeT, PDType_kSizeT,
			PDType_kFloatPointer, PDType_kBoolPointer, PDType_kFloat,
			PDType_kUInt16T, PDType_kCsplineCoeffPointer, PDType_kFloatPointer,
			PDType_kFloatPointer, PDType_kBoolPointer, PDType_kFloatPointer,
			PDType_kSizeTPointer, PDType_kFitStatusPointer};
		break;
		case ApiName_kLSQFitSinusoidFloat:
		num_params = 15;
		param_name = {"context", "num_nwave",
			"nwave", "num_data", "data",
			"mask", "clip_threshold_sigma", "num_fitting_max",
			"num_coeff", "coeff", "best_fit",
			"residual", "final_mask",
			"rms", "lsqfit_status"};
		param_dtype = {PDType_kContextPointer, PDType_kSizeT,
			PDType_kSizeTPointer, PDType_kSizeT, PDType_kFloatPointer,
			PDType_kBoolPointer, PDType_kFloat, PDType_kUInt16T,
			PDType_kSizeT, PDType_kDoublePointer, PDType_kFloatPointer,
			PDType_kFloatPointer, PDType_kBoolPointer,
			PDType_kFloatPointer, PDType_kFitStatusPointer};
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
		tc->title = ApiNameStr[api_name] + " : " + title_suffix;
		tc->desc = test_description;
		tc->category = test_category;
		vector<string> names = Split(param_name, ',');
		for (size_t i = 0; i < names.size(); ++i) {
			tc->GetParamAttrByName(names[i]).value_type = param_vtype;
		}
		tc->expect_status = expect_status;
	};
	switch (api_name) {
	case ApiName_kCreateLSQFitContextPolynomialFloat:
		add_test_case("lsqfit_typeVP", "chebyshev", TCat_kOK, "lsqfit_type",
				PVType_kIR, Ret_kOK);
		add_test_case("orderLB", "=0", TCat_kOK, "order", PVType_kLB, Ret_kOK);
		//add_test_case("orderUB", "=UINT16_MAX", TCat_kOK, "order", PVType_kUB, Ret_kNM);
		add_test_case("orderLBnum_dataLB", "=0,1", TCat_kOK, "order",
				PVType_kLB, Ret_kOK);
		add_test_case("orderUBnum_dataUB", "=UINT16_MAX,SIZE_MAX", TCat_kOK,
				"order", PVType_kUB, Ret_kNM);
		add_test_case("num_dataTL", "=3", TCat_kNG, "num_data", PVType_kTL,
				Ret_kIA);
		add_test_case("num_dataLB", "=4", TCat_kOK, "num_data", PVType_kLB,
				Ret_kOK);
		add_test_case("num_dataUB", "=SIZE_MAX", TCat_kOK, "num_data",
				PVType_kUB, Ret_kNM);
		add_test_case("contextNULL", "", TCat_kNG, "context", PVType_kNULL,
				Ret_kIA);
		break;
	case ApiName_kCreateLSQFitContextCubicSplineFloat:
		add_test_case("npieceLB", "=1", TCat_kOK, "npiece", PVType_kLB,
				Ret_kOK);
		//add_test_case("npieceUB", "=UINT16_MAX", TCat_kNG, "npiece", PVType_kUB, Ret_kNM);
		add_test_case("npieceLBnum_dataLB", "=1,4", TCat_kOK, "npiece",
				PVType_kLB, Ret_kOK);
		add_test_case("npieceUBnum_dataUB", "=UINT16_MAX,SIZE_MAX", TCat_kNG,
				"npiece", PVType_kUB, Ret_kNM);
		add_test_case("num_dataTL", "=3", TCat_kNG, "num_data", PVType_kTL,
				Ret_kIA);
		add_test_case("num_dataLB", "=4", TCat_kOK, "num_data", PVType_kLB,
				Ret_kOK);
		add_test_case("num_dataUB", "=SIZE_MAX", TCat_kNG, "num_data",
				PVType_kUB, Ret_kNM);
		add_test_case("contextNULL", "", TCat_kNG, "context", PVType_kNULL,
				Ret_kIA);
		break;
	case ApiName_kCreateLSQFitContextSinusoidFloat:
		add_test_case("nwaveLB", "=0", TCat_kOK, "nwave", PVType_kLB, Ret_kOK);
		//add_test_case("nwaveUB", "=UINT16_MAX", TCat_kNG, "nwave", PVType_kUB, Ret_kNM);
		add_test_case("nwaveLBnum_dataLB", "=0,1", TCat_kOK, "nwave",
				PVType_kLB, Ret_kOK);
		add_test_case("nwaveUBnum_dataUB", "=UINT16_MAX,SIZE_MAX", TCat_kNG,
				"nwave", PVType_kUB, Ret_kNM);
		add_test_case("num_dataTL", "=6", TCat_kNG, "num_data", PVType_kTL,
				Ret_kIA);
		add_test_case("num_dataLB", "=7", TCat_kOK, "num_data", PVType_kLB,
				Ret_kOK);
		add_test_case("num_dataUB", "=SIZE_MAX", TCat_kNG, "num_data",
				PVType_kUB, Ret_kNM);
		add_test_case("contextNULL", "", TCat_kNG, "context", PVType_kNULL,
				Ret_kIA);
		break;
	case ApiName_kDestroyLSQFitContextFloat:
		//new_add_test_case(Ret_kIA, "", "context", PVType_kNULL);
		add_test_case("contextNULL", "", TCat_kNG, "context", PVType_kNULL,
				Ret_kIA);
		break;
	case ApiName_kGetNumberOfCoefficientsFloat:
		add_test_case("contextVP", "chebyshev", TCat_kOK, "context", PVType_kVP,
				Ret_kOK);
		add_test_case("contextVP", "cspline", TCat_kOK, "context", PVType_kVP,
				Ret_kOK);
		add_test_case("contextOR", "sinusoid", TCat_kNG, "context", PVType_kOR,
				Ret_kIA);
		add_test_case("contextNULL", "", TCat_kNG, "context", PVType_kNULL,
				Ret_kIA);
		add_test_case("orderTL", "=0 for cspline", TCat_kNG, "order",
				PVType_kTL, Ret_kIA);
		add_test_case("orderLB", "=1 for cspline", TCat_kOK, "order",
				PVType_kLB, Ret_kOK);
		add_test_case("orderLB", "=0", TCat_kOK, "order", PVType_kLB, Ret_kOK);
		add_test_case("orderUB", "=3", TCat_kOK, "order", PVType_kUB, Ret_kOK);
		add_test_case("orderTG", "=4", TCat_kNG, "order", PVType_kTG, Ret_kIA);
		add_test_case("num_coeffNULL", "", TCat_kNG, "num_coeff", PVType_kNULL,
				Ret_kIA);
		break;
	case ApiName_kLSQFitPolynomialFloat:
		add_test_case("contextVP", "chebyshev", TCat_kOK, "context", PVType_kVP,
				Ret_kOK);
		add_test_case("contextOR", "cspline", TCat_kNG, "context", PVType_kOR,
				Ret_kIA);
		add_test_case("contextOR", "sinusoid", TCat_kNG, "context", PVType_kOR,
				Ret_kIA);
		add_test_case("contextNULL", "", TCat_kNG, "context", PVType_kNULL,
				Ret_kIA);
		add_test_case("orderLB", "order=0", TCat_kOK, "order", PVType_kLB,
				Ret_kOK);
		add_test_case("orderTG", "order=4", TCat_kNG, "order", PVType_kTG,
				Ret_kIA);
		//new_add_test_case(Ret_kIA, "14", "num_data", PVType_kTL);
		add_test_case("num_dataTL", "num_data=14", TCat_kNG, "num_data",
				PVType_kTL, Ret_kIA);
		add_test_case("num_dataTG", "num_data=16", TCat_kNG, "num_data",
				PVType_kTG, Ret_kIA);
		add_test_case("num_dataLB", "num_data=2", TCat_kOK, "num_data",
				PVType_kLB, Ret_kOK);
		add_test_case("dataNULL", "", TCat_kNG, "data", PVType_kNULL, Ret_kIA);
		add_test_case("dataNA", "", TCat_kNG, "data", PVType_kNAL, Ret_kIA);
		add_test_case("maskNULL", "", TCat_kNG, "mask", PVType_kNULL, Ret_kIA);
		add_test_case("maskNA", "", TCat_kNG, "mask", PVType_kNAL, Ret_kIA);
		add_test_case("effdataLB", "minimum required effective data", TCat_kOK,
				"mask", PVType_kVP, Ret_kOK);
		add_test_case("effdataTL", "less than minimum required effective data",
				TCat_kNG, "mask", PVType_kVP, Ret_kNG);
		add_test_case("noeffdata", "no effective data", TCat_kNG, "mask",
				PVType_kVP, Ret_kNG);
		add_test_case("clip_threshold_sigmaOR", "clip_threshold_sigma<0",
				TCat_kNG, "clip_threshold_sigma", PVType_kOR, Ret_kIA);
		add_test_case("clip_threshold_sigmaOR", "clip_threshold_sigma=0",
				TCat_kNG, "clip_threshold_sigma", PVType_kOR, Ret_kIA);
		add_test_case("clip_threshold_sigmaUB", "clip_threshold_sigma=FLT_MAX",
				TCat_kOK, "clip_threshold_sigma", PVType_kUB, Ret_kOK);
		add_test_case("clip_threshold_sigmaTL",
				"clip_threshold_sigma is positive but too small", TCat_kNG,
				"clip_threshold_sigma", PVType_kTL, Ret_kNG);
		add_test_case("num_fitting_maxLB", "", TCat_kOK, "num_fitting_max",
				PVType_kLB, Ret_kOK);
		add_test_case("num_fitting_maxUB", "", TCat_kOK, "num_fitting_max",
				PVType_kUB, Ret_kOK);
		add_test_case("num_coeffTL", "", TCat_kNG, "num_coeff", PVType_kTL,
				Ret_kIA);
		add_test_case("num_coeffLB", "", TCat_kOK, "num_coeff", PVType_kLB,
				Ret_kOK);
		add_test_case("num_coeffUB", "", TCat_kOK, "num_coeff", PVType_kUB,
				Ret_kOK);
		add_test_case("num_coeffTG", "", TCat_kNG, "num_coeff", PVType_kTG,
				Ret_kIA);
		add_test_case("coeffNULL", "", TCat_kOK, "coeff", PVType_kNULL,
				Ret_kOK);
		add_test_case("best_fitNULL", "", TCat_kOK, "best_fit", PVType_kNULL,
				Ret_kOK);
		add_test_case("residualNULL", "", TCat_kOK, "residual", PVType_kNULL,
				Ret_kOK);
		add_test_case("coeffNULLbest_fitNULL", "", TCat_kOK, "coeff,best_fit",
				PVType_kNULL, Ret_kOK);
		add_test_case("coeffNULLresidualNULL", "", TCat_kOK, "coeff,residual",
				PVType_kNULL, Ret_kOK);
		add_test_case("best_fitNULLresidualNULL", "", TCat_kOK,
				"best_fit,residual", PVType_kNULL, Ret_kOK);
		add_test_case("coeffNULLbest_fitNULLresidualNULL", "", TCat_kOK,
				"coeff,best_fit,residual", PVType_kNULL, Ret_kOK);
		add_test_case("coeffNA", "", TCat_kNG, "coeff", PVType_kNAL, Ret_kIA);
		add_test_case("best_fitNA", "", TCat_kNG, "best_fit", PVType_kNAL,
				Ret_kIA);
		add_test_case("best_fitVP", "best_fit=data", TCat_kOK, "best_fit",
				PVType_kVP, Ret_kOK);
		add_test_case("residualNA", "", TCat_kNG, "residual", PVType_kNAL,
				Ret_kIA);
		add_test_case("residualVP", "residual=data", TCat_kOK, "residual",
				PVType_kVP, Ret_kOK);
		add_test_case("final_maskNULL", "", TCat_kNG, "final_mask",
				PVType_kNULL, Ret_kIA);
		add_test_case("final_maskNA", "", TCat_kNG, "final_mask", PVType_kNAL,
				Ret_kIA);
		add_test_case("final_maskVP", "final_mask=mask", TCat_kOK, "final_mask",
				PVType_kVP, Ret_kOK);
		add_test_case("rmsNULL", "", TCat_kNG, "rms", PVType_kNULL, Ret_kIA);
		add_test_case("lsqfit_statusNULL", "", TCat_kNG, "lsqfit_status",
				PVType_kNULL, Ret_kIA);
		add_test_case("performance_order", "order=99", TCat_kPF, "order",
				PVType_kIR, Ret_kOK);
		add_test_case("performance_num_data", "num_data=10000", TCat_kPF,
				"order", PVType_kIR, Ret_kOK);
		add_test_case("performance_num_fitting_max", "num_fitting_max=100",
				TCat_kPF, "order", PVType_kIR, Ret_kOK);
		break;
	case ApiName_kLSQFitCubicSplineFloat:
		add_test_case("contextOR", "poly", TCat_kNG, "context", PVType_kOR,
				Ret_kIA);
		add_test_case("contextOR", "chebyshev", TCat_kNG, "context", PVType_kOR,
				Ret_kIA);
		add_test_case("contextOR", "sinusoid", TCat_kNG, "context", PVType_kOR,
				Ret_kIA);
		add_test_case("contextNULL", "", TCat_kNG, "context", PVType_kNULL,
				Ret_kIA);
		add_test_case("num_piecesLB", "=1", TCat_kOK, "num_pieces", PVType_kLB,
				Ret_kOK);
		add_test_case("num_piecesLBnum_dataLB", "=1,4", TCat_kOK, "num_pieces",
				PVType_kLB, Ret_kOK);
		add_test_case("num_piecesTG", "=13", TCat_kNG, "num_pieces", PVType_kTG,
				Ret_kIA);
		add_test_case("num_dataTL", "=14", TCat_kNG, "num_data", PVType_kTL,
				Ret_kIA);
		add_test_case("num_dataTG", "=16", TCat_kNG, "num_data", PVType_kTG,
				Ret_kIA);
		add_test_case("dataNULL", "", TCat_kNG, "data", PVType_kNULL, Ret_kIA);
		add_test_case("dataNA", "", TCat_kNG, "data", PVType_kNAL, Ret_kIA);
		add_test_case("maskNULL", "", TCat_kNG, "mask", PVType_kNULL, Ret_kIA);
		add_test_case("maskNA", "", TCat_kNG, "mask", PVType_kNAL, Ret_kIA);
		add_test_case("effdataLB", "minimum required effective data", TCat_kOK,
				"mask", PVType_kVP, Ret_kOK);
		add_test_case("effdataTL", "less than minimum required effective data",
				TCat_kNG, "mask", PVType_kVP, Ret_kNG);
		add_test_case("noeffdata", "no effective data", TCat_kNG, "mask",
				PVType_kVP, Ret_kNG);
		add_test_case("clip_threshold_sigmaOR", "clip_threshold_sigma<0",
				TCat_kNG, "clip_threshold_sigma", PVType_kOR, Ret_kIA);
		add_test_case("clip_threshold_sigmaOR", "clip_threshold_sigma=0",
				TCat_kNG, "clip_threshold_sigma", PVType_kOR, Ret_kIA);
		add_test_case("clip_threshold_sigmaUB", "clip_threshold_sigma=FLT_MAX",
				TCat_kOK, "clip_threshold_sigma", PVType_kUB, Ret_kOK);
		add_test_case("clip_threshold_sigmaTL",
				"clip_threshold_sigma is positive but too small", TCat_kNG,
				"clip_threshold_sigma", PVType_kTL, Ret_kNG);
		add_test_case("num_fitting_maxLB", "", TCat_kOK, "num_fitting_max",
				PVType_kLB, Ret_kOK);
		add_test_case("num_fitting_maxUB", "", TCat_kOK, "num_fitting_max",
				PVType_kUB, Ret_kOK);
		add_test_case("coeffNULL", "", TCat_kOK, "coeff", PVType_kNULL,
				Ret_kOK);
		add_test_case("best_fitNULL", "", TCat_kOK, "best_fit", PVType_kNULL,
				Ret_kOK);
		add_test_case("residualNULL", "", TCat_kOK, "residual", PVType_kNULL,
				Ret_kOK);
		add_test_case("coeffNULLbest_fitNULL", "", TCat_kOK, "coeff,best_fit",
				PVType_kNULL, Ret_kOK);
		add_test_case("coeffNULLresidualNULL", "", TCat_kOK, "coeff,residual",
				PVType_kNULL, Ret_kOK);
		add_test_case("best_fitNULLresidualNULL", "", TCat_kOK,
				"best_fit,residual", PVType_kNULL, Ret_kOK);
		add_test_case("coeffNULLbest_fitNULLresidualNULL", "", TCat_kOK,
				"coeff,best_fit,residual", PVType_kNULL, Ret_kOK);
		add_test_case("coeffNA", "", TCat_kNG, "coeff", PVType_kNAL, Ret_kIA);
		add_test_case("best_fitNA", "", TCat_kNG, "best_fit", PVType_kNAL,
				Ret_kIA);
		add_test_case("best_fitVP", "best_fit=data", TCat_kOK, "best_fit",
				PVType_kVP, Ret_kOK);
		add_test_case("residualNA", "", TCat_kNG, "residual", PVType_kNAL,
				Ret_kIA);
		add_test_case("residualVP", "residual=data", TCat_kOK, "residual",
				PVType_kVP, Ret_kOK);
		add_test_case("final_maskNULL", "", TCat_kNG, "final_mask",
				PVType_kNULL, Ret_kIA);
		add_test_case("final_maskNA", "", TCat_kNG, "final_mask", PVType_kNAL,
				Ret_kIA);
		add_test_case("final_maskVP", "final_mask=mask", TCat_kOK, "final_mask",
				PVType_kVP, Ret_kOK);
		add_test_case("rmsNULL", "", TCat_kNG, "rms", PVType_kNULL, Ret_kIA);
		add_test_case("lsqfit_statusNULL", "", TCat_kNG, "lsqfit_status",
				PVType_kNULL, Ret_kIA);
		add_test_case("performance_num_pieces", "num_pieces=99", TCat_kPF,
				"num_pieces", PVType_kIR, Ret_kOK);
		add_test_case("performance_num_data", "num_data=10000", TCat_kPF,
				"num_pieces", PVType_kIR, Ret_kOK);
		add_test_case("performance_num_fitting_max", "num_fitting_max=100",
				TCat_kPF, "num_pieces", PVType_kIR, Ret_kOK);
		break;
	case ApiName_kLSQFitSinusoidFloat:
		add_test_case("contextOR", "poly", TCat_kNG, "context", PVType_kOR,
				Ret_kIA);
		add_test_case("contextOR", "chebyshev", TCat_kNG, "context", PVType_kOR,
				Ret_kIA);
		add_test_case("contextOR", "cspline", TCat_kNG, "context", PVType_kOR,
				Ret_kIA);
		add_test_case("contextNULL", "", TCat_kNG, "context", PVType_kNULL,
				Ret_kIA);
		add_test_case("num_nwaveTL", "=0", TCat_kNG, "num_nwave", PVType_kTL,
				Ret_kIA);
		add_test_case("num_nwaveLB", "=1", TCat_kOK, "num_nwave", PVType_kLB,
				Ret_kOK);
		add_test_case("nwaveNULL", "", TCat_kNG, "nwave", PVType_kNULL,
				Ret_kIA);
		add_test_case("nwaveDUP", "0,0,1,2", TCat_kNG, "nwave", PVType_kDUP,
				Ret_kIA);
		add_test_case("nwaveDUP", "0,1,1,2", TCat_kNG, "nwave", PVType_kDUP,
				Ret_kIA);
		add_test_case("nwaveDUP", "0,1,2,2", TCat_kNG, "nwave", PVType_kDUP,
				Ret_kIA);
		add_test_case("nwaveIVO", "0,2,1", TCat_kNG, "nwave", PVType_kIVO,
				Ret_kIA);
		add_test_case("nwaveOR", "0,1,2,3", TCat_kNG, "nwave", PVType_kOR,
				Ret_kIA);
		add_test_case("nwavePartial", "0", TCat_kOK, "nwave", PVType_kIR,
				Ret_kOK);
		add_test_case("nwavePartial", "1", TCat_kOK, "nwave", PVType_kIR,
				Ret_kOK);
		add_test_case("nwavePartial", "2", TCat_kOK, "nwave", PVType_kIR,
				Ret_kOK);
		add_test_case("nwavePartial", "0,1", TCat_kOK, "nwave", PVType_kIR,
				Ret_kOK);
		add_test_case("nwavePartial", "0,2", TCat_kOK, "nwave", PVType_kIR,
				Ret_kOK);
		add_test_case("nwavePartial", "1,2", TCat_kOK, "nwave", PVType_kIR,
				Ret_kOK);
		add_test_case("num_dataTL", "num_data=14", TCat_kNG, "num_data",
				PVType_kTL, Ret_kIA);
		add_test_case("num_dataTG", "num_data=16", TCat_kNG, "num_data",
				PVType_kTG, Ret_kIA);
		add_test_case("num_dataLB", "num_data=2", TCat_kOK, "num_data",
				PVType_kLB, Ret_kOK);
		add_test_case("dataNULL", "", TCat_kNG, "data", PVType_kNULL, Ret_kIA);
		add_test_case("dataNA", "", TCat_kNG, "data", PVType_kNAL, Ret_kIA);
		add_test_case("maskNULL", "", TCat_kNG, "mask", PVType_kNULL, Ret_kIA);
		add_test_case("maskNA", "", TCat_kNG, "mask", PVType_kNAL, Ret_kIA);
		add_test_case("effdataLB",
				"minimum required effective data ; 0_in_nwave", TCat_kOK,
				"mask", PVType_kVP, Ret_kOK);
		add_test_case("effdataLB",
				"minimum required effective data ; 0_not_in_nwave", TCat_kOK,
				"mask", PVType_kVP, Ret_kOK);
		add_test_case("effdataTL",
				"less than minimum required effective data ; 0_in_nwave",
				TCat_kNG, "mask", PVType_kVP, Ret_kNG);
		add_test_case("effdataTL",
				"less than minimum required effective data ; 0_not_in_nwave",
				TCat_kNG, "mask", PVType_kVP, Ret_kNG);
		add_test_case("noeffdata", "no effective data", TCat_kNG, "mask",
				PVType_kVP, Ret_kNG);
		add_test_case("clip_threshold_sigmaOR", "clip_threshold_sigma<0",
				TCat_kNG, "clip_threshold_sigma", PVType_kOR, Ret_kIA);
		add_test_case("clip_threshold_sigmaOR", "clip_threshold_sigma=0",
				TCat_kNG, "clip_threshold_sigma", PVType_kOR, Ret_kIA);
		add_test_case("clip_threshold_sigmaUB", "clip_threshold_sigma=FLT_MAX",
				TCat_kOK, "clip_threshold_sigma", PVType_kUB, Ret_kOK);
		add_test_case("clip_threshold_sigmaTL",
				"clip_threshold_sigma is positive but too small", TCat_kNG,
				"clip_threshold_sigma", PVType_kTL, Ret_kNG);
		add_test_case("num_fitting_maxLB", "", TCat_kOK, "num_fitting_max",
				PVType_kLB, Ret_kOK);
		add_test_case("num_fitting_maxUB", "", TCat_kOK, "num_fitting_max",
				PVType_kUB, Ret_kOK);
		add_test_case("num_coeffTL", "", TCat_kNG, "num_coeff", PVType_kTL,
				Ret_kIA);
		add_test_case("num_coeffLB", "", TCat_kOK, "num_coeff", PVType_kLB,
				Ret_kOK);
		add_test_case("num_coeffUB", "", TCat_kOK, "num_coeff", PVType_kUB,
				Ret_kOK);
		add_test_case("num_coeffTG", "", TCat_kNG, "num_coeff", PVType_kTG,
				Ret_kIA);
		add_test_case("coeffNULL", "", TCat_kOK, "coeff", PVType_kNULL,
				Ret_kOK);
		add_test_case("best_fitNULL", "", TCat_kOK, "best_fit", PVType_kNULL,
				Ret_kOK);
		add_test_case("residualNULL", "", TCat_kOK, "residual", PVType_kNULL,
				Ret_kOK);
		add_test_case("coeffNULLbest_fitNULL", "", TCat_kOK, "coeff,best_fit",
				PVType_kNULL, Ret_kOK);
		add_test_case("coeffNULLresidualNULL", "", TCat_kOK, "coeff,residual",
				PVType_kNULL, Ret_kOK);
		add_test_case("best_fitNULLresidualNULL", "", TCat_kOK,
				"best_fit,residual", PVType_kNULL, Ret_kOK);
		add_test_case("coeffNULLbest_fitNULLresidualNULL", "", TCat_kOK,
				"coeff,best_fit,residual", PVType_kNULL, Ret_kOK);
		add_test_case("coeffNA", "", TCat_kNG, "coeff", PVType_kNAL, Ret_kIA);
		add_test_case("best_fitNA", "", TCat_kNG, "best_fit", PVType_kNAL,
				Ret_kIA);
		add_test_case("best_fitVP", "best_fit=data", TCat_kOK, "best_fit",
				PVType_kVP, Ret_kOK);
		add_test_case("residualNA", "", TCat_kNG, "residual", PVType_kNAL,
				Ret_kIA);
		add_test_case("residualVP", "residual=data", TCat_kOK, "residual",
				PVType_kVP, Ret_kOK);
		add_test_case("final_maskNULL", "", TCat_kNG, "final_mask",
				PVType_kNULL, Ret_kIA);
		add_test_case("final_maskNA", "", TCat_kNG, "final_mask", PVType_kNAL,
				Ret_kIA);
		add_test_case("final_maskVP", "final_mask=mask", TCat_kOK, "final_mask",
				PVType_kVP, Ret_kOK);
		add_test_case("rmsNULL", "", TCat_kNG, "rms", PVType_kNULL, Ret_kIA);
		add_test_case("lsqfit_statusNULL", "", TCat_kNG, "lsqfit_status",
				PVType_kNULL, Ret_kIA);
		add_test_case("performance_num_nwave", "num_nwave=51", TCat_kPF,
				"num_nwave", PVType_kIR, Ret_kOK);
		add_test_case("performance_num_data", "num_data=10000", TCat_kPF,
				"num_nwave", PVType_kIR, Ret_kOK);
		add_test_case("performance_num_fitting_max", "num_fitting_max=100",
				TCat_kPF, "num_nwave", PVType_kIR, Ret_kOK);
		break;
	default:
		assert(false);
		break;
	}
	return test_cases;
}

void AllocateAligned(size_t const length, double (**data)[4],
		void **data_storage) {
	*data = nullptr;
	unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(double) * length * 4, data));
	assert(LIBSAKURA_SYMBOL(IsAligned)(*data));
	*data_storage = storage_for_data.release();
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

void AllocateNotAligned(size_t const length, double (**data)[4],
		void **data_storage) {
	*data = nullptr;
	unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(double) * (length * 4 + 1), data));
	// ----------------------------------------------------
	// Note: in order to get unaligned address (*data),
	// ++(*data);
	// keeps it 64bit-aligned in case T is double.
	// Therefore, do as follows to shift just 1byte to
	// get unaligned address except for pointers to char.
	char *tmp = reinterpret_cast<char *>(*data);
	*data = reinterpret_cast<double (*)[4]>(++tmp);
	// ----------------------------------------------------
	assert(!LIBSAKURA_SYMBOL(IsAligned)(*data));
	*data_storage = storage_for_data.release();
}
template<typename T>
void AllocateNotAligned(size_t const length, T **data, void **data_storage) {
	*data = nullptr;
	unique_ptr<void, LIBSAKURA_PREFIX::Memory> storage_for_data(
			LIBSAKURA_PREFIX::Memory::AlignedAllocateOrException(
					sizeof(T) * (length + 1), data));
	// ----------------------------------------------------
	// Note: in order to get unaligned address (*data),
	// ++(*data);
	// keeps it 64bit-aligned in case T is double.
	// Therefore, do as follows to shift just 1byte to
	// get unaligned address except for pointers to char.
	char *tmp = reinterpret_cast<char *>(*data);
	*data = reinterpret_cast<T *>(++tmp);
	// ----------------------------------------------------
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
	map<string, LIBSAKURA_SYMBOL(LSQFitType)> ftyp;
	map<string, size_t *> sptr;
	map<string, float *> fptr;
	map<string, double *> dptr;
	map<string, bool *> bptr;
	map<string, LIBSAKURA_SYMBOL(LSQFitContextFloat) *> ctxt;
	map<string, LIBSAKURA_SYMBOL(LSQFitContextFloat) **> cptr;
	map<string, double (*)[4]> coef;
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
		DoSetNullPointer(name, sptr);
		DoSetNullPointer(name, fptr);
		DoSetNullPointer(name, dptr);
		DoSetNullPointer(name, bptr);
		DoSetNullPointer(name, ctxt);
		DoSetNullPointer(name, coef);
		DoSetNullPointer(name, fsta);
	}
	void AllocateAlignedSizeT(string const &name, size_t const length) {
		Free(name);
		AllocateAligned(length, &sptr[name], &sto[name]);
	}
	void AllocateAlignedSizeT(string const &name, string const &length_name) {
		Free(name);
		AllocateAligned(size[length_name], &sptr[name], &sto[name]);
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
	void AllocateSizeT(string const &name) {
		Free(name);
		Allocate(1, &sptr[name]);
	}
	void AllocateFloat(string const &name) {
		Free(name);
		Allocate(1, &fptr[name]);
	}
	void AllocateStatus(string const &name) {
		Free(name);
		Allocate(1, &fsta[name]);
	}
	void AllocateAlignedCsplineCoeff(string const &name,
			string const &length_name) {
		Free(name);
		AllocateAligned(size[length_name], &coef[name], &sto[name]);
	}
	void AllocateNotAlignedCsplineCoeff(string const &name,
			string const &length_name) {
		Free(name);
		AllocateNotAligned(size[length_name], &coef[name], &sto[name]);
	}
	void AllocateContexts(string const &name, string const &defaultvalue_name) {
		LIBSAKURA_SYMBOL(Status) create_status;
		uint16_t order = 3;
		if (ui16.count("order") == 1) {
			order = ui16["order"];
		}
		uint16_t npiece = 3;
		if (ui16.count("npiece") == 1) {
			npiece = ui16["npiece"];
		}
		uint16_t nwave = 2;
		if (ui16.count("nwave") == 1) {
			nwave = ui16["nwave"];
		}
		size_t num_data = 15;
		if (size.count("num_data") == 1) {
			num_data = size["num_data"];
		}

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
		LIBSAKURA_SYMBOL(CreateLSQFitContextCubicSplineFloat)(npiece, num_data,
				&context_cspline);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["cspline"] = context_cspline;

		//cspline(num_data=4,num_pieces=1)
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_csplinend4np1 = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextCubicSplineFloat)(1, 4,
				&context_csplinend4np1);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["csplinend4np1"] = context_csplinend4np1;

		//cspline(num_data=102,num_pieces=99)
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_csplinend102np99 = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextCubicSplineFloat)(99, 102,
				&context_csplinend102np99);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["csplinend102np99"] = context_csplinend102np99;

		//cspline(num_data=10000)
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_csplinend10000 = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextCubicSplineFloat)(npiece, 10000,
				&context_csplinend10000);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["csplinend10000"] = context_csplinend10000;

		//sinusoid
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_sinusoid = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextSinusoidFloat)(nwave, num_data,
				&context_sinusoid);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["sinusoid"] = context_sinusoid;

		//sinusoid(num_data=2,nwave=0)
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_sinusoidnd2nw0 = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextSinusoidFloat)(0, 2,
				&context_sinusoidnd2nw0);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["sinusoidnd2nw0"] = context_sinusoidnd2nw0;

		//sinusoid(nwave=3)
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_sinusoidnw3 = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextSinusoidFloat)(3, num_data,
				&context_sinusoidnw3);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["sinusoidnw3"] = context_sinusoidnw3;

		//sinusoid(num_data=102,nwave=50)
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_sinusoidnd102nw50 =
				nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextSinusoidFloat)(50, 102,
				&context_sinusoidnd102nw50);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["sinusoidnd102nw50"] = context_sinusoidnd102nw50;

		//sinusoid(num_data=10000)
		LIBSAKURA_SYMBOL(LSQFitContextFloat) *context_sinusoidnd10000 = nullptr;
		create_status =
		LIBSAKURA_SYMBOL(CreateLSQFitContextSinusoidFloat)(nwave, 10000,
				&context_sinusoidnd10000);
		EXPECT_EQ(Ret_kOK, create_status);
		ctxt["sinusoidnd10000"] = context_sinusoidnd10000;

		ctxt[name] = ctxt[defaultvalue_name];
		cptr[name] = &ctxt[name];
	}
	void Free(string const &name) {
		if (sto.count(name) == 1) {
			if (sto[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(sto[name]);
				sto[name] = nullptr;
			}
		} else if (sptr.count(name) == 1) {
			if (sptr[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(sptr[name]);
				sptr[name] = nullptr;
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
		} else if (coef.count(name) == 1) {
			if (coef[name] != nullptr) {
				LIBSAKURA_PREFIX::Memory::Free(coef[name]);
				coef[name] = nullptr;
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
	case ApiName_kCreateLSQFitContextPolynomialFloat:
		ps.ftyp["lsqfit_type"] = LIBSAKURA_SYMBOL(LSQFitType_kPolynomial);
		ps.ui16["order"] = 3;
		ps.size["num_data"] = 15;
		ps.AllocateContexts("context", "poly");
		break;
	case ApiName_kCreateLSQFitContextCubicSplineFloat:
		ps.ui16["npiece"] = 3;
		ps.size["num_data"] = 15;
		ps.AllocateContexts("context", "cspline");
		break;
	case ApiName_kCreateLSQFitContextSinusoidFloat:
		ps.ui16["nwave"] = 3;
		ps.size["num_data"] = 15;
		ps.AllocateContexts("context", "sinusoid");
		break;
	case ApiName_kDestroyLSQFitContextFloat:
		ps.AllocateContexts("context", "poly");
		break;
	case ApiName_kGetNumberOfCoefficientsFloat:
		ps.AllocateContexts("context", "poly");
		ps.ui16["order"] = 2;
		ps.AllocateSizeT("num_coeff");
		break;
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
	case ApiName_kLSQFitCubicSplineFloat:
		ps.size["num_pieces"] = 3;
		ps.size["num_data"] = 15;
		ps.AllocateContexts("context", "cspline");
		ps.AllocateAlignedFloat("data", "num_data");
		ps.AllocateAlignedBool("mask", "num_data");
		for (size_t i = 0; i < ps.size["num_data"]; ++i) {
			ps.fptr["data"][i] = 10.0;
			ps.bptr["mask"][i] = true;
		}
		ps.fval["clip_threshold_sigma"] = 3.0;
		ps.ui16["num_fitting_max"] = 2;
		ps.AllocateAlignedCsplineCoeff("coeff", "num_pieces");
		ps.AllocateAlignedFloat("best_fit", "num_data");
		ps.AllocateAlignedFloat("residual", "num_data");
		ps.AllocateAlignedBool("final_mask", "num_data");
		ps.AllocateFloat("rms");
		ps.AllocateAlignedSizeT("boundary", ps.size["num_pieces"] + 1);
		ps.AllocateStatus("lsqfit_status");
		break;
	case ApiName_kLSQFitSinusoidFloat:
		ps.size["num_nwave"] = 3;
		ps.AllocateAlignedSizeT("nwave", "num_nwave");
		for (size_t i = 0; i < ps.size["num_nwave"]; ++i) {
			ps.sptr["nwave"][i] = i;
		}
		ps.AllocateContexts("context", "sinusoid");
		ps.size["num_data"] = 15;
		ps.AllocateAlignedFloat("data", "num_data");
		ps.AllocateAlignedBool("mask", "num_data");
		for (size_t i = 0; i < ps.size["num_data"]; ++i) {
			ps.fptr["data"][i] = 10.0;
			ps.bptr["mask"][i] = true;
		}
		ps.fval["clip_threshold_sigma"] = 3.0;
		ps.ui16["num_fitting_max"] = 2;
		ps.size["num_coeff"] = ps.size["num_nwave"] * 2 - 1;
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

template<typename T> bool HasKey(T const &m, string const &key) {
	return (m.count(key) > 0);
}

//interpret test cases, setup all parameter values to be used for testing
void Prologue(TestCase &tc, TestCase const &default_tc, ParamSet &ps) {
	auto title_has = [&](string const &s) {return HasString(tc.title, s);};
	auto desc_has = [&](string const &s) {return HasString(tc.desc, s);};
	cout << "    " << tc.title << (tc.desc != "" ? (" (" + tc.desc + ") ") : "")
			<< flush << endl;
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
				ps.size["num_data"] = ps.size["num_coeff"];
				ps.ctxt["context"] = ps.ctxt["polynd100o99"];
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
				ps.size["num_data"] = 10000;
				ps.ctxt["context"] = ps.ctxt["polynd10000"];
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
		} else if (tc.api_name == ApiName_kLSQFitCubicSplineFloat) {
			if (desc_has("num_pieces=99")) {
				tc.num_repeat = 40;
				ps.size["num_pieces"] = 99;
				ps.size["num_data"] = ps.size["num_pieces"] + 3;
				ps.ctxt["context"] = ps.ctxt["csplinend102np99"];
				ps.AllocateAlignedFloat("data", "num_data");
				ps.AllocateAlignedBool("mask", "num_data");
				for (size_t i = 0; i < ps.size["num_data"]; ++i) {
					ps.fptr["data"][i] = 10.0f;
					ps.bptr["mask"][i] = true;
				}
				ps.ui16["num_fitting_max"] = 1;
				ps.AllocateAlignedCsplineCoeff("coeff", "num_pieces");
				ps.AllocateAlignedFloat("best_fit", "num_data");
				ps.AllocateAlignedFloat("residual", "num_data");
				ps.AllocateAlignedBool("final_mask", "num_data");
				ps.AllocateAlignedSizeT("boundary", ps.size["num_pieces"] + 1);
			} else if (desc_has("num_data=10000")) {
				tc.num_repeat = 300;
				ps.size["num_data"] = 10000;
				ps.ctxt["context"] = ps.ctxt["csplinend10000"];
				ps.AllocateAlignedFloat("data", "num_data");
				ps.AllocateAlignedBool("mask", "num_data");
				for (size_t i = 0; i < ps.size["num_data"]; ++i) {
					ps.fptr["data"][i] = 10.0f;
					ps.bptr["mask"][i] = true;
				}
				ps.ui16["num_fitting_max"] = 1;
				ps.AllocateAlignedCsplineCoeff("coeff", "num_pieces");
				ps.AllocateAlignedFloat("best_fit", "num_data");
				ps.AllocateAlignedFloat("residual", "num_data");
				ps.AllocateAlignedBool("final_mask", "num_data");
				ps.AllocateAlignedSizeT("boundary", ps.size["num_pieces"] + 1);
			} else if (desc_has("num_fitting_max=100")) {
				tc.num_repeat = 12000;
				ps.fval["clip_threshold_sigma"] = FLT_MAX;
				ps.ui16["num_fitting_max"] = 100;
				ps.fptr["data"][0] = 0.0;
			}
		} else if (tc.api_name == ApiName_kLSQFitSinusoidFloat) {
			if (desc_has("num_nwave=51")) {
				tc.num_repeat = 40;
				ps.size["num_nwave"] = 51;
				ps.AllocateAlignedSizeT("nwave", "num_nwave");
				for (size_t i = 0; i < ps.size["num_nwave"]; ++i) {
					ps.sptr["nwave"][i] = i;
				}
				ps.size["num_coeff"] = ps.size["num_nwave"] * 2 - 1;
				ps.size["num_data"] = ps.size["num_coeff"] + 1;
				ps.ctxt["context"] = ps.ctxt["sinusoidnd102nw50"];
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
				tc.num_repeat = 400;
				ps.size["num_data"] = 10000;
				ps.ctxt["context"] = ps.ctxt["sinusoidnd10000"];
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
			} else if (desc_has("num_fitting_max=100")) {
				tc.num_repeat = 15000;
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
			if ((tc.api_name == ApiName_kCreateLSQFitContextPolynomialFloat)
					|| (tc.api_name
							== ApiName_kCreateLSQFitContextCubicSplineFloat)
					|| (tc.api_name == ApiName_kCreateLSQFitContextSinusoidFloat)) {
				ps.cptr[name] = nullptr;
			} else {
				ps.SetNullPointer(name);
			}
		} else if (vtype == PVType_kNAL) {
			if ((name == "data") || (name == "best_fit")
					|| (name == "residual")) {
				ps.AllocateNotAlignedFloat(name, "num_data");
			} else if (name == "coeff") {
				if (tc.api_name == ApiName_kLSQFitCubicSplineFloat) {
					ps.AllocateNotAlignedCsplineCoeff(name, "num_pieces");
				} else {
					ps.AllocateNotAlignedDouble(name, "num_coeff");
				}
			} else if ((name == "mask") || (name == "final_mask")) {
				ps.AllocateNotAlignedBool(name, "num_data");
			}
		} else {
			if (name == "lsqfit_type") {
				ps.ftyp[name] = LIBSAKURA_SYMBOL(LSQFitType_kChebyshev);
			} else if (name == "context") {
				if ((tc.api_name == ApiName_kCreateLSQFitContextPolynomialFloat)
						|| (tc.api_name
								== ApiName_kCreateLSQFitContextCubicSplineFloat)
						|| (tc.api_name
								== ApiName_kCreateLSQFitContextSinusoidFloat)) {
					assert(dtype == PDType_kContextPointerPointer);
				} else {
					assert(dtype == PDType_kContextPointer);
					for (size_t j = 0; j < LsqFuncTypeStr.size(); ++j) {
						if (desc_has(LsqFuncTypeStr[j])) {
							ps.ctxt[name] = ps.ctxt[LsqFuncTypeStr[j]];
							break;
						}
					}
				}
			} else if (name == "order") {
				if (tc.api_name
						== ApiName_kCreateLSQFitContextPolynomialFloat) {
					if (vtype == PVType_kUB) {
						ps.ui16[name] = UINT16_MAX;
						ps.size["num_data"] = ps.ui16[name] + 1;
						if (title_has("num_dataUB")) {
							ps.size["num_data"] = SIZE_MAX;
						}
					} else if (vtype == PVType_kLB) {
						ps.ui16[name] = 0;
						if (title_has("num_dataLB")) {
							ps.size["num_data"] = ps.ui16[name] + 1;
						}
					}
				} else if (tc.api_name
						== ApiName_kGetNumberOfCoefficientsFloat) {
					if (vtype == PVType_kTL) {
						assert(desc_has("cspline"));
						ps.ui16[name] = 0;
					} else if (vtype == PVType_kLB) {
						ps.ui16[name] = desc_has("cspline") ? 1 : 0;
					} else if (vtype == PVType_kUB) {
						ps.ui16[name] = ps.ctxt["context"]->lsqfit_param;
					} else if (vtype == PVType_kTG) {
						ps.ui16[name] = ps.ctxt["context"]->lsqfit_param + 1;
					}
				} else {
					if (vtype == PVType_kLB) {
						ps.ui16[name] = 0;
						ps.size["num_coeff"] = ps.ui16[name] + 1;
						ps.AllocateAlignedDouble("coeff", "num_coeff");
					} else if (vtype == PVType_kTG) {
						ps.ui16[name] = 4;
						ps.size["num_coeff"] = ps.ui16[name] + 1;
						ps.AllocateAlignedDouble("coeff", "num_coeff");
					}
				}
			} else if (name == "npiece") {
				if (vtype == PVType_kUB) {
					ps.ui16[name] = UINT16_MAX;
					ps.size["num_data"] = ps.ui16[name] + 3;
					if (title_has("num_dataUB")) {
						ps.size["num_data"] = SIZE_MAX;
					}
				} else if (vtype == PVType_kLB) {
					ps.ui16[name] = 1;
					if (title_has("num_dataLB")) {
						ps.size["num_data"] = ps.ui16[name] + 3;
					}
				}
			} else if (name == "num_pieces") {
				if (vtype == PVType_kUB) {
					ps.size["num_data"] = SIZE_MAX;
					ps.size[name] = ps.size["num_data"] - 3;
				} else if (vtype == PVType_kLB) {
					ps.size[name] = 1;
					if (title_has("num_dataLB")) {
						ps.size["num_data"] = ps.size[name] + 3;
						ps.ctxt["context"] = ps.ctxt["csplinend4np1"];
					}
				} else if (vtype == PVType_kTG) {
					ps.size[name] = ps.size["num_data"] - 2;
				}
			} else if (name == "num_nwave") {
				if (vtype == PVType_kTL) {
					ps.size[name] = 0;
				} else if (vtype == PVType_kLB) {
					ps.size[name] = 1;
					ps.AllocateAlignedSizeT("nwave", name);
					ps.sptr["nwave"][0] = 0;
					ps.size["num_coeff"] = 1;
					ps.AllocateAlignedDouble("coeff", "num_coeff");
				}
			} else if (name == "nwave") {
				if (tc.api_name == ApiName_kCreateLSQFitContextSinusoidFloat) {
					if (vtype == PVType_kUB) {
						ps.ui16[name] = UINT16_MAX;
						ps.size["num_data"] = ps.ui16[name] * 2 + 2;
						if (title_has("num_dataUB")) {
							ps.size["num_data"] = SIZE_MAX;
						}
					} else if (vtype == PVType_kLB) {
						ps.ui16[name] = 0;
						if (title_has("num_dataLB")) {
							ps.size["num_data"] = ps.ui16[name] * 2 + 2;
						}
					}
				} else if ((tc.api_name == ApiName_kLSQFitSinusoidFloat)
						&& (title_has("nwaveDUP") || (title_has("nwaveIVO"))
								|| (title_has("nwaveOR"))
								|| (title_has("nwavePartial")))) {
					vector<string> nw = Split(tc.desc, ',');
					ps.size["num_nwave"] = nw.size();
					ps.AllocateAlignedSizeT(name, "num_nwave");
					bool has_zero = false;
					for (size_t j = 0; j < ps.size["num_nwave"]; ++j) {
						ps.sptr[name][j] = stoi(nw[j]);
						if (ps.sptr[name][j] == 0) {
							has_zero = true;
						}
					}
					ps.size["num_coeff"] = ps.size["num_nwave"] * 2
							- (has_zero ? 1 : 0);
					ps.AllocateAlignedDouble("coeff", "num_coeff");
				}
			} else if (name == "num_data") {
				if (tc.api_name
						== ApiName_kCreateLSQFitContextPolynomialFloat) {
					if (vtype == PVType_kUB) {
						ps.size[name] = SIZE_MAX;
					} else if (vtype == PVType_kLB) {
						ps.size[name] = ps.ui16["order"] + 1;
					} else if (vtype == PVType_kTL) {
						ps.size[name] = ps.ui16["order"];
					}
				} else if (tc.api_name
						== ApiName_kCreateLSQFitContextCubicSplineFloat) {
					if (vtype == PVType_kUB) {
						ps.size[name] = SIZE_MAX;
					} else if (vtype == PVType_kLB) {
						ps.size[name] = ps.ui16["npiece"] + 3;
					} else if (vtype == PVType_kTL) {
						ps.size[name] = ps.ui16["npiece"];
					}
				} else if (tc.api_name
						== ApiName_kCreateLSQFitContextSinusoidFloat) {
					if (vtype == PVType_kUB) {
						ps.size[name] = SIZE_MAX;
					} else if (vtype == PVType_kLB) {
						ps.size[name] = ps.ui16["nwave"] * 2 + 2;
					} else if (vtype == PVType_kTL) {
						ps.size[name] = ps.ui16["nwave"];
					}
				} else {
					if (vtype == PVType_kTL) {
						ps.size[name]--;
					} else if (vtype == PVType_kTG) {
						ps.size[name]++;
					} else if (vtype == PVType_kLB) {
						if (tc.api_name == ApiName_kLSQFitPolynomialFloat) {
							ps.ui16["order"] = 1;
							ps.size["num_coeff"] = ps.ui16["order"] + 1;
							ps.size[name] = ps.size["num_coeff"];
							ps.ctxt["context"] = ps.ctxt["polynd2o1"];
						} else if (tc.api_name
								== ApiName_kLSQFitSinusoidFloat) {
							ps.size["num_nwave"] = 1;
							ps.AllocateAlignedSizeT("nwave",
									ps.size["num_nwave"]);
							ps.sptr["nwave"][0] = 0;
							ps.size["num_coeff"] = 1;
							ps.size[name] = ps.size["num_coeff"] + 1;
							ps.ctxt["context"] = ps.ctxt["sinusoidnd2nw0"];
						}
					}
				}
			} else if (name == "mask") {
				if (title_has("effdataLB")) {
					size_t idx = 0;
					if (tc.api_name == ApiName_kLSQFitPolynomialFloat) {
						idx = ps.ui16["order"] + 1;
					} else if (tc.api_name == ApiName_kLSQFitCubicSplineFloat) {
						idx = ps.size["num_pieces"] + 3;
					} else if (tc.api_name == ApiName_kLSQFitSinusoidFloat) {
						ps.size["num_nwave"] = 2;
						ps.AllocateAlignedSizeT("nwave", ps.size["num_nwave"]);
						if (desc_has("0_in_nwave")) {
							ps.sptr["nwave"][0] = 0;
							ps.sptr["nwave"][1] = 1;
							ps.size["num_coeff"] = ps.size["num_nwave"] * 2 - 1;
						} else {
							ps.sptr["nwave"][0] = 1;
							ps.sptr["nwave"][1] = 2;
							ps.size["num_coeff"] = ps.size["num_nwave"] * 2;
						}
						idx = ps.size["num_coeff"];
					}
					for (size_t j = idx; j < ps.size["num_data"]; ++j) {
						ps.bptr[name][j] = false;
					}
				} else if (title_has("effdataTL")) {
					size_t idx = 0;
					if (tc.api_name == ApiName_kLSQFitPolynomialFloat) {
						idx = ps.ui16["order"];
					} else if (tc.api_name == ApiName_kLSQFitCubicSplineFloat) {
						idx = ps.size["num_pieces"] + 2;
					} else if (tc.api_name == ApiName_kLSQFitSinusoidFloat) {
						ps.size["num_nwave"] = 2;
						ps.AllocateAlignedSizeT("nwave", ps.size["num_nwave"]);
						if (desc_has("0_in_nwave")) {
							ps.sptr["nwave"][0] = 0;
							ps.sptr["nwave"][1] = 1;
							ps.size["num_coeff"] = ps.size["num_nwave"] * 2 - 1;
						} else {
							ps.sptr["nwave"][0] = 1;
							ps.sptr["nwave"][1] = 2;
							ps.size["num_coeff"] = ps.size["num_nwave"] * 2;
						}
						idx = ps.size["num_coeff"] - 1;
					}
					for (size_t j = idx; j < ps.size["num_data"]; ++j) {
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
				if (tc.api_name == ApiName_kLSQFitPolynomialFloat) {
					auto set_params = [&](uint16_t const order, size_t const ncoeff){
						ps.ui16["order"] = order;
						ps.size[name] = ncoeff;
						ps.AllocateAlignedDouble("coeff", name);
					};
					if (vtype == PVType_kTL) {
						set_params(2, 2);
					} else if (vtype == PVType_kLB) {
						set_params(2, 3);
					} else if (vtype == PVType_kUB) {
						set_params(3, 4);
					} else if (vtype == PVType_kTG) {
						set_params(3, 5);
					}
				} else if (tc.api_name == ApiName_kLSQFitSinusoidFloat) {
					auto set_params = [&](uint16_t const nnwave, size_t const ncoeff){
						ps.size["num_nwave"] = nnwave;
						ps.AllocateAlignedSizeT("nwave", "num_nwave");
						for (size_t j = 0; j < ps.size["num_nwave"]; ++j) {
							ps.sptr["nwave"][j] = j;
						}
						ps.size[name] = ncoeff;
						ps.AllocateAlignedDouble("coeff", name);
						ps.ctxt["context"] = ps.ctxt["sinusoidnw3"];
					};
					if (vtype == PVType_kTL) {
						set_params(3, 4);
					} else if (vtype == PVType_kLB) {
						set_params(3, 5);
					} else if (vtype == PVType_kUB) {
						set_params(4, 7);
					} else if (vtype == PVType_kTG) {
						set_params(4, 8);
					}
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
	case ApiName_kCreateLSQFitContextPolynomialFloat:
		run_status = LIBSAKURA_SYMBOL(CreateLSQFitContextPolynomialFloat)(
				ps.ftyp[tc.param[0].name], ps.ui16[tc.param[1].name],
				ps.size[tc.param[2].name], ps.cptr[tc.param[3].name]);
		break;
	case ApiName_kCreateLSQFitContextCubicSplineFloat:
		run_status = LIBSAKURA_SYMBOL(CreateLSQFitContextCubicSplineFloat)(
				ps.ui16[tc.param[0].name], ps.size[tc.param[1].name],
				ps.cptr[tc.param[2].name]);
		break;
	case ApiName_kCreateLSQFitContextSinusoidFloat:
		run_status = LIBSAKURA_SYMBOL(CreateLSQFitContextSinusoidFloat)(
				ps.ui16[tc.param[0].name], ps.size[tc.param[1].name],
				ps.cptr[tc.param[2].name]);
		break;
	case ApiName_kDestroyLSQFitContextFloat:
		run_status = LIBSAKURA_SYMBOL(DestroyLSQFitContextFloat)(
				ps.ctxt[tc.param[0].name]);
		if (run_status == LIBSAKURA_SYMBOL(Status_kOK)) {
			ps.ctxt["poly"] = nullptr;
		}
		break;
	case ApiName_kGetNumberOfCoefficientsFloat:
		run_status = LIBSAKURA_SYMBOL(GetNumberOfCoefficientsFloat)(
				ps.ctxt[tc.param[0].name], ps.ui16[tc.param[1].name],
				ps.sptr[tc.param[2].name]);
		break;
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
	case ApiName_kLSQFitCubicSplineFloat:
		run_status = LIBSAKURA_SYMBOL(LSQFitCubicSplineFloat)(
				ps.ctxt[tc.param[0].name], ps.size[tc.param[1].name],
				ps.size[tc.param[2].name], ps.fptr[tc.param[3].name],
				ps.bptr[tc.param[4].name], ps.fval[tc.param[5].name],
				ps.ui16[tc.param[6].name], ps.coef[tc.param[7].name],
				ps.fptr[tc.param[8].name], ps.fptr[tc.param[9].name],
				ps.bptr[tc.param[10].name], ps.fptr[tc.param[11].name],
				ps.sptr[tc.param[12].name], ps.fsta[tc.param[13].name]);
		break;
	case ApiName_kLSQFitSinusoidFloat:
		//LIBSAKURA_SYMBOL(LSQFitContextFloat) const *context,
		// size_t num_nwave,
		// size_t const nwave[/*num_nwave*/],
		// size_t num_data,
		// float const data[/*num_data*/],

		// bool const mask[/*num_data*/],
		// float clip_threshold_sigma,
		// uint16_t num_fitting_max,
		// size_t num_coeff,
		// double coeff[/*num_coeff*/],

		// float best_fit[/*num_data*/],
		// float residual[/*num_data*/],
		// bool final_mask[/*num_data*/],
		// float *rms,
		// LIBSAKURA_SYMBOL(LSQFitStatus) *lsqfit_status

		run_status = LIBSAKURA_SYMBOL(LSQFitSinusoidFloat)(
				ps.ctxt[tc.param[0].name], ps.size[tc.param[1].name],
				ps.sptr[tc.param[2].name], ps.size[tc.param[3].name],
				ps.fptr[tc.param[4].name], ps.bptr[tc.param[5].name],
				ps.fval[tc.param[6].name], ps.ui16[tc.param[7].name],
				ps.size[tc.param[8].name], ps.dptr[tc.param[9].name],
				ps.fptr[tc.param[10].name], ps.fptr[tc.param[11].name],
				ps.bptr[tc.param[12].name], ps.fptr[tc.param[13].name],
				ps.fsta[tc.param[14].name]);
		break;
	default:
		assert(false);
		break;
	}
}

void CheckStatus(TestCase const &tc, LIBSAKURA_SYMBOL(Status) const &status) {
	EXPECT_EQ(tc.expect_status, status);
}

void CheckValues(TestCase const &tc, ParamSet &ps) {
	if (!HasString(tc.title, "coeffNULL")
			&& (HasKey(ps.dptr, "coeff") || (HasKey(ps.coef, "coeff")))
			&& (HasKey(ps.size, "num_coeff") || (HasKey(ps.size, "num_pieces")))) {
		if (tc.api_name == ApiName_kLSQFitCubicSplineFloat) {
			for (size_t i = 0; i < ps.size["num_pieces"]; ++i) {
				cout << "        coeff[" << i << "]: ";
				for (size_t j = 0; j < 4; ++j) {
					if (1 <= j) {
						cout << ", ";
					}
					cout << ps.coef["coeff"][i][j];
				}
				cout << endl;
			}
		} else {
			cout << "        coeff: ";
			for (size_t i = 0; i < ps.size["num_coeff"]; ++i) {
				if (1 <= i) {
					cout << ", ";
				}
				cout << ps.dptr["coeff"][i];
			}
			cout << endl;
		}
	}
	if (!HasString(tc.title, "best_fitNULL") && (HasKey(ps.fptr, "best_fit"))
			&& (HasKey(ps.size, "num_data"))) {
		if (tc.api_name == ApiName_kLSQFitCubicSplineFloat) {
			cout << "        best_fit: ";
			for (size_t i = 0; i < ps.size["num_data"]; ++i) {
				if (1 <= i) {
					cout << ", ";
				}
				cout << ps.fptr["best_fit"][i];
			}
			cout << endl;
		} else {
			cout << "        best_fit: ";
			for (size_t i = 0; i < ps.size["num_data"]; ++i) {
				if (1 <= i) {
					cout << ", ";
				}
				cout << ps.fptr["best_fit"][i];
			}
			cout << endl;
		}
	}
	if (!HasString(tc.title, "residualNULL") && (HasKey(ps.fptr, "residual"))
			&& (HasKey(ps.size, "num_data"))) {
		if (tc.api_name == ApiName_kLSQFitCubicSplineFloat) {
			cout << "        residual: ";
			for (size_t i = 0; i < ps.size["num_data"]; ++i) {
				if (1 <= i) {
					cout << ", ";
				}
				cout << ps.fptr["residual"][i];
			}
			cout << endl;
		} else {
			cout << "        residual: ";
			for (size_t i = 0; i < ps.size["num_data"]; ++i) {
				if (1 <= i) {
					cout << ", ";
				}
				cout << ps.fptr["residual"][i];
			}
			cout << endl;
		}
	}
}

void Execute(TestCase const &tc, ParamSet &ps) {
	if (tc.category != TCat_kPF) {
		assert(tc.num_repeat == 1);
	}
	double time_elapsed = 0.0;
	LIBSAKURA_SYMBOL (Status) run_status;
	for (size_t i = 0; i < tc.num_repeat; ++i) {
		double time_start = GetCurrentTime();
		RunApi(tc, ps, run_status);
		double time_end = GetCurrentTime();
		time_elapsed += (time_end - time_start);
	}
	if (tc.category == TCat_kPF) {
		cout << setprecision(5) << "#x# benchmark Lsq_" << tc.title << " "
				<< time_elapsed << endl;
	}
	CheckStatus(tc, run_status);
	if ((tc.category == TCat_kOK) || (tc.category == TCat_kPF)) {
		//CheckValues(tc, ps);
	}
}

template<typename T> void FreeParamSet(ParamSet &ps, T &ptr) {
	for (auto x : ptr) {
		ps.Free(x.first);
	}
}

void Epilogue(ParamSet &ps) {
	FreeParamSet(ps, ps.sto);
	FreeParamSet(ps, ps.sptr);
	FreeParamSet(ps, ps.fptr);
	FreeParamSet(ps, ps.dptr);
	FreeParamSet(ps, ps.bptr);
	FreeParamSet(ps, ps.ctxt);
	FreeParamSet(ps, ps.coef);
	FreeParamSet(ps, ps.fsta);
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

//Tests for CreateLSQFitContextPolynomialFloat()
TEST_F(Lsq, CreateLSQFitContextPolynomialFloat) {
	RunTest(ApiName_kCreateLSQFitContextPolynomialFloat);
}
//Tests for CreateLSQFitContextCubicSplineFloat()
TEST_F(Lsq, CreateLSQFitContextCubicSplineFloat) {
	RunTest(ApiName_kCreateLSQFitContextCubicSplineFloat);
}
//Tests for CreateLSQFitContextSinusoidFloat()
TEST_F(Lsq, CreateLSQFitContextSinusoidFloat) {
	RunTest(ApiName_kCreateLSQFitContextSinusoidFloat);
}
//Tests for DestroyLSQFitContextFloat()
TEST_F(Lsq, DestroyLSQFitContextFloat) {
	RunTest(ApiName_kDestroyLSQFitContextFloat);
}
//Tests for GetNumberOfCoefficientsFloat()
TEST_F(Lsq, GetNumberOfCoefficientsFloat) {
	RunTest(ApiName_kGetNumberOfCoefficientsFloat);
}
//Tests for LSQFitPolynomialFloat()
TEST_F(Lsq, LSQFitPolynomialFloat) {
	RunTest(ApiName_kLSQFitPolynomialFloat);
}
//Tests for LSQFitCubicSplineFloat()
TEST_F(Lsq, LSQFitCubicSplineFloat) {
	RunTest(ApiName_kLSQFitCubicSplineFloat);
}
//Tests for LSQFitSinusoidFloat()
TEST_F(Lsq, LSQFitSinusoidFloat) {
	RunTest(ApiName_kLSQFitSinusoidFloat);
}

//*******************************************************
//Test variadic template handling
//*******************************************************
/*
 static map<string, PVType> params_test;
 static bool has_plength = false;
 void AddTestCase(LIBSAKURA_SYMBOL(Status) &ret_status, string const &desc,
 string const &name, PVType const &pvtype) {
 if (!has_plength) {
 cout << "***sizeof...tail = " << 0 << flush << endl;
 has_plength = true;
 }
 cout << "returned_status = " << flush;
 string status = "undefined";
 if (ret_status == Ret_kOK) {
 status = "OK";
 } else if (ret_status == Ret_kNG) {
 status = "NG";
 } else if (ret_status == Ret_kIA) {
 status = "IA";
 } else if (ret_status == Ret_kNM) {
 status = "NM";
 }
 cout << status << flush << endl;
 cout << "description = " << desc << flush << endl;
 params_test[name] = pvtype;
 }
 template<typename ... Types>
 void AddTestCase(LIBSAKURA_SYMBOL(Status) &ret_status, string const &desc,
 string const &name, PVType const &pvtype, Types ... tail) {
 if (!has_plength) {
 cout << "***sizeof...tail = " << sizeof...(tail) << flush << endl;
 has_plength = true;
 }
 params_test[name] = pvtype;
 AddTestCase(ret_status, desc, tail...);
 }
 void ShowParamsTest() {
 cout << "specify " << params_test.size() << " parameter" << (params_test.size() > 1 ? "s" : "") << flush << endl;
 for (auto x : params_test) {
 cout << "{param name = " << x.first << ", type = ";
 PVType pvtype = x.second;
 string vtype = "undefined";
 if (pvtype == PVType_kNULL) {
 vtype = "NULLPointer";
 } else if (pvtype == PVType_kNAL) {
 vtype = "NotAligned";
 } else if (pvtype == PVType_kLB) {
 vtype = "LowerBoundary";
 } else if (pvtype == PVType_kIR) {
 vtype = "InRange";
 } else if (pvtype == PVType_kVP) {
 vtype = "ValidPointer";
 }
 cout << vtype << "}" << flush << endl;
 }
 }
 TEST_F(Lsq, AddTestCase) {
 cout << "----------" << flush << endl;
 params_test.clear();
 has_plength = false;
 AddTestCase(Ret_kOK, "desc1", "data", PVType_kNULL);
 ShowParamsTest();
 cout << "----------" << flush << endl;
 params_test.clear();
 has_plength = false;
 AddTestCase(Ret_kNG, "desc2", "mask", PVType_kNAL, "num_data", PVType_kIR);
 ShowParamsTest();
 cout << "----------" << flush << endl;
 params_test.clear();
 has_plength = false;
 AddTestCase(Ret_kIA, "desc3", "coeff", PVType_kVP, "num_coeff", PVType_kLB,
 "final_mask", PVType_kNULL);
 ShowParamsTest();
 cout << "----------" << flush << endl;
 }
 */
