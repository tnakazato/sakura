#ifndef _SAKURA_E2E_OPTION_PARSER_H_
#define _SAKURA_E2E_OPTION_PARSER_H_

#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <utility>
#include <locale>
#include <map>

#include <libsakura/sakura.h>

#include "config_file_reader.h"

namespace {
struct CalibrationOptions {
	std::string sky_table;
	std::string tsys_table;
	unsigned int tsys_ifno;
};

struct FlaggingOptions {
	size_t edge_channels;
	float clipping_threshold;
	bool do_clipping;
};

struct BaselineOptions {
	sakura_BaselineType baseline_type;
	uint16_t order;
	float clipping_threshold;
	uint16_t num_fitting_max;
	std::vector<size_t> line_mask;
};

struct SmoothingOptions {
	sakura_Convolve1DKernelType kernel_type;
	size_t kernel_width;
	bool use_fft;
};

struct E2EOptions {
	std::string config_file;
	std::string input_file;
	std::string output_file;
	bool serialize;
	unsigned int max_threads;
	unsigned int ifno;
	CalibrationOptions calibration;
	FlaggingOptions flagging;
	BaselineOptions baseline;
	SmoothingOptions smoothing;
};

class OptionParser {
public:
	static void Parse(std::string const input_file, E2EOptions *option_list) {
		OptionMap options;
		ConfigFileReader::read(input_file, &options);
		// parse config file name
		option_list->config_file = GetValue(options, GetKey("config"));
		ParseE2e(options, option_list);
		ParseCalibration(options, &(option_list->calibration));
		ParseFlagging(options, &(option_list->flagging));
		ParseBaseline(options, &(option_list->baseline));
		ParseSmoothing(options, &(option_list->smoothing));
	}
	static std::string GetSummary(E2EOptions const option_list) {
		std::ostringstream oss;
		oss << "config file (" << option_list.config_file << ") summary:\n";
		oss << "\tinput filename=" << option_list.input_file
				<< "\n\toutput filename=" << option_list.output_file
				<< "\n\tspw=" << option_list.ifno << "\n\tserialize="
				<< option_list.serialize << "\n\tmax_threads="
				<< option_list.max_threads << "\n";
		oss << "\tsky filename=" << option_list.calibration.sky_table
				<< "\n\ttsys filename=" << option_list.calibration.tsys_table
				<< "\n\ttsys_spw=" << option_list.calibration.tsys_ifno << "\n";
		oss << "\tedge channels=" << option_list.flagging.edge_channels
				<< "\n\tclipping threshold=";
		if (option_list.flagging.do_clipping) {
			oss << option_list.flagging.clipping_threshold << "\n";
		} else {
			oss << "\n";
		}
		char separator = '[';
		oss << "\tline mask=";
		for (auto i = option_list.baseline.line_mask.begin();
				i != option_list.baseline.line_mask.end(); ++i) {
			oss << separator << *i;
			separator = ',';
		}
		oss << "]\n";
		oss << "\tbaseline function=" << option_list.baseline.baseline_type
				<< "\n";
		oss << "\torder=" << option_list.baseline.order << "\n";
		oss << "\tclipping_threshold=" << option_list.baseline.clipping_threshold << "\n";
		oss << "\tnum_fitting_max=" << option_list.baseline.num_fitting_max << "\n";
		oss << "\tkernel type=" << option_list.smoothing.kernel_type << "\n";
		oss << "\tkernel width=" << option_list.smoothing.kernel_width << "\n";
		oss << "\tuse_fft=" << option_list.smoothing.use_fft << "\n";
		return oss.str();
	}
private:
	static void ParseE2e(OptionMap const options, E2EOptions *option_list) {
		// parse input file name
		option_list->input_file = GetValue(options, GetKey("input"));

		// parse output file name
		option_list->output_file = GetValue(options, GetKey("output"),
				((option_list->input_file.rfind("/")
						== option_list->input_file.length() - 1) ?
						option_list->input_file.substr(0,
								option_list->input_file.length() - 1) :
						option_list->input_file) + "_out");

		// parse IFNO
		option_list->ifno = std::atoi(GetValue(options, GetKey("spw")).c_str());
		option_list->serialize =
				"TRUE" == ToUpper(GetValue(options, GetKey("serialize"))) ?
						true : false;
		option_list->max_threads = std::max(1,
				std::atoi(GetValue(options, GetKey("max_threads")).c_str()));
	}

	static void ParseCalibration(OptionMap const options,
			CalibrationOptions *option_list) {
		std::string const category = "calibration";

		// parse sky table name
		option_list->sky_table = GetValue(options, GetKey(category, "sky"));

		// parse tsys table name
		option_list->tsys_table = GetValue(options, GetKey(category, "tsys"));

		// parse IFNO for Tsys
		option_list->tsys_ifno = std::atoi(
				GetValue(options, GetKey(category, "tsys_spw")).c_str());
	}

	static void ParseFlagging(OptionMap const options,
			FlaggingOptions *option_list) {
		std::string const category = "flagging";

		// parse edge channels
		std::string s = GetValue(options, GetKey(category, "edge"));
		option_list->edge_channels = (s.size() > 0) ? std::atoi(s.c_str()) : 0;

		// parse clipping threshold
		s = GetValue(options, GetKey(category, "clipping_threshold"));
		if (s.size() > 0) {
			option_list->do_clipping = true;
			option_list->clipping_threshold = std::atof(s.c_str());
		} else {
			// no clipping by default
			option_list->do_clipping = false;
			option_list->clipping_threshold = 0.0;
		}
	}

	static void ParseBaseline(OptionMap const options,
			BaselineOptions *option_list) {
		std::string const category = "baseline";

		// parse baseline type
		option_list->baseline_type = StringToBaselineType(
				GetValue(options, GetKey(category, "fitfunc")));

		// parse order
		std::string s = GetValue(options, GetKey(category, "order"));
		option_list->order = (s.size() > 0) ? std::atoi(s.c_str()) : 0;

		// parse clipping threshold
		s = GetValue(options, GetKey(category, "clipping_threshold"));
		option_list->clipping_threshold = (s.size() > 0) ? std::atof(s.c_str()) : 0.0;

		// parse num_fitting_max
		s = GetValue(options, GetKey(category, "num_fitting_max"));
		option_list->num_fitting_max = (s.size() > 0) ? std::atoi(s.c_str()) : 1;

		// parse line mask
		s = GetValue(options, GetKey(category, "mask"));
		if (s.size() > 0) {
			ToVector(s, &(option_list->line_mask));
		} else {
			option_list->line_mask.clear();
		}
	}
	static void ParseSmoothing(OptionMap const options,
			SmoothingOptions *option_list) {
		std::string const category = "smoothing";

		// parse kernel type
		option_list->kernel_type = StringToKernelType(
				GetValue(options, GetKey(category, "kernel")));

		// parse kernel width
		std::string s = GetValue(options, GetKey(category, "kernel_width"));
		if (s.size() > 0) {
			option_list->kernel_width = std::atoi(s.c_str());
		} else {
			option_list->kernel_width = 5; // default is 5
		}

		// parse use_fft
		s = ToUpper(GetValue(options, GetKey(category, "use_fft")));
		option_list->use_fft = (s == "FALSE") ? false : true;
	}
	static std::string GetKey(std::string const category,
			std::string const attribute = "") {
		if (category.size() > 0) {
			std::string s = prefix + delimiter + category;
			if (attribute.size() > 0) {
				s += delimiter + attribute;
			}
			return s;
		} else {
			return "";
		}
	}
	static std::string GetValue(OptionMap const options, std::string const key,
			std::string const default_value = "") {
		auto entry = options.find(key);
		if (entry != options.end()) {
			return entry->second;
		} else {
			return default_value;
		}
	}
	static void ToVector(std::string const s, std::vector<size_t> *v) {
		// convert [0,0] or [[0,0],[1,1]] type string to vector
		char const separator = ',';
		char const left_bracket = '[';
		char const right_bracket = ']';

		// input string is empty
		if (s.size() == 0) {
			v->clear();
			return;
		}

		// first character must be '[' and the last character must be ']'
		if (s[0] == left_bracket && s[s.size() - 1] == right_bracket) {
			std::string substring = Trim(s.substr(1, s.size() - 2));
			if (substring[0] == left_bracket) {
				v->clear();
				// multiple entries
				size_t pos1 = substring.find_first_of(left_bracket);
				size_t pos2 = substring.find_first_of(right_bracket, pos1);
				while (pos1 != std::string::npos && pos2 != std::string::npos) {
					std::vector<size_t> temporary_vector;
					ToVector(substring.substr(pos1, pos2 - pos1 + 1),
							&temporary_vector);
					v->insert(v->end(), temporary_vector.begin(),
							temporary_vector.end());
					pos1 = substring.find_first_of(left_bracket, pos2);
					pos2 = substring.find_first_of(right_bracket, pos1);
				}
			} else {
				// single entry
				size_t pos = substring.find(separator);
				v->resize(2);
				(*v)[0] = std::atoi(Trim(substring.substr(0, pos)).c_str());
				(*v)[1] = std::atoi(Trim(substring.substr(pos + 1)).c_str());
				if ((*v)[0] > (*v)[1]) {
					std::swap((*v)[0], (*v)[1]);
				}
				return;
			}
		} else {
			v->clear();
			return;
		}
	}
	static inline sakura_Convolve1DKernelType StringToKernelType(
			std::string const type) {
		std::string type_upper = ToUpper(type);
		std::map<std::string, sakura_Convolve1DKernelType> typemap;
		typemap["GAUSSIAN"] = sakura_Convolve1DKernelType_kGaussian;
		typemap["BOXCAR"] = sakura_Convolve1DKernelType_kBoxcar;
		typemap["HANNING"] = sakura_Convolve1DKernelType_kHanning;
		typemap["HAMMING"] = sakura_Convolve1DKernelType_kHamming;
		typemap["DEFAULT"] = sakura_Convolve1DKernelType_kGaussian;
		std::map<std::string, sakura_Convolve1DKernelType>::iterator entry =
				typemap.find(type_upper);
		if (entry != typemap.end()) {
			return entry->second;
		} else {
			return typemap["DEFAULT"];
		}
	}
	static inline sakura_BaselineType StringToBaselineType(
			std::string const type) {
		std::string type_upper = ToUpper(type);
		std::map<std::string, sakura_BaselineType> typemap;
		typemap["POLYNOMIAL"] = sakura_BaselineType_kPolynomial;
		typemap["CHEBYSHEV"] = sakura_BaselineType_kChebyshev;
		typemap["CUBICSPLINE"] = sakura_BaselineType_kCubicSpline;
		typemap["SINUSOID"] = sakura_BaselineType_kSinusoid;
		typemap["DEFAULT"] = sakura_BaselineType_kChebyshev;
		std::map<std::string, sakura_BaselineType>::iterator entry =
				typemap.find(type_upper);
		if (entry != typemap.end()) {
			return entry->second;
		} else {
			return typemap["DEFAULT"];
		}
	}
	static inline std::string ToUpper(std::string const str) {
		std::locale loc;
		std::string out(str);
		for (std::string::iterator i = out.begin(); i != out.end(); ++i) {
			*i = std::toupper(*i, loc);
		}
		return out;
	}
	static std::string const prefix;
	static char const delimiter;
};

std::string const OptionParser::prefix = "sakura_e2e";
char const OptionParser::delimiter = '.';

} /* anonymous namespace */

#endif /* _SAKURA_E2E_OPTION_PARSER_H_ */
