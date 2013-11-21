#ifndef _SAKURA_E2E_OPTION_PARSER_H_
#define _SAKURA_E2E_OPTION_PARSER_H_

#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <utility>

#include "config_file_reader.h"

namespace {
struct CalibrationOptions {
	std::string sky_table;
	std::string tsys_table;
	unsigned int tsys_ifno;
};

struct FlaggingOptions {
	int edge_channels;
	float clipping_threshold;
	bool do_clipping;
};

struct BaselineOptions {
	std::vector<uint64_t> line_mask;
};

struct SmoothingOptions {
	std::string kernel_type;
	size_t kernel_width;
};

struct E2EOptions {
	std::string config_file;
	std::string input_file;
	std::string output_file;
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
				<< "\n\tspw=" << option_list.ifno << "\n";
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
		oss << "\tkernel type=" << option_list.smoothing.kernel_type << "\n";
		oss << "\tkernel width=" << option_list.smoothing.kernel_width << "\n";
		return oss.str();
	}
private:
	static void ParseE2e(OptionMap const options, E2EOptions *option_list) {
		// parse input file name
		option_list->input_file = GetValue(options, GetKey("input"));

		// parse output file name
		option_list->output_file = GetValue(options, GetKey("output"),
				((option_list->input_file.rfind("/") == option_list->input_file.length() - 1) ?
						option_list->input_file.substr(0, option_list->input_file.length() - 1) :
						option_list->input_file) + "_out");

		// parse IFNO
		option_list->ifno = std::atoi(GetValue(options, GetKey("spw")).c_str());
	}

	static void ParseCalibration(OptionMap const options, CalibrationOptions *option_list) {
		std::string const category = "calibration";

		// parse sky table name
		option_list->sky_table = GetValue(options, GetKey(category, "sky"));

		// parse tsys table name
		option_list->tsys_table = GetValue(options, GetKey(category, "tsys"));

		// parse IFNO for Tsys
		option_list->tsys_ifno = std::atoi(
				GetValue(options, GetKey(category, "tsys_spw")).c_str());
	}

	static void ParseFlagging(OptionMap const options, FlaggingOptions *option_list) {
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

	static void ParseBaseline(OptionMap const options, BaselineOptions *option_list) {
		std::string const category = "baseline";

		// parse line mask
		std::string s = GetValue(options, GetKey(category, "mask"));
		if (s.size() > 0) {
			ToVector(s, &(option_list->line_mask));
		} else {
			option_list->line_mask.clear();
		}
	}
	static void ParseSmoothing(OptionMap const options, SmoothingOptions *option_list) {
		std::string const category = "smoothing";

		// parse kernel type
		option_list->kernel_type = GetValue(options, GetKey(category, "kernel"));

		// parse kernel width
		std::string s = GetValue(options, GetKey(category, "kernel_width"));
		if (s.size() > 0) {
			option_list->kernel_width = std::atoi(s.c_str());
		}
		else {
			option_list->kernel_width = 5; // default is 5
		}
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
	static void ToVector(std::string const s, std::vector<uint64_t> *v) {
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
					std::vector<uint64_t> temporary_vector;
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
	static std::string const prefix;
	static char const delimiter;
};

std::string const OptionParser::prefix = "sakura_e2e";
char const OptionParser::delimiter = '.';

} /* anonymous namespace */

#endif /* _SAKURA_E2E_OPTION_PARSER_H_ */
