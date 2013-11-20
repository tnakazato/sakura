#ifndef _SAKURA_E2E_OPTION_PARSER_H_
#define _SAKURA_E2E_OPTION_PARSER_H_

#include <string>
#include <vector>
#include <cstdlib>

#include "config_file_reader.h"

namespace {

class OptionParser {
public:
	static void ParseE2e(OptionList const options, std::string *input_file,
			std::string *output_file, unsigned int *ifno) {
		// parse input file name
		*input_file = GetValue(options, GetKey("input"));

		// parse output file name
		*output_file = GetValue(options, GetKey("output"),
				((input_file->rfind("/") == input_file->length() - 1) ?
						input_file->substr(0, input_file->length() - 1) :
						*input_file) + "_out");

		// parse IFNO
		*ifno = std::atoi(GetValue(options, GetKey("spw")).c_str());
	}

	static void ParseCalibration(OptionList const options,
			std::string *sky_table, std::string *tsys_table, unsigned int *tsys_ifno) {
		std::string const category = "calibration";

		// parse sky table name
		*sky_table = GetValue(options, GetKey(category, "sky"));

		// parse tsys table name
		*tsys_table = GetValue(options, GetKey(category, "tsys"));

		// parse IFNO for Tsys
		*tsys_ifno = std::atoi(GetValue(options, GetKey(category, "tsys_spw")).c_str());
	}

	static void ParseFlagging(OptionList const options, int *edge_channels,
			float *clipping_threshold, bool *do_clipping) {
		std::string const category = "flagging";

		// parse edge channels
		std::string s = GetValue(options, GetKey(category, "edge"));
		*edge_channels = (s.size() > 0) ? std::atoi(s.c_str()) : 0;

		// parse clipping threshold
		s = GetValue(options, GetKey(category, "clipping_threshold"));
		if (s.size() > 0) {
			*do_clipping = true;
			*clipping_threshold = std::atof(s.c_str());
		}
		else {
			*do_clipping = false;
			*clipping_threshold = 0.0;
		}
	}

	static void ParseBaseline(OptionList const options, std::vector<uint64_t> *line_mask) {
		std::string const category = "baseline";

		// parse line mask
		std::string s = GetValue(options, GetKey(category, "mask"));
		if (s.size() > 0) {
			ToVector(s, line_mask);
		}
		else {
			line_mask->clear();
		}
	}
private:
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
	static std::string GetValue(OptionList const options, std::string const key,
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
		if (s[0] == left_bracket && s[s.size()-1] == right_bracket) {
			std::string substring = Trim(s.substr(1,s.size()-2));
			if (substring[0] == left_bracket) {
				v->clear();
				// multiple entries
				size_t pos1 = substring.find_first_of(left_bracket);
				size_t pos2 = substring.find_first_of(right_bracket, pos1);
				while (pos1 != std::string::npos && pos2 != std::string::npos) {
					std::vector<uint64_t> temporary_vector;
					ToVector(substring.substr(pos1,pos2-pos1+1), &temporary_vector);
					v->insert(v->end(), temporary_vector.begin(), temporary_vector.end());
					pos1 = substring.find_first_of(left_bracket, pos2);
					pos2 = substring.find_first_of(right_bracket, pos1);
				}
			}
			else {
				// single entry
				size_t pos = substring.find(separator);
				v->resize(2);
				(*v)[0] = std::atoi(Trim(substring.substr(0,pos)).c_str());
				(*v)[1] = std::atoi(Trim(substring.substr(pos+1)).c_str());
				if ((*v)[0] > (*v)[1]) {
					uint64_t tmp = (*v)[0];
					(*v)[0] = (*v)[1];
					(*v)[1] = tmp;
				}
				return;
			}
		}
		else {
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
