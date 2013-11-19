#ifndef _SAKURA_E2E_CONFIG_FILE_READER_H_
#define _SAKURA_E2E_CONFIG_FILE_READER_H_

#include <string>
#include <fstream>
#include <map>

namespace {

typedef std::map<std::string, std::string> OptionList;

std::string Trim(std::string const s) {
  std::string trimmed_string = "";
  std::string const trim_char_list(" \t\n");
  std::string::size_type left = s.find_first_not_of(trim_char_list);
  if (left != std::string::npos) {
    std::string::size_type right = s.find_last_not_of(trim_char_list);
    trimmed_string = s.substr(left,right-left+1);
  }
  return trimmed_string;
}

class ConfigFileReader {
public:
	static void read(std::string const input_file, OptionList *options) {
		// clear options
		options->clear();

		std::ifstream f(input_file);
		if (!f.is_open())
			return;

		std::string s, key, value;
		while (!f.eof()) {
			std::getline(f, s);
			if (s.size() > 0 && s.at(s.find_first_not_of(" \t")) != '#') {
				SplitKeyAndValue(s, key, value);
				(*options)[key] = value;
			}
		}
	}

	static void SplitKeyAndValue(std::string const &s, std::string &key, std::string &value) {
	  std::string::size_type pos = s.find("=");
	  // TODO: support inline comment
	  if (pos == std::string::npos) {
	    key = "";
	    value = "";
	  }
	  else {
	    key = Trim(s.substr(0,pos));
	    value = Trim(s.substr(pos+1));
	  }
	}

	static int SplitKey(std::string const key, std::vector<std::string> &key_list) {
	  int count = 0;
	  for (std::string::const_iterator i = key.begin(); i != key.end(); ++i) {
	    count += (*i == '.') ? 1 : 0;
	  }
	  std::string::size_type start = 0;
	  std::string::size_type end = std::string::npos;
	  for (int i = 0; i < count; ++i) {
	    end = key.find('.',start);
	    key_list.push_back(key.substr(start,end-start));
	    start = end + 1;
	  }
	  std::string s = key.substr(start);
	  if (s.size() > 0) {
	    key_list.push_back(s);
	  }
	  return static_cast<int>(key_list.size());
	}
};

} /* anonymous namespace */

#endif /* _SAKURA_E2E_CONFIG_FILE_READER_H_ */
