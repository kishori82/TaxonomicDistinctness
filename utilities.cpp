#include "utilities.h"

using namespace std;

std::string trim(const std::string& str, const std::string& whitespace = " \t") {
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}


vector<string> split(const string& input, const char delim) {
    vector<std::string> result;
    std::istringstream ss(input);
    std::string token;

    while (std::getline(ss, token, delim)) {
       result.push_back(token);
    }
    return result;
}


std::string taxon_concat(const std::vector<std::string> &taxons) {
    std::string taxon;
    for (auto it = taxons.begin(); it != taxons.end(); it++) {
        taxon += *it + ";";
    }
    return taxon;
}

