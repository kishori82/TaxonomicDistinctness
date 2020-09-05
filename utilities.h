#ifndef __OPTIMUS_UTILITES__
#define __OPTIMUS_UTILITES__
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unordered_map>
#include <regex>



#include "tax_distinct_types.h"

using namespace std;

void read_options(int, char **, INPUT_OPTIONS &); 

std::string taxon_concat(const std::vector<std::string> &taxons) ;

vector<string> split(const string& input, const char delim);
std::string trim(const std::string& str, const std::string& whitespace);

#endif
