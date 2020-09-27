#ifndef __METATAXA_UTILITES__
#define __METATAXA_UTILITES__
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <regex>
#include <string>

#include "metataxa_types.h"

using namespace std;

/**

*/
std::string taxon_concat(const std::vector<std::string> &taxons) ;

/**


*/
vector<string> split(const string& input, const char delim);

/**

*/
std::string trim(const std::string& str, const std::string& whitespace);

// functor for getting sum of previous result and square of current element

template<typename T>
struct square
{
    T operator()(const T& Left, const T& Right) const
    {   
        return (Left + Right*Right);
    }
};

#endif
