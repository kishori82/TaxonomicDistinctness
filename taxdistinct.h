#ifndef __TAX_DISTINCT__
#define __TAX_DISTINCT__

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
#include "utilities.h"

using namespace std;

NODE *read_tax_file(const string &tax_file_name);

void print_tree(NODE *node) ;
int sum_tree_count(NODE *node);

int subtree_count(NODE *node);

double compute_delta_star(NODE *node); 
#endif
