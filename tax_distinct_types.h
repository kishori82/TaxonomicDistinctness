#ifndef __TAX_DISTINCT_TYPES__
#define __TAX_DISTINCT_TYPES__

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
using namespace std;

typedef struct _input_options {
  std::string tax_file_name;
} INPUT_OPTIONS;

typedef struct node_ NODE;   
typedef std::unordered_map <std::string, NODE *> CHILD_NODES;

typedef struct node_ {
    CHILD_NODES children;   
    int count;
    int subtree_count;
    std::string name;

    node_() {
      count = 0;
      subtree_count = 0;
    }
} NODE;







#endif
