#ifndef __METATAXA_TYPES__
#define __METATAXA_TYPES__

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

typedef struct node_ NODE;   
typedef std::unordered_map <std::string, NODE *> CHILD_NODES;

typedef std::unordered_map<std::string, std::vector<std::string> *> STRING_STRING_VECTOR_MAP;
typedef std::unordered_map<std::string, std::string> STRING_STRING_MAP;

// Structure for storing the information related to a node
// in the tree contstructed from the data
typedef struct node_ {
    // map of taxon name to NODE * for the tree
    CHILD_NODES children;   

    // number of counts for the taxa
    int count;

    // sum total of the counts in the subtree at the node
    int subtree_count;

    // depth in the tree
    int depth;
   
    // name of the taxon
    std::string name;

    // initialize the counts to 0
    node_() {
      count = 0;
      subtree_count = 0;
      depth = 0;
    }
} NODE;







#endif
