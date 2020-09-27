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
#include <numeric>
#include <iomanip>      // std::setw


#include "metataxa_types.h"
#include "utilities.h"

using namespace std;

// Type of delta
typedef enum _DeltaType {DELTA_PLUS, DELTA_STAR, DELTA} DeltaType;

/**
 * @brief read_tax_file  read the tha functional and taxonomy table 
 *
 * @param tax_file_name reads the functional and taxonomic file 
 * @return NODE for the root of a tree
*/
NODE *read_tax_file(const string &tax_file_name);

void print_tree(NODE *node) ;
int sum_tree_count(NODE *node);

int subtree_count(NODE *node);
/**
 * @brief compute_delta_star This function computes the Delta* for the tax tree
 *   rooted at the NODE* node
 *
 * @detail The basic idea is toe compute the number of paths that crosses an 
 *   edge when we consider a pair of taxons. Suppose we have a an edge e in the 
 *   tree such that the path between two taxa A and B goes through e (n->m). The the total 
 *   number of pairs that though e is simply (tot -m->count)*m->count.  Then if we want
 *   to count the number of paths that goes through e for all pairs of taxons then we 
 *   simply compute the sum-product as (tot- m->count)*m->count. 
 *    
 * @param node        the node of a tree    
 * @param deltatype   the type of Delta    
 * @param use_wtd     use WTD    
 *
 * @return delta* for the tree rooted at node    
*/
double compute_delta(NODE *node, DeltaType deltatype, bool use_wtd); 

/**
 * @brief compute_tax_distance computes the three taxonomic distance from the
 *   file as outputted by Metapathways
 *
 * @param tax_file_name  The taxonomic tax file
 * @param use_wtd  Flag to turn on Weighted Taxonomic Distance (WTD)
 *
*/
void compute_tax_distance(const std::string &tax_file_name);

/**
 * @brief compute_pwy_tax_distance  computes the three taxonomic distance from the
 *   file as outputted by Metapathways
 *
 * @param tax_file_name      The taxonomic tax file
 * @param pwy_orf_file_name  The pwy to list of  orfs in a row
 *
*/
void compute_pwy_tax_distance(const std::string &tax_file_name, 
  const std::string &pwy_orf_fliename);
#endif
