#include "options.h"
#include "utilities.h"
#include "taxdistinct.h"
#include <chrono>

#include <time.h>

using namespace std;
using namespace std::chrono; 


/* Flag set by ‘--verbose’. */

int main (int argc, char **argv)
{
    time_t timer;

  INPUT_OPTIONS options;

  read_options(argc, argv, options);

  auto start = high_resolution_clock::now(); 
  //compute_tax_distance(options.tax_file_name);
  compute_pwy_tax_distance(options.tax_file_name, options.pwy_file_name);
  
  //std::cout << "Subtree count at root: " << subtree_count(tree) << endl; 
  //std::cout << "Delta* " << compute_delta_star(tree) << endl; 

  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 
  std::cout << "Time taken by function: " << duration.count()/1000000 << " seconds" << endl; 

  //std::cout << "Subtree count at root: " << subtree_count(tree) << endl; 
  //std::cout << "Total tree sum : " << sum_tree_count(tree) << std::endl;

  return 0;

}

