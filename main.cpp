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

  NODE *tree = read_tax_file(options.tax_file_name) ;

  print_tree(tree);

  auto start = high_resolution_clock::now(); 
  sum_tree_count(tree);
  
  std::cout << "Subtree count at root: " << subtree_count(tree) << endl; 
  std::cout << "Delta* " << compute_delta_star(tree) << endl; 

  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 
  std::cout << "Time taken by function: " << duration.count() << " microseconds" << endl; 

  std::cout << "Subtree count at root: " << subtree_count(tree) << endl; 
  std::cout << "Total tree sum : " << sum_tree_count(tree) << std::endl;

  return 0;

}

