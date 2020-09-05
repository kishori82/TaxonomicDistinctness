#include "utilities.h"
#include <cinttypes>
using namespace std;

void insert_node(NODE **node, const vector<std::string> &tax_hierarchy, unsigned int depth) {
    NODE *n;

    if (*node == nullptr) {
      *node  = new NODE;
      (*node)->name = tax_hierarchy[depth];
    }
    if (tax_hierarchy.size() == 1 ) (*node)->count++;

    if (tax_hierarchy.size() > depth + 1 ) {
       depth = depth + 1;
       if ((*node)->children.find(tax_hierarchy[depth]) == (*node)->children.end() ) {
         n = new NODE;
         n->name = tax_hierarchy[depth];
         (*node)->children.insert( {tax_hierarchy[depth], n }) ;
            //   std::cout <<  "\t" << n->name << " : " << depth << std::endl;
       } else {
         n = (*node)->children.at(tax_hierarchy[depth]);
       }

       if( tax_hierarchy.size() == depth + 1 ) n->count++;

       //std::cout << "\t" << depth << std::endl;
       insert_node(&n, tax_hierarchy, depth) ;
    }
}

void _sum_tree_count_(NODE *node, int &sum) {
   sum = sum + node->count;
   for (auto it = node->children.begin(); it != node->children.end(); it++) {
      _sum_tree_count_(it->second, sum);
   }
}

int sum_tree_count(NODE *node) {
   int sum = 0;
   _sum_tree_count_(node, sum);
   return sum;
}


int subtree_count(NODE *node) {
   int sum = node->count;

   for (auto it = node->children.begin(); it != node->children.end(); it++) {
     sum = sum + subtree_count(it->second);
   }
   node->subtree_count = sum;

   return sum;
}


int _count_edge_crossings_(NODE *node, const int tot_count) {
   int sum = 0;
   for (auto it = node->children.begin(); it != node->children.end(); it++) {
     sum = (tot_count - it->second->subtree_count)* it->second->subtree_count + _count_edge_crossings_(it->second, tot_count);
   }

   return sum;
}


void _get_count_vector(NODE *node, vector<int> &count_vector) {

   if (node->count > 0) {
      count_vector.push_back(node->count);
   }

   for (auto it = node->children.begin(); it != node->children.end(); it++) {
      _get_count_vector(it->second, count_vector);
   }
}


double compute_delta_star(NODE *node) {
   int tot_count = node->subtree_count;
   int prod_count = _count_edge_crossings_(node, tot_count);

   vector<int> count_vector;
   _get_count_vector(node, count_vector);

   int sum_product = 0;
   for (uint32_t i = 0;  i < count_vector.size(); i++) {
      for (uint32_t j = i+1;  j < count_vector.size(); j++) {
         sum_product +=  count_vector[i]*count_vector[j];
      }
   }
   std::cout << "vector size " << count_vector.size() << std::endl;
   

   std::cout << "prod_count " << prod_count << "  " << sum_product << std::endl;
   return static_cast<double>(prod_count)/static_cast<double>(sum_product);
} 

void _print_tree_(NODE *node, std::vector<std::string> &taxons) {

   taxons.push_back(node->name);
   if (node->count >= 0) {
   //   std::cout << "\t<" << taxon_concat(taxons) << ">\t" << node->count << std::endl;
     // std::cout <<  taxon <<  std::endl;
   }

   for (auto it = node->children.begin(); it != node->children.end(); it++) {
      _print_tree_(it->second, taxons);
   }
   taxons.pop_back();
}


void print_tree(NODE *node) {
   std::vector<std::string> taxons;

   _print_tree_(node, taxons);
}



NODE *read_tax_file(const string &tax_file_name) {
   fstream newfile;

   regex root_regex_patt("root");
   regex count_regex_patt("(\\([0-9]+\\)$)");
   smatch sm;

   NODE *tree = nullptr;

   vector<std::string> taxstring;
   newfile.open(tax_file_name, ios::in); //open a file to perform read operation using file object
   int k = 0;
   if (newfile.is_open()){   //checking whether the file is open
      string tp;
      while (getline(newfile, tp)){ //read data from file object and put it into string.
       //  cout << tp << "\n"; //print the data of the string
         taxstring = split(tp, '\t');
         if (taxstring.size() >= 10) {
              regex_search(taxstring[8], sm,  root_regex_patt);
              if (sm.size() > 0) { 

                 std::string taxon = regex_replace(taxstring[8], count_regex_patt, "");
                 if(sm.size() >= 0) {
                    /*
                     for(auto x: split(trim(taxon, " "), ';') ) {
                         std::cout << '\t' << "<" <<  x << ">" << std::endl;
                     }
                    */
                   // std::cout << taxon << std::endl;
                    k++;
                    insert_node(&tree, split(trim(taxon, " "),';'), 0);
                 }
              }
         }
      }
      newfile.close(); //close the file object.
   }
   std::cout << "No of nodes inserted " << k << std::endl;
   return tree;
} 

