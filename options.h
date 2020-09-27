#ifndef __OPTIONS__
#define __OPTIONS__
#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>


typedef struct _input_options {
  std::string tax_file_name;
  std::string pwy_file_name;
} INPUT_OPTIONS;

void read_options(int, char **, INPUT_OPTIONS &); 

#endif
