#include "options.h"
using namespace std;

/**
  @brief read_options Reads the command line options

  @param argc same argc as main, i.e., argument counts
  @param argv same argv as main, i.e., the array of character arrays
  @param options intput options structure 
*/
void read_options(int argc, char **argv, INPUT_OPTIONS &options) {
  int c;
  int i;

  int verbose_flag ;
  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"verbose", no_argument,  &verbose_flag, false},
   /* These options don’t set a flag.
   We distinguish them by their indices. */
    {"tax-file",    required_argument, 0, 'f'},
    {"pwy-file",    required_argument, 0, 'p'},
    {0, 0, 0, 0}
  };

  const char *help_messages[] = {
    "verbose messages", 
    "taxonomic file",
    "pathway file" 
  };
  int num_flags = 3;


  while (true) {
      /* getopt_long stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "f:p:h", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c) {
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0) {
            printf ("option %s", long_options[option_index].name);
            if (optarg) printf (" with arg %s", optarg);
            printf ("\n");
           }
           break;
        case 'f':
          options.tax_file_name = std::string(optarg); 
          break;
        case 'p':
          options.pwy_file_name = std::string(optarg); 
          break;
        case '?':
        case 'h':
          i = 0;
          printf("Usage: %s [options] \n", argv[0]);
          for (i = 0; i < num_flags; i++) {
            printf("\t--%-20s  %-25s  %-35s\n", long_options[i].name, 
              long_options[i].has_arg==no_argument? "no argument" : "required_argument", 
              help_messages[i]
            );
          }
          exit(0);
          /* getopt_long already printed an error message. */
          return;
          break;
        default:
          abort();
        }
    }

  if (options.tax_file_name == std::string(""))  {
     printf("Must specify at least one taxonomic file\n");
     exit(0);
  }

  if (verbose_flag) {
       std::cout << "Tax files" <<  options.tax_file_name <<  std::endl;
  }
}
