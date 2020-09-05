#include "utilities.h"

using namespace std;

std::string trim(const std::string& str, const std::string& whitespace = " \t") {
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}


vector<string> split(const string& input, const char delim) {
    vector<std::string> result;
    std::istringstream ss(input);
    std::string token;

    while(std::getline(ss, token, delim)) {
       result.push_back(token);
    }
    return result;
}


std::string taxon_concat(const std::vector<std::string> &taxons) {
    std::string taxon;
    for(auto it = taxons.begin(); it != taxons.end(); it++) {
        taxon += *it + ";";
    }
    return taxon;
}


void read_options(int argc, char **argv, INPUT_OPTIONS &options) {
  int c;
  int i;

  int verbose_flag ;
  while (true) {
      static struct option long_options[] =
        {
          /* These options set a flag. */
          {"verbose", no_argument,  &verbose_flag, false},
          /* These options don’t set a flag.
             We distinguish them by their indices. */
          {"tax-file",    required_argument, 0, 'f'},
          {0, 0, 0, 0}
        };

      const char *help_messages[] = {
           "verbose messages", 
           "taxonomic file" 
      };


      /* getopt_long stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "f:", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c)
        {
        case 0:
          /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf ("option %s", long_options[option_index].name);
            if (optarg)
                printf (" with arg %s", optarg);
            printf ("\n");
            break;
        case 'f':
            options.tax_file_name = std::string(optarg); 
            break;
        case '?':
        case 'h':
          i = 0;
          printf("Usage: %s [options] \n", argv[0]);
          while(long_options[i].name != 0) {
             printf("\t--%-20s  %-25s  %-35s\n", long_options[i].name, \
                     long_options[i].has_arg==no_argument? "no argument" : "required_argument", \
                     help_messages[i]
                   );
             i = i + 1;
          }
          /* getopt_long already printed an error message. */
          return;
          break;
        default:
          abort ();
        }
    }

  if( options.tax_file_name == std::string("") )  {
     printf("Must specify at least one taxonomic file\n");
     exit(0);
  }

  if (verbose_flag) {
       std::cout << "Tax files" <<  options.tax_file_name <<  std::endl;
  }

}
