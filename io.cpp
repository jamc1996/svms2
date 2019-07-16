#include "io.h"

static int max_line_len;

int readline(FILE *input, char **line)
/* Function to read lines from file */
{
  int len;
  max_line_len = 1024;
  *line = (char*)realloc(*line,sizeof(char)*max_line_len);
  if(fgets(*line,max_line_len,input) == NULL)
  {
    return 0;
  }

  while(strrchr(*line,'\n') == NULL)
  {
    max_line_len *= 2;
    *line = (char *) realloc(*line, max_line_len);
    len = (int) strlen(*line);
    if (fgets(*line+len,max_line_len-len,input) == NULL) {
      break;
    }
  }
  return 1;
}


int parse_arguments(int argc, char *argv[], char** filename, struct svm_args *parameters)
/* Function to parse command line arguments with getopt */
{
  int c;

  // Default values set:
  parameters->type = 0;
  parameters->kernel = 1;
  parameters->degree = 1;
  parameters->verbose = false;
  parameters->C = 1;

  while ((c = getopt( argc, argv, "f:t:c:d:vh")) != -1){
    switch (c) {
      case 'f':
        *filename = optarg;
        break;
      case 't':
        parameters->type = atoi(optarg);
        break;
      case 'c':
        parameters->C = atoi(optarg);
        break;
      case 'd':
        parameters->degree = atoi(optarg);
        break;
      case 'v':
        parameters->verbose = true;
        break;
      case 'h':
        std::cout << "I was supposed to put in a help message here." << '\n';
        break;
    }
  }

  if (*filename == NULL) {
    printf("io.cpp: parse_arguments: no input file selected.\n");
    exit(1);
  }
  std::cout << "all ok at end of parsing" << '\n';
  return 0;
}
