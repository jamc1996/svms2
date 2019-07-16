#include "io.h"

void read_file(char* const filename, struct denseData* ds){
  FILE *fp = fopen(filename, "r");
  count_entries(fp, ds);
  ds->data1d = (double*)malloc(sizeof(double)*ds->nInstances*ds->nFeatures);
  ds->data = (double**)malloc(sizeof(double*)*ds->nInstances);
  ds->instanceLabels = (char**)malloc(sizeof(char*)*ds->nInstances);
  ds->featureLabels = (char**)malloc(sizeof(char*)*ds->nFeatures);
  ds->y = (double*)malloc(sizeof(double)*ds->nInstances);

  char* line = NULL;
  char* endptr;
  if (false) {
    for (int i = 0; i < ds->nFeatures; i++) {
      ds->featureLabels[i] = strtok(NULL, " \t");
  //    std::cout << ds->featureLabels[i] << '\n';
    }
  }


  for (int i = 0; i < ds->nInstances; i++) {
    ds->data[i] = &ds->data1d[i*ds->nFeatures];
    readline(fp,&line);
    ds->instanceLabels[i] = strtok(line, " \t");
    if (ds->instanceLabels[i] == NULL || *(ds->instanceLabels[i]) == '\n') {
      fprintf(stderr, "main.cpp: read_file(): bad read at %d\n",i );
      exit(1);
    }
    //std::cout << ds->instanceLabels[i] << '\n';

    for (int j = 0; j < ds->nFeatures; j++) {
      char* p = strtok(NULL, " \t");
      if (p == NULL || *p == '\n') {
        fprintf(stderr, "Oh dear\n" );
        exit(1);
      }
      ds->data[i][j] = strtod(p, &endptr);
      //std::cout << "\t" << ds->data[i][j];
    }

    char* p = strtok(NULL, " \t");
    ds->y[i] = strtod(p, &endptr);
    //std::cout << '\n';
  }
}

void count_entries(FILE *input, struct denseData* ds)
{
  ds->nInstances = 0;
  ds->nFeatures = -1;
  char* line = NULL;


  // Find size of dataset:
  while (readline(input, &line)) {
    if (ds->nInstances==0) {
      char *p = strtok(line," \t");
      while (true) {
        p  = strtok(NULL, " \t");
        if (p == NULL || *p == '\n') {
          break;
        }
        ds->nFeatures++;
      }
    }
    ds->nInstances++;
  }
  rewind(input);
}

int readline(FILE *input, char **line)
/* Function to read lines from file */
{
  int len;
  int max_line_len = 1024;
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
  return 0;
}
