#include "io.h"

void read_file(char* const filename, struct denseData* ds){
  FILE *fp = fopen(filename, "r");
  count_entries(fp, ds);

  ds->data1d = (double*)malloc(sizeof(double)*ds->nInstances*ds->nFeatures);
  ds->data = (double**)malloc(sizeof(double*)*ds->nInstances);
  ds->instanceLabels = (char**)malloc(sizeof(char*)*ds->nInstances);
  ds->featureLabels = (char**)malloc(sizeof(char*)*ds->nFeatures);
  ds->y = (double*)malloc(sizeof(double)*ds->nInstances);
  double* temp= (double*)malloc(sizeof(double)*ds->nFeatures);

  char* line = NULL;
  char* endptr;
  if (0) {
    for (int i = 0; i < ds->nFeatures; i++) {
      ds->featureLabels[i] = strtok(NULL, " \t");
    }
  }
  for (int i = 0; i < ds->nInstances; i++) {
    ds->data[i] = &ds->data1d[i*ds->nFeatures];
  }
  int r = 0;
  int q = 0;
  for (int i = 0; i < ds->nInstances; i++) {
    readline(fp,&line);
    ds->instanceLabels[i] = strtok(line, " \t");
    if (ds->instanceLabels[i] == NULL || *(ds->instanceLabels[i]) == '\n') {
      fprintf(stderr, "main.cpp: read_file(): bad read at %d\n",i );
      exit(1);
    }

    for (int j = 0; j < ds->nFeatures; j++) {
      char* p = strtok(NULL, " \t");
      if (p == NULL || *p == '\n') {
        fprintf(stderr, "Oh dear\n" );
        exit(1);
      }
      temp[j] = strtod(p, &endptr);
//      ds->data[i][j] = strtod(p, &endptr);
    }
    char* p = strtok(NULL, " \t");
    if (atoi(p) == 1) {
      for (int j = 0; j < ds->nFeatures; j++) {
        ds->data[r][j] = temp[j];
      }
      ds->y[r] = 1;

      r++;
    }else if (atoi(p) == -1){
      for (int j = 0; j < ds->nFeatures; j++) {
        ds->data[ds->nPos+q][j] = temp[j];
      }
      ds->y[ds->nPos+q] = -1;

      q++;

    }
  }

  free(temp);
}

void count_entries(FILE *input, struct denseData* ds)
{
  ds->nInstances = 0;
  ds->nFeatures = -1;
  ds->nPos = 0;
  ds->nNeg = 0;
  char* line = NULL;


  // Find size of dataset:
  while (readline(input, &line)) {
    if (ds->nFeatures==-1) {
      char *p = strtok(line," \t");
      while (1) {
        p  = strtok(NULL, " \t");
        if (p == NULL || *p == '\n') {
          break;
        }
        ds->nFeatures++;
      }
      rewind(input);
      continue;
    }
    char* p = strtok(line," \t");
    for (int i = 0; i < ds->nFeatures+1; i++) {
      p = strtok(NULL, " \t");
    }
    int num = atoi(p);
    if (num == 1) {
      ds->nPos++;
    }else if(num == -1){
      ds->nNeg++;
    }else{
      fprintf(stderr, "invalid classes (should be 1 or -1)\n");
      exit(1);
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
  parameters->verbose = 0;
  parameters->C = 1;
  parameters->test = 0;

  while ((c = getopt( argc, argv, "f:kt:c:d:vh")) != -1){
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
      case 'k':
        parameters->test = 1;
        break;
      case 'd':
        parameters->degree = atoi(optarg);
        break;
      case 'v':
        parameters->verbose = 1;
        break;
      case 'h':
        printf("I was supposed to put in a help message here.\n");
        break;
    }
  }

  if (*filename == NULL) {
    printf("io.cpp: parse_arguments: no input file selected.\n");
    exit(1);
  }
  return 0;
}

void preprocess(struct denseData *ds)
{
  double* means = (double*)calloc(ds->nFeatures,sizeof(double));
  double* stdDev = (double*)calloc(ds->nFeatures,sizeof(double));

  calcMeans(means, ds);
  calcStdDev(stdDev,means,ds);
  normalise(means,stdDev,ds);
  free(means);
  free(stdDev);
}


void calcMeans(double *mean, struct denseData *ds)
{
  for (int i = 0; i < ds->nInstances; i++) {
    for (int j = 0; j < ds->nFeatures; j++) {
      mean[j] += ds->data[i][j];
    }
  }
  for (int i = 0; i < ds->nFeatures; i++) {
    mean[i]/=(double)(ds->nInstances);
  }
}

void normalise(double* mean, double* stdDev, struct denseData* ds)
{
  for (int i = 0; i < ds->nInstances; i++) {
    for (int j = 0; j < ds->nFeatures; j++) {
      ds->data[i][j]-=mean[j];
      ds->data[i][j]/=stdDev[j];
    }
  }
}

void calcStdDev(double* stdDev, double* mean, struct denseData *ds)
{
  for (int i = 0; i < ds->nInstances; i++) {
    for (int j = 0; j < ds->nFeatures; j++) {
      stdDev[j]+=(ds->data[i][j]-mean[j])*(ds->data[i][j]-mean[j]);
    }
  }
  for (int i = 0; i < ds->nFeatures; i++) {
    stdDev[i] = sqrt(stdDev[i]/((double)(ds->nInstances)-1.0));
  }
}
