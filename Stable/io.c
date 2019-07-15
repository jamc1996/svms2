#include "io.h"

void read_file(char* filename, struct denseData* ds){
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

void saveTrainedModel(struct Fullproblem *fp, struct denseData *ds, char *filename){
  FILE *file = fopen(filename, "w");
  fprintf(file, "%d\n",ds->nFeatures );

  int count = 0;
  for (int i = 0; i < fp->n; i++) {
    if (fp->alpha[i] > 0.0) {
      count++;
    }
  }
  int * active = malloc(sizeof(int)*count);
  int j = 0;
  for (int i = 0; i < fp->n; i++) {
    if (fp->alpha[i] > 0.0) {
      active[j] = i;
      j++;
    }
  }
  double *h = malloc(sizeof(double)*fp->n*count);
  double **H = malloc(sizeof(double*)*fp->n);
  for (int i = 0; i < fp->n; i++) {
    H[i] = &h[i*count];
  }

  for (int i = 0; i < fp->n; i++) {
    for (int j = 0; j < count; j++) {
      H[i][j] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        H[i][j] += ds->data[i][k]*ds->data[active[j]][k];
      }
      H[i][j] *= ds->y[active[j]]*ds->y[i];
    }
  }
  int r = -1;
  for (int i = 0; i < count; i++) {
    if (fp->alpha[active[i]] < fp->C*0.99) {
      r = active[i];
      break;
    }
  }
  if (r<0) {
    exit(77);
  }

  double b = 1.0;
  for (int i = 0; i < count; i++) {
    b -= H[r][i]*fp->alpha[active[i]];
  }
  b *= ds->y[r];

  double *w = malloc(sizeof(double)*ds->nFeatures);
  for (int i = 0; i < ds->nFeatures; i++) {
    w[i] = 0.0;
    for (int j = 0; j < count; j++) {
      w[i] += fp->alpha[active[j]]*ds->data[active[j]][i]*ds->y[active[j]];
    }
  }



  fprintf(file, "%lf\n",b );
  for (int i = 0; i < ds->nFeatures; i++) {
    fprintf(file, "%lf\n",w[i] );
  }
  free(active);
  free(h);
  free(H);
  free(w);
  fclose(file);
}

void testSavedModel(struct denseData *ds, char* fn)
{
  FILE *fp = fopen(fn, "r");
  int k;
  double b;
  int res = fscanf(fp, "%d",&k);
  if (res == EOF) {
    exit(2222);
  }
  res = fscanf(fp, "%lf", &b);
  double *w = malloc(sizeof(double)*k);
  for (size_t i = 0; i < k; i++) {
    res = fscanf(fp, "%lf",&w[i]);
  }
  int wrong = 0;
  for (int i = 0; i < ds->nInstances; i++) {
    double res = b;
    for (int j = 0; j < ds->nFeatures; j++) {
      res += w[j]*ds->data[i][j];
    }
    if (ds->y[i]*res < 0.0) {
      printf("%sres[%d] = %.3lf%s\n",RED,i,res,RESET );
      wrong++;
    }
    else{
      printf("%sres[%d] = %.3lf%s\n",GRN,i,res,RESET );
    }
  }
  int right = ds->nInstances - wrong;
  double pct = 100.0*((double)right/(double)ds->nInstances);
  printf("%d correct classifications out of %d. %.2lf%% correct.\n",right,ds->nInstances,pct );

  free(w);
  fclose(fp);

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
  parameters->modelfile = NULL;
  parameters->save = 0;
  parameters->savename = NULL;
  parameters->Gamma = 1;
  while ((c = getopt( argc, argv, "f:k:t:c:d:vhs:g:")) != -1){
    switch (c) {
      case 'f':
        *filename = optarg;
        break;
      case 't':
        parameters->type = atoi(optarg);
        break;
      case 'c':
        parameters->C = atof(optarg);
        break;
      case 'k':
        parameters->test = 1;
        parameters->modelfile = optarg;
        break;
      case 'd':
        parameters->degree = atoi(optarg);
        break;
      case 's':
        parameters->save = 1;
        parameters->savename = optarg;
        break;
      case 'v':
        parameters->verbose = 1;
        break;
      case 'g':
        parameters->Gamma = atof(optarg);
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


void cleanData( struct denseData *ds){
  for (int i = 0; i < ds->nInstances - 1; i++) {
    printf("%d\n",i );
    for (int j = i; j < ds->nInstances; j++) {
      int flag = 1;
      double check;
      double previous = ds->data[i][0]/ds->data[j][0];
      for (int k = 1; k < ds->nFeatures; k++) {
        check = ds->data[i][k]/ds->data[j][k];
        printf("%lf and %lf\n",previous,check );
        if (fabs(check - previous) > 0.0001 ) {
          flag = 0;
          break;
        }
        previous = check;
      }
      if (flag == 1) {
        printf("%d is broken\n",j );
      }
    }
  }
}


void freeDenseData(struct denseData *ds)
/* Function to free dynamically allocated memory in dense data set struct. */
{
  free(ds->data);
  free(ds->data1d);
  free(ds->y);
  free(ds->instanceLabels);
  free(ds->featureLabels);
}
