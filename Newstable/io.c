#include "io.h"

void read_file(char* filename, struct denseData* ds){
  FILE *fp = fopen(filename, "r");
  if (fp == NULL){
    fprintf(stderr, "gert: io.c: read_file() - file %s not found\n",filename );
    exit(1);
  }
  count_entries(fp, ds);

  ds->data1d = (double*)malloc(sizeof(double)*ds->nInstances*ds->nFeatures);
  ds->data = (double**)malloc(sizeof(double*)*ds->nInstances);
  ds->instanceLabels = (char**)malloc(sizeof(char*)*ds->nInstances);
  ds->featureLabels = (char**)malloc(sizeof(char*)*ds->nFeatures);

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
      r++;
    }
    else if (atoi(p) == -1){
      for (int j = 0; j < ds->nFeatures; j++) {
        ds->data[ds->nPos+q][j] = temp[j];
      }
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

  int counter = 0;
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
    counter++;
    int num = atoi(p);
    if (num == 1) {
      ds->nPos++;
    }else if(num == -1){
      ds->nNeg++;
    }else{
      fprintf(stderr, "invalid classes (should be 1 or -1, %d found %d)\n",num,counter);
      exit(1);
    }
    ds->nInstances++;
  }

  rewind(input);
}

void saveTrainedModel(struct Fullproblem *fp, struct denseData *ds, char *filename, struct svm_args *params)
{
  FILE *file = fopen(filename, "w");
  fprintf(file, "%d\n",params->kernel );

  fprintf(file, "%d\n",ds->nFeatures );

  int count = 0;
  for (int i = 0; i < fp->n; i++) {
    if (fp->alpha[i] > 0.0) {
      count++;
    }
  }
  int *active = malloc(sizeof(int)*count);
  int j = 0;
  for (int i = 0; i < fp->n; i++) {
    if (fp->alpha[i] > 0.0) {
      active[j] = i;
      j++;
    }
  }
  double *h = malloc(sizeof(double)*fp->n*count);
  double **H = malloc(sizeof(double*)*fp->n);
  int r = -1;

  for (int i = 0; i < fp->n; i++) {
    H[i] = &h[i*count];
  }

  if (params->kernel == LINEAR) {
    for (int i = 0; i < fp->n; i++) {
      for (int j = 0; j < count; j++) {
        H[i][j] = 0.0;
        for (int k = 0; k < ds->nFeatures; k++) {
          H[i][j] += ds->data[i][k]*ds->data[active[j]][k];
        }
				if( (active[j] < ds->nPos) ^ (i < ds->nPos)  ){
        	H[i][j] = -H[i][j];
      	}
    	}
		}

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
			printf("b is %lf\n",b);
      b -= H[r][i]*fp->alpha[active[i]];
    }
			printf("b is %lf\n",b);
		if(r >= ds->nPos){
	    b = -b;
		}
    double *w = malloc(sizeof(double)*ds->nFeatures);

    for (int i = 0; i < ds->nFeatures; i++) {
      w[i] = 0.0;
      for (int j = 0; j < count; j++) {
				if(active[j] < ds->nPos){
	        w[i] += fp->alpha[active[j]]*ds->data[active[j]][i];
				}
				else{
	        w[i] -= fp->alpha[active[j]]*ds->data[active[j]][i];
				}
      }
    }



    fprintf(file, "%lf\n",b );
    for (int i = 0; i < ds->nFeatures; i++) {
      fprintf(file, "%lf\n",w[i] );
    }
    free(w);

  }
  else if(params->kernel == POLYNOMIAL){
    fprintf(file, "%d\n", count );
    for (int i = 0; i < fp->n; i++) {
      for (int j = 0; j < count; j++) {
        H[i][j] = 0.0;
        for (int k = 0; k < ds->nFeatures; k++) {
          H[i][j] += ds->data[i][k]*ds->data[active[j]][k];
        }
        H[i][j] = pow(H[i][j]+params->Gamma, params->degree);
        if(   (active[j] < ds->nPos )  ^  ( i < ds->nPos)  ){
        	H[i][j] = -H[i][j];
      	}
      }
    }
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
		if(r>= ds->nPos){
	    b = -b;
		}  
	  fprintf(file, "%lf\n",b );

    for (int i = 0; i < count; i++) {
      for (int j = 0; j < ds->nFeatures; j++) {
        fprintf(file, "%lf\n",ds->data[active[i]][j] );
      }
    }

    for (int i = 0; i < count; i++) {
			if(active[i] < ds->nPos){
	      fprintf(file, "%lf\n",fp->alpha[active[i]] );
  	  }
			else{
	      fprintf(file, "%lf\n",-fp->alpha[active[i]] );
			}
  	}
	}
  else if(params->kernel == EXPONENTIAL)
  {
    fprintf(file, "%d\n", count );
    double x,y;
    for (int i = 0; i < fp->n; i++) {
      for (int j = 0; j < count; j++) {
        y = 0.0;
        for (int k = 0; k < ds->nFeatures; k++) {
          x = ds->data[i][k] - ds->data[active[j]][k];
          y -= x*x;
        }
        y *= params->Gamma;
        H[i][j] = exp(y);
				if(  (active[j] < ds->nPos  )   ^   (i < ds->nPos)  ){
        	H[i][j] = -H[i][j];
      	}
      }
    }
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
		if(r>=ds->nPos){
	    b = -b;
		}
    fprintf(file, "%lf\n",b );

    for (int i = 0; i < count; i++) {
      for (int j = 0; j < ds->nFeatures; j++) {
        fprintf(file, "%lf\n",ds->data[active[i]][j] );
      }
    }

    for (int i = 0; i < count; i++) {
			if(active[i] < ds->nPos){
	      fprintf(file, "%lf\n",fp->alpha[active[i]] );
  	  }
			else{
	      fprintf(file, "%lf\n",-fp->alpha[active[i]] );
  	  }
  	}
	}
  free(active);
  free(h);
  free(H);
  fclose(file);
}

void testSavedModel(struct denseData *ds, char* fn, struct svm_args *params)
{
  FILE *fp = fopen(fn, "r");
  int k;
  double b;
  int kernel;

  int res = fscanf(fp, "%d",&kernel);
  res = fscanf(fp, "%d",&k);

  if (k!= ds->nFeatures) {
    fprintf(stderr, "io.c: \n" );
    exit(1);
  }
  int wrong = 0;

  if (kernel == LINEAR) {
    res = fscanf(fp, "%lf", &b);

    double *w = malloc(sizeof(double)*k);
    for (size_t i = 0; i < k; i++) {
      res = fscanf(fp, "%lf",&w[i]);
    }
    double value;
    for (int i = 0; i < ds->nInstances; i++) {
      value = b;
      for (int j = 0; j < ds->nFeatures; j++) {
        value += w[j]*ds->data[i][j];
      }
			if(i<ds->nPos){
	      if (value < 0.0) {
  	      printf("%sres[%d] = %.3lf%s\n",RED,i,value,RESET );
  	      wrong++;
  	    }
  	    else{
  	      printf("%sres[%d] = %.3lf%s\n",GRN,i,value,RESET );
  	    }
    	}
			else{
	      if (value > 0.0) {
  	      printf("%sres[%d] = %.3lf%s\n",RED,i,value,RESET );
  	      wrong++;
  	    }
  	    else{
  	      printf("%sres[%d] = %.3lf%s\n",GRN,i,value,RESET );
  	    }
    	}
		}
    free(w);
  }
  else if(params->kernel == POLYNOMIAL)  {
    int count;
    res = fscanf(fp, "%d", &count);
    res = fscanf(fp, "%lf", &b);
    double value;
    double *alphaY = malloc(sizeof(double)*count);
    double *x = malloc(sizeof(double)*count*k);
    double **X = malloc(sizeof(double*)*count);
    for (int i = 0; i < count; i++) {
      X[i] = &x[i*k];
    }
    for (int i = 0; i < count; i++) {
      for (int j = 0; j < k; j++) {
        res = fscanf(fp, "%lf", &X[i][j]);
      }
    }
    for (int i = 0; i < count; i++) {
      res = fscanf(fp, "%lf", &alphaY[i]);
    }
    double contrib = 0;
    for (int i = 0; i < ds->nInstances; i++) {
      value = b;
      for (int j = 0; j < count; j++) {
        contrib = 0;
        for (int k = 0; k < ds->nFeatures; k++) {
          contrib += X[j][k]*ds->data[i][k];
        }
        contrib = pow(contrib + params->Gamma, params->degree);
        value += alphaY[j]*contrib;
      }
			if(i<ds->nPos){
	      if (value < 0.0) {
  	      printf("%sres[%d] = %.3lf%s\n",RED,i,value,RESET );
  	      wrong++;
  	    }
  	    else{
  	      printf("%sres[%d] = %.3lf%s\n",GRN,i,value,RESET );
  	    }
    	}
			else{
	      if (value > 0.0) {
  	      printf("%sres[%d] = %.3lf%s\n",RED,i,value,RESET );
  	      wrong++;
  	    }
  	    else{
  	      printf("%sres[%d] = %.3lf%s\n",GRN,i,value,RESET );
  	    }
    	}
    }
  }
  else if(params->kernel == EXPONENTIAL)  {
    int count;
    res = fscanf(fp, "%d", &count);
    res = fscanf(fp, "%lf", &b);
    double value;
    double *alphaY = malloc(sizeof(double)*count);
    double *x = malloc(sizeof(double)*count*k);
    double **X = malloc(sizeof(double*)*count);
    for (int i = 0; i < count; i++) {
      X[i] = &x[i*k];
    }
    for (int i = 0; i < count; i++) {
      for (int j = 0; j < k; j++) {
        res = fscanf(fp, "%lf", &X[i][j]);
      }
    }
    for (int i = 0; i < count; i++) {
      res = fscanf(fp, "%lf", &alphaY[i]);
    }
    double contrib = 0;
    double y;
    for (int i = 0; i < ds->nInstances; i++) {
      value = b;
      for (int j = 0; j < count; j++) {
        contrib = 0;
        for (int k = 0; k < ds->nFeatures; k++) {
          y = X[j][k] - ds->data[i][k];
          contrib -= y*y;
        }
        contrib *= params->Gamma;
        value += alphaY[j]*exp(contrib);
      }
			if(i<ds->nPos){
	      if (value < 0.0) {
  	      printf("%sres[%d] = %.3lf%s\n",RED,i,value,RESET );
  	      wrong++;
  	    }
  	    else{
  	      printf("%sres[%d] = %.3lf%s\n",GRN,i,value,RESET );
  	    }
    	}
			else{
	      if (value > 0.0) {
  	      printf("%sres[%d] = %.3lf%s\n",RED,i,value,RESET );
  	      wrong++;
  	    }
  	    else{
  	      printf("%sres[%d] = %.3lf%s\n",GRN,i,value,RESET );
  	    }
    	}
    }
  }

  int right = ds->nInstances - wrong;
  double pct = 100.0*((double)right/(double)ds->nInstances);
  printf("%d correct classifications out of %d. %.2lf%% correct.\n",right,ds->nInstances,pct );

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
  parameters->kernel = 0;
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
      case 'k':
        parameters->kernel = atoi(optarg);
        break;
      case 'c':
        parameters->C = atof(optarg);
        break;
      case 't':
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
  free(ds->instanceLabels);
  free(ds->featureLabels);
}
