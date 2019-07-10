#include "svm.h"
#include "io.h"
#include "fullproblem.h"
#include "subproblem.h"

#include <stdio.h>

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

int main(int argc, char *argv[]) {
  char* filename = NULL;
  struct svm_args parameters;
  struct denseData ds;
  struct Fullproblem fp;
  struct Projected sp;

  // Input processed:
  parse_arguments(argc, argv, &filename, &parameters);
  read_file(filename, &ds);
  //preprocess(&ds);

  //double* bigH = malloc(sizeof(double)*ds.nInstances*ds.nInstances);
  //double** fullBigH = malloc(sizeof(double*)*ds.nInstances);
  for (int i = 0; i < ds.nInstances; i++) {
    //fullBigH[i] = &bigH[i*ds.nInstances];
    for (int j = 0; j < ds.nInstances; j++) {
      //fullBigH[i][j] = 0.0;
      for (int k = 0; k < ds.nFeatures; k++) {
        //fullBigH[i][j] += ds.data[i][k]*ds.data[j][k];
      }
      //fullBigH[i][j] *= ds.y[j]*ds.y[i];
    }
  }
  double* check = malloc(sizeof(double)*ds.nInstances);
  //cleanData(&ds);

  // Projected problem size chosen temporarily
  int p = 4;

  if(parameters.test){
    testSavedModel(&ds, parameters.modelfile);
    return 0;
  }

  // Full problem allocated and filled in, all alpha = 0.0 all gradF = 1.0:
  alloc_prob(&fp, &ds, p);
  init_prob(&fp, &ds);



  // Subproblem allocated:
  alloc_subprob(&sp, p);


  // We loop until no negative entries in beta:
  int k = 1;
  int max_iters = 10000000;
  int itt = 0;
  int n = 0;

  while(k){
  //  break;
    // H matrix columns re-set and subproblem changed
    //printf("\n\nitt = %d\n\n\n",itt );
    if(itt%100000 == 0){
      printf("itt = %d\n",itt );
    }
    setH(&fp, &ds);
    init_subprob(&sp, &fp, &ds);

    //  congjugate gradient algorithm
    //  n = 0 if algorithm completes
    //  n != 0 if algorithm interrupt
    n = cg(&sp, &fp);


    updateAlphaR(&fp, &sp);
    calcYTR(&sp, &fp);
    calculateBeta(&fp, &sp, &ds);

    for (int i = 0; i < fp.n; i++) {
      //check[i] = 1.0;
      for (int j = 0; j < fp.n; j++) {
        //check[i] -= fullBigH[i][j]*fp.alpha[j];
      }
    }
    for (int i = 0; i < fp.n; i++) {
      if (fabs(fp.gradF[i] - check[i]) > 0.00001 ) {
        //printf("n is %d/%d\n",n,fp.p );
        //printf("here %lf\n",fabs(fp.gradF[i] - check[i]) );
        //exit(22);
      }
    }



    if (n==0) {
      printf("Converged!\n" );

      int add = 2;
      int *temp = malloc(sizeof(int)*add);
      int *temp2 = malloc(sizeof(int)*add);
      add = findWorstest(&fp,add,temp,temp2);
      if (add == 0){
        break;
      }
      changeP(&fp, &sp, add);
      reinitprob(&fp, &sp, add, temp, temp2);

      free(temp);
      free(temp2);
    }

    if (n) {
      // BCs broken, fix one at a time for the moment
      k = singleswap(&ds, &fp, &sp, n);
      if (k < 0) {
        shrinkSize(&fp, &sp, k+fp.p);
      }
      n = checkfpConstraints(&fp);
    }

    itt++;
    if(itt == max_iters){
      printf("Reached max iters (%d)!!!!!\n\n\n",itt );
      break;
    }
  }



  if (parameters.save) {
    saveTrainedModel(&fp, &ds, parameters.savename);
  }

  //Memory freed
  freeDenseData(&ds);
  freeFullproblem(&fp);
  freeSubProblem(&sp);
  return 0;
}
