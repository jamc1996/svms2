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

  double* bigH = malloc(sizeof(double)*ds.nInstances*ds.nInstances);
  double** fullBigH = malloc(sizeof(double*)*ds.nInstances);
  for (int i = 0; i < ds.nInstances; i++) {
    fullBigH[i] = &bigH[i*ds.nInstances];
    for (int j = 0; j < ds.nInstances; j++) {
      fullBigH[i][j] = 0.0;
      for (int k = 0; k < ds.nFeatures; k++) {
        fullBigH[i][j] += ds.data[i][k]*ds.data[j][k];
      }
      fullBigH[i][j] *= ds.y[j]*ds.y[i];
    }
  }
  double* check = malloc(sizeof(double)*ds.nInstances);
  //cleanData(&ds);

  // Projected problem size chosen temporarily
  int p = 4;
  char* fname = "TrainedStates.txt";

  if(parameters.test){
    testSavedModel(&ds, fname);
    return 0;
  }

  // Full problem allocated and filled in, all alpha = 0.0 all gradF = 1.0:
  alloc_prob(&fp, &ds, p);
  init_prob(&fp, &ds);



  // Subproblem allocated:
  alloc_subprob(&sp, p, &fp, &ds);


  // We loop until no negative entries in beta:
  int k = 1;
  int max_iters = 1000;
  int itt = 0;
  int n = 0;

  while(k){
  //  break;
    // H matrix columns re-set and subproblem changed
    printf("\n\nitt = %d\n\n\n",itt );
    setH(&fp, &ds);
    init_subprob(&sp, &fp, &ds);



    //  congjugate gradient algorithm
    //  n = 0 if algorithm completes
    //  n != 0 if algorithm interrupt
    n = cg(&sp, &fp);
    printf("n is %d\n",n );


    updateAlphaR(&fp, &sp);


    for (int i = 0; i < fp.n; i++) {
      check[i] = 1.0;
      for (int j = 0; j < fp.n; j++) {
        check[i] -= fullBigH[i][j]*fp.alpha[j];
      }
    }
    for (int i = 0; i < fp.n; i++) {
      printf("check is %lf\n",fp.gradF[i]-check[i] );
    }
    for (int i = 0; i < fp.n; i++) {
      printf("alpha[%d] == %lf (%lf)\n",i,fp.alpha[i], fp.gradF[i]-check[i] );
    }
    for (int i = 0; i < fp.p; i++) {
      //printf("active is %d\n",fp.active[i] );
    }
    for (int i = 0; i < fp.q; i++) {
      //printf("INactive is %d\n",fp.inactive[i] );
    }

    for (int i = 0; i < fp.n; i++) {
      printf("check is %lf\n",fp.gradF[i]-check[i] );
      if (fabs(fp.gradF[i] - check[i]) > 0.00001 ) {
        printf("here\n" );
        exit(22);
      }
    }

    calcYTR(&sp, &fp);
    printf("ytr %lf\n",sp.ytr/(double)(fp.p) );
    double b = 1.0;
    for (int i = 0; i < fp.p; i++) {
      printf("h is %lf and alpha is %lf b i s %lf\n",sp.H[0][i],fp.alpha[fp.active[i]],b );
      b -= sp.H[0][i]*fp.alpha[fp.active[i]];
    }
    printf("gradF is %lf %lf %lf\n",fp.gradF[fp.active[0]],fp.gradF[fp.active[1]],fp.gradF[fp.active[2]] );
    b*=ds.y[fp.active[0]];
    printf("b = %lf\n",b );
  //  exit(0);
    calculateBeta(&fp, &sp, &ds);
    for (int i = 0; i < fp.q; i++) {
      if (fp.beta[i] < 0.0) {
        printf("beta[%d](%d) = %lf\n",i,fp.inactive[i],fp.beta[i] );
      }
    }
    if (n==0) {
      printf("Converged!\n" );
      for (int i = 0; i < sp.p; i++) {
        printf("%lf\n",fp.gradF[fp.active[i]] );
      }
      printf("ytr %lf\n",sp.ytr/(double)(fp.p) );

      int add = 2;
      int *temp = malloc(sizeof(int)*add);
      int *temp2 = malloc(sizeof(int)*add);
      add = findWorstest(&fp,add,temp,temp2);
      printf("adding %d\n",add );
      if (add == 0){
        printf("Add == 0\n" );
        break;
      }
      changeP(&fp, &sp, add);
      reinitprob(&fp, &sp, add, temp, temp2);

      free(temp);
      free(temp2);
    }
    int prqe = 0;

    while (n) {
      // While BCs broken
      k = singleswap(&ds, &fp, &sp, n);
      printf("k is %d\n",k );
      if (n >= fp.p) {
        printf("Exceeding %d with %d\n",fp.p,n );
        //exit(5);
      }
      for (int i = 0; i < fp.q; i++) {
      //  printf("%lf\n",fp.beta[i] );
      }
      if (k < 0) {
        printf("shrinking\n" );
        shrinkSize(&fp, &sp, k+fp.p);
        break;
      }
      n = checkfpConstraints(&fp);
      printf("now n is %d\n", n);
      prqe++;
      break;
      if (prqe == fp.n) {
        return 1;
      }
    }

    itt++;
    if(itt == max_iters){
      printf("Reached max iters (%d)!!!!!\n\n\n",itt );
      break;
    }
  }

  saveTrainedModel(&fp, &ds, fname);
  testSavedModel(&ds, fname);


  //Memory freed
  freeDenseData(&ds);
  freeFullproblem(&fp);
  freeSubProblem(&sp);
  return 0;
}
