#include "svm.h"
#include "io.h"
#include "fullproblem.h"
#include "subproblem.h"
#include "kernels.h"
#include "linked.h"

#include <sys/time.h>
#include <stdio.h>
#include <math.h>

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
  struct denseData ds;
  struct Fullproblem fp;
  struct Projected sp;

  struct timeval start, trainStart, trainEnd, end;


  gettimeofday(&start, 0);


  // Input processed:
  parse_arguments(argc, argv, &filename, &parameters);
  read_file(filename, &ds);
  //preprocess(&ds);




  //cleanData(&ds);

  // Projected problem size chosen temporarily
  int p = 6;
  if(parameters.test){
    testSavedModel(&ds, parameters.modelfile, &parameters);
    return 0;
  }

  gettimeofday(&trainStart, 0);
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

  setH(&fp, &ds, &parameters);

  while(k){
    // H matrix columns re-set and subproblem changed
    if(itt%10000 == 0){
      printf("itt = %d\n",itt );
    }

    init_subprob(&sp, &fp, &ds, &parameters, 1);

    //  congjugate gradient algorithm
    //  if algorithm completes n == 0
    //  if algorithm interrupt n != 0
    n = cg(&sp, &fp);


    updateAlphaR(&fp, &sp);
    calcYTR(&sp, &fp);
    calculateBeta(&fp, &sp, &ds);

    if (n==0) {
      printf("Converged! (itt = %d) %d\n", itt, fp.p );

      int add = 2;
      int *temp = malloc(sizeof(int)*add);
      int *temp2 = malloc(sizeof(int)*add);
      add = findWorstest(&fp,add,temp,temp2);
      if (add == 0){
        break;
      }

      changeP(&fp, &sp, add);
      reinitprob(&ds, &fp, &sp, add, temp, temp2);

      free(temp);
      free(temp2);
    }

    if (n) {
      // BCs broken, fix one at a time for the moment
      k = singleswap(&ds, &fp, &sp, n, &parameters);
      if (k < 0) {
        shrinkSize(&fp, &sp, k+fp.p);
      }
      else{
        n = checkfpConstraints(&fp);
      }
    }

    itt++;
    if(itt == max_iters){
      printf("Reached max iters (%d)!!!!!\n\n\n",itt );
      break;
    }

  }
  gettimeofday(&trainEnd, 0);


  for (int i = 0; i < fp.n; i++) {
    printf("%lf\n",fp.alpha[i] );
  }

  if (parameters.save) {
    saveTrainedModel(&fp, &ds, parameters.savename, &parameters);
  }

  gettimeofday(&end, 0);

  long totalelapsed = (end.tv_sec-start.tv_sec)*1000000 + end.tv_usec-start.tv_usec;
  long trainelapsed = (trainEnd.tv_sec-trainStart.tv_sec)*1000000 + trainEnd.tv_usec-trainStart.tv_usec;

  printf("Training Complete\n" );
  printf("Total Time spent: %ld micro seconds\n",totalelapsed );
  printf("Time spent training: %ld micro seconds\n",trainelapsed );

  //Memory freed
  freeDenseData(&ds);
  freeFullproblem(&fp);
  freeSubProblem(&sp);
  return 0;
}
