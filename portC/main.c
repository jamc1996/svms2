#include "svm.h"
#include "io.h"
#include "fullproblem.h"
#include "subproblem.h"

#include <stdio.h>

void freeDenseData(struct denseData *ds);
void  freeFullproblem(struct Fullproblem *fp);
void   freeSubProblem( struct Projected* sp);

int main(int argc, char *argv[]) {
  char* filename = NULL;
  struct svm_args parameters;
  struct denseData ds;
  struct Fullproblem fp;
  struct Projected sp;

  // Input processed:
  parse_arguments(argc, argv, &filename, &parameters);
  read_file(filename, &ds);
  preprocess(&ds);

  // Projected problem size chosen temporarily
  int p = 4;


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

    // H matrix columns fixed and subproblem changed
    setH(&fp, &ds);
    init_subprob(&sp, &fp, &ds);

    //  congjugate gradient algorithm
    //  n = 0 if algorithm completes
    //  n != 0 if algorithm interrupt
    n = cg(&sp, &fp);

    updateAlphaR(&fp, &sp);
    calcYTR(&sp, &fp);
    calculateBeta(&fp, &sp, &ds);

    printf("n is %d\n",n );
    for (int i = 0; i < fp.p; i++) {
      //printf("alpha[%d] = %lf\n",fp.active[i],fp.alpha[fp.active[i]] );
    }

    if (n==0) {
      printf("Converged!\n" );
      for (int i = 0; i < fp.n; i++) {
        //printf("alpha[%d] = %lf\n",i,fp.alpha[i] );
      }
      k = swapMostNegative(&fp);
    }

    while (n) {
      // While BCs broken
      k = singleswap(&ds, &fp, &sp, n);
      n = checkfpConstraints(&fp);
      printf("n is %d and alpha is %lf\n",n,fp.alpha[fp.active[0]]);
    }

    printf("active is %d %d %d %d\n",fp.active[0],fp.active[1],fp.active[2],fp.active[3] );

    for (int i = 0; i < fp.n; i++) {
    //  printf("alpha[%d] = %lf\n",i,fp.alpha[i] );
    }
    for (int i = 0; i < fp.q; i++) {
    //  printf("beta[%d] = % lf\n",i,fp.beta[i] );
    }
    itt++;
    if(itt == max_iters){
      break;
    }

  }
  for (int i = 0; i < fp.n; i++) {
    printf("alpha[%d] = %lf\n",i,fp.alpha[i] );
  }
  for (int i = 0; i < fp.q; i++) {
    printf("beta[%d] = % lf\n",i,fp.beta[i] );
  }
  if (k==0) {
    printf("k goes to 0.\n");
  }

  //Memory freed
  freeDenseData(&ds);
  freeFullproblem(&fp);
  freeSubProblem(&sp);

  return 0;
}



void   freeSubProblem( struct Projected* sp)
/* Function to free dynamically allocated memory in subproblem stuct. */
{
  free(sp->alphaHat);
  free(sp->yHat);
  free(sp->rHat);
  free(sp->H);

  free(sp->gamma);
  free(sp->rho);
  free(sp->Hrho);

  free(sp->h);
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

void  freeFullproblem(struct Fullproblem *fp)
/* Function to free dynamically allocated memory in Fullproblem struct */
{
  free(fp->alpha);
  free(fp->beta);
  free(fp->gradF);

  free(fp->active);
  free(fp->inactive);

  free(fp->partialH);
  free(fp->h);
}
