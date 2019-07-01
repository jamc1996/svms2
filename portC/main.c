#include "svm.h"
#include "io.h"
#include "fullproblem.h"
#include "subproblem.h"

#include <stdio.h>

void freeDenseData(struct denseData *ds);
void  freeFullproblem(struct Fullproblem *fp);
void   freeSubProblem( struct Projected* sp);
void changeP( struct Fullproblem *fp, struct Projected *sp, int add);
int findWorstest(struct Fullproblem *fp , int add, int* temp, int* temp2);
void reinitprob( struct Fullproblem *fp, struct Projected *sp, int add, int* temp, int* temp2);

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
  int max_iters = 30;
  int itt = 0;
  int n = 0;

  while(k){
    // H matrix columns re-set and subproblem changed
    printf("back we are\n");
    for (int i = 0; i < sp.p; i++) {
      printf("active[%d] = %d\n",i,fp.active[i] );
    }
    setH(&fp, &ds);
    printf("possible probs\n");

    init_subprob(&sp, &fp, &ds);
    printf("possible probs\n");

    //  congjugate gradient algorithm
    //  n = 0 if algorithm completes
    //  n != 0 if algorithm interrupt
    printf("\n\n\n\n" );
    for (int i = 0; i < fp.p; i++) {
      printf("alpha[%d] = %lf\n",fp.active[i],fp.alpha[fp.active[i]] );
    }
    n = cg(&sp, &fp);

    updateAlphaR(&fp, &sp);
    calcYTR(&sp, &fp);
    calculateBeta(&fp, &sp, &ds);


    if (n==0) {
      printf("Converged!\n" );
      for (int i = 0; i < fp.q; i++) {
      printf("beta[%d](%d) = % lf\n",i,fp.inactive[i],fp.beta[i] );
      }
      int add = 2;
      int *temp = malloc(sizeof(int)*add);
      int *temp2 = malloc(sizeof(int)*add);
      add = findWorstest(&fp,add,temp,temp2);
      if (add == 0){
        break;
      }
      for (int i = 0; i < fp.p; i++) {
        printf("active[%d] = %d\n",i, fp.active[i] );
      }
      printf("\n" );
      for (int i = 0; i < fp.q; i++) {
        printf("inactive[%d] = %d\n",i, fp.inactive[i] );
      }
      changeP(&fp, &sp, add);
      reinitprob(&fp, &sp, add, temp, temp2);
      for (int i = 0; i < fp.p; i++) {
        printf("active[%d] = %d\n",i, fp.active[i] );
      }
      printf("\n" );
      for (int i = 0; i < fp.q; i++) {
        printf("inactive[%d] = %d\n",i, fp.inactive[i] );
      }
      free(temp);
      free(temp2);
      //break;
    }
    int prqe = 0;

    while (n) {
      // While BCs broken
      k = singleswap(&ds, &fp, &sp, n);
      n = checkfpConstraints(&fp);
      prqe++;
      if (prqe == fp.n) {
        return 1;
      }
      printf("n is %d and alpha is %lf\n",n,fp.alpha[fp.active[0]]);
    }

    printf("active is %d %d %d %d\n",fp.active[0],fp.active[1],fp.active[2],fp.active[3] );

    itt++;
    if(itt == max_iters){
      break;
    }
  }

  for (int i = 0; i < fp.n; i++) {
    printf("alpha[%d] = %lf\n",i,fp.alpha[i] );
  }
  for (int i = 0; i < fp.q; i++) {
  printf("beta[%d](%d) = % lf\n",i,fp.inactive[i],fp.beta[i] );
  }
  if (k==0) {
    //printf("k goes to 0.\n");
  }

  //Memory freed
  freeDenseData(&ds);
  freeFullproblem(&fp);
  freeSubProblem(&sp);

  return 0;
}

void changeP( struct Fullproblem *fp, struct Projected *sp, int add)
{
  printf("p is %d and q is %d\n",fp->p, fp->q );
  fp->p += add;
  fp->q -= add;
  printf("p is %d and q is %d\n",fp->p, fp->q );
  fp->active = realloc(fp->active,sizeof(int)*fp->p);
  fp->inactive = realloc(fp->inactive,sizeof(int)*fp->q);
  fp->beta = realloc(fp->beta,sizeof(double)*fp->q);

  fp->h = realloc(fp->h,sizeof(double)*fp->q*fp->p);
  fp->partialH = realloc(fp->partialH,sizeof(double*)*fp->q);
  for (int i = 0; i < fp->q; i++) {
    fp->partialH[i] = &fp->h[i*fp->p];
  }

  printf("reallocing\n" );

  sp->p += add;

  sp->alphaHat = realloc(sp->alphaHat,sizeof(double)*sp->p);
  printf("reallocing\n" );
  sp->yHat = realloc(sp->yHat,sizeof(double)*sp->p);

  sp->yHat = realloc(sp->yHat,sizeof(double)*sp->p);
  printf("reallocing\n" );

  sp->rHat = realloc(sp->rHat,sizeof(double)*sp->p);

  printf("reallocing\n" );

  sp->gamma = realloc(sp->gamma,sizeof(double)*sp->p);
  sp->rho = realloc(sp->rho,sizeof(double)*sp->p);
  sp->Hrho = realloc(sp->Hrho,sizeof(double)*sp->p);

  printf("reallocing\n" );

  sp->H = realloc(sp->H,sizeof(double*)*sp->p);
  sp->h = realloc(sp->h,sizeof(double)*((sp->p*(sp->p+1))/2));

  printf("reallocing\n" );

  int j = 0;
  for (int i = 0; i < sp->p; i++) {
    sp->H[i] = &(sp->h[j]);
    j+=(sp->p-i-1);
  }
}

int findWorstest(struct Fullproblem *fp , int add, int* temp, int* temp2)
{
  double *betaVal = malloc(sizeof(double)*add);
  for (int i = 0; i < add; i++) {
    temp2[i] = fp->inactive[fp->q - 2 + i];
    betaVal[i] = DBL_MAX;
  }
  for (int i = 0; i < fp->q; i++)
  {
    for (int j = 0; j < fp->p; j++)
    {
      if (fp->beta[i]<betaVal[j])
      {
        for (int k = add - 1; k > j ; k--)
        {
          temp[k] = temp[k-1];
          betaVal[k] = betaVal[k-1];
        }
        temp[j] = fp->inactive[i];
        betaVal[j] = fp->beta[i];
        break;
      }
    }
  }

  for (int i = 0; i < add; i++) {
    if (betaVal[i] > 0) {
      return i;
    }
  }

  printf("temp1 is %d and %d\n", temp[0],temp[1]);
  printf("temp2 is %d and %d\n", temp2[0],temp2[1]);

  return add;
}


void reinitprob( struct Fullproblem *fp, struct Projected *sp, int add, int* temp, int* temp2)
{
  for (int i = 0; i < add; i++) {
    fp->active[fp->p - add + i] = temp[i];
  }

  int k = 0;
  for (int i = 0; i < fp->q; i++) {
    for (int j = 0; j < add; j++) {
      if (fp->inactive[i] == temp[j]) {
        fp->inactive[i] = temp2[k];
        k++;
        if (k == add) {
          return;
        }
      }
    }
  }
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
