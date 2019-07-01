#include "fullproblem.h"

void alloc_prob(struct Fullproblem *prob, struct denseData *ds, int p)
/* Function to allocate space necessary for a full problem of size n,
 * that will be projected to size p.  */
{
  prob->n = ds->nInstances;
  prob->p = p;
  prob->q = prob->n-prob->p;
  prob->C = 1000.0;
  prob->alpha = (double*)calloc(prob->n, sizeof(double) );
  prob->gradF = (double*)malloc(sizeof(double)*prob->n);
  prob->active = (int*)malloc(sizeof(int)*prob->p);
  prob->inactive = (int*)malloc(sizeof(int)*prob->q);
  prob->beta = (double*)malloc(sizeof(double)*(prob->q));
  prob->h = (double*)malloc(sizeof(double)*prob->q*prob->p);
  prob->partialH = (double**)malloc(sizeof(double*)*prob->q);
  for (int i = 0; i < prob->q; i++) {
    prob->partialH[i] = &prob->h[i*prob->p];
  }

}

void init_prob(struct Fullproblem *prob, struct denseData *ds)
/* Function to initialize values for full problem of size n,
 * that will be projected to size p.  */
{
  for (int i = 0; i < prob->n; i++) {
    prob->gradF[i] = 1.0;
  }

  for (int i = 0; i < prob->p/2; i++) {
    prob->active[i] = i;
  }

  for (int i = 0; i < prob->p/2; i++) {
    prob->active[i+prob->p/2] = ds->nPos + i;
  }

  if (prob->p%2 != 0) {
    prob->active[prob->p-1] = ds->nPos+ prob->p/2;
    for (int i = prob->p/2; i < ds->nPos; i++) {
      prob->inactive[i-(prob->p/2)] = i;
    }
    for (int i = 1+(prob->p/2); i < ds->nNeg; i++) {
      prob->inactive[ds->nPos-(prob->p)+i] = ds->nPos + i;
    }
    fprintf(stderr, "fullproblem.cpp: init_prob(): Not yet working for odd\n");
    exit(1);
  }
  else{
    for (int i = prob->p/2; i < ds->nPos; i++) {
      prob->inactive[i-(prob->p/2)] = i;
    }
    for (int i = prob->p/2; i < ds->nNeg; i++) {
      prob->inactive[ds->nPos-(prob->p)+i] = ds->nPos + i;
    }
  }

}

void setH(struct Fullproblem *prob, struct denseData *ds)
/*  Function to update the values of the matrix partialH
 *  TO DO -> ALLOW ONLY UPDATE THE SWAPPED OUT VALUES - or not?   */
{
  for (int i = 0; i < prob->q; i++) {
    for (int j = 0; j < prob->p; j++) {
      prob->partialH[i][j] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        prob->partialH[i][j] += ds->data[prob->inactive[i]][k]*ds->data[prob->active[j]][k];
      }
      prob->partialH[i][j] *= ds->y[prob->inactive[i]]*ds->y[prob->active[j]];
    }
  }
}

void updateAlphaR(struct Fullproblem *fp, struct Projected *sp)
/* Function to update the values of the alpha, r vectors after a cg sweep*/
{
  // Use rho as temporary place to store P*alpha:
  for (int i = 0; i < sp->p; i++) {
    printf("alphaHat = %lf\n",sp->alphaHat[i] );
    sp->rho[i] = sp->alphaHat[i];
    for (int j = 0; j < sp->p; j++) {
      sp->rho[i] -= ((sp->yHat[i]*sp->yHat[j])/((double)(sp->p)))*sp->alphaHat[j];
    }
    printf("rho[%d] = %lf\n",i,sp->rho[i] );
  }

  // Alpha of each active point is updated:
  for (int i = 0; i < sp->p; i++) {
    fp->alpha[fp->active[i]] += sp->rho[i];
  }

  // gradF of each inactive point is updated:
  for (int i = 0; i < fp-> q; i++) {
    for (int j = 0; j < sp->p; j++) {
      fp->gradF[fp->inactive[i]] -= fp->partialH[i][j] * sp->rho[j];
    }
  }

  // gradF of each active point is updated:
  for (int i = 0; i < fp->p; i++) {
    for (int j = 0; j < i; j++) {
      fp->gradF[fp->active[i]] -= sp->H[j][i]*sp->rho[j];
    }
    for (int j = i; j < fp->p; j++) {
      fp->gradF[fp->active[i]] -= sp->H[i][j]*sp->rho[j];
    }
  }
}

void calculateBeta(struct Fullproblem *fp, struct Projected *sp, struct denseData *ds)
/* Function to calculate the beta vector - tracks how suboptimal the values for*/
{
  for (int i = 0; i < fp->q; i++) {
    if (fp->alpha[fp->inactive[i]] < 0.0000000001) {
      fp->beta[i] =  - fp->gradF[fp->inactive[i]] + ((ds->y[fp->inactive[i]]*sp->ytr)/((double)sp->p));
    }
    else if (fp->alpha[fp->inactive[i]] >= sp->C - 0.001) {
      fp->beta[i] =  fp->gradF[fp->inactive[i]] - ds->y[fp->inactive[i]]*sp->ytr/((double)sp->p);
    }
  }
}

void findWorst(int *worst, int* target, int* change, int *n, struct denseData *ds, struct Fullproblem *fp)
{
  double tester = DBL_MAX;

  if (*n > 0) {
    (*change) = 1;
    *n -= fp->p;
  }
  else {
    (*change) = -1;
    *n+=fp->p;
  }
  *target = ds->y[fp->active[*n]]*(*change);

  for (int i = 0; i < fp->q; i++) {
    if( ds->y[fp->inactive[i]] == *target )  {
      if (fp->beta[i] < tester) {
        if (fp->alpha[fp->inactive[i]] < fp->C){
          *worst = i;
          tester = fp->beta[i];
        }
      }
    }
    if( ds->y[fp->inactive[i]] != *target )  {
      if (fp->beta[i] < tester) {
        if (fp->alpha[fp->inactive[i]] > 0.0){
          *worst = i;
          tester = fp->beta[i];
        }
      }
    }
  }
  if (tester > 0.0) {
    *worst = -1;
  }
}

int singleswap(struct denseData *ds, struct Fullproblem *fp, struct Projected *sp, int n)
{
  int worst = -1;
  int target, change=1;


  findWorst(&worst,&target,&change,&n, ds, fp);

  printf("worst is %d\n",worst );
  double diff;
  if (change<0) {
    diff = fp->alpha[fp->active[n]];
  }
  else{
    diff = fp->C - fp->alpha[fp->active[n]];
  }

  if(worst < 0)
  {
    for (int i = 0; i < fp->p; i++) {
      printf("fp->alpha[i] = %lf\n",fp->alpha[i]);
    }
    printf("\n" );
    for (int i = 0; i < sp->p; i++) {
      if (i != n) {
        if (ds->y[fp->active[i]] == target){
          fp->alpha[fp->active[i]] += diff/(double)(fp->p-1);
        }
        else{
          fp->alpha[fp->active[i]] -= diff/(double)(fp->p-1);
        }
      }
    }
    if (change < 0) {
      fp->alpha[fp->active[n]] = 0;
    }
    else{
      fp->alpha[fp->active[n]] = fp->C;
    }
    for (int i = 0; i < fp->p; i++) {
      printf("fp->alpha[i] = %lf\n",fp->alpha[i]);
    }
    return n - fp->p;
  }


  int temp = fp->active[n];

  if (change < 0)
  {
    if (ds->y[fp->inactive[worst]] == target) {
      adjustGradF(fp, ds, sp, n, worst, change, 1);
      fp->alpha[fp->inactive[worst]] -= fp->alpha[ fp->active[n] ];
      fp->alpha[fp->active[n]] = 0.0;//sp->C ;
      fp->active[n] = fp->inactive[worst];
      fp->inactive[worst] = temp;
      fp->beta[worst] = DBL_MAX-1.0;
    }
    else {
      adjustGradF(fp, ds, sp, n, worst, change, 0);
      fp->alpha[fp->inactive[worst]] += fp->alpha[fp->active[n]];
      fp->alpha[fp->active[n]] = 0.0;//sp->C ;
      fp->active[n] = fp->inactive[worst];
      fp->inactive[worst] = temp;
      fp->beta[worst] = DBL_MAX;
    }
  }
  else
  {
    if (ds->y[fp->inactive[worst]] == target) {
      adjustGradF(fp, ds, sp, n, worst, change, 1);
      fp->alpha[fp->inactive[worst]] -= sp->C - fp->alpha[fp->active[n]];
      fp->alpha[fp->active[n]] = sp->C;
      fp->active[n] = fp->inactive[worst];
      fp->inactive[worst] = temp;
      fp->beta[worst] = DBL_MAX;
    }
    else{
      adjustGradF(fp, ds, sp, n, worst, change, 0);
      //printf("alpha is %lf  and  target is %d and ina is %lf\n\n\n",fp->alpha[fp->active[n]], target,ds->y[fp->inactive[worst]]);
      fp->alpha[fp->inactive[worst]] += sp->C - fp->alpha[fp->active[n]];
      fp->alpha[fp->active[n]] = sp->C;
      fp->active[n] = fp->inactive[worst];
      fp->inactive[worst] = temp;
      fp->beta[worst] = DBL_MAX;
    }
  }

  if(fp->alpha[fp->active[n]] < -0.0001 || fp->alpha[fp->active[n]] > sp->C + 0.0001){
    //swapWithException(fp, ds, n, worst, );
  }

  return 1;
}

int checkfpConstraints(struct Fullproblem *fp)
{
  printf("\n" );
  for (int i = 0; i < fp->p; i++) {
    //printf("alpha is %lf for %d\n",fp->alpha[fp->active[i]],fp->active[i] );
    if(fp->alpha[fp->active[i]]>fp->C){
      //printf("greater than C fp[%d]-> = %lf \n",fp->active[i],fp->alpha[fp->active[i]] );
      return i+fp->p;
    }
    else if(fp->alpha[fp->active[i]] < 0.0){
      //printf("less than 0 fp[%d]-> = %lf \n",fp->active[i],fp->alpha[fp->active[i]] );
      return i-fp->p;
    }
  }
  return 0;
}

void adjustGradF(struct Fullproblem *fp, struct denseData *ds, struct Projected *sp, int n, int worst, int signal, int target)
{
  printf("n is %d and worst is %d and signal is %d targ %d\n", n, fp->inactive[worst], signal, target);
  // Update based on change of H matrix:
  if (signal == -1) {
    printf("fp == %lf\n",fp->alpha[n] );
    for (int i = 0; i < fp->q; i++) {
      fp->gradF[ fp->inactive[i] ] += fp->partialH[i][n]*fp->alpha[fp->active[n]] ;
    }
    for (int i = 0; i < fp->p; i++) {
      if(i<n){
        fp->gradF[ fp->active[i] ] += sp->H[i][n]*fp->alpha[fp->active[n]];
      }
      else{
        fp->gradF[ fp->active[i] ] += sp->H[n][i]*fp->alpha[fp->active[n]];
      }
    }
  }
  else
  {
    for (int i = 0; i < fp->q; i++) {
      fp->gradF[ fp->inactive[i] ] += fp->partialH[i][n]*(fp->alpha[fp->active[n]] - fp->C) ;
    }
    for (int i = 0; i < fp->p; i++) {
      if(i<n){
        fp->gradF[ fp->active[i] ] += sp->H[i][n]*(fp->alpha[fp->active[n]] - fp->C) ;
      }
      else{
        fp->gradF[ fp->active[i] ] += sp->H[n][i]*(fp->alpha[fp->active[n]] - fp->C) ;
      }
    }
  }

  // Update based on change of
  double* temp = (double*) malloc(sizeof(double)*fp->n);

  for (int j = 0; j < fp->n; j++) {
    temp[j] = 0.0;
    for (int k = 0; k < ds->nFeatures; k++) {
      temp[j]+=ds->data[fp->inactive[worst]][k]*ds->data[j][k];
    }
    temp[j]*=ds->y[fp->inactive[worst]]*ds->y[j];
  }

  if (signal == -1) {
    if (target) {
      for (int i = 0; i < fp->n; i++) {
        fp->gradF[i] += temp[i]*fp->alpha[fp->active[n]] ;
      }
    }
    else{
      for (int i = 0; i < fp->n; i++) {
        fp->gradF[i] -= temp[i]*fp->alpha[fp->active[n]] ;
      }
    }
  }
  else
  {
    if (target) {
      for (int i = 0; i < fp->n; i++) {
        fp->gradF[i] += temp[i]*(fp->C - fp->alpha[fp->active[n]]) ;
      }
    }
    else{
      for (int i = 0; i < fp->n; i++) {
        fp->gradF[i] -= temp[i]*(fp->C - fp->alpha[fp->active[n]]) ;
      }
    }
  }
  free(temp);
}

int swapMostNegative(struct Fullproblem *fp)
/*  If the cg has fully converged we swap all p values out for a new set
 *  based on how negative their beta value is.    */
{
  int* location = (int*)malloc(sizeof(int)*fp->p);
  int* index = (int*)malloc(sizeof(int)*fp->p);
  double* betaVal = (double*)malloc(sizeof(double)*fp->p);
  for (int i = 0; i < fp->p; i++)
  {
    location[i] = -1;
    betaVal[i] = DBL_MAX;
  }
  for (int i = 0; i < fp->q; i++)
  {
    for (int j = 0; j < fp->p; j++)
    {
      if (fp->beta[i]<betaVal[j])
      {
        for (int k = (fp->p) - 1; k > j ; k--)
        {
          location[k] = location[k-1];
          betaVal[k] = betaVal[k-1];
        }
        location[j] = i;
        betaVal[j] = fp->beta[i];
        break;
      }
    }
  }
  if (betaVal[0]>=0.0)
  {
    free(location);
    free(index);
    free(betaVal);
    return 0;
  }

  for (int i = 0; i < fp->p; i++)
  {
    index[i] = fp->inactive[location[i]];
    fp->inactive[location[i]] = fp->active[i];
    fp->active[i] = index[i];
  }


  free( location );
  free( index );
  free( betaVal );

  return 1;
}
