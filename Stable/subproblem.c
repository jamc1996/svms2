#include "subproblem.h"


void alloc_subprob(struct Projected *sp, int p)
/*  Function to allocate space to solve the projected problem of size p.
 */
{
  sp->p = p;
  sp->C = 1000.0;

  sp->alphaHat = malloc(sizeof(double)*p);
  sp->yHat = malloc(sizeof(double)*p);
  sp->rHat = malloc(sizeof(double)*p);
  sp->gamma = malloc(sizeof(double)*p);
  sp->rho = malloc(sizeof(double)*p);
  sp->Hrho = malloc(sizeof(double)*p);

  // H symmetric so can save space:
  sp->H = malloc(sizeof(double*)*p);
  sp->h = malloc(sizeof(double)*((p*(p+1))/2));
  int j = 0;
  for (int i = 0; i < p; i++) {
    sp->H[i] = &(sp->h[j]);
    j+=(p-i-1);
  }
}

void init_subprob(struct Projected *sp, struct Fullproblem *fp, struct denseData *ds)
/* Function to initialize the values of the subproblem struct, given information
 * in a dataset and fp struct with active/inactive vectors initialized.

 *  Required: Everything allocated, active and inactive correct.
 *

 */
{
  for (int i = 0; i < sp->p; i++) {
    sp->yHat[i] = ds->y[fp->active[i]];
    sp->alphaHat[i] = 0.0;    // This is the change from original
    sp->rHat[i] = fp->gradF[fp->active[i]];
  }


  for (int i = 0; i < sp->p; i++) {
    for (int j = i; j < sp->p; j++) {
      sp->H[i][j] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        sp->H[i][j] += ds->data[fp->active[i]][k]*ds->data[fp->active[j]][k];
      }
      sp->H[i][j] *= ds->y[ fp->active[i] ] * ds->y[ fp->active[j] ];
    }
  }
}

int cgSLS(struct Projected *sp, struct Fullproblem *fp)
/* Maybe a more stable build. */
{
  double* x = malloc(sizeof(double)*sp->p);
  double* y = malloc(sizeof(double)*sp->p);
  double* g = malloc(sizeof(double)*sp->p);
  double* gA = malloc(sizeof(double)*sp->p);
  double* p = malloc(sizeof(double)*sp->p);
  for (int i = 0; i < sp->p; i++) {
    x[i] = 0.0;
    y[i] = 0.0;
    g[i] = -sp->rHat[i];
    for (int j = 0; j < sp->p; j++) {
      g[i] += sp->yHat[j]*sp->rHat[j]/(double)sp->p;
    }
    gA[i] = 0.0;
    for (int j = 0; j < sp->p; j++) {
      if(i<j){
        gA[i] += sp->H[i][j]*g[j];
      }
      else{
        gA[i] += sp->H[j][i]*g[j];
      }
    }
  }
  for (int i = 0; i < sp->p; i++) {
    for (int j = 0; j < sp->p; j++) {
      if(i<j){
        gA[i] += sp->H[i][j]*g[j];
      }
      else{
        gA[i] += sp->H[j][i]*g[j];
      }
    }
    p[i] = gA[i];
  }
  double AxySq =0.0;
  double gAsQ =0.0;
  double lambda = 0.0, Alambda = 0.0;
  while(AxySq + gAsQ > 0.001)
  {
    lambda = 8;
    Alambda = 9;
  }
  return 77;
}

int cg(struct Projected *sp, struct Fullproblem *fp)
/* Conjugate gradient method to solve projected subproblem. */
{
  double lambda, mu;
  init_error(sp);
  double rSq = inner_prod(sp->gamma,sp->gamma,sp->p);

  double newRSQ;
  int problem = 0;

  while (rSq > 0.0000001) {
    calc_Hrho(sp);

    if (fabs(inner_prod(sp->Hrho, sp->rho, sp->p)) < 0.00000000000000000000000000001) {
      printf("%lf\n",inner_prod(sp->Hrho, sp->rho, sp->p) );
      for (int j = 0; j < fp->p; j++) {
        printf("act = %d\n",fp->active[j] );
      }
      exit(250);
    }
    lambda = rSq/inner_prod(sp->Hrho, sp->rho, sp->p);
    linearOp(sp->alphaHat, sp->rho, lambda, sp->p);

    problem = checkConstraints(sp, fp);

    if(problem){
      if (problem >= sp->p*2) {
        printf("Oh no! (%d) %d\n",sp->p,fp->active[problem - (sp->p*2)] );
        linearOp(sp->alphaHat, sp->rho, -lambda, sp->p);
        return problem;
      }
      else if( problem < -sp->p*2){
        exit(79);
        //linearOp(sp->alphaHat, sp->rho, -lambda, sp->p);
        //return problem + sp->p;
      }
      return problem;
    }

    updateGamma(sp, lambda);

    newRSQ = inner_prod(sp->gamma, sp->gamma, sp->p);

    mu = newRSQ/rSq;

    linearOp2(sp->rho, sp->gamma, mu, sp->p);

    rSq = newRSQ;

  }

  return 0;
}

void calcYTR(struct Projected *sp, struct Fullproblem *fp)
/* Function to calculate the innter product of the projected*/
{
  sp->ytr = 0.0;
  for (int i = 0; i < sp->p; i++) {
    sp->ytr += sp->yHat[i]*fp->gradF[fp->active[i]];
  }
}

void linearOp2(double* vecOut, double* vecIn, double a, int p)
{
  for (int i = 0; i < p; i++) {
    vecOut[i] *= a;
    vecOut[i] += vecIn[i];
  }
}

void updateGamma(struct Projected *sp, double lambda)
{
  for (int i = 0; i < sp->p; i++) {
    sp->Hrho[i]*=lambda;
  }
  for (int i = 0; i < sp->p; i++) {
    sp->gamma[i] -= sp->Hrho[i];
    for (int j = 0; j < sp->p; j++) {
      sp->gamma[i] += sp->yHat[i]*sp->yHat[j]*sp->Hrho[j]/((double)sp->p);
    }
  }
}


void linearOp(double* vecOut, double* vecIn, double a, int p)
/* Function to add constant times */
{
  for (int i = 0; i < p; i++) {
    vecOut[i] += a*vecIn[i];
  }
}

void calc_Hrho(struct Projected *sp)
/* Function to multiply rho by H */
{
  for (int i = 0; i < sp->p; i++) {
    sp->Hrho[i] = 0.0;
    for (int j = 0; j < i; j++) {
      sp->Hrho[i] += sp->H[j][i]*sp->rho[j];
    }
    for (int j = i; j < sp->p; j++) {
      sp->Hrho[i] += sp->H[i][j]*sp->rho[j];
    }
  }
}



int checkConstraints(struct Projected* sp, struct Fullproblem *fp)
/* Function to check if any constrainst have been violated by the cg process. */
{
  int flag = -1;
  double* temp = (double*)malloc(sizeof(double)*sp->p);
  constraint_projection(temp, sp->alphaHat, sp->yHat, sp->p);
  for (int i = 0; i < sp->p; i++) {
    temp[i] += fp->alpha[ fp->active[i] ];
  }
  for (int i = 0; i < sp->p; i++) {
    if(temp[i]>2*sp->C){
      flag = i;
      for (int j = 0; j < sp->p; j++) {
        if (temp[i] < temp[j]) {
          flag = j;
        }
      }
      free(temp);
      return flag+sp->p+sp->p;
    }
    //printf("temp[%d] = %lf\n",i,temp[i] );
  }

  for (int i = 0; i < sp->p; i++) {
    if(temp[i]>sp->C){
      if (temp[i]>sp->C*2) {
        //free(temp);
        //return i+(2*sp->p);
      }
      free(temp);
      return i+sp->p;
    }
    else if(temp[i]<0.0){
      if (temp[i]<sp->C) {
        //free(temp);
        //return i+(2*sp->p);
      }
      free(temp);
      return i-sp->p;
    }
  }
  free(temp);
  return 0;
}

void init_error(struct Projected* sp)
/* Function to initialize the gamma and rho vectors for cg. */
{
  constraint_projection(sp->gamma, sp->rHat, sp->yHat, sp->p);
  copy_vector(sp->rho, sp->gamma, sp->p);
}

void copy_vector(double* a, double* b, int p)
{
  for (int i = 0; i < p; i++) {
    a[i] = b[i];
  }
}

void constraint_projection(double* vecOut, double* vecIn, double* y, int p)
{
  for (int i = 0; i < p; i++) {
    vecOut[i] = vecIn[i];
    for (int j = 0; j < p; j++) {
      vecOut[i] -= vecIn[j]*(y[i]*y[j]/((double)p));
    }
  }
}

double inner_prod(double *a, double *b, int p)
{
  double val = 0.0;
  for (int i = 0; i < p; i++) {
    val+=a[i]*b[i];
  }
  return val;
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
