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
  sp->Hgamma = malloc(sizeof(double)*p);
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
/* Conjugate gradient method to solve projected subproblem. */
{
  double lambda, Hlambda, mu;
  init_error(sp);
  double rSq = inner_prod(sp->gamma,sp->gamma,sp->p);

  double newRSQ;
  int problem = 0;
  double divisor;

  while (rSq > 0.000000001) {
    calc_Hrho(sp);

    divisor = inner_prod(sp->Hrho, sp->rho, sp->p);
    if (fabs(divisor) < 0.00000000000000000000000000001) {
      printf("%lf\n",inner_prod(sp->Hrho, sp->rho, sp->p) );
      for (int j = 0; j < fp->p; j++) {
        printf("act = %d\n",fp->active[j] );
      }
      exit(250);
    }

    lambda = inner_prod(sp->gamma, sp->rho, sp->p)/divisor;
    Hlambda = inner_prod(sp->Hgamma, sp->rho, sp->p)/divisor;

    linearOp(sp->alphaHat, sp->rho, lambda, sp->p);

    problem = checkConstraints(sp, fp);

    if(problem){
      if (problem >= sp->p*2) {
        linearOp(sp->alphaHat, sp->rho, -lambda, sp->p);
        return problem;
      }
      else if( problem < -sp->p){
        exit(79);
        //linearOp(sp->alphaHat, sp->rho, -lambda, sp->p);
        //return problem + sp->p;
      }
      return problem;
    }

    linearOp(sp->gamma, sp->Hrho, lambda, sp->p);
    linearOp(sp->Hgamma, sp->Hrho, Hlambda, sp->p);

    newRSQ = inner_prod(sp->gamma, sp->gamma, sp->p);

    mu = inner_prod(sp->Hgamma, sp->Hrho, sp->p)/divisor;

    linearOp2(sp->rho, sp->gamma, mu, sp->p);

    rSq = newRSQ;
  }
  return 0;
}

void linearOp2(double* vecOut, double* vecIn, double a, int p)
{
  for (int i = 0; i < p; i++) {
    vecOut[i] *= -a;
    vecOut[i] += vecIn[i];
  }
}

void linearOp(double* vecOut, double* vecIn, double a, int p)
/* Function to add constant times */
{
  for (int i = 0; i < p; i++) {
    vecOut[i] -= a*vecIn[i];
  }
}

void calcYTR(struct Projected *sp, struct Fullproblem *fp)
/* Function to calculate the innter product of the projected*/
{
  sp->ytr = 0.0;
  for (int i = 0; i < sp->p; i++) {
    sp->ytr += sp->yHat[i]*fp->gradF[fp->active[i]];
  }
}



void updateGamma(struct Projected *sp, double lambda, double Hlambda)
{

  for (int i = 0; i < sp->p; i++) {
    sp->gamma[i] -= sp->Hrho[i]*lambda;
    for (int j = 0; j < sp->p; j++) {
      sp->gamma[i] += lambda*sp->yHat[i]*sp->yHat[j]*sp->Hrho[j]/((double)sp->p);
    }
  }

  for (int i = 0; i < sp->p; i++) {
    sp->gamma[i] -= sp->Hrho[i]*Hlambda;
    for (int j = 0; j < sp->p; j++) {
      sp->gamma[i] += Hlambda*sp->yHat[i]*sp->yHat[j]*sp->Hrho[j]/((double)sp->p);
    }
  }
}




void calc_Hrho(struct Projected *sp)
/* Function to multiply rho by H */
{
  double *temp = malloc(sizeof(double)*sp->p);
  for (int i = 0; i < sp->p; i++) {
    temp[i] = 0.0;
    for (int j = 0; j < i; j++) {
      temp[i] += sp->H[j][i]*sp->rho[j];
    }
    for (int j = i; j < sp->p; j++) {
      temp[i] += sp->H[i][j]*sp->rho[j];
    }
  }
  for (int i = 0; i < sp->p; i++) {
    sp->Hrho[i] = temp[i];
    for (int j = 0; j < sp->p; j++) {
      sp->Hrho[i] -= sp->yHat[i]*sp->yHat[j]*temp[j]/(double)sp->p;
    }
  }
  free(temp);
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
  Neg_constraint_projection(sp->gamma, sp->rHat, sp->yHat, sp->p);

  for (int i = 0; i < sp->p; i++) {
    sp->Hgamma[i] = 0.0;
    for (int j = 0; j < i; j++) {
      sp->Hgamma[i]+= sp->H[j][i]*sp->gamma[j];
    }
    for (int j = i; j < sp->p; j++) {
      sp->Hgamma[i]+= sp->H[i][j]*sp->gamma[j];
    }
  }

  for (int i = 0; i < sp->p; i++) {
    sp->rho[i] = sp->Hgamma[i];
    for (int j = 0; j < sp->p; j++) {
      sp->rho[i] -= sp->Hgamma[j]/(double)sp->p;
    }
  }

  copy_vector(sp->Hgamma, sp->rho, sp->p);
}

void copy_vector(double* a, double* b, int p)
{
  for (int i = 0; i < p; i++) {
    a[i] = b[i];
  }
}

void Neg_constraint_projection(double* vecOut, double* vecIn, double* y, int p)
{
  for (int i = 0; i < p; i++) {
    vecOut[i] = -vecIn[i];
    for (int j = 0; j < p; j++) {
      vecOut[i] += vecIn[j]*(y[i]*y[j]/((double)p));
    }
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
