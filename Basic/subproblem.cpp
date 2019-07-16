#include "subproblem.h"
#include <iostream>

Subproblem::Subproblem(double** supH, double* supY, int sp, int* sP)
{
  h = NULL;
  H = NULL;
  y = NULL;
  p = sp;
  P = (int*)malloc(sizeof(int)*p);
  changeOut = (bool*)malloc(sizeof(bool)*p);
  for (int i = 0; i < p; i++) {
    changeOut[i] = false;
  }
  c = 3.0;
  for (int i = 0; i < p; i++) {
    P[i] = sP[i];
  }
  alloc_space(p);
  project_down(supH, supY, p, P);
  alpha = (double*)malloc(sizeof(double)*p);
  gradF = (double*)malloc(sizeof(double)*p);
  for (int i = 0; i < p; i++) {
    alpha[i] = 0.0;
    gradF[i] = 1.0;
  }
}


void Subproblem::alloc_space(int p)
{
  h = (double*)(realloc(h,sizeof(double)*((p*(p+1))/2)));
  H = (double**)realloc(H,sizeof(double*)*(p));
  y = (double*)realloc(y,sizeof(double*)*p);
  int j = 0;
  for (int i = 0; i < p; i++) {
    H[i] = &h[j];
    j+=(p-i-1);
    //std::cout << "i = " << i << " j == " << j << '\n';
  }
}

void Subproblem::project_down(double** supH, double* supY, int p, int* P)
{
  for (int i = 0; i < p; i++) {
    for (int j = i; j < p; j++) {
      //std::cout << P[i] << P[j]<< '\n';
      H[i][j] = supH [P[i]][P[j]];
      y[i] = supY[P[i]];
      //std::cout << "H[" << i << "][" << j << "]"<< H[i][j] << '\t';
    }
  }

  for (int i = 0; i < p; i++) {
    for (int j = i; j < p; j++) {
      std::cout << "H[" << i << "][" << j << "]"<< H[i][j] << '\t';
    }
    std::cout << '\n';
  }
}

void Subproblem::constraint_projection(double* vecOut, double* vecIn)
{
  for (int i = 0; i < p; i++) {
    vecOut[i] = vecIn[i];
    for (int j = 0; j < p; j++) {
      vecOut[i] -= vecIn[j]*(y[i]*y[j]/(double)p);
    }
  }
}


void Subproblem::HtimesVec(double* vecOut, double* vecIn)
{
  for (int i = 0; i < p; i++) {
    vecOut[i] = 0;
    for (int j = 0; j < i; j++) {
      vecOut[i] += vecIn[j]*H[j][i];
    }
    for (int j = i; j < p; j++) {
      vecOut[i] += vecIn[j]*H[i][j];
    }
  }
}

void Subproblem::solve()
{
  double* err = (double*)malloc(sizeof(double)*p);
  double* perr = (double*)malloc(sizeof(double)*p);
  double* Aperr = (double*)malloc(sizeof(double)*p);
  cg(err,perr,Aperr);
  for (int i = 0; i < p; i++) {
    std::cout << alpha[i] << '\n';
  }
}

void Subproblem::cg(double* err, double* perr, double* Aperr)
{
  double alp, beta;
  init_error(err, perr);
  double rSq = inner_prod(err,err);
  double newRSQ;
  int flag;
  for (int i = 0; i < p; i++) {
    std::cout << perr[i] << '\n';
  }
  while (rSq>0.001) {
    php(Aperr, perr);
    for (int i = 0; i < p; i++) {
      std::cout << "A = " << Aperr[i] << '\n';
    }
    alp = inner_prod(Aperr, perr);
    linOp(alpha, perr, alp);
    linOp(err, Aperr, -alp);
    newRSQ = inner_prod(err, err);
    beta = newRSQ/rSq;
    linOp2(perr,err,beta);
    flag = 0;
    std::cout << alpha[0] << '\n';
    for (int i = 0; i < p; i++) {
      if (alpha[i]>=c) {
        flag = 1;
        alpha[i] = c;
        changeOut[i] = true;
      }
      else if (alpha[i]<=0.0)  {
        changeOut[i] = true;
        flag = 1;
        alpha[i] = 0.0;
      }
    }
    if (flag == 1) {
      break;
    }
  }
}

void Subproblem::linOp(double* err, double* Aperr, double alp)
{
  for (int i = 0; i < p; i++) {
    err[i] += alp*Aperr[i];
  }
}

void Subproblem::linOp2(double* err, double* Aperr, double alp)
{
  for (int i = 0; i < p; i++) {
    err[i] *= alp;
    err[i] += Aperr[i];
  }
}


void Subproblem::php(double* Aperr, double* perr)
{
  double* temp1 = (double*)malloc(sizeof(double)*p);
  constraint_projection(Aperr, perr);
  HtimesVec(temp1,Aperr);
  constraint_projection(Aperr,temp1);
  free(temp1);
}



double Subproblem::inner_prod(double* a1, double* a2)
{
  double res = 0;
  for (int i = 0; i < p; i++) {
    res+=a1[i]*a2[i];
  }
  return res;
}

void Subproblem::init_error(double* err, double* perr)
{
  double* temp1 = (double*)malloc(sizeof(double)*p);
  double* temp2 = (double*)malloc(sizeof(double)*p);
  constraint_projection(err, gradF);

  constraint_projection(temp1, alpha);
  HtimesVec(temp2, temp1);
  constraint_projection(temp1, temp2);
  for (int i = 0; i < p; i++) {
    err[i]-= temp1[i];
    perr[i] = err[i];
  }
  free(temp1);
  free(temp2);
}
