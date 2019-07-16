#include <iostream>
#include "svm.h"
#include "io.h"
#include "subproblem.h"
#include "fullproblem.h"

void inner_prod(double* res, double* x, double* y,int n);
void project_H(double** H, int p, int* P, double* h_hat, double** H_hat);

int main(int argc, char *argv[]) {
  char* filename = NULL;
  struct svm_args parameters;
  struct denseData ds;

  parse_arguments(argc, argv, &filename, &parameters);
  read_file(filename, &ds);

  Fullproblem states(ds);
  //states.preprocess();
  int p = 8;
  int P[8];
  for (int i = 0; i < p; i++) {
    P[i] = i;
  }
  Subproblem mini(states.H, states.y, p, P);

  mini.solve();



  /*for (int i = 0; i < ds.nInstances; i++) {
    for (int j = 0; j < ds.nFeatures; j++) {
      std::cout << ds.data[i][j] << "\t";
    }
    std::cout << ds.y[i] << '\n';
  }*/


  return 0;
}

void reallocH(double* h_hat, double** H_hat, int p)
{
  h_hat = (double*)(realloc(h_hat,sizeof(double)*((p*p+1)/2)));
  H_hat = (double**)realloc(H_hat,sizeof(double*)*(p));
  int j=0;
  for (int i = 0; i < p; i++) {
    j+=i;
    H_hat[i] = &h_hat[j];
  }
}

void project_H(double** H, int p, int* P, double* h_hat, double** H_hat)
{
  for (int i = 0; i < p; i++) {
    for (int j = p-1; j >= i; j--) {
      H_hat[i][j] = H [P[i]][P[j]];
    }
  }
}





void p_times_vec(double* alphaHat, double alpha, int p, int m, denseData ds){
  for (int i = 0; i < p; i++) {
    //alpha += h;
  }
}
