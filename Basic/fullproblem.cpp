#include "fullproblem.h"


Fullproblem::Fullproblem(struct denseData ds)
{
  nInstances = ds.nInstances;
  nFeatures = ds.nFeatures;
  data = ds.data;
  data1d = ds.data1d;
  y = ds.y;
  instanceLabels = ds.instanceLabels;
  featureLabels = ds.featureLabels;
  preprocess();
  alpha = (double*)malloc(sizeof(double)*nInstances);
  gradF = (double*)malloc(sizeof(double)*nInstances);
  for (int i = 0; i < nInstances; i++) {
    alpha[i] = 0.0;
    gradF[i] = 1.0;
  }
  fillH();

}

void Fullproblem::preprocess()
{
  double* means = (double*)calloc(nFeatures,sizeof(double));
  double* stdDev = (double*)calloc(nFeatures,sizeof(double));

  calcMeans(means);
  calcStdDev(stdDev,means);
  normalise(means,stdDev);
  for (int i = 0; i < nFeatures; i++) {
    std::cout << means[i] << '\n';
    std::cout << stdDev[i] << '\n';
  }
  free(means);
  free(stdDev);
}

void Fullproblem::normalise(double* mean, double* stdDev)
{
  for (int i = 0; i < nInstances; i++) {
    for (int j = 0; j < nFeatures; j++) {
      data[i][j]-=mean[j];
      data[i][j]/=stdDev[j];
    }
  }
}

void Fullproblem::calcMeans(double *mean)
{
  for (int i = 0; i < nInstances; i++) {
    for (int j = 0; j < nFeatures; j++) {
      mean[j] += data[i][j];
    }
  }
  for (int i = 0; i < nFeatures; i++) {
    mean[i]/=(double)nInstances;
  }
}

void Fullproblem::calcStdDev(double* stdDev, double* mean)
{
  for (int i = 0; i < nInstances; i++) {
    for (int j = 0; j < nFeatures; j++) {
      stdDev[j]+=(data[i][j]-mean[j])*(data[i][j]-mean[j]);
    }
  }
  for (int i = 0; i < nFeatures; i++) {
    stdDev[i] = sqrt(stdDev[i]/((double)nInstances-1.0));
  }
}

double Fullproblem::getHij(int i, int j){
  return H[i][j];
}

void Fullproblem::fillH()
{
  h = (double*)malloc((sizeof(double)*(nInstances+1)*nInstances)/2);
  H = (double**)malloc(sizeof(double*)*nInstances);
  int j = 0;
  for (int i = 0; i < nInstances; i++) {
    H[i] = &h[j];
    j+=(nInstances-i);
    j--;
  }

  for (int i = 0; i < nInstances; i++) {
    for (int j = i; j < nInstances; j++) {
      inner_prod(i,j);
    }
  }
}


void Fullproblem::inner_prod(int p, int q){
  H[p][q] = 0;
  for (int i = 0; i < nFeatures; i++) {
    H[p][q]+=data[p][i]*data[q][i];
  }
  if (p!=q) {
    H[p][q]*=y[p]*y[q];
  }

}
