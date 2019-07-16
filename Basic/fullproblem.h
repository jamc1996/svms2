#ifndef FULLPROBLEM_H
#define FULLPROBLEM_H

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "svm.h"


class Fullproblem
{
private:
  int nInstances;
  int nFeatures;
  double** data;
  double* data1d;
  char** instanceLabels;
  char** featureLabels;


public:
  double* h;
  double** H;
  double* y;
  double* alpha;
  double* gradF;
  Fullproblem(struct denseData ds);
  void fillH();
  void inner_prod(int p, int q);
  double getHij(int i, int j);
  void calcStdDev(double* stdDev, double* mean);
  void calcMeans(double *mean);
  void preprocess();
  void normalise(double* mean, double* stdDev);

};


#endif
