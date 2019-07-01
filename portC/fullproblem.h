#ifndef FULLPROBLEM_H
#define FULLPROBLEM_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include "svm.h"

void alloc_prob(struct Fullproblem *prob, struct denseData *ds, int p);
void init_prob(struct Fullproblem *prob, struct denseData *ds);
void updateAlphaR(struct Fullproblem *fp, struct Projected *sp);
void calculateBeta(struct Fullproblem *fp, struct Projected *sp, struct denseData *ds);
int swapMostNegative(struct Fullproblem *fp);
void setH(struct Fullproblem *prob, struct denseData *ds);
int singleswap(struct denseData *ds, struct Fullproblem *fp, struct Projected *sp, int n);
void adjustGradF(struct Fullproblem *fp, struct denseData *ds, int n, int worst, int signal, int target);
int checkfpConstraints(struct Fullproblem *fp);
void findWorst(int *worst, int* target, int* change, int *n, struct denseData *ds, struct Fullproblem *fp);

#endif
