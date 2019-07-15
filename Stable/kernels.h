#ifndef KERNELS_H
#define KERNELS_H

#include <stdlib.h>
#include <math.h>


#include "svm.h"


int setH(struct Fullproblem *prob, struct denseData *ds, struct svm_args *params);
int updateSubH(struct Fullproblem *fp, struct Projected *sp, struct denseData *ds, struct svm_args *params);
void partialHupdate(struct Fullproblem *fp, struct Projected *sp, struct denseData *ds, struct svm_args *params, int n, int worst);


#endif
