#include "svm.h"
#include "io.h"
#include "fullproblem.h"
#include "subproblem.h"
#include "kernels.h"
#include "linked.h"
#include "algorithm.h"

#include <sys/time.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

void tradeInfo(struct receiveData *rd, struct denseData *ds, struct yDenseData *nds, struct Fullproblem *fp, struct Fullproblem *newfp, MPI_Comm mini_comm, int commSz, int commID, int myid);
