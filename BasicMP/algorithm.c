#include "algorithm.h"

void run_Yserial_problem(struct yDenseData *ds, struct Fullproblem *fp, struct Projected *sp){
  int k = 1;
  int max_iters = 10000000;
  int itt = 0;
  int n = 0;

  while(k){
			
   	// H matrix columns re-set and subproblem changed
   	if(itt%10000 == 0){
   	  printf("itt = %d\n",itt );
   	}

   	Yinit_subprob(sp, fp, ds, &parameters, 1);

   	//  congjugate gradient algorithm
   	//  if algorithm completes n == 0
   	//  if algorithm interrupt n != 0
   	n = cg(sp, fp);

   	updateAlphaR(fp, sp);
   	calcYTR(sp, fp);
   	YcalculateBeta(fp, sp, ds);

   	if (n==0) {
     	printf("Converged! (itt = %d) %d %d\n", itt, fp->p, fp->n );

     	int add = 2;
	    int *temp = malloc(sizeof(int)*add);
	    int *temp2 = malloc(sizeof(int)*add);
	    add = findWorstest(fp,add,temp,temp2);

	    if (add == 0){
	      break;
	    }

	    changeP(fp, sp, add);
	    Yreinitprob(ds, fp, sp, add, temp, temp2);
	    free(temp);
	    free(temp2);
	  }
	
	  if (n) {
	    // BCs broken, fix one at a time for the moment
	    k = Ysingleswap(ds, fp, sp, n, &parameters);
	    if (k < 0) {
	      shrinkSize(fp, sp, k+fp->p);
	    }
	    else{
	      n = checkfpConstraints(fp);
	    }
	  }
									
	  itt++;
	  if(itt == max_iters){
	    printf("Reached max iters (%d)!!!!!\n\n\n",itt );
	    break;
	  }										
	}
}

void run_serial_problem(struct denseData *ds, struct Fullproblem *fp, struct Projected *sp){
  int k = 1;
  int max_iters = 10000000;
  int itt = 0;
  int n = 0;

  while(k){
			
   	// H matrix columns re-set and subproblem changed
   	if(itt%10000 == 0){
   	  printf("itt = %d\n",itt );
   	}

   	init_subprob(sp, fp, ds, &parameters, 1);

   	//  congjugate gradient algorithm
   	//  if algorithm completes n == 0
   	//  if algorithm interrupt n != 0
   	n = cg(sp, fp);

   	updateAlphaR(fp, sp);
   	calcYTR(sp, fp);
   	calculateBeta(fp, sp, ds);

   	if (n==0) {
     	printf("Converged! (itt = %d) %d %d\n", itt, fp->p, fp->n );

     	int add = 2;
	    int *temp = malloc(sizeof(int)*add);
	    int *temp2 = malloc(sizeof(int)*add);
		printf("finding add\n");
	    add = findWorstest(fp,add,temp,temp2);

	    if (add == 0){
	      break;
	    }

	    changeP(fp, sp, add);
	    reinitprob(ds, fp, sp, add, temp, temp2);
		for(int i=0;i<fp->p; i++){
			printf("fp.active[%d] = %d\n",i,fp->active[i]);
		}
		for(int i=0;i<fp->q; i++){
			printf("fp.inactive[%d] = %d\n",i,fp->inactive[i]);
		}

	    free(temp);
	    free(temp2);
	  }
	
	  if (n) {
	    // BCs broken, fix one at a time for the moment
	    k = singleswap(ds, fp, sp, n, &parameters);
	    if (k < 0) {
	      shrinkSize(fp, sp, k+fp->p);
	    }
	    else{
	      n = checkfpConstraints(fp);
	    }
	  }
									
	  itt++;
	  if(itt == max_iters){
	    printf("Reached max iters (%d)!!!!!\n\n\n",itt );
	    break;
	  }										
	}
}

void rootCalcW(struct receiveData *rd, struct yDenseData *nds, struct Fullproblem *nfp){
	for(int i =0; i<nds->nInstances; i++){
		if(nfp->alpha[i] > 0.0){
			for(int j=0; j<nds->nFeatures; j++){
				rd->w[j] += (nds->y[i])*nfp->alpha[i]*nds->data[i][j];
			}
		}
	}
}


void calcW(struct receiveData *rd){
	for(int i =0; i<rd->total; i++){
		if(rd->alpha[i] > 0.0){
			for(int j=0; j<rd->nFeatures; j++){
				if(i < rd->nPos){
					rd->w[j] += rd->alpha[i]*rd->data[i][j];
				}else{
					rd->w[j] -= rd->alpha[i]*rd->data[i][j];
				}
			}
		}
	}
}

void ReceiveCalcBeta(struct Fullproblem *fp, struct receiveData *rd, struct denseData *ds){
	for(int i=0; i<fp->q; i++){
		fp->beta[i] = rd->ytr;
		for(int j=0; j<ds->nFeatures; j++){
			fp->beta[i] += rd->w[j]*ds->data[fp->inactive[i]][j];
		}
		printf("beta[%d] = %lf\n",i,fp->beta[i]);
	}
}

void find_n_worst(int *temp, int n, struct Fullproblem *fp){
  double *betaVal = malloc(sizeof(double)*n);
  for (int i = 0; i < n; i++) {
    betaVal[i] = DBL_MAX;
  }
  for (int i = 0; i < fp->q; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (fp->beta[i] < betaVal[j])
      {
        for (int k = n - 1; k > j ; k--)
        {
          temp[k] = temp[k-1];
          betaVal[k] = betaVal[k-1];
        }
        temp[j] = fp->inactive[i];
        betaVal[j] = fp->beta[i];
        break;
      }
    }
  }


  for (int i = 0; i < n; i++) {
    if (betaVal[i] > 0) {
      n = i;
    }
  }

  return n;
}





