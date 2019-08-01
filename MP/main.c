#include "svm.h"
#include "io.h"
#include "fullproblem.h"
#include "subproblem.h"
#include "kernels.h"
#include "linked.h"

#include <sys/time.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

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
	    add = findWorstest(fp,add,temp,temp2);

	    if (add == 0){
	      break;
	    }

	    changeP(fp, sp, add);
	    reinitprob(ds, fp, sp, add, temp, temp2);
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

void tradeInfo(struct receiveData *rd, struct denseData *ds, struct denseData *nds, struct Fullproblem *fp, struct Fullproblem *newfp, MPI_Comm mini_comm, int commSz, int commID, int myid);

int main(int argc, char *argv[]) {
  char* filename = NULL;
  struct denseData ds;
  struct Fullproblem fp;
  struct Projected sp;
	struct receiveData rd;
//  struct timeval start, trainStart, trainEnd, end;

  //gettimeofday(&start, 0);

	int nprocs = 1, myid = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // Input processed:
  parse_arguments(argc, argv, &filename, &parameters);
  read_file(filename, &ds, nprocs, myid);
  //preprocess(&ds);

  //cleanData(&ds);

  // Projected problem size chosen temporarily
  int p = 6;
  if(parameters.test){
    testSavedModel( &ds, parameters.modelfile, &parameters);
    return 0;
  }

  //gettimeofday(&trainStart, 0);
  // Full problem allocated and filled in, all alpha = 0.0 all gradF = 1.0:

  alloc_prob(&fp, &ds, p);
  init_prob(&fp, &ds);

  // Subproblem allocated:
  alloc_subprob(&sp, p);

	int nProcGroups = nprocs;
	int nStep = 0;
	int colour = myid;
	int groupID, groupSz;
	MPI_Comm mini_comm;

  // We loop until no negative entries in beta:
	run_serial_problem(&ds, &fp, &sp); 
	if(nProcGroups < nprocs){
		//MPI_Comm_free(&mini_comm);
	}
	nProcGroups/=2;
	colour/=2;
	MPI_Comm_split(MPI_COMM_WORLD, colour, myid, &mini_comm);
		
	MPI_Comm_rank(mini_comm, &groupID);
	MPI_Comm_size(mini_comm, &groupSz);
 	struct Fullproblem nfp;
	struct denseData nds;
	struct Projected nsp;
	tradeInfo(&rd, &ds, &nds, &fp, &nfp, mini_comm, groupSz, groupID, myid);
	MPI_Barrier(MPI_COMM_WORLD);

	if(myid == 0){
			printf("ndsF is %d\n",nds.nFeatures);
			printf("nds ==> nInst %d procInst %d nPos %d nNeg %d procPost %d procNeg %d\n", nds.nInstances,nds.procInstances, nds.nPos, nds.nNeg, nds.procPos, nds.procNeg);
		
			printf("nfp ==> n %d, p %d, q %d C %lf\n",nfp.n,nfp.p,nfp.q,nfp.C);		
		for(int i=0; i<nfp.p; i++){
			printf("alpha[%d] = %lf ----- gradF[%d] = %lf\n",i,nfp.alpha[i],i, nfp.gradF[i]);
			printf("active[%d] = %d\n",i,nfp.active[i]);
		}
		for(int i=0; i<nds.nInstances; i++){
			printf("%d\t\t",i);
			for(int j=0; j<nds.nFeatures; j++){
				printf("%lf\t",nds.data[i][j]);
			}
			printf("\n");
		}
	  alloc_subprob(&nsp, nfp.p);

	}

	MPI_Win dataWin, alphaWin, ytrWin;

	MPI_Barrier(MPI_COMM_WORLD);
	if(myid == 0){
		MPI_Win_create(nds.data1d, nds.nInstances*nds.nFeatures*sizeof(double), sizeof(double), MPI_INFO_NULL, mini_comm, &dataWin);
		MPI_Win_create(nfp.alpha, nfp.n*sizeof(double), sizeof(double), MPI_INFO_NULL, mini_comm, &alphaWin);
		MPI_Win_create(&(nsp.ytr), sizeof(double), sizeof(double), MPI_INFO_NULL, mini_comm, &ytrWin);
	}else{
		MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, mini_comm, &dataWin);
		MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, mini_comm, &alphaWin);
		MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, mini_comm, &ytrWin);
	}

	MPI_Barrier(mini_comm);

	MPI_Win_fence(MPI_MODE_NOPRECEDE, dataWin);
	MPI_Win_fence(MPI_MODE_NOPRECEDE, alphaWin);
	MPI_Win_fence(MPI_MODE_NOPRECEDE, ytrWin);

	if(groupID != 0){
		printf(" count is HEHAIOEHATOHFRT %d\n",rd.total*ds.nFeatures);
		MPI_Get(rd.data1d, rd.total*ds.nFeatures, MPI_DOUBLE, 0, 0 ,rd.total*nds.nFeatures, MPI_DOUBLE, dataWin);
		printf(" got your get\n");
	}

printf("%d is this fence?\n",groupID);

	MPI_Win_fence(0, dataWin);
		printf("%d is this fence?\n",groupID);
	if(groupID == 0){
		run_serial_problem(&nds, &nfp, &nsp);
		for(int i=0; i<fp.n; i++){
			printf("%lf\n",nfp.alpha[i]);
		}
		printf("myid is %d and ytr is %lf\n",myid,nfp.alpha[0]);
	}

	MPI_Barrier(mini_comm);
	printf("myid is %d \n",myid);

//  	if (parameters.save) {
    	//saveTrainedModel( &nfp, &nds, parameters.savename, &parameters);
  //	}

	if(groupID !=0 ){
		printf(" count is HEHAIOEHATOHFRT %d\n",rd.total);
		MPI_Get(rd.alpha, rd.total, MPI_DOUBLE, 0, 0 ,rd.total, MPI_DOUBLE, alphaWin);
		printf(" cgettin\n");
		MPI_Get(&(rd.ytr), 1, MPI_DOUBLE, 0, 0 ,1, MPI_DOUBLE, ytrWin);

	}

	MPI_Win_fence(0, alphaWin);
	MPI_Win_fence(0, ytrWin);
	if(myid == 1){
	printf("WHy?\n");
	printf("myid is %d and ytr is %lf\n",myid,rd.alpha[0]);
  ///gettimeofday(&trainEnd, 0);
}

  //gettimeofday(&end, 0);

  //long totalelapsed = (end.tv_sec-start.tv_sec)*1000000 + end.tv_usec-start.tv_usec;
  //long trainelapsed = (trainEnd.tv_sec-trainStart.tv_sec)*1000000 + trainEnd.tv_usec-trainStart.tv_usec;

  printf("Training Complete\n" );
  //printf("Total Time spent: %ld micro seconds\n",totalelapsed );
  //printf("Time spent training: %ld micro seconds\n",trainelapsed );

  //Memory freed
	printf("free dens %d\n",myid);
  freeDenseData(&ds);
	printf("free full %d\n", myid);
  freeFullproblem(&fp);
	printf("free sub %d\n",myid);
  freeSubProblem(&sp);

	MPI_Win_fence(MPI_MODE_NOSUCCEED, dataWin);
	MPI_Win_fence(MPI_MODE_NOSUCCEED, alphaWin);
	MPI_Win_fence(MPI_MODE_NOSUCCEED, ytrWin);

	MPI_Win_free(&dataWin);
	MPI_Win_free(&alphaWin);
	MPI_Win_free(&ytrWin);
	MPI_Finalize();
  return 0;
}



void tradeInfo(struct receiveData *rd, struct denseData *ds, struct denseData *nds, struct Fullproblem *fp, struct Fullproblem *nfp, MPI_Comm mini_comm, int commSz, int commID, int myid)
{
	int my_missed = 0;
	int *my_missedInds;	
	int *myIndex;
	for(int i=0; i<fp->n; i++){
		if(fabs(fp->alpha[i]-fp->C) < 0.00005){
			my_missed++;
		}
	}
		my_missedInds = malloc(sizeof(int)*my_missed);
	if(my_missed > 0){
	

		int my_missedMark = ds->procInstances;
		int j;
		for(int i=0; i<fp->n; i++){
			if(fabs(fp->alpha[i]-fp->C) < 0.00005){
				if(i >= ds->procPos && my_missedMark == ds->procInstances){
					my_missedMark = j;
				}
				my_missedInds[j] = i;
				j++;
			}
		}
	}
	int my_p = fp->p;
		myIndex = malloc(sizeof(int)*fp->p);
	if(my_p > 0){

		for(int i=0; i<fp->p; i++){
			myIndex[i] = fp->active[i];
		}
	}
	double* alp;
	int tMissed = my_missed;
	int tP = my_p;
	printf("myid is %d and tmissed %d and tp %d\n",commID, tMissed, tP);
fflush(stdout);
MPI_Barrier(mini_comm);
	MPI_Allreduce( &(my_missed), &tMissed, 1, MPI_INT, MPI_SUM, mini_comm);
	MPI_Allreduce( &(my_p), &tP, 1, MPI_INT, MPI_SUM, mini_comm);
	rd->total = tMissed+tP;

printf("Reduction ok %d\n",commID);
	int *otherP = malloc(sizeof(int)*commSz);
	otherP[commID] = fp->p + my_missed;
	if(commID == 0){
		if(tMissed ==0){
			nfp->inactive = NULL;
			nfp->beta = NULL;
		}
		nds->nFeatures = ds->nFeatures;
		nds->nInstances = tMissed + tP;

		nds->data1d = malloc(sizeof(double)*nds->nFeatures*nds->nInstances);
		nds->data = malloc(sizeof(double*)*nds->nInstances);
		for(int i=0; i< nds->nInstances; i++){
			nds->data[i] = &nds->data1d[i*ds->nFeatures];
		}
		nfp->alpha = malloc(sizeof(double)*nds->nInstances);
		nfp->gradF = malloc(sizeof(double)*nds->nInstances);
		nfp->n = nds->nInstances;
		nfp->p = tP;
		nfp->q = tMissed;
		nfp->C = 2*fp->C;
		nfp->active = malloc(sizeof(int)*nfp->p);
		nfp->inactive = malloc(sizeof(int)*nfp->q);
		nfp->beta = malloc(sizeof(double)*nfp->q);
printf("%d\n",nfp->n);
	}
	else{
		alp = malloc(sizeof(double)*(fp->p+my_missed));
		for(int i=0; i<fp->p; i++){
			if(fp->active[i] < ds->procPos){
				alp[i] = fp->alpha[fp->active[i]];
			}
			else{
				alp[i] = -fp->alpha[fp->active[i]];
			}
		}
		for(int i=fp->p; i<fp->p+my_missed; i++){
			if(my_missedInds[i - fp->p] < ds->procPos){
				alp[i] = fp->C;
			}
			else{
				alp[i] = -fp->C;
			}
		}
	}


for(int i=0; i<nfp->n; i++){

}
	int q = otherP[commID];
	for(int i = 1; i<commSz; i++){
		if(i == commID){
			MPI_Send(&otherP[i], 1, MPI_INT, 0, i, mini_comm);
			printf("%lf is alp\n",alp[0]);
			MPI_Send(alp, otherP[i], MPI_DOUBLE, 0, i, mini_comm);
		}else if(commID == 0){
			MPI_Recv(&otherP[i] , 1 , MPI_INT, i, i, mini_comm, MPI_STATUS_IGNORE);
			printf("q is %d otherp is %d\n",q,otherP[i]);
			MPI_Recv(&(nfp->gradF[q]) , otherP[i], MPI_DOUBLE, i, i, mini_comm, MPI_STATUS_IGNORE);
			printf("after rec %lf\n",nfp->gradF[q]);
			q+=otherP[i];
		}
	}
	if(commID == 0){
		nds->procPos = 0;

		for(int i=0; i<fp->p; i++){
			nfp->gradF[i] = -fp->alpha[fp->active[i]];
			if(fp->active[i] < ds->procPos){
				nfp->gradF[i] = -nfp->gradF[i];
			}
		}
	printf("ok %d\n",my_missed);
		for(int i=0; i<my_missed; i++){
	printf("i = %d\n",fp->p+i);
			nfp->gradF[fp->p + i] = fp->C;
	printf("ok\n");
			if(my_missedInds[i] < ds->procPos){
	printf("ok1\n");
				nfp->gradF[fp->p + i] = -nfp->gradF[fp->p + i];
	printf("ok2\n");
			}
	printf("ok3\n");
		}

		for(int i = 0; i< nfp->n; i++){

			if(nfp->gradF[i] > 0){
				nds->procPos++;
			}
			printf("%d is procPos and i is %d %lf\n",nds->procPos,i,nfp->gradF[i]);
		}
		int tot = 0, pos = 0, neg =0;

		for(int id = 0; id< commSz; id++){
			printf("%d id \n",id);
			printf("%d is opid\n",otherP[id]);
			for(int j=0; j<otherP[id]; j++){
				printf("j is %d\n",j);
				if(id == 0){
					if(nfp->gradF[tot] > 0){
						nfp->alpha[pos] = nfp->gradF[tot];
						printf("alpah %d = %lf\n",pos,nfp->alpha[pos]);
						for(int k =0; k< ds->nFeatures; k++){
							nds->data[pos][k] = ds->data[fp->active[j]][k];
						}
						pos++;
					}else{
						nfp->alpha[nds->procPos + neg] = -nfp->gradF[tot];
						printf("alpah %d = %lf\n",nds->procPos + neg,nfp->alpha[nds->procPos + neg]);
						for(int k =0; k< ds->nFeatures; k++){
							nds->data[nds->procPos + neg][k] = ds->data[fp->active[j]][k];
						}
						neg++;
					}
					tot++;
				}else{
						printf("tot = %d\n",tot);
					if(nfp->gradF[tot] > 0){
						nfp->alpha[pos] = nfp->gradF[tot];
						printf("%d is pos\n",pos);
						MPI_Recv((nds->data[pos]), ds->nFeatures, MPI_DOUBLE, id, j, mini_comm, MPI_STATUS_IGNORE);
						pos++;
					}else{
						nfp->alpha[nds->procPos + neg] = -nfp->gradF[tot];
						printf("%d is neg\n",nds->procPos + neg);
						MPI_Recv((nds->data[nds->procPos + neg]), ds->nFeatures, MPI_DOUBLE, id, j, mini_comm, MPI_STATUS_IGNORE);
						neg++;
					}
					tot++;
				}
			}
		}
		printf("ok\n");
	}else{
		printf("myid %d got here\n",commID);
		for(int j=0; j<fp->p; j++){
			MPI_Send(ds->data[fp->active[j]], ds->nFeatures, MPI_DOUBLE, 0, j, mini_comm );
		}
		for(int j = fp->p ; j<otherP[commID] ; j++ ){
			MPI_Send(ds->data[my_missedInds[j - fp->p]], ds->nFeatures, MPI_DOUBLE, 0, j, mini_comm);
		}
	}


	if(commID == 0){
		int act=0, inact = 0;
		for(int i=0; i<nds->nInstances; i++){
			if(fabs(nfp->alpha[i] - nfp->C)< 0.0005){
				nfp->inactive[inact] = i;
				inact++;
			}
			else{
				nfp->active[act] = i;
				act++;
			}
		}
	}

	if(commID ==0){
		for(int i=0; i<nds->nInstances; i++){
			for(int j=0; j<nds->nFeatures; j++){
				printf("%lf\t",nds->data[i][j]);
			}
			printf("\n");
		}
	}
	if(commID == 0){
		double* fulH = malloc(sizeof(double)*ds->nInstances*ds->nInstances);
		double** hH = malloc(sizeof(double*)*ds->nInstances);

		for(int i=0; i<ds->nInstances ; i++){
			hH[i] = &fulH[i*ds->nInstances];
		}
	
		for(int i=0; i<nds->nInstances ; i++){
			for(int j=0; j<nds->nInstances ; j++){
				hH[i][j] = 0.0;
				for(int k=0; k<nds->nInstances ; k++){
					hH[i][j] += nds->data[i][k]*nds->data[j][k];
				}
				if( ( i < nds->procPos) ^ (j <nds->procPos) ){
					hH[i][j] = -hH[i][j];
				}
			}
		}

		for(int i=0; i<nds->nInstances ; i++){
			nfp->gradF[i] = 1.0;			
			for(int j=0; j<nds->nInstances ; j++){
				nfp->gradF[i] -= hH[i][j]*nfp->alpha[j];
			}
			printf("%lf\n",nfp->alpha[i]);
			printf("%lf\n",nfp->gradF[i]);
		}
		nds->procInstances = nds->nInstances;
		nds->nPos = nds->procPos;
		nds->nNeg = nds->nInstances - nds->nPos;
		nds->procNeg = nds->nInstances - nds->nPos;
	}

	if(myid == 0){
		nfp->partialH = Init_Empty_List();
	  for (int i = 0; i < nfp->p; i++) {
  	  nfp->partialH = append(nds,nfp->partialH,nfp->active[i]);
		}
	}else{
		printf("%d is total\n",ds->nFeatures*rd->total);
		rd->data1d = malloc(sizeof(double)*ds->nFeatures*rd->total);
		rd->data = malloc(sizeof(double*)*rd->total);
		for(int i=0; i<rd->total; i++){
			rd->data[i] = &(rd->data1d[i*ds->nFeatures]);
		}

		rd->alpha = malloc(sizeof(double)*rd->total);

	}

	free(my_missedInds);	
	free(myIndex);
}






















