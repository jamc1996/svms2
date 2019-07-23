#include "svm.h"
#include "io.h"
#include "fullproblem.h"
#include "subproblem.h"
#include "kernels.h"
#include "linked.h"

#include <sys/time.h>
#include <stdio.h>
#include <math.h>

#include <mpi.h>

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

void tradeInfo(struct denseData *ds, struct Fullproblem *fp, MPI_Comm mini_comm, int commSz, int commID);

int main(int argc, char *argv[]) {
  char* filename = NULL;
  struct denseData ds;
  struct Fullproblem fp;
  struct Projected sp;

  struct timeval start, trainStart, trainEnd, end;


  gettimeofday(&start, 0);

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
  int p = 4;
  if(parameters.test){
    testSavedModel(&ds, parameters.modelfile, &parameters);
    return 0;
  }

  gettimeofday(&trainStart, 0);
  // Full problem allocated and filled in, all alpha = 0.0 all gradF = 1.0:

	MPI_Barrier(MPI_COMM_WORLD);


  alloc_prob(&fp, &ds, p);

  init_prob(&fp, &ds);


  // Subproblem allocated:
  alloc_subprob(&sp, p);


  // We loop until no negative entries in beta:
  int k = 1;
  int max_iters = 10000000;
  int itt = 0;
  int n = 0;
	int nProcGroups = nprocs;
	int nStep = 0;
	int colour = myid;
	int groupID, groupSz;
	MPI_Comm mini_comm;
	while(nProcGroups > 0){
  	while(k){
			
    	// H matrix columns re-set and subproblem changed
    	if(itt%10000 == 0){
    	  printf("itt = %d\n",itt );
    	}

    	init_subprob(&sp, &fp, &ds, &parameters, 1);

    	//  congjugate gradient algorithm
    	//  if algorithm completes n == 0
    	//  if algorithm interrupt n != 0
    	n = cg(&sp, &fp);


    	updateAlphaR(&fp, &sp);
    	calcYTR(&sp, &fp);
    	calculateBeta(&fp, &sp, &ds);

    	if (n==0) {
      	printf("Converged! (itt = %d) %d\n", itt, fp.p );

      	int add = 2;
	      int *temp = malloc(sizeof(int)*add);
	      int *temp2 = malloc(sizeof(int)*add);
	      add = findWorstest(&fp,add,temp,temp2);

	      if (add == 0){
	        break;
	      }

	      changeP(&fp, &sp, add);
	      reinitprob(&ds, &fp, &sp, add, temp, temp2);
	      free(temp);
	      free(temp2);
	    }
	
	    if (n) {
	      // BCs broken, fix one at a time for the moment
	      k = singleswap(&ds, &fp, &sp, n, &parameters);
	      if (k < 0) {
	        shrinkSize(&fp, &sp, k+fp.p);
	      }
	      else{
	        n = checkfpConstraints(&fp);
	      }
	    }
												
	    itt++;
	    if(itt == max_iters){
	      printf("Reached max iters (%d)!!!!!\n\n\n",itt );
	      break;
	    }
										
	  }
		if(nProcGroups < nprocs){
			MPI_Comm_free(&mini_comm);
		}
		nProcGroups/=2;
		colour/=2;
		MPI_Comm_split(MPI_COMM_WORLD, colour, myid, &mini_comm);
		
		MPI_Comm_rank(mini_comm, &groupID);
		MPI_Comm_size(mini_comm, &groupSz);
		printf("myid is %d and mygroupID is %d\n",myid,groupID);
		tradeInfo(&ds, &fp, mini_comm, groupSz, groupID);

		nProcGroups/=2;

		MPI_Barrier(MPI_COMM_WORLD);
	}

  gettimeofday(&trainEnd, 0);
	MPI_Finalize();
	return 0;

  if (parameters.save) {
    saveTrainedModel( &fp, &ds, parameters.savename, &parameters);
  }

  gettimeofday(&end, 0);

  long totalelapsed = (end.tv_sec-start.tv_sec)*1000000 + end.tv_usec-start.tv_usec;
  long trainelapsed = (trainEnd.tv_sec-trainStart.tv_sec)*1000000 + trainEnd.tv_usec-trainStart.tv_usec;

  printf("Training Complete\n" );
  printf("Total Time spent: %ld micro seconds\n",totalelapsed );
  printf("Time spent training: %ld micro seconds\n",trainelapsed );

  //Memory freed
  freeDenseData(&ds);
  freeFullproblem(&fp);
  freeSubProblem(&sp);
  return 0;
}

void tradeInfo(struct denseData *ds, struct Fullproblem *fp, MPI_Comm mini_comm, int commSz, int commID)
{
	int p[commSz];

	int missed[commSz];
	int *missedInds;
	int *missedRecv;
	int missedTotal=0;;

	int j=0;
	int pTot =0;


	p[commID] = fp->p;	
	missed[commID] = 0;	
	for(int i=0; i<fp->n; i++){
		if(fabs(fp->alpha[i]-fp->C) < 0.00005){
			missed[commID]++;
		}
	}
	printf("%d is\n",missed[commID]);
	missedInds = malloc(sizeof(int)*missed[commID]);
	for(int i=0; i<fp->n; i++){
		if(fabs(fp->alpha[i]-fp->C) < 0.00005){
			missedInds[j]=ds->global[i];
		}
	}

	for(int i=0; i<commSz; i++){
		MPI_Bcast(&p[i],1,MPI_INT,i,mini_comm);
		MPI_Bcast(&missed[i],1,MPI_INT,i,mini_comm);
		if(i!=commID) {
			pTot+=p[i];
			missedTotal += missed[i];
		}
	}
	missedRecv = malloc(sizeof(int)*missedTotal);
	int k=0;
	for(int i=0; i<commSz; i++){
		MPI_Bcast(&missedIds[k],nMisC[i],MPI_INT,i,mini_comm);
		if(i!=commID){
			k+=missed[i];
		}
	}


	printf("myGrid is %d and commSz is %d\n",commID,commSz);
	return;
	for(int i=0; i<commSz; i++){
		printf("p is %d\n",p[i]);
		printf("nMisC is %d\n",nMisC[i]);
	}


	int *pI = malloc(sizeof(int)*pTot);
	double *alpO = malloc(sizeof(double)*pTot);

	k = 0;
	for(int i=0; i<commSz; i++){
		MPI_Bcast(&pI[k],fp->active[i],MPI_INT,i,mini_comm);
		if(i!=commID){
			k+=p[i];
		}
	}
	for(int i=0; i<commSz; i++){
//		MPI_Bcast(&pI[k],p[i],MPI_,i,mini_comm);
		if(i!=commID){
	//		k+=p[i];
		}
	}
}























