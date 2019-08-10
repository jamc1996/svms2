#include "svm.h"
#include "io.h"
#include "fullproblem.h"
#include "subproblem.h"
#include "kernels.h"
#include "linked.h"
#include "algorithm.h"
#include "managetrade.h"

#include <sys/time.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>



int main(int argc, char *argv[]) {
  char* filename = NULL;
  struct denseData ds;
  struct Fullproblem fp;
  struct Projected sp;
	struct receiveData rd;

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
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Run\n");

	if(nProcGroups < nprocs){
		MPI_Comm_free(&mini_comm);
	}
	nProcGroups/=2;
	colour/=2;
	MPI_Comm_split(MPI_COMM_WORLD, colour, myid, &mini_comm);
	MPI_Comm_rank(mini_comm, &groupID);
	MPI_Comm_size(mini_comm, &groupSz);

 	struct Fullproblem nfp;
	struct yDenseData nds;
	struct Projected nsp;

	printf("Trading\n");
	tradeInfo(&rd, &ds, &nds, &fp, &nfp, mini_comm, groupSz, groupID, myid);
	printf("Traded\n");

	if(myid == 0){
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



	MPI_Win_fence(MPI_MODE_NOPRECEDE, dataWin);

	if(groupID != 0){
		MPI_Get(rd.data1d, rd.total*ds.nFeatures, MPI_DOUBLE, 0, 0 , 66, MPI_DOUBLE, dataWin);
	}

	MPI_Win_fence(0, dataWin);



	if(groupID == 0){
		run_Yserial_problem(&nds, &nfp, &nsp);
	}

	MPI_Barrier(mini_comm);
	MPI_Win_fence(0, dataWin);
	MPI_Win_fence(MPI_MODE_NOPRECEDE, alphaWin);
	MPI_Win_fence(MPI_MODE_NOPRECEDE, ytrWin);

	if(groupID !=0 ){
		MPI_Get(rd.alpha, rd.total, MPI_DOUBLE, 0, 0 ,rd.total, MPI_DOUBLE, alphaWin);
		MPI_Get(&(rd.ytr), 1, MPI_DOUBLE, 0, 0 ,1, MPI_DOUBLE, ytrWin);
	}

	MPI_Win_fence(0, alphaWin);
	MPI_Win_fence(0, ytrWin);

	if(myid == 1){
		calcW(&rd);
	}else{
		rootCalcW(&rd, &nds, &nfp);
		rd.ytr = nsp.ytr;
	}

	ReceiveCalcBeta(&fp, &rd, &ds);
	int* temp = malloc(sizeof(int)*5);
	find_n_worst(temp, 5, &fp);
	for(int i=0; i<5; i++){
		printf("%d\n",temp[i]);
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

	MPI_Barrier(MPI_COMM_WORLD);
	printf("myid is %d and I'm finalizing\n",myid);

	MPI_Win_fence(MPI_MODE_NOSUCCEED, dataWin);
	MPI_Win_fence(MPI_MODE_NOSUCCEED, alphaWin);
	MPI_Win_fence(MPI_MODE_NOSUCCEED, ytrWin);

	MPI_Win_free(&dataWin);
	MPI_Win_free(&alphaWin);
	MPI_Win_free(&ytrWin);
	MPI_Finalize();
  return 0;
}


























