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

	tradeInfo(&rd, &ds, &nds, &fp, &nfp, mini_comm, groupSz, groupID, myid);

	if(myid == 0){
	  alloc_subprob(&nsp, nfp.p);
	}

	MPI_Win dataWin, alphaWin, ytrWin, gradWin, yWin;
	MPI_Barrier(MPI_COMM_WORLD);

	if(myid == 0){
		MPI_Win_create(nds.data1d, 5*nds.nInstances*nds.nFeatures*sizeof(double), sizeof(double), MPI_INFO_NULL, mini_comm, &dataWin);
		MPI_Win_create(nfp.alpha, 5*nfp.n*sizeof(double), sizeof(double), MPI_INFO_NULL, mini_comm, &alphaWin);
		MPI_Win_create(nfp.gradF, 5*nfp.n*sizeof(double), sizeof(double), MPI_INFO_NULL, mini_comm, &gradWin);
		MPI_Win_create(nds.y, 5*nds.nInstances*sizeof(int), sizeof(int), MPI_INFO_NULL, mini_comm, &yWin);
		MPI_Win_create(&(nsp.ytr), sizeof(double), sizeof(double), MPI_INFO_NULL, mini_comm, &ytrWin);
	}else{
		MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, mini_comm, &dataWin);
		MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, mini_comm, &alphaWin);
		MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, mini_comm, &gradWin);
		MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, mini_comm, &yWin);
		MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, mini_comm, &ytrWin);
	}
	

	MPI_Win_fence(MPI_MODE_NOPRECEDE, dataWin);

	if(groupID != 0){
		MPI_Get(rd.data1d, rd.total*ds.nFeatures, MPI_DOUBLE, 0, 0 , rd.total*ds.nFeatures, MPI_DOUBLE, dataWin);
	}
	MPI_Win_fence(0, dataWin);
	
	
	if(groupID == 0){
		for (int i =0; i<nds.nInstances; i++){
			printf("i = %d y[i] = %d gradf =  %lf,\n\n alpha = %lf\n\n",i, nds.y[i], nfp.gradF[i], nfp.alpha[i]);
			for(int j=0; j<nds.nFeatures; j++){
				printf("%lf\t",nds.data[i][j] );
			}
			printf("\n");
		}
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
	int sending = 5;
	int* temp = malloc(sizeof(int)*sending);
	find_n_worst(temp, sending, &fp);
	MPI_Win_fence(MPI_MODE_NOPRECEDE, yWin);
	MPI_Win_fence(MPI_MODE_NOPRECEDE, gradWin);
	
	for (int i=0; i < sending; i++){
		MPI_Put(ds.data[temp[i]], ds.nFeatures, MPI_DOUBLE, 0, ds.nFeatures*(rd.total+(myid*sending)+i), ds.nFeatures, MPI_DOUBLE, dataWin);
		MPI_Put(&fp.gradF[temp[i]], 1, MPI_DOUBLE, 0, rd.total+(myid*sending)+i, 1, MPI_DOUBLE, gradWin);
		int y_send = -1;
		if(temp[i] < ds.procPos){
			y_send = 1;
		}
		MPI_Put(&y_send, 1, MPI_INT, 0, rd.total+(myid*sending)+i, 1, MPI_INT, yWin);
	}
	
	if(myid == 0){
		nfp.inactive = realloc(nfp.inactive, sizeof(int)*(nfp.q+(nprocs*sending))  );
		nfp.beta = realloc(nfp.beta,sizeof(double)*(nfp.q+(nprocs*sending)));
		for( int i=0; i<nprocs*sending; i++){
			nfp.alpha[i+nds.nInstances] = 0.0;
			nfp.inactive[i+nfp.q] = i + nds.nInstances;
		}
		nfp.n += (nprocs*sending);
		nfp.q += (nprocs*sending);
		printf("%d\n",nfp.p);
		Cell *temp = nfp.partialH.head;
		while(temp!= NULL){
			printf("lable is %d\n",temp->label);
			temp = temp->next;
		}

		printf("That's all\n");
		temp = nfp.partialH.head;
		while(temp != NULL){
			printf("reallocing with label is %d\n",temp->label);
			temp->line = realloc(temp->line, sizeof(double)*(nfp.n));
			temp = temp->next;
		}
		printf("should be ok now\n");
	}
	
	MPI_Win_fence(0, yWin);
	MPI_Win_fence(0, gradWin);
	MPI_Win_fence(0, dataWin);
	MPI_Barrier(MPI_COMM_WORLD);

	rd.total += nprocs*sending;
	if(myid == 0){
		Cell* temp = nfp.partialH.head;
		while (temp!=NULL){
			for (int i=nds.nInstances; i<nfp.n; i++){
				temp->line[i] = 0.0;
				for(int j = 0; j<nds.nFeatures; j++){
					temp->line[i] += nds.data[temp->label][j]*nds.data[i][j];
				}
				temp->line[i] *= nds.y[temp->label]*nds.y[i];
			}
			temp = temp->next;
		}

		nds.nInstances += nprocs*sending;
		for (int i =0; i<nds.nInstances; i++){
			printf("i = %d y[i] = %d gradf =  %lf, \n\nalpha = %lf\n\n",i, nds.y[i], nfp.gradF[i], nfp.alpha[i]);
			for(int j=0; j<nds.nFeatures; j++){
				printf("%lf\t",nds.data[i][j] );
			}
			printf("\n");
		}
}


	if(groupID == 0){
		//run_Yserial_problem(&nds, &nfp, &nsp);
	}




  freeDenseData(&ds);
  freeFullproblem(&fp);
  freeSubProblem(&sp);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Win_fence(MPI_MODE_NOSUCCEED, dataWin);
	MPI_Win_fence(MPI_MODE_NOSUCCEED, alphaWin);
	MPI_Win_fence(MPI_MODE_NOSUCCEED, ytrWin);
	MPI_Win_fence(MPI_MODE_NOSUCCEED, gradWin);
	MPI_Win_fence(MPI_MODE_NOSUCCEED, yWin);

	MPI_Win_free(&yWin);
	MPI_Win_free(&gradWin);
	MPI_Win_free(&dataWin);
	MPI_Win_free(&alphaWin);
	MPI_Win_free(&ytrWin);
	MPI_Finalize();
  return 0;
}


























