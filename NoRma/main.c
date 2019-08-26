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

void run_parallel_algorithm( struct receiveData *rd, struct denseData *ds, struct Fullproblem *fp, struct yDenseData *nds, struct Fullproblem *nfp, struct Projected *nsp, int myid, int nprocs, MPI_Comm Comm, int sending );



void update_root_nfp(struct yDenseData *nds, struct Fullproblem *nfp, int nprocs, int global_n ){
	nfp->inactive = realloc(nfp->inactive, sizeof(int)*(nfp->q+global_n)  );
	nfp->beta = realloc(nfp->beta,sizeof(double)*(nfp->q+global_n));
	for( int i=0; i<global_n; i++){
		nfp->alpha[i+nds->nInstances] = 0.0;
		nfp->inactive[i+nfp->q] = i + nds->nInstances;
	}
	nfp->n += global_n;
	nfp->q += global_n;
	Cell *temp = nfp->partialH.head;

	temp = nfp->partialH.head;
	while(temp != NULL){
		temp->line = realloc(temp->line, sizeof(double)*(nfp->n));
		temp = temp->next;
	}
}

void updatePartialH(struct yDenseData *nds, struct Fullproblem *nfp, int global_n){
	Cell* temp = nfp->partialH.head;
	while (temp!=NULL){
		for (int i=nds->nInstances; i<nfp->n; i++){
			temp->line[i] = 0.0;
			for(int j = 0; j<nds->nFeatures; j++){
				temp->line[i] += nds->data[temp->label][j]*nds->data[i][j];
			}
			temp->line[i] *= nds->y[temp->label]*nds->y[i];
		}
		temp = temp->next;
	}

	nds->nInstances += global_n;
}




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
	int colour = myid;
	int groupID, groupSz;
	MPI_Comm mini_comm;

  // We loop until no negative entries in beta:
//	run_serial_problem(&ds, &fp, &sp); 
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

	tradeInfo(&rd, &ds, &nds, &fp, &nfp, nprocs, myid, MPI_COMM_WORLD);

	if(myid == 0){
	  alloc_subprob(&nsp, nfp.p);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	int sending = 5;	
	run_parallel_algorithm( &rd, &ds, &fp, &nds, &nfp, &nsp, myid, nprocs, MPI_COMM_WORLD, sending );



	freeDenseData(&ds);
	freeFullproblem(&fp);
	freeSubProblem(&sp);
	MPI_Finalize();

	return 0;
}

void run_parallel_algorithm(struct receiveData *rd, struct denseData *ds, struct Fullproblem *fp, struct yDenseData *nds, struct Fullproblem *nfp, struct Projected *nsp, int myid, int nprocs, MPI_Comm Comm, int sending )
{
	int *temp = malloc(sizeof(int)*sending);
	int *local_n = malloc(sizeof(int)*nprocs);
	int itt = 0;
	while ( 1 ) 
	{
		
		itt++;
		if(myid != 0)
		{
			MPI_Recv(rd->data1d, rd->total*ds->nFeatures, MPI_DOUBLE, 0, myid, Comm, MPI_STATUS_IGNORE);
			MPI_Recv(rd->y, rd->total, MPI_INT, 0, myid, Comm, MPI_STATUS_IGNORE);
		}
		if(myid == 0){
			for(int i=1; i<nprocs; i++){
				MPI_Send(nds->data1d, nds->nInstances*nds->nFeatures, MPI_DOUBLE, i, i, Comm);
				MPI_Send(nds->y, nds->nInstances, MPI_INT, i, i, Comm);
			}
			run_Yserial_problem(nds, nfp, nsp);
		}

		if(myid !=0 ){
			MPI_Recv(rd->alpha, rd->total, MPI_DOUBLE, 0, myid, Comm, MPI_STATUS_IGNORE);
			MPI_Recv(&(rd->ytr), 1, MPI_DOUBLE, 0, myid, Comm, MPI_STATUS_IGNORE);
		}else{
			for(int i=0; i< nprocs; i++){
				MPI_Send(nfp->alpha, nfp->n, MPI_DOUBLE, i, i, Comm);
				MPI_Send(&(nsp->ytr), 1, MPI_DOUBLE, i, i, Comm);
			}
		}
	
		if(myid == 1){
			calcW(rd);
		}else{
			rootCalcW(rd, nds, nfp);
			rd->ytr = nsp->ytr;
		}

		ReceiveCalcBeta(fp, rd, ds);

		local_n[myid] = find_n_worst(temp, sending, fp);
		int global_n = 0;
		for(int i=0; i<nprocs; i++){
			MPI_Bcast(&local_n[i], 1, MPI_INT, i, Comm);
			global_n += local_n[i];
		}
		int local_start = 0;
		for(int i=0; i<myid; i++){
			local_start += local_n[i];
		}
		if(global_n == 0){
			break;
		}

		if(myid != 0){		
			for (int i=0; i < local_n[myid]; i++){
				MPI_Send(ds->data[temp[i]], ds->nFeatures, MPI_DOUBLE, 0, myid, Comm);
				MPI_Send(&(fp->gradF[temp[i]]), 1, MPI_DOUBLE, 0, myid, Comm);
				int y_send = -1;
				if(temp[i] < ds->procPos){
					y_send = 1;
				}
				MPI_Send(&(y_send), 1, MPI_INT, 0, myid, Comm);
			}
		}else{
			for(int i=0; i< local_n[0]; i++){
				for(int j=0; j<nds->nFeatures; j++){
					nds->data[rd->total+i][j]=ds->data[temp[i]][j];
				}
				nfp->gradF[rd->total+i] = fp->gradF[temp[i]];
				if(temp[i] < ds->procPos){
					nds->y[rd->total+i] = 1;
				}else{
					nds->y[rd->total+i] = -1;
				}
			}
			local_start+= local_n[0];
			for (int i=1 ; i<nprocs; i++){
				for (int j=0; j<local_n[i]; j++){
					MPI_Recv(nds->data[rd->total+local_start+j], ds->nFeatures, MPI_DOUBLE, i, i, Comm, MPI_STATUS_IGNORE);
					MPI_Recv(&(nfp->gradF[rd->total+local_start+j]), 1, MPI_DOUBLE, i, i, Comm, MPI_STATUS_IGNORE);
					MPI_Recv(&(nds->y[rd->total+local_start+j]), 1, MPI_DOUBLE, i, i, Comm, MPI_STATUS_IGNORE);
			}
				local_start+= local_n[i];
			}
		}
		MPI_Barrier(Comm);
		if(myid == 0){
			update_root_nfp(nds, nfp, nprocs,global_n);
		}
	
		rd->total += global_n;

		if(myid == 0){
			updatePartialH(nds, nfp, global_n);
		}

	}


	free(temp);
	free(local_n);
}























