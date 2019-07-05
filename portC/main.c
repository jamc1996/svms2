#include "svm.h"
#include "io.h"
#include "fullproblem.h"
#include "subproblem.h"

#include <stdio.h>

void freeDenseData(struct denseData *ds);
void  freeFullproblem(struct Fullproblem *fp);
void   freeSubProblem( struct Projected* sp);
void changeP( struct Fullproblem *fp, struct Projected *sp, int add);
int findWorstest(struct Fullproblem *fp , int add, int* temp, int* temp2);
void reinitprob( struct Fullproblem *fp, struct Projected *sp, int add, int* temp, int* temp2);
void shrinkSize( struct Fullproblem *fp, struct Projected *sp, int k);
void cleanData( struct denseData *ds){
  for (int i = 0; i < ds->nInstances - 1; i++) {
    printf("%d\n",i );
    for (int j = i; j < ds->nInstances; j++) {
      int flag = 1;
      double check;
      double previous = ds->data[i][0]/ds->data[j][0];
      for (int k = 1; k < ds->nFeatures; k++) {
        check = ds->data[i][k]/ds->data[j][k];
        printf("%lf and %lf\n",previous,check );
        if (fabs(check - previous) > 0.0001 ) {
          flag = 0;
          break;
        }
        previous = check;
      }
      if (flag == 1) {
        printf("%d is broken\n",j );
      }
    }
  }
}
int main(int argc, char *argv[]) {
  char* filename = NULL;
  struct svm_args parameters;
  struct denseData ds;
  struct Fullproblem fp;
  struct Projected sp;

  // Input processed:
  parse_arguments(argc, argv, &filename, &parameters);
  read_file(filename, &ds);
  preprocess(&ds);

  double* bigH = malloc(sizeof(double)*ds.nInstances*ds.nInstances);
  double** fullBigH = malloc(sizeof(double*)*ds.nInstances);
  for (int i = 0; i < ds.nInstances; i++) {
    fullBigH[i] = &bigH[i*ds.nInstances];
    for (int j = 0; j < ds.nInstances; j++) {
      fullBigH[i][j] = 0.0;
      for (int k = 0; k < ds.nFeatures; k++) {
        fullBigH[i][j] += ds.data[i][k]*ds.data[j][k];
      }
      fullBigH[i][j] *= ds.y[j]*ds.y[i];
    }
  }
  double* check = malloc(sizeof(double)*ds.nInstances);
  //cleanData(&ds);

  // Projected problem size chosen temporarily
  int p = 4;

  if(parameters.test){
    p = 3;
    alloc_prob(&fp, &ds, p);
    init_prob(&fp, &ds);
    alloc_subprob(&sp, p, &fp, &ds);
    int j = 0;
    for (int i = 0; i < fp.n; i++) {
      if (i != 0 && i != 25 && i != 29) {
        fp.inactive[j] = i;
        j++;
      }
    }
    fp.active[0] = 0;
    fp.active[1] = 25;
    fp.active[2] = 29;
    fp.alpha[0] = 1.992153;
    fp.alpha[25] = 1.640305;
    fp.alpha[29] = 0.351848;

    setH(&fp, &ds);
    init_subprob(&sp, &fp, &ds);
    for (int i = 0; i < fp.n; i++) {
      printf("alpha[%d] = %lf\n",i,fp.alpha[i] );
    }

    double b = 1.0;
    double w[6] = {0};
    for (int i = 0; i < fp.p; i++) {
      b -= sp.H[0][i]*fp.alpha[fp.active[i]];
    }
    b*=ds.y[fp.active[0]];
    for (int i = 0; i < ds.nFeatures; i++) {
      w[i] = 0.0;
      for (int j = 0; j < fp.p; j++) {
        w[i] += fp.alpha[fp.active[j]]*ds.data[fp.active[j]][i]*ds.y[fp.active[j]];
      }
      printf("w[%d] = %lf\n",i,w[i] );
    }
    printf("b = %lf\n",b );
    for (int i = 0; i < fp.n; i++) {
      double res = b;
      for (int j = 0; j < ds.nFeatures; j++) {
        res += w[j]*ds.data[i][j];
      }
      printf("res[%d] = %lf\n",i,res );
    }
    return 0;
  }

  // Full problem allocated and filled in, all alpha = 0.0 all gradF = 1.0:
  alloc_prob(&fp, &ds, p);
  init_prob(&fp, &ds);



  // Subproblem allocated:
  alloc_subprob(&sp, p, &fp, &ds);


  // We loop until no negative entries in beta:
  int k = 1;
  int max_iters = 40;
  int itt = 0;
  int n = 0;

  while(k){
  //  break;
    // H matrix columns re-set and subproblem changed
    printf("i = %d\n",itt );
    setH(&fp, &ds);
    init_subprob(&sp, &fp, &ds);



    //  congjugate gradient algorithm
    //  n = 0 if algorithm completes
    //  n != 0 if algorithm interrupt
    n = cg(&sp, &fp);

    updateAlphaR(&fp, &sp);
    for (int i = 0; i < fp.n; i++) {
      check[i] = 1.0;
      for (int j = 0; j < fp.n; j++) {
        check[i] -= fullBigH[i][j]*fp.alpha[j];
      }
    }
    for (int i = 0; i < fp.p; i++) {
      printf("active is %d\n",fp.active[i] );
    }
    for (int i = 0; i < fp.q; i++) {
      printf("INactive is %d\n",fp.inactive[i] );
    }
    for (int i = 0; i < fp.n; i++) {
      printf("%lf\n",fp.gradF[i]-check[i] );
    }
    calcYTR(&sp, &fp);
    printf("total ytr %lf\n",sp.ytr );
    printf("ytr %lf\n",sp.ytr/(double)(fp.p) );
    double b = 1.0;
    for (int i = 0; i < fp.p; i++) {
      printf("h is %lf and alpha is %lf b i s %lf\n",sp.H[0][i],fp.alpha[fp.active[i]],b );
      b -= sp.H[0][i]*fp.alpha[fp.active[i]];
    }
    printf("gradF is %lf %lf %lf\n",fp.gradF[fp.active[0]],fp.gradF[fp.active[1]],fp.gradF[fp.active[2]] );
    b*=ds.y[fp.active[0]];
    printf("b = %lf\n",b );
  //  exit(0);
    calculateBeta(&fp, &sp, &ds);
    for (int i = 0; i < sp.p; i++) {
      printf("After: active[%d] = %d, alpha[] = %lf\n",i,fp.active[i],fp.alpha[fp.active[i]] );
    }
    for (int i = 0; i < fp.q; i++) {
      if (fp.beta[i] < 0.0) {
        printf("beta[%d] = %lf\n",fp.inactive[i],fp.beta[i] );
      }
    }
    if (n==0) {
      printf("Converged!\n" );
      for (int i = 0; i < sp.p; i++) {
        printf("%lf\n",fp.gradF[fp.active[i]] );
      }
      printf("ytr %lf\n",sp.ytr/(double)(fp.p) );

      int add = 2;
      int *temp = malloc(sizeof(int)*add);
      int *temp2 = malloc(sizeof(int)*add);
      add = findWorstest(&fp,add,temp,temp2);
      printf("adding %d\n",add );
      if (add == 0){
        printf("Add == 0\n" );
        break;
      }
      changeP(&fp, &sp, add);
      reinitprob(&fp, &sp, add, temp, temp2);

      free(temp);
      free(temp2);
    }
    int prqe = 0;

    while (n) {
      // While BCs broken
      k = singleswap(&ds, &fp, &sp, n);
      printf("k is %d\n",k );
      for (int i = 0; i < fp.q; i++) {
        printf("%lf\n",fp.beta[i] );
      }
      if (k < 0) {
        printf("shrinking\n" );
        shrinkSize(&fp, &sp, k+fp.p);
        break;
      }
      n = checkfpConstraints(&fp);
      prqe++;
      break;
      if (prqe == fp.n) {
        return 1;
      }
    }

    itt++;
    printf("%d\n",itt );
    if(itt == max_iters){
      printf("Reached max iters!!!!!\n\n\n" );
      break;
    }
  }

  for (int i = 0; i < fp.n; i++) {
    printf("alpha[%d] = %lf\n",i,fp.alpha[i] );
  }
  for (int i = 0; i < fp.q; i++) {
  printf("beta[%d](%d) = % lf\n",i,fp.inactive[i],fp.beta[i] );
  }
  if (k==0) {
    //printf("k goes to 0.\n");
  }

  //Memory freed
  freeDenseData(&ds);
  freeFullproblem(&fp);
  freeSubProblem(&sp);

  return 0;
}

void shrinkSize( struct Fullproblem *fp, struct Projected *sp, int k)
{
  int temp = fp->active[k];
  for (int i = k; i < fp->p - 1; i++) {
    fp->active[i] = fp->active[i+1];
  }
  printf("p is %d and q is %d\n",fp->p, fp->q );
  fp->p--;
  fp->q++;
  printf("p is %d and q is %d\n",fp->p, fp->q );
  fp->active = realloc(fp->active,sizeof(int)*fp->p);
  fp->inactive = realloc(fp->inactive,sizeof(int)*fp->q);
  fp->inactive[fp->q-1] = temp;
  fp->beta = realloc(fp->beta,sizeof(double)*fp->q);

  fp->h = realloc(fp->h,sizeof(double)*fp->q*fp->p);
  fp->partialH = realloc(fp->partialH,sizeof(double*)*fp->q);
  for (int i = 0; i < fp->q; i++) {
    fp->partialH[i] = &fp->h[i*fp->p];
  }

  sp->p--;

  sp->alphaHat = realloc(sp->alphaHat,sizeof(double)*sp->p);
  sp->yHat = realloc(sp->yHat,sizeof(double)*sp->p);
  sp->rHat = realloc(sp->rHat,sizeof(double)*sp->p);

  sp->gamma = realloc(sp->gamma,sizeof(double)*sp->p);
  sp->rho = realloc(sp->rho,sizeof(double)*sp->p);
  sp->Hrho = realloc(sp->Hrho,sizeof(double)*sp->p);

  sp->H = realloc(sp->H,sizeof(double*)*sp->p);
  sp->h = realloc(sp->h,sizeof(double)*((sp->p*(sp->p+1))/2));

  int j = 0;
  for (int i = 0; i < sp->p; i++) {
    sp->H[i] = &(sp->h[j]);
    j+=(sp->p-i-1);
  }
}

void changeP( struct Fullproblem *fp, struct Projected *sp, int add)
{
  fp->p += add;
  fp->q -= add;
  sp->p += add;

  fp->active = realloc(fp->active,sizeof(int)*fp->p);
  fp->inactive = realloc(fp->inactive,sizeof(int)*fp->q);
  fp->beta = realloc(fp->beta,sizeof(double)*fp->q);
  fp->h = realloc(fp->h,sizeof(double)*fp->q*fp->p);
  fp->partialH = realloc(fp->partialH,sizeof(double*)*fp->q);
  for (int i = 0; i < fp->q; i++) {
    fp->partialH[i] = &fp->h[i*fp->p];
  }

  sp->alphaHat = realloc(sp->alphaHat,sizeof(double)*sp->p);
  sp->yHat = realloc(sp->yHat,sizeof(double)*sp->p);
  sp->rHat = realloc(sp->rHat,sizeof(double)*sp->p);
  sp->gamma = realloc(sp->gamma,sizeof(double)*sp->p);
  sp->rho = realloc(sp->rho,sizeof(double)*sp->p);
  sp->Hrho = realloc(sp->Hrho,sizeof(double)*sp->p);
  sp->H = realloc(sp->H,sizeof(double*)*sp->p);
  sp->h = realloc(sp->h,sizeof(double)*((sp->p*(sp->p+1))/2));

  int j = 0;
  for (int i = 0; i < sp->p; i++) {
    sp->H[i] = &(sp->h[j]);
    j+=(sp->p-i-1);
  }
}

int findWorstest(struct Fullproblem *fp , int add, int* temp, int* temp2)
{
  double *betaVal = malloc(sizeof(double)*add);
  for (int i = 0; i < add; i++) {
    betaVal[i] = DBL_MAX;
  }
  for (int i = 0; i < fp->q; i++)
  {
    printf("beta is %lf\n",fp->beta[i] );
    for (int j = 0; j < add; j++)
    {
      if (fp->beta[i] < betaVal[j])
      {
        printf("ok %d %lf\n", i, betaVal[j]);
        for (int k = add - 1; k > j ; k--)
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
  int flag;
  for (int i = 0; i < add; i++) {
    flag = 0;
    for (int j = 0; j < add; j++) {
      if (fp->inactive[fp->q - 2 + i] == temp[i]) {
        flag = 1;
      }
    }
    if (flag == 0) {
      temp2[i] = fp->inactive[fp->q - 2 + i];
//      betaVal[i] = DBL_MAX;
    }
  }

  for (int i = 0; i < add; i++) {
    if (betaVal[i] > 0) {
      printf("%lf\n",betaVal[i] );
      return i;
    }
  }
  printf("temp are %d and %d\n",temp[0],temp[1] );
  printf("temp2 are %d and %d\n",temp2[0],temp2[1] );

  return add;
}


void reinitprob( struct Fullproblem *fp, struct Projected *sp, int add, int* temp, int* temp2)
{
  for (int i = 0; i < add; i++) {
    fp->active[fp->p - add + i] = temp[i];
  }

  int k = 0;
  for (int i = 0; i < fp->q; i++) {
    for (int j = 0; j < add; j++) {
      if (fp->inactive[i] == temp[j]) {
        fp->inactive[i] = temp2[k];
        k++;
        if (k == add) {
          return;
        }
      }
    }
  }
}



void   freeSubProblem( struct Projected* sp)
/* Function to free dynamically allocated memory in subproblem stuct. */
{
  free(sp->alphaHat);
  free(sp->yHat);
  free(sp->rHat);
  free(sp->H);

  free(sp->gamma);
  free(sp->rho);
  free(sp->Hrho);

  free(sp->h);
}

void freeDenseData(struct denseData *ds)
/* Function to free dynamically allocated memory in dense data set struct. */
{
  free(ds->data);
  free(ds->data1d);
  free(ds->y);
  free(ds->instanceLabels);
  free(ds->featureLabels);
}

void  freeFullproblem(struct Fullproblem *fp)
/* Function to free dynamically allocated memory in Fullproblem struct */
{
  free(fp->alpha);
  free(fp->beta);
  free(fp->gradF);

  free(fp->active);
  free(fp->inactive);

  free(fp->partialH);
  free(fp->h);
}
