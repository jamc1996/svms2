#ifndef SVM_H
#define SVM_H

struct svm_args
{
	int type;
	int kernel;
	int degree;
	double C;	/* for C_SVC, EPSILON_SVR and NU_SVR */
	int verbose;
	int test;
	char* modelfile;
	int save;
	char* savename;
} clargs;


struct denseData{
  int nInstances;
  int nFeatures;
	int nPos;
	int nNeg;
  double** data;
  double* data1d;
	double* y;
  char** instanceLabels;
  char** featureLabels;
};

struct CGSLSProj{
	// Re-inited each time:
	double* alphaHat;
	double* yHat;
	double* rHat;
	double** H;

	// Will be set during CG iterations
	double* gamma;
	double* rho;
	double* Hrho;

	// Changed independent (p size, ytr calc) or unchanging
	int p;
	double C;
	double* h;
	double ytr;
};

struct Projected{
	// Re-inited each time:
	double* alphaHat;
	double* yHat;
	double* rHat;
	double** H;

	// Will be set during CG iterations
	double* gamma;
	double* Hgamma;
	double* rho;
	double* Hrho;

	// Changed independent (p size, ytr calc) or unchanging
	int p;
	double C;
	double* h;
	double ytr;
}Subproblem;

struct Fullproblem{
	int n;
	int p;
	int q;

	double C;

	double* alpha;
	double* beta;
	double* gradF;

	int* active;
	int* inactive;

	double** partialH;
	double* h;
} fullprob;

#endif