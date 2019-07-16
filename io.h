#ifndef IO_H
#define IO_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <iostream>

struct svm_args
{
	int type;
	int kernel;
	int degree;
	double C;	/* for C_SVC, EPSILON_SVR and NU_SVR */
	bool verbose;

	/* Parameters for later */

//	double gamma;	/* for poly/rbf/sigmoid */
//	double coef0;	/* for poly/sigmoid */
//	double cache_size; /* in MB */
//	double eps;	/* stopping criteria */
//	int nr_weight;		/* for C_SVC */
//	int *weight_label;	/* for C_SVC */
//	double* weight;		/* for C_SVC */
//	double nu;	/* for NU_SVC, ONE_CLASS, and NU_SVR */
//	double p;	/* for EPSILON_SVR */
//	int shrinking;	/* use the shrinking heuristics */
//	int probability; /* do probability estimates */

};


int parse_arguments(int argc, char *argv[], char** filename, struct svm_args *parameters);

int readline(FILE *input, char **line);


#endif
