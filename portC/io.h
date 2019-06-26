#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include "svm.h"


void read_file(char* const filename, struct denseData* ds);
int readline(FILE *input, char **line);
void count_entries(FILE *input, struct denseData* ds);
int parse_arguments(int argc, char *argv[], char** filename, struct svm_args *parameters);

void preprocess(struct denseData *ds);
void calcMeans(double *mean, struct denseData *ds);

void normalise(double* mean, double* stdDev, struct denseData* ds);
void calcStdDev(double* stdDev, double* mean, struct denseData *ds);

#endif
