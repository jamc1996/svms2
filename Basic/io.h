#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>

#include "svm.h"


void read_file(char* const filename, struct denseData* ds);
int readline(FILE *input, char **line);
void count_entries(FILE *input, struct denseData* ds);
int parse_arguments(int argc, char *argv[], char** filename, struct svm_args *parameters);


#endif
