#include <stdio.h>
#include <string>
#include <unistd.h>
#include <stdlib.h>


int parse_arguments(int argc, char *argv[], char** filename);
void read_file(char* const filename);


int main(int argc, char *argv[]) {
  char* filename = NULL;

  parse_arguments(argc, argv, &filename);
  read_file(filename);

  return 0;
}

int parse_arguments(int argc, char *argv[], char** filename){
  int c;

  while ((c = getopt( argc, argv, "f:")) != -1){
    switch (c) {
      case 'f':
        *filename = optarg;
    }
  }
  if (*filename == NULL) {
    printf("main.cpp: parse_arguments: no input file selected.\n");
    exit(1);
  }

  return 0;
}


void read_file(char* const filename) {
  FILE *fp = fopen(filename, "r");

  if (fp == NULL) {
    fprintf(stderr, "main.cpp: read_file: file '%s' not found\n",filename);
    exit(1);
  }
}
