#ifndef SVM_H
#define SVM_H


struct instance{
	int index;
	double value;
};

struct dataset {
  int l;
	double *y;
	struct instance *x_space;
	struct instance **x;
};

#endif
