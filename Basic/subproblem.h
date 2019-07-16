#ifndef SUBPROBLEM_H
#define SUBPROBLEM_H

#include <stdlib.h>
//#include "fullproblem.h"

class Subproblem
{
private:
  double* h;
  double** H;
  double*y;
  int p;
  int* P;
  double* gradF;
  double* alpha;
  double c;
  bool* changeOut;


public:
  Subproblem(double** supH, double* supY, int p, int* P);
  void alloc_space(int p);
  void project_down(double** supH, double* supY, int p, int* P);
  void constraint_projection(double* vecOut, double* vecIn);
  void HtimesVec(double* vecOut, double* vecIn);
  double inner_prod(double* a1, double* a2);
  void init_error(double* err, double* perr);
  void linOp2(double* err, double* Aperr, double alp);
  void linOp(double* err, double* Aperr, double alp);
  void php(double* Aperr, double* perr);
  void solve();

  void cg(double* err, double* perr, double* Aperr);



};

#endif
