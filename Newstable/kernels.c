#include "kernels.h"


int setH(struct Fullproblem *prob, struct denseData *ds, struct svm_args *params)
/*  Function to update the values of the matrix partialH
 *  TO DO -> ALLOW ONLY UPDATE THE SWAPPED OUT VALUES - or not?   */
{
  Cell *temp = prob->partialH.head;
  int i = 0;
  if (params->kernel == LINEAR) {
    while(temp != NULL) {
      for (int j = 0; j < prob->n; j++) {
        temp->line[j] = 0.0;
        for (int k = 0; k < ds->nFeatures; k++) {
          temp->line[j] += ds->data[j][k]*ds->data[prob->active[i]][k];
        }
        temp->line[j] *= ds->y[prob->active[i]]*ds->y[j];
      }
      temp = temp->next;
      i++;
    }
  }
  else if (params->kernel == POLYNOMIAL) {
    while(temp != NULL) {
      for (int j = 0; j < prob->n; j++) {
        temp->line[j] = 0.0;
        for (int k = 0; k < ds->nFeatures; k++) {
          temp->line[j] += ds->data[j][k]*ds->data[prob->active[i]][k];
        }
        temp->line[j] = pow(temp->line[j]+params->Gamma, params->degree);
        temp->line[j] *= ds->y[prob->active[i]]*ds->y[j];
      }
      i++;
      temp = temp->next;
    }
  }
  else if (params->kernel == EXPONENTIAL) {
    double x;
    double y;
    while (temp != NULL) {
      for (int j = 0; j < prob->p; j++) {
        y = 0.0;
        for (int k = 0; k < ds->nFeatures; k++) {
          x = (ds->data[prob->active[i]][k]-ds->data[j][k]);
          y -= x*x;
        }
        y *= params->Gamma;
        temp->line[j] = exp(y)*ds->y[prob->active[i]]*ds->y[j];
      }
      i++;
      temp = temp->next;
    }
  }
  else {
    return 1;
  }
  return 0;
}

int updateSubH(struct Fullproblem *fp, struct Projected *sp, struct denseData *ds, struct svm_args *params)
{
  Cell* temp = fp->partialH.head;
  int i = 0;
  while ( temp != NULL ) {
    for (int j = i; j < fp->p ; j++) {
      sp->H[i][j] = temp->line[fp->active[j]];
    }
    i++;
    temp = temp->next;
  }

  return 0;
}

void appendUpdate(struct denseData *ds, double *line, int n)
{
  if (parameters.kernel == LINEAR) {
    for (int i = 0; i < ds->nInstances; i++) {
      line[i] = 0.0;
      for (int j = 0; j < ds->nFeatures; j++) {
        line[i] += ds->data[i][j]*ds->data[n][j];
      }
      line[i]*= ds->y[i]*ds->y[n];
    }
  }
  if (parameters.kernel == POLYNOMIAL) {
    for (int i = 0; i < ds->nInstances; i++) {
      line[i] = 0.0;
      for (int j = 0; j < ds->nFeatures; j++) {
        line[i] += ds->data[i][j]*ds->data[n][j];
      }
      line[i] = pow(line[i]+parameters.Gamma, parameters.degree);
      line[i]*= ds->y[i]*ds->y[n];
    }
  }
  else if (parameters.kernel == EXPONENTIAL) {
    double x,y;
    for (int i = 0; i < ds->nInstances; i++) {
      y = 0.0;
      for (int j = 0; j < ds->nFeatures; j++) {
        x = ds->data[i][j] - ds->data[n][j];
        y -= x*x;
      }
      y *= (parameters.Gamma);
      line[i] = exp(y)*ds->y[i]*ds->y[n];
    }
  }
}


void partialHupdate(struct Fullproblem *fp, struct Projected *sp, struct denseData *ds, struct svm_args *params, int n, int worst)
{
  double *nline = findListLineSetLabel(fp->partialH, fp->active[n],fp->inactive[worst]);
  if (params->kernel == LINEAR) {
    for (int j = 0; j < fp->n; j++) {
      nline[j] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        nline[j]+=ds->data[fp->inactive[worst]][k]*ds->data[j][k];
      }
      nline[j]*=ds->y[fp->inactive[worst]]*ds->y[j];
    }
  }
  else if (params->kernel == POLYNOMIAL) {
    for (int j = 0; j < fp->n; j++) {
      nline[j] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        nline[j]+=ds->data[fp->inactive[worst]][k]*ds->data[j][k];
      }
      nline[j] = pow(nline[j]+params->Gamma, params->degree);
      nline[j] *= ds->y[fp->inactive[worst]]*ds->y[j];
    }
  }
  else if (params->kernel == EXPONENTIAL) {
    double x,y;
    for (int j = 0; j < fp->n; j++) {
      y = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        x = ds->data[fp->inactive[worst]][k] - ds->data[j][k];
        y -= x*x;
      }
      y *= (params->Gamma);
      nline[j] = exp(y)*ds->y[fp->inactive[worst]]*ds->y[j];
    }
  }

  for (int i = 0; i < n; i++) {
    sp->H[i][n] = nline[fp->active[i]];
  }
  sp->H[n][n] = nline[fp->inactive[worst]];

  for (int i = n+1; i < sp->p; i++) {
    sp->H[n][i] = nline[fp->active[i]];
  }
}
