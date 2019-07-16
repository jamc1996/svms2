#include "kernels.h"


int setH(struct Fullproblem *prob, struct denseData *ds, struct svm_args *params)
/*  Function to update the values of the matrix partialH
 *  TO DO -> ALLOW ONLY UPDATE THE SWAPPED OUT VALUES - or not?   */
{
  if (params->kernel == LINEAR) {
    for (int i = 0; i < prob->q; i++) {
      for (int j = 0; j < prob->p; j++) {
        prob->partialH[i][j] = 0.0;
        for (int k = 0; k < ds->nFeatures; k++) {
          prob->partialH[i][j] += ds->data[prob->inactive[i]][k]*ds->data[prob->active[j]][k];
        }
        prob->partialH[i][j] *= ds->y[prob->inactive[i]]*ds->y[prob->active[j]];
      }
    }
  }
  else if (params->kernel == POLYNOMIAL) {
    for (int i = 0; i < prob->q; i++) {
      for (int j = 0; j < prob->p; j++) {
        prob->partialH[i][j] = 0.0;
        for (int k = 0; k < ds->nFeatures; k++) {
          prob->partialH[i][j] += ds->data[prob->inactive[i]][k]*ds->data[prob->active[j]][k];
        }
        prob->partialH[i][j] = pow(prob->partialH[i][j]+params->Gamma, params->degree);
        prob->partialH[i][j] *= ds->y[prob->inactive[i]]*ds->y[prob->active[j]];
      }
    }
  }
  else if (params->kernel == EXPONENTIAL) {
    double x;
    double y;
    for (int i = 0; i < prob->q; i++) {
      for (int j = 0; j < prob->p; j++) {
        y = 0.0;
        for (int k = 0; k < ds->nFeatures; k++) {
          x = (ds->data[prob->inactive[i]][k]-ds->data[prob->active[j]][k]);
          y -= x*x;
        }
        y *= params->Gamma;
        prob->partialH[i][j] = exp(y)*ds->y[prob->inactive[i]]*ds->y[prob->active[j]];
      }
    }
  }
  else {
    return 1;
  }
  return 0;
}

int updateSubH(struct Fullproblem *fp, struct Projected *sp, struct denseData *ds, struct svm_args *params)
{
  if (params->kernel == LINEAR) {
    for (int i = 0; i < sp->p; i++) {
      for (int j = i; j < sp->p; j++) {
        sp->H[i][j] = 0.0;
        for (int k = 0; k < ds->nFeatures; k++) {
          sp->H[i][j] += ds->data[fp->active[i]][k]*ds->data[fp->active[j]][k];
        }
        sp->H[i][j] *= ds->y[ fp->active[i] ] * ds->y[ fp->active[j] ];
      }
    }
  }
  else if (params->kernel == POLYNOMIAL) {
    for (int i = 0; i < fp->p; i++) {
      for (int j = i; j < fp->p; j++) {
        sp->H[i][j] = 0.0;
        for (int k = 0; k < ds->nFeatures; k++) {
          sp->H[i][j] += ds->data[fp->active[i]][k]*ds->data[fp->active[j]][k];
        }
        sp->H[i][j] = pow(sp->H[i][j]+params->Gamma, params->degree);
        sp->H[i][j] *= ds->y[fp->active[i]]*ds->y[fp->active[j]];
      }
    }
  }
  else if (params->kernel == EXPONENTIAL) {
    double x;
    double y;
    for (int i = 0; i < fp->p; i++) {
      for (int j = i; j < fp->p; j++) {
        y = 0.0;
        for (int k = 0; k < ds->nFeatures; k++) {
          x = (ds->data[fp->active[i]][k]- ds->data[fp->active[j]][k]);
          y -= x*x;
        }
        y *= params->Gamma;
        sp->H[i][j] = exp(y)*ds->y[fp->active[i]]*ds->y[fp->active[j]];
      }
    }
  }
  else{
    return 1;
  }
  return 0;
}


void partialHupdate(struct Fullproblem *fp, struct Projected *sp, struct denseData *ds, struct svm_args *params, int n, int worst)
{
  if (params->kernel == LINEAR) {
    for (int j = 0; j < fp->q; j++) {
      if (j == worst) {
        for (int k = 0; k < n; k++) {
          fp->partialH[j][k] = sp->H[k][n];
        }
        for (int k = n+1; k < fp->p; k++) {
          fp->partialH[j][k] = sp->H[n][k];
        }
        continue;
      }
      fp->partialH[j][n] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        fp->partialH[j][n]+=ds->data[fp->inactive[worst]][k]*ds->data[fp->inactive[j]][k];
      }
      fp->partialH[j][n]*=ds->y[fp->inactive[worst]]*ds->y[fp->inactive[j]];
    }
    for (int i = 0; i < n; i++) {
      sp->H[i][n] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        sp->H[i][n]+=ds->data[fp->inactive[worst]][k]*ds->data[fp->active[i]][k];
      }
      sp->H[i][n]*=ds->y[fp->inactive[worst]]*ds->y[fp->active[i]];
    }

    sp->H[n][n] = 0.0;
    for (int k = 0; k < ds->nFeatures; k++) {
      sp->H[n][n]+=ds->data[fp->inactive[worst]][k]*ds->data[fp->inactive[worst]][k];
    }

    for (int i = n+1; i < sp->p; i++) {
      sp->H[n][i] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        sp->H[n][i]+=ds->data[fp->inactive[worst]][k]*ds->data[fp->active[i]][k];
      }
      sp->H[n][i]*=ds->y[fp->inactive[worst]]*ds->y[fp->active[i]];
    }
  }
  else if (params->kernel == POLYNOMIAL) {
    for (int j = 0; j < fp->q; j++) {
      if (j == worst) {
        for (int k = 0; k < n; k++) {
          fp->partialH[j][k] = sp->H[k][n];
        }
        for (int k = n+1; k < fp->p; k++) {
          fp->partialH[j][k] = sp->H[n][k];
        }
        continue;
      }
      fp->partialH[j][n] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        fp->partialH[j][n]+=ds->data[fp->inactive[worst]][k]*ds->data[fp->inactive[j]][k];
      }
      fp->partialH[j][n] = pow(fp->partialH[j][n]+params->Gamma, params->degree);
      fp->partialH[j][n] *= ds->y[fp->inactive[worst]]*ds->y[fp->inactive[j]];
    }

    for (int i = 0; i < n; i++) {
      sp->H[i][n] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        sp->H[i][n]+=ds->data[fp->inactive[worst]][k]*ds->data[fp->active[i]][k];
      }
      sp->H[i][n] = pow(sp->H[i][n]+params->Gamma, params->degree);
      sp->H[i][n]*=ds->y[fp->inactive[worst]]*ds->y[fp->active[i]];
    }

    sp->H[n][n] = 0.0;
    for (int k = 0; k < ds->nFeatures; k++) {
      sp->H[n][n]+=ds->data[fp->inactive[worst]][k]*ds->data[fp->inactive[worst]][k];
    }
    sp->H[n][n] = pow(sp->H[n][n]+params->Gamma, params->degree);

    for (int i = n+1; i < sp->p; i++) {
      sp->H[n][i] = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        sp->H[n][i]+=ds->data[fp->inactive[worst]][k]*ds->data[fp->active[i]][k];
      }
      sp->H[n][i] = pow(sp->H[n][i]+params->Gamma, params->degree);
      sp->H[n][i]*=ds->y[fp->inactive[worst]]*ds->y[fp->active[i]];
    }
  }
  else if (params->kernel == EXPONENTIAL) {
    double x,y;
    for (int j = 0; j < fp->q; j++) {
      if (j == worst) {
        for (int k = 0; k < n; k++) {
          fp->partialH[j][k] = sp->H[k][n];
        }
        for (int k = n+1; k < fp->p; k++) {
          fp->partialH[j][k] = sp->H[n][k];
        }
        continue;
      }
      y = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        x = ds->data[fp->inactive[worst]][k] - ds->data[fp->inactive[j]][k];
        y -= x*x;
      }
      y *= (params->Gamma);
      fp->partialH[j][n] = exp(y)*ds->y[fp->inactive[worst]]*ds->y[fp->inactive[j]];
    }

    for (int i = 0; i < n; i++) {
      y = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        x = ds->data[fp->inactive[worst]][k] - ds->data[fp->active[i]][k];
        y -= x*x;
      }
      y *= (params->Gamma);
      sp->H[i][n] = exp(y)*ds->y[fp->inactive[worst]]*ds->y[fp->active[i]];
    }

    y = 0.0;
    for (int k = 0; k < ds->nFeatures; k++) {
      x = ds->data[fp->inactive[worst]][k] - ds->data[fp->inactive[worst]][k];
      y -= x*x;
    }
    y *= params->Gamma;
    sp->H[n][n] = exp(y);

    for (int i = n+1; i < sp->p; i++) {
      y = 0.0;
      for (int k = 0; k < ds->nFeatures; k++) {
        x = ds->data[fp->inactive[worst]][k] - ds->data[fp->active[i]][k];
        y += x*x;
      }
      y *= params->Gamma;
      sp->H[n][i] = exp(y)*ds->y[fp->inactive[worst]]*ds->y[fp->active[i]];
    }

  }
}
