#ifndef RMAPLM_PSEUDO_H
#define RMAPLM_PSEUDO_H

void compute_pseudoweights(double *resids, double *weights,int rows, int cols, int psi_code,double psi_k);
void compute_pseudoSE(double *resids, double *pseudoSE,int rows, int cols, int psi_code,double psi_k);

#endif
