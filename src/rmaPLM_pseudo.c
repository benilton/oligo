/*********************************************************************
 **
 ** file: rmaPLM_pseudo.c
 **
 ** Aim: When rma is fit as a PLM compute pseudo weights and
 **      standard errors..
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Sept 14, 2003
 **
 **
 ** Note that the weights and weights given here in no way reflect
 ** what might really be happening with the median polish. Instead
 ** these are given to allow usage of some of the quality assessment
 ** functionality.
 **
 **
 ** History
 ** Sept 14, 2003 - Initial version
 ** Sept 16, 2003 - estimate residual SE better
 ** Oct  9, 2003 - adjust Pseudo variance matrix a little.
 **                slightly more documentation
 ** March 1, 2006 - change comments to ansi c style
 **
 *********************************************************************/


#include <math.h>
#include "psi_fns.h"
#include "rlm.h"

#include "rmaPLM_pseudo.h"

/*********************************************************************
 **
 ** void compute_pseudoweights(double *resids, double *weights,int rows, 
 **                             int cols, int psi_code,double psi_k)
 **
 ** This function computes pseudo weights using a specified M estimation
 ** function.
 **
 *********************************************************************/



void compute_pseudoweights(double *resids, double *weights,int rows, int cols, int psi_code,double psi_k){

  int i,j;
  int n = rows*cols;
  double MAD;
 
  pt2psi PsiFn = PsiFunc(psi_code);
  
  MAD = med_abs(resids,n)/0.6745;

  for (i = 0; i < rows; i++){
    for (j = 0; j < cols; j++){
      weights[j*rows + i] = PsiFn(resids[j*rows + i]/MAD, psi_k, 0);
    }
  }
}



/*********************************************************************
 **
 ** void compute_pseudoSE(double *resids, double *pseudoSE,int rows, 
 **                         int cols, int psi_code,double psi_k)
 **
 ** This function computes pseudo SE using a specified M estimation
 ** function.
 **
 *********************************************************************/


void compute_pseudoSE(double *resids, double *pseudoSE,int rows, int cols, int psi_code,double psi_k){


  int i,j;
  int n = rows*cols;
  double MAD;
  double residSE = 0.0;


  double sum_weights=0.0;
  
  pt2psi PsiFn = PsiFunc(psi_code);
  

   
  MAD = med_abs(resids,n)/0.6745;

  for (i=0; i < rows; i++){
    for (j=0; j < cols; j++){
      residSE+= PsiFn(resids[j*rows + i]/MAD, psi_k, 0)*resids[j*rows+i]*resids[j*rows+i];
    }
  }
  
  residSE = sqrt(residSE/(double)(n - (rows + cols -1)));
  
  
  for (i = 0; i < rows; i++){
    sum_weights = 0.0;
    for (j = 0; j < cols; j++){
      sum_weights +=  PsiFn(resids[j*rows + i]/MAD, psi_k, 0); /*   *PsiFn(resids[j*rows + i]/MAD, psi_k, 0); */
    }
    pseudoSE[i] = residSE/sqrt(sum_weights);
  }

  for (j = 0; j < cols; j++){
    sum_weights = 0.0;
    for (i = 0; i < rows; i++){
      sum_weights +=  PsiFn(resids[j*rows + i]/MAD, psi_k, 0); /*    *PsiFn(resids[j*rows + i]/MAD, psi_k, 0); */
    }
    pseudoSE[j + rows] = residSE/sqrt(sum_weights);
  }

}
