/*********************************************************************
 **
 ** file: rlm_threestep.c
 **
 ** Aim: connect robust linear model fitting routine with threestep framework.
 **
 ** Copyright (C) 2003 Ben Bolstad 
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Jan 13, 2003
 **
 ** Last modified: Jan 13, 2003
 **
 ** Manipulate Affy data into a form that can be fitted using a
 ** robust linear model.
 **
 ** History
 **
 ** Jan 13, 2003 - Initial version
 ** Feb 24, 2003 - fix a unused variable compiler warning, make model_matrix_test static 
 **                should probably be removed at some point.
 ** Jun 4, 2003  - add new parameter to rlm_fit call.
 ** Jul 23, 2003 - SE parameter added for threestep call (and implemented)
 **                commented out model_matrix_test
 ** Sep 13, 2003 - altered call to rlm_compute_se to reflect new structure
 ** Oct 5, 2003 - added summary_param, allow the user set parameters to control
 **               tuning constant, and M estimation function
 ** Oct 10, 2003 - a PLM threestep version added.
 ** May 26, 2004 - changed the rlm fitting algorithm
 ** Aug 05, 2004 - changed order in output to fit new version
 ** March 1, 2006 - changes all comments to ansi c standard
 **
 ********************************************************************/

#include "rlm_threestep.h"
#include "rlm.h"
#include "rlm_se.h"
#include "psi_fns.h"

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>



/***********************************************************************
 **
 ** void rlm_threestep(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** double *data - Probe intensity matrix
 ** int rows - number of rows in matrix *data (probes)
 ** int cols - number of cols in matrix *data (chips)
 ** int *cur_rows - indicies of rows corresponding to current probeset
 ** double *results - already allocated location to store expression measures (cols length)
 ** int nprobes - number of probes in current probeset.
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix
 **      corresponding to a single probeset, compute tukey biweight expression measure.
 **      Note that data is a probes by chips matrix, we will fit a robust linear
 **      model using rlm() to fit a model with probe and chip effects.
 **
 ** note that expression will be based on a model
 **  
 ** log2(B(PM)) = alpha_i + beta_j + epsilon_ij
 **
 ** alpha_i will be probe effect.
 ** beta_j will be chip effect.
 **
 ** expression values will be beta_j. alpha_i will be constrained to sum to zero
 ** 
 **
 **********************************************************************/

void rlm_threestep(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param){

  int i, j, row,curcol;
  
  int n = (nprobes*cols);
  int p = (cols+(nprobes-1));
  
  double *Y = Calloc(n,double);
  double *X = Calloc(n*p,double);
  double *out_beta=Calloc(p,double);
  double *out_se_estimates=Calloc(p,double);
  double *out_resids=Calloc(n,double);
  double *out_weights=Calloc(n,double);
  double *residSE = Calloc(2,double);


  /* log2 transform and create Y vector */
 
  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      Y[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);
    }
  }

  /* 
     now make an X matrix. First columns for probe effect then chip effect columns 
     Calloc puts everything to zero, so need change only non zero elements.
     
  */

  for (row =0; row < nprobes*cols; row++){
    curcol = row%nprobes;
    
    if (curcol == nprobes -1){
      for (j=0; j < (nprobes-1); j++){
	X[j*n + row] = -1.0;
      }
    } else {
      X[curcol*n + row] = 1.0;
    }
  }

 /* now do chip effects */

  for (row =0; row < nprobes*cols; row++){
    curcol = row/nprobes;         /* integer division */
    
    X[(curcol+(nprobes-1))*n + row] = 1.0;
    
  }


  /*  THE OLD WAY USING IRLS 
  rlm_fit(X,Y, n, p, out_beta, out_resids, out_weights,PsiFunc(summary_param->psi_method),summary_param->psi_k,20,0);  */

   rlm_fit_anova(Y, nprobes, cols, out_beta, out_resids, out_weights,PsiFunc(summary_param->psi_method),summary_param->psi_k,20,0);

  rlm_compute_se(X,Y, n, p, out_beta, out_resids, out_weights, out_se_estimates,NULL, residSE, 1, PsiFunc(summary_param->psi_method),summary_param->psi_k);
  

  for (i=0; i< cols; i++){
    results[i] = out_beta[i];    
    resultsSE[i] =  out_se_estimates[i];            
  }



  Free(out_se_estimates);
  Free(out_beta);
  Free(out_resids);
  Free(out_weights);
  Free(residSE);
  Free(X);
  Free(Y);
}


void rlm_threestep_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param){

  int i, j, row,curcol;
  
  int n = (nprobes*cols);
  int p = (cols+(nprobes-1));
  
  double *Y = Calloc(n,double);
  double *X = Calloc(n*p,double);
  double *out_beta=Calloc(p,double);
  double *out_se_estimates=Calloc(p,double);

  double *out_weights=Calloc(n,double);
  double *residSE = Calloc(2,double);


  /* log2 transform and create Y vector */
 
  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      Y[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);
    }
  }

  /* 
     now make an X matrix. First columns for probe effect then chip effect columns 
     Calloc puts everything to zero, so need change only non zero elements.
     
  */
  
  for (row =0; row < nprobes*cols; row++){
    curcol = row%nprobes;
    
    if (curcol == nprobes -1){
      for (j=0; j < (nprobes-1); j++){
	X[j*n + row] = -1.0;
      }
    } else {
      X[curcol*n + row] = 1.0;
    }
  }

  /* now do chip effects */

  for (row =0; row < nprobes*cols; row++){
    curcol = row/nprobes;         /*integer division */
    
    X[(curcol+(nprobes-1))*n + row] = 1.0;
    
  }

  /* THE OLD WAY USING IRLS 
  rlm_fit(X,Y, n, p, out_beta, residuals, out_weights,PsiFunc(summary_param->psi_method),summary_param->psi_k,20,0); */

  rlm_fit_anova(Y, nprobes, cols, out_beta, residuals, out_weights,PsiFunc(summary_param->psi_method),summary_param->psi_k,20,0);
  rlm_compute_se(X,Y, n, p, out_beta, residuals, out_weights, out_se_estimates,NULL, residSE, 1, PsiFunc(summary_param->psi_method),summary_param->psi_k);
  

  for (i=0; i< cols; i++){
    results[i] =  out_beta[i];               
    resultsSE[i] = out_se_estimates[i];      
  }



  Free(out_se_estimates);
  Free(out_beta);
  Free(out_weights);
  Free(residSE);
  Free(X);
  Free(Y);
}



/*
  .C("model_matrix_test",as.integer(2),as.integer(3),as.integer(3))

*/




/*
static void model_matrix_test(int *rows_in,int *cols_in, int *nprobes_in){

  int i, j, row,curcol;
  
  int cols = *cols_in;
  int nprobes = *nprobes_in;



  double *X = Calloc(nprobes*cols*((nprobes-1)+cols),double);
  


 for (row =0; row < nprobes*cols; row++){
    curcol = row%nprobes;
    
    if (curcol == nprobes -1){
      for (j=0; j < (nprobes-1); j++){
	X[j*(nprobes*cols) + row] = -1.0;
      }
    } else {
      X[curcol*(nprobes*cols) + row] = 1.0;
    }
  }

  // now do chip effects 

  for (row =0; row < nprobes*cols; row++){
    curcol = row/nprobes;         // *integer division 
    
    X[(curcol+(nprobes-1))*(nprobes*cols) + row] = 1.0;
    
  }
  
  for (i=0; i < nprobes*cols; i++){
    for (j =0; j < (nprobes-1)+cols; j++){
      printf("%f ",X[j*nprobes*cols + i]);
    }
    printf("\n");
  }


  free(X);
  
}
*/
