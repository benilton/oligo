/************************************************************************
 **
 ** file: medianpolish.c
 **
 ** Copyright (C) 2002-2007 Ben Bolstad
 **
 ** created by: B. M. Bolstad   <bolstad@stat.berkeley.edu>
 ** created on: Jan 7, 2003 (but based on code dating back as far as June 2002)
 **
 ** last modified: Jan 7, 2003
 **
 ** License: GPL V2 or later (same as the rest of the Affy package)
 **
 ** Median polish summary measure (used in the RMA expression measure)
 **
 **
 ** History
 **
 ** Jan 7, 2003 - Initial version to fit into the three-step framework.
 ** Jan 13, 2003 - move median() into threestep_common.c
 ** Feb 24, 2003 - make maxiter get used.
 ** Jul 23, 2003 - add ability to accept SE parameter
 ** Sept 13, 2003 - introduced medianpolishPLM which returns 
 **                 most of what is required by the fitting
 **                 algorithm
 ** Oct 05, 2003 - added in summary_param
 ** Apr 5, 2004 - change malloc/free to Calloc/Free
 ** Nov 13, 2006 - make median calls to median_nocopy
 **
 ************************************************************************/


#include "PLM_medianpolish.h"
#include "rma_common.h"
#include "threestep_common.h"

#include <medianpolish.h>

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "preprocessCore_summarization_stubs.h"






/*************************************************************************************
 **
 ** void median_polish(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** double *data - a data matrix of dimension rows by cols (the entire PM matrix)
 ** int rows, cols - rows and columns dimensions of matrix
 ** int cur_rows - vector of length nprobes containg row indicies of *data matrix which apply for a 
 **                particular probeset
 ** double *results - a vector of length cols already allocated. on output contains expression values
 ** int nprobes - number of probes in current probeset.
 **
 ** a function to do median polish expression summary.
 **
 *************************************************************************************/

void median_polish_threestep(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param){

  int i,j;
 
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
  median_polish_no_copy(z, nprobes, cols, results, resultsSE);
    
  Free(z); 
}

void median_polishPLM(double *data, int rows, int cols, int *cur_rows, double *probe_param, double *chip_param, double *intercept_param, int nprobes, double *residuals){

  int i,j;

  double t = 0.0;
  
  double *r = Calloc(nprobes,double);
  double *c = Calloc(cols,double);
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
 
  median_polish_fit_no_copy(z, nprobes, cols, r, c, &t);

  for (i=0; i < nprobes; i++){
    probe_param[i] = r[i];
  }

  
  for (j=0; j < cols; j++){
    chip_param[j] =  c[j];
  }

  intercept_param[0] = t;

  for (j =0; j < cols; j++){
    for (i=0; i < nprobes; i++){
      residuals[j*nprobes +i] = z[j*nprobes +i];
    }
  }


  
  Free(r);
  Free(c);
  Free(z); 
}





void median_polish_threestep_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param){

  int i,j;

  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
 
  median_polish_no_copy(z, nprobes, cols, results, resultsSE);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      residuals[j*nprobes+i] = z[j*nprobes + i];
    }
  } 

  Free(z); 
}



