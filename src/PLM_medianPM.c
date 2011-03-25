/************************************************************************
 **
 ** PLM_medianPM.c   (was medianPM.c)
 **
 ** Copyright (C) 2002-2007 Ben Bolstad
 **
 ** created by: B. M. Bolstad   <bmb@bmbolstad.com>  
 ** created on: Feb 6, 2003  (but based on earlier work from Nov 2002)
 **
 ** last modified: Feb 6, 2003
 **
 ** License: GPL V2 or later (same as the rest of the Affy package)
 **
 ** General discussion
 **
 ** Implement log2 median pm summarization.
 **
 ** Feb 6, 2003 - Initial version of this summarization method
 ** Feb 24, 2003 - remove unused int i from LogMedian()
 ** Jul 23, 2003 - add a parameter for storing SE (not yet implemented)
 ** Oct 10, 2003 - PLM version of threestep
 **
 ************************************************************************/


#include "PLM_medianPM.h"
#include "threestep_common.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <log_median.h>


void LogMedianPM_threestep(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param){

  LogMedian(data, rows, cols, cur_rows, results, nprobes, resultsSE);
}


void LogMedianPM_threestep_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param){
  int i,j;
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = data[j*rows + cur_rows[i]];  
    }
  } 
  
  logmedian_no_copy(z, nprobes, cols, results, resultsSE);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      residuals[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0) - results[j];  
    }
  } 
    


  Free(z);
}
