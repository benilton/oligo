/************************************************************************
 **
 ** PLM_median_logPM.c
 **
 ** created by: B. M. Bolstad   <bmb@bmbolstad.com>
 ** created on: Feb 6, 2003  (but based on earlier work from Nov 2002)
 **
 ** Copyright (C) 2003-2007   Ben Bolstad
 **
 ** last modified: Feb 6, 2003
 **
 ** License: GPL V2 or later (same as the rest of the Affy package)
 **
 ** General discussion
 **
 ** Implement median log2 pm summarization.
 **
 ** Feb 6, 2003 - Initial version of this summarization method
 ** Feb 24, 2003 - Remove unused variable in i from MedianLog
 ** Jul 23, 2003 - add SE parameter (but not yet implemented)
 ** Oct 10, 2003 - added PLM version
 ** Sep 10, 2007 - renamed this file to PLM_median_logPM.c from median_logPM.c
 **
 ************************************************************************/

#include "threestep_common.h"
#include "PLM_median_logPM.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <median_log.h>



/***************************************************************************
 **
 ** double MedianLogPM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute log2 Median expression measure. 
 **      Note that data is a probes by chips matrix.
 **
 ** double *data - Probe intensity matrix
 ** int rows - number of rows in matrix *data (probes)
 ** int cols - number of cols in matrix *data (chips)
 ** int *cur_rows - indicies of rows corresponding to current probeset
 ** double *results - already allocated location to store expression measures (cols length)
 ** int nprobes - number of probes in current probeset.
 **
 ***************************************************************************/

void MedianLogPM_threestep(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE,  summary_plist *summary_param){
  MedianLog(data, rows, cols, cur_rows, results, nprobes, resultsSE);
}



/***************************************************************************
 **
 ** double MedianLogPM_noSE(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute log2 Median expression measure. 
 **      Note that data is a probes by chips matrix.
 **
 ** double *data - Probe intensity matrix
 ** int rows - number of rows in matrix *data (probes)
 ** int cols - number of cols in matrix *data (chips)
 ** int *cur_rows - indicies of rows corresponding to current probeset
 ** double *results - already allocated location to store expression measures (cols length)
 ** int nprobes - number of probes in current probeset.
 **
 ***************************************************************************/

void MedianLogPM_noSE(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes){
  MedianLog_noSE(data, rows, cols, cur_rows, results, nprobes);
}





void MedianLogPM_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals,  summary_plist *summary_param){
  int i,j;
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = data[j*rows + cur_rows[i]];  
    }
  } 


  medianlog_no_copy(z, nprobes, cols, results, resultsSE);
  
  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      residuals[j*nprobes + i] = z[j*nprobes + i] - results[j];  
    }
  } 

  
  Free(z);
}
