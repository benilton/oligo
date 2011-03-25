/************************************************************************
 **
 ** PLM_log_avg.c
 **
 ** Copyright (C) 2003 Ben Bolstad
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
 ** Implement log2 average pm summarization, with or without normalization
 **
 ** Feb 6, 2002 - Initial version of this summarization method
 ** Jul 23, 2003 - parameter for storing SE added (not yet implemented)
 ** Oct 5, 2003 - method of adding parameters.
 ** May 26, 2007 - clean up code. Core functionality has been moved to preprocessCore
 ** Sept 9, 2007 - rename this file PLM_log_avg.c from log_avg.c
 **
 ************************************************************************/

#include "PLM_log_avg.h"

#include <log_avg.h>

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



void LogAverage_threestep(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param){
  
  LogAverage(data,rows,cols,cur_rows,results,nprobes,resultsSE);

}


void LogAverage_threestep_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param){
  int i,j;

  LogAverage(data,rows,cols,cur_rows,results,nprobes,resultsSE);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      residuals[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0) - results[j] ;  
    }
  }
}
