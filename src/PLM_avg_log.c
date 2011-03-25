/************************************************************************
 **
 ** PLM_avg_log.c
 **
 ** created by: B. M. Bolstad   <bmb@bmbolstad.com>
 ** created on: Jan 7, 2003  (but based on earlier work from Nov 2002)
 **
 ** Copyright (C) 2002-2007 Ben Bolstad
 **
 ** last modified: Jan 7, 2003
 **
 ** License: GPL V2 or later (same as the rest of the Affy package)
 **
 ** General discussion
 **
 ** Implement avgerage log2 pm summarization, with or without normalization
 **
 ** Nov 2, 2002 - modify so that it will work efficently with affy2
 ** Jan 3, 2003 - Clean up commenting, prepare for integration in AffyExtensions
 ** Jan 7, 2003 - Make function standalone, to prepare for later combination into
 **               a more general framework.
 ** Jul 23, 2003 - add parameter for computing SE and SE implemented
 ** Oct 5, 2003 - add output_param
 ** Oct 10, 2003 - added threestepPLM version of this summary.
 ** May 26, 2007 - clean up. most of the core functionality has been moved
 **                to preprocessCore
 ** Sept 3, 2007 - renamed this file from avg_log.c to PLM_avg_log.c
 **
 ************************************************************************/

#include "PLM_avg_log.h"

/* the preprocessCore one */
#include <preprocessCore_summarization_stubs.c>

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/***************************************************************************
 **
 ** double AverageLog(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute average log2 expression measure. 
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

void AverageLog_threestep(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param){


  AverageLog(data,rows,cols,cur_rows,results,nprobes,resultsSE);


}



void AverageLog_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param){
  int i,j;
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = data[j*rows + cur_rows[i]];  
    }
  } 
  averagelog_no_copy(z, nprobes, cols, results, resultsSE);

  for (j =0; j < cols; j++){
    for (i=0; i < nprobes; i++){
      residuals[j*nprobes + i] = z[j*nprobes + i] - results[j] ;
    }
  }

  Free(z);
}
