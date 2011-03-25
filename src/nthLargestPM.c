/************************************************************************
 **
 ** nthLargestPM.c
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** created by: B. M. Bolstad   <bolstad@stat.berkeley.edu>
 ** created on: Feb 6, 2003  (but based on earlier work from Nov 2002)
 **
 ** last modified: Feb 6, 2003
 **
 ** License: GPL V2 or later (same as the rest of the Affy package)
 **
 ** General discussion
 **
 ** Implement log2  nth largest pm summarization.
 **
 ** NOTE: currently hardwired to give log2 of second largest PM, at
 ** some point this will be generalized. If there is anly one PM
 ** in a probeset then that will be returned
 **
 ** Feb 6, 2003 - Initial version of this summarization method
 ** Oct 5, 2003 - added summary_param
 ** Apr 5, 2004 - all malloc/free are now Calloc/Free
 **
 ************************************************************************/


#include "threestep_common.h"
#include "rma_common.h"
#include "nthLargestPM.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/***************************************************************************
 **
 ** double LogNthLargest(double *x, int length)
 **
 ** double *x - a vector of PM intensities 
 ** int length - length of *x
 **
 ** take the log2 of the median of PM intensities.
 **
 ***************************************************************************/

double LogNthLargest(double *x, int length,int n){
  int i;
  double nthLargest;
  double *buffer = Calloc(length,double);

  for (i = 0; i < length; i++)
    buffer[i] = x[i];

  qsort(buffer,length,sizeof(double), (int(*)(const void*, const void*))sort_double);

  if (length == 1){
    nthLargest = buffer[0];
  } else {
    nthLargest = buffer[length-n];
  }
  nthLargest = log(nthLargest)/log(2.0);

  Free(buffer);
  return (nthLargest);    
}

/***************************************************************************
 **
 ** double LogNthLargestPM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute log2 NthLargest expression measure. 
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

void LogNthLargestPM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param){
  int i,j;
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = data[j*rows + cur_rows[i]];  
    }
  } 
  
  for (j=0; j < cols; j++){
    results[j] = LogNthLargest(&z[j*nprobes],nprobes,2);
    resultsSE[j] = R_NaReal;
  }
  Free(z);
}




void LogNthLargestPM_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param){
  int i,j;
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = data[j*rows + cur_rows[i]];  
    }
  } 
  
  for (j=0; j < cols; j++){
    results[j] = LogNthLargest(&z[j*nprobes],nprobes,2);
    resultsSE[j] = R_NaReal;
  }

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      residuals[j*nprobes + i] = log(z[j*nprobes + i])/log(2.0)  - results[j];  
    }
  } 
  
  Free(z);
}
