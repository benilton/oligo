/************************************************************
 **
 ** file: qnorm_probeset.c
 **
 ** aim: a probeset specific quantile normalization
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** written by: B. M. Bolstad  <bolstad@stat.berkeley.edu>
 **
 ** written: May 19, 2003
 **
 ** This implements an extension to the quantile normalization 
 ** method. In particular we will  normalize some summary of 
 ** the probeset (perhaps the mean or median) and then adjust the
 ** probes accordingly. Working in the log scale (base 2) will be supported.
 ** The goal here will be to keep stricter parallism than the standard 
 ** quantile normalization
 **
 ** History:
 **
 ** May 19, 2003 - Initial verison
 ** May 29, 2003 - Further implementation, added parameter
 **                for using median
 ** Aug 25, 2003 - log and now non log scales are handled
 **
 ************************************************************/

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <math.h>

#include <median_log.h>
#include <avg_log.h>
#include <log_median.h>
#include <log_avg.h>

#include "preprocessCore_normalization_stubs.h"


/*************************************************************
 **
 ** void AdjustProbes(double *data, int rows, int cols, int *cur_rows, 
 **                   double *results, double *results_original, int nprobes)
 **
 ** double *data
 ** int rows, cols;
 ** int *cur_rows - row indices in PM matrix
 ** double *results - should be normalized summary
 ** double *results_original - the unnormalized summary
 ** int nprobes - number of probes in current set (that is being adjusted)
 ** int current_probeset - index of current probeset
 **
 **
 ************************************************************/


void AdjustProbes(double *data, int rows, int cols, int *cur_rows, double *results, double *results_original, int nprobes, int n_probesets, int current_probeset,int uselogs){
  int i,j;
  double adj;
  double *z = Calloc(nprobes*cols,double);
 
  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);
    }
  }
   
  for (j=0; j < cols; j++){
    if (uselogs){
      adj = results_original[j*n_probesets + current_probeset] - results[j*n_probesets + current_probeset];
    } else {
      adj = log(results_original[j*n_probesets + current_probeset])/log(2.0) - log(results[j*n_probesets + current_probeset])/log(2.0);
    }
    for (i =0; i < nprobes; i++){
      data[j*rows  + cur_rows[i]] = pow(2.0,z[j*nprobes + i]  - adj) ;
    }     
  }
  
  Free(z);
}

/*************************************************************
 **
 ** void qnorm_probeset(double *data, int rows, int cols, 
 **                int n_probesets, const char **ProbeNames)
 **
 ** double *data - a matrix of probe intensities.
 ** int rows, cols - dimensions of *data
 ** int n_probesets - number of probesets
 ** const char **ProbeNames - names for each probe.
 **
 **
 ** Perform a quantile normalization, but do it at the 
 ** probeset level. To do this we take the mean or median, of the probeset
 ** normalize these, then adjust probes accordingly.
 **
 **
 **
 **
 ************************************************************/

void qnorm_probeset_c(double *data, int rows, int cols,int n_probesets, const char **ProbeNames, int usemedian, int uselog){
  int j = 0;
  int i = 0;
  int k = 0;
  const char *first;
  int first_ind;
  int max_nrows = 1000;
  /* buffers of size 1000 should be enough. */

  int *cur_rows=Calloc(max_nrows,int);
  int nprobes=0;

  double *cur_exprs = Calloc(cols,double);
  double *cur_SE = Calloc(cols,double);

  double *results = Calloc(cols*n_probesets,double);
  double *results_original = Calloc(cols*n_probesets,double);


  /* Compute the summary measure that will normalized */
  
  first = ProbeNames[0];
  first_ind = 0;
  i =0;
  nprobes = 1;
  for (j = 1; j < rows; j++){
    if ((strcmp(first,ProbeNames[j]) != 0) | (j == (rows -1))){
      if (j == (rows -1)){
	nprobes++;
       	for (k = 0; k < nprobes; k++){
	  if (k >= max_nrows){
	    max_nrows = 2*max_nrows;
	    cur_rows = Realloc(cur_rows, max_nrows, int);
	  }
	  cur_rows[k] = (j+1 - nprobes)+k; 
	}
      } else {
	for (k = 0; k < nprobes; k++){
	  if (k >= max_nrows){
	    max_nrows = 2*max_nrows;
	    cur_rows = Realloc(cur_rows, max_nrows, int);
	  }
	  cur_rows[k] = (j - nprobes)+k; 
	}
      }
      if (uselog){
	if (usemedian){
	  MedianLog_noSE(data, rows, cols, cur_rows, cur_exprs, nprobes);
	} else {
	  AverageLog_noSE(data, rows, cols, cur_rows, cur_exprs, nprobes);
	}
      } else {
	if (usemedian){
	  LogMedian(data, rows, cols, cur_rows, cur_exprs, nprobes,cur_SE);
	} else {
	  LogAverage(data, rows, cols, cur_rows, cur_exprs, nprobes,cur_SE);
	}
      }
      for (k =0; k < cols; k++){
	if (uselog){
	  results[k*n_probesets + i] = cur_exprs[k];
	  results_original[k*n_probesets + i] = cur_exprs[k];
	} else {
	  results[k*n_probesets + i] = pow(2.0,cur_exprs[k]);
	  results_original[k*n_probesets + i] = pow(2.0,cur_exprs[k]);
	}
      } 
      i++;
      first = ProbeNames[j];
      first_ind = j;
      nprobes = 0;
    }
    nprobes++;
  }

  /* Now normalize the summary measure */

  qnorm_c(results, &n_probesets, &cols);

  
  /* Now adjust based upon the normalization to the summary measure */

  first = ProbeNames[0];
  first_ind = 0;
  i =0;
  nprobes = 1;
  for (j = 1; j < rows; j++){
    if ((strcmp(first,ProbeNames[j]) != 0) | (j == (rows -1))){
      if (j == (rows -1)){
	nprobes++;
       	for (k = 0; k < nprobes; k++){
	  if (k >= max_nrows){
	    max_nrows = 2*max_nrows;
	    cur_rows = Realloc(cur_rows, max_nrows, int);
	  }
	  cur_rows[k] = (j+1 - nprobes)+k; 
	}
      } else {
	for (k = 0; k < nprobes; k++){
	  if (k >= max_nrows){
	    max_nrows = 2*max_nrows;
	    cur_rows = Realloc(cur_rows, max_nrows, int);
	  }
	  cur_rows[k] = (j - nprobes)+k; 
	}
      }
      
      AdjustProbes(data, rows, cols, cur_rows, results, results_original, nprobes,n_probesets, i,uselog);
      
      i++;
      first = ProbeNames[j];
      first_ind = j;
      nprobes = 0;
    }
    nprobes++;
  }



  Free(results_original);
  Free(results);
  Free(cur_exprs);
  Free(cur_SE);
  Free(cur_rows);

}

void qnorm_probeset_R(double *data, int *rows, int *cols,int *n_probesets, const char **ProbeNames, int *usemedian, int *uselog){

  qnorm_probeset_c(data, *rows, *cols, *n_probesets, ProbeNames, *usemedian, *uselog);

}
