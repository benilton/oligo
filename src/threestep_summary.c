/************************************************************************
 **
 ** file: threestep_summary.c
 ** 
 ** aim: provide a summarization framework for three step computation of expression 
 **
 ** Copyright (C) 2003-2005 Ben Bolstad
 **
 ** created by: B. M. Bolstad  <bolstad@stat.berkeley.edu>
 ** created on: Jan 7, 2002 (based on code extending back as far as June 2002)
 **
 ** last modified: Jan 7, 2003
 **
 ** License: GPL V2 or later (same as the rest of the Affy package)
 **
 ** General discussion
 **
 ** Implement three step summarization.
 **
 ** June, 2002 - RMA code
 ** Oct - Reformulated to deal with changes in R package affy structure
 ** Nov, 2002 - a tukey biweight and a avglog implementation branched off code base
 ** Jan 7, 2003 - Prepare for integration into AffyExtensions version 0.4. A more general framework
 ** Feb 24, 2003 - eliminate a compiler warning (for some compilers) unused variable OLDPM
 ** Apr 4, 2003 - make the number of rows allowed in a probeset dynamic.
 ** Jul 23, 2003 - standard errors should now be returned by threestep summarization routines
 ** May 11, 2004 - fix a small memory leak
 ** Mar 13, 2005 - change loop
 **
 ************************************************************************/

#include "threestep_summary.h"
#include "threestep_summary_methods.h"
#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/************************************************************************************
 **
 ** void do_summary(double *PM, char **ProbeNames, int *rows, int *cols, double *results, 
 **                 char **outnames, int nps, void (* SummaryMeth)(double*, int, int, int *,double *, int))
 **
 ** double *PM - matrix of dimension rows by cols (probes by chips) should already be 
 **              normalized and background corrected.
 ** char **ProbeNames - Probeset names, one for each probe.
 ** int *rows, *cols - dimensions of matrix
 ** double *results - space already allocated to store chip expression estimates
 ** char **outnames - space already allocated to store probeset names
 ** int nps - number of probesets
 ** void (*SummaryMeth)(double*, int, int, int *,double *, int, double *) - pointer to the summarization function
 ** double *resultsSE - space already allocated to store SE of chip expression estimates
 **
 ** PM should be background corrected and normalized
 ** assumed that Probes are sorted, by ProbeNames, so that we can just look at 
 ** consecutive rows in PM matrix when doing summarized
 **
 ** each item is then used to create a matrix that is summarized to give
 ** expression estimates.
 **
 ************************************************************************************/

void do_3summary(double *PM, const char **ProbeNames, int *rows, int *cols, double *results, char **outNames, int nps,void (* SummaryMeth)(double*, int, int, int *,double *, int, double *, summary_plist *),double *resultsSE, summary_plist *summary_param){
  int j = 0;
  int i = 0;
  int k = 0;
  int size;
  const char *first;
  int first_ind;
  int max_nrows = 1000;
  /* buffers of size 1000 should be enough. */

  int *cur_rows=Calloc(max_nrows,int);
  int nprobes=0;

  double *cur_exprs = Calloc(*cols,double);
  double *cur_se = Calloc(*cols,double);
  /* double *OLDPM = NULL;*/

  first = ProbeNames[0];
  first_ind = 0;
  i = 0;     /* indexes current probeset */
  j = 0;    /* indexes current row in PM matrix */
  k = 0;    /* indexes current probe in probeset */
  
  while ( j < *rows){
    if (strcmp(first,ProbeNames[j]) == 0){
      if (k >= max_nrows){
	max_nrows = 2*max_nrows;
	cur_rows = Realloc(cur_rows, max_nrows, int);
      }
      cur_rows[k] = j;
      k++;
      j++;
      
    } else {
      nprobes = k;
      SummaryMeth(PM, *rows, *cols, cur_rows, cur_exprs, nprobes,cur_se,summary_param);
      for (k =0; k < *cols; k++){
	results[k*nps + i] = cur_exprs[k];
	resultsSE[k*nps + i] = cur_se[k];
      } 
      size = strlen(first);
      outNames[i] = Calloc(size+1,char);
      strcpy(outNames[i],first);
      i++;
      first = ProbeNames[j];
      k = 0;
    }
  }
  nprobes = k;
  SummaryMeth(PM, *rows, *cols, cur_rows, cur_exprs, nprobes,cur_se,summary_param);
  for (k =0; k < *cols; k++){
    results[k*nps + i] = cur_exprs[k];
    resultsSE[k*nps + i] = cur_se[k];
  } 
  size = strlen(first);
  outNames[i] = Calloc(size+1,char);
  strcpy(outNames[i],first);
  
  
  Free(cur_exprs);
  Free(cur_se);
  Free(cur_rows);
}
