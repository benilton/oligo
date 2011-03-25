/************************************************************************
 **
 ** file: PLM_biweight.c
 **
 ** Copyright (C) 2003-2007 Ben Bolstad
 ** 
 ** aim: implement the tukey biweight - one step method of summarizing a probeset 
 **
 ** created by: B. M. Bolstad   <bolstad@stat.berkeley.edu>
 ** created on: Jan 7, 2003 (But based on a file mas5.c created in Nov 2002)
 **
 ** last modified: Jan 7, 2003
 **
 ** License: GPL V2 or later (same as the rest of the Affy package)
 **
 ** General discussion
 **
 ** Implement Tukey Biweight Summarization method.
 **
 **
 ** Nov, 2002 - Initial versions
 ** Jan 2, 2003 - Clean up commenting, prepare for integration into AffyExtensions version 0.4
 ** Jan 7, 2003 - make the code a standalone file, data structure manipulation will be handled 
 **               elsewhere.
 ** Jul 23, 2003 - SE parameter added and implemented
 ** Oct 10, 2003 - added in PLM version
 ** Apr 5, 2004 - Change mallocs to Callocs
 ** May 26, 2007 - clean code. Core routines are now in preprocessCore
 ** Sept 9, 2007 - renames this file PLM_biweight.c from biweight.c
 **
 ************************************************************************/

/*#include "threestep_common.h" */

#include "PLM_biweight.h"
#include <biweight.h>
#include "rma_common.h"
#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>





/**********************************************************************************
 **
 ** void tukeybiweight(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute tukey biweight expression measure. 
 **      Note that data is a probes by chips matrix, apply tukeys biweight to columns
 **
 ** double *data - Probe intensity matrix
 ** int rows - number of rows in matrix *data (probes)
 ** int cols - number of cols in matrix *data (chips)
 ** int *cur_rows - indicies of rows corresponding to current probeset
 ** double *results - already allocated location to store expression measures (cols length)
 ** int nprobes - number of probes in current probeset.
 **
 ***********************************************************************************/ 

void TukeyBiweight_threestep(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param){
  
  TukeyBiweight(data, rows, cols, cur_rows, results, nprobes, resultsSE);
}


void TukeyBiweight_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param){
  int i,j;

  TukeyBiweight(data, rows, cols, cur_rows, results, nprobes, resultsSE);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      residuals[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0) - results[j];  
    }
  } 
  
}
