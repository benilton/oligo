/************************************************************************
 **
 ** file: rma.c
 **
 ** Copyright (C) 2002 - 2003   B. M. Bolstad
 **
 ** created by: B. M. Bolstad   <bolstad@stat.berkeley.edu>
 ** created on: June 26, 2002
 **
 ** last modified: January 6, 2003
 ** 
 ** last modified: Apr 4, 2003
 **
 ** License: GPL V2 or later (same as the rest of the Affy package)
 **
 ** version 1.1 - Initial release to affy package
 **
 ** Version History (LEADING UP TO AND INCLUDING AFFYEXTENSIONS)
 ** 0.1 - Initial version released on July 14, 2002. Implements median
 **       polish RMA method with 
 ** 0.2 - background implemented in c with the density estimation still carried
 **       out by the R function density()
 ** 0.25 - correct background implementation, version 0.2 is broken.
 **        background is implemented in rma_background.c
 ** 0.30 - Have a copy and none copy path. ie we can either work inplace or on
 **        duplicates. the purpose of this is to reduce memory overhea. For 
 **        someone with an interest only in expression estimates this should not be a problem
 ** 
 ** Version History (AFTER INCLUSION INTO AFFY PACKAGE)
 ** 1.1 - Initial inclusion into Affy package, heavy modification to how PM data structure
 **       dealt with.
 **
 ** OLD COMMENTS
 **
 ** a c language implementation of the RMA method as given in the RMA.R file I 
 ** received from Rafael at an earlier point, but assume already had background 
 ** correction to PM's at somepoint (perhaps in the c code) bg will be written in later.
 ** Possibly another background method will be inserted in at this stage. <-- COMMENT DEPRECIATED
 ** 
 ** Note that the normalization code that is used in this algorithm is updated
 ** from that in the affy version 1.1.1 (there will be slight differences in the 
 ** expression values because of this), there is also slight differences in the 
 ** ordering of the results.  <-- THIS COMMENT IS DEPRECIATED. Quantile
 ** normalization updates will happen in the bioconductor cvs.
 **
 ** Ideally and at some later point a more modular approach that can be called 
 ** in a better manner from R than this will be written. This is a quick and 
 ** dirty approach to get something that will run acceptably.in terms of memory 
 ** and speed. From a software engineering viewpoint expect the code to be poorly 
 ** constructed.  <-- SOMEWHAT DEPRECIATED. some work should be done to
 ** clean things up. The user will generally only be dealing with the R 
 ** interface.
 **
 ** Input to the function should be processed from within R
 **
 ** NEW COMMENTS
 **
 ** This is the main c function for implementing the RMA method
 ** it provides c interfaces to be called from R.
 **
 ** Specific Modification History
 **
 ** Note that the qnorm code here will not be the development tree
 ** LG: what do you mean ?
 **
 ** BMB: legacy comment, from when this code was outside affy, in AffyExtensions
 **      and before that as raw c code that was floating around.
 **
 ** Specific Modification History
 ** Nov 2, 2002 - modify so that it will work efficently with affy2
 ** Nov 3, 2002 - More modifications, remove cruft from old version
 ** Nov 4, 2002 - testing, check docs etc
 ** Nov 10,2002 - remove pesky debug printf()
 ** Dec 5, 2002 - add ability to turn background off
 ** Dec 31, 2002 - add ability to change to type 2 background
 ** Jan 2, 2003 - clean up old/incorrect documentation/comments
 **
 ** Dec 26, 2002 - '//' is not a valid way to comment out (and some C compilers complain about it)
 **                (Laurent)
 ** Jan 6, 2003 - fix merging. Note "//" is valid according to the language standards (http://anubis.dkuug.dk/jtc1/sc22/open/n2794/n2794.txt)
 ** Feb 6, 2003 - change some printfs to Rprintfs this will allow the windows users to see some
 **               verbage when running rma
 ** Feb 25, 2003 - try to reduce or eliminate compiler warnings (from gcc -Wall) 
 ** Apr 4, 2003 - fix up so that the number of probes in a probeset is allowed to be more dynamic
 ** Dec 9, 2003 - fix a bug in do_RMA (max_nrows in Calloc)
 ** Mar 6, 2004 - all mallocs/frees are now Calloc/Frees. Removed
 **               the function R_median_polish
 ** Jul 27, 2004 - fix a small memory leak
 ** Aug 4, 2004 - move the "Background correcting" message. 
 ** Nov 8, 2004 - change how things are structured in do_RMA()
 **
 ************************************************************************/


/* #include "rma_structures.h" */
#include "rma_common.h"
#include "rma_background2.h"
#include "qnorm.h" 

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void do_RMA(double *PM, char **ProbeNames, int *rows, int * cols,double *results,char **outNames,int nps);


/**************************************************************************
 **
 ** double median(double *x, int length)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the median of *x
 **
 *************************************************************************/

double  median(double *x, int length){
  int i;
  int half;
  double med;
  double *buffer = Calloc(length,double);
  
  for (i = 0; i < length; i++)
    buffer[i] = x[i];
  
  qsort(buffer,length,sizeof(double), (int(*)(const void*, const void*))sort_double);
  half = (length + 1)/2;
  if (length % 2 == 1){
    med = buffer[half - 1];
  } else {
    med = (buffer[half] + buffer[half-1])/2.0;
  }
  
  Free(buffer);
  return med;
}

/*******************************************************************************
 **
 ** double sum_abs(double *z, int rows, int cols)
 **
 ** double *z - matrix of doubles
 ** int rows - dimension of matrix
 ** int cols - dimension of matrix
 **
 ** returns the sum of the absolute values of elements of the matrix *z
 **
 ******************************************************************************/

double sum_abs(double *z, int rows, int cols){
 
  int i, j;
  double sum = 0.0;

  for (i=0; i < rows; i++)
    for (j=0; j < cols; j++)
      sum+=fabs(z[j*rows+i]);

  return sum;
}

/********************************************************************************
 **
 ** void get_row_median(double *z, double *rdelta, int rows, int cols)
 **
 ** double *z - matrix of dimension  rows*cols
 ** double *rdelta - on output will contain row medians (vector of length rows)
 ** int rows, cols - dimesion of matrix
 **
 ** get the row medians of a matrix 
 **
 ********************************************************************************/

void get_row_median(double *z, double *rdelta, int rows, int cols){
  int i,j;
  double *buffer = Calloc(cols,double);

  for (i = 0; i < rows; i++){ 
    for (j = 0; j < cols; j++){
      buffer[j] = z[j*rows + i];
    }
    rdelta[i] = median(buffer,cols);
  }
  
  Free(buffer);
}

/********************************************************************************
 **
 ** void get_col_median(double *z, double *cdelta, int rows, int cols)
 **
 ** double *z - matrix of dimension  rows*cols
 ** double *cdelta - on output will contain col medians (vector of length cols)
 ** int rows, cols - dimesion of matrix
 **
 ** get the col medians of a matrix 
 **
 ********************************************************************************/

void get_col_median(double *z, double *cdelta, int rows, int cols){
  
  int i, j;
  
  double *buffer = Calloc(rows,double);
  for (j = 0; j < cols; j++){
    for (i = 0; i < rows; i++){  
      buffer[i] = z[j*rows + i];
    }
    cdelta[j] = median(buffer,rows);
  }
  
  Free(buffer);

}

/***********************************************************************************
 **
 ** void subtract_by_row(double *z, double *rdelta, int rows, int cols)
 ** 
 ** double *z - matrix of dimension rows by cols
 ** double *rdelta - vector of length rows
 ** int rows, cols dimensions of matrix
 **
 ** subtract the elements of *rdelta off each row of *z
 **
 ***********************************************************************************/

void subtract_by_row(double *z, double *rdelta, int rows, int cols){
  
  int i,j;

  for (i = 0; i < rows; i++){
    for (j = 0; j < cols; j++){
      z[j*rows +i]-= rdelta[i];
    }
  }
}


/***********************************************************************************
 **
 ** void subtract_by_col(double *z, double *cdelta, int rows, int cols)
 ** 
 ** double *z - matrix of dimension rows by cols
 ** double *cdelta - vector of length rows
 ** int rows, cols dimensions of matrix
 **
 ** subtract the elements of *cdelta off each col of *z
 **
 ***********************************************************************************/

void subtract_by_col(double *z, double *cdelta, int rows, int cols){
  
  int i,j;
  for (j = 0; j < cols; j++){
    for (i = 0; i < rows; i++){
      z[j*rows +i]-= cdelta[j];
    }
  }

}

/***********************************************************************************
 **
 ** void rmod(double *r, double *rdelta, int rows)
 ** 
 ** double *r - vector of length rows
 ** double *rdelta - vector of length rows
 ** int rows, cols dimensions of matrix
 **
 ** add elementwise *rdelta to *r
 **
 ***********************************************************************************/


void rmod(double *r, double *rdelta, int rows){
  int i;

  for (i = 0; i < rows; i++){
    r[i]= r[i] + rdelta[i];
  }
}

/***********************************************************************************
 **
 ** void cmod(double *c, double *cdelta, int cols)
 ** 
 ** double *c - vector of length rows
 ** double *cdelta - vector of length rows
 ** int cols length of vector
 **
 ** add elementwise *cdelta to *c
 **
 ***********************************************************************************/

void cmod(double *c, double *cdelta, int cols){
  int j;

  for (j = 0; j < cols; j++){
    c[j]= c[j] + cdelta[j];
  }
}


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
 ** a function to do median polish.
 **
 *************************************************************************************/

void median_polish(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes){

  int i,j,iter;
  int maxiter = 10;
  double eps=0.01;
  double oldsum = 0.0,newsum = 0.0;
  double t = 0.0;
  double delta;
  double *rdelta = Calloc(nprobes,double);
  double *cdelta = Calloc(cols,double);
  
  double *r = Calloc(nprobes,double);
  double *c = Calloc(cols,double);
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
  
  for (iter = 1; iter <= maxiter; iter++){
    get_row_median(z,rdelta,nprobes,cols);
    subtract_by_row(z,rdelta,nprobes,cols);
    rmod(r,rdelta,nprobes);
    delta = median(c,cols);
    for (j = 0; j < cols; j++){
      c[j] = c[j] - delta;
    }
    t = t + delta;
    get_col_median(z,cdelta,nprobes,cols);
    subtract_by_col(z,cdelta,nprobes,cols);
    cmod(c,cdelta,cols);
    delta = median(r,nprobes);
    for (i =0; i < nprobes; i ++){
      r[i] = r[i] - delta;
    }
    t = t+delta;
    newsum = sum_abs(z,nprobes,cols);
    if (newsum == 0.0 || fabs(1.0 - oldsum/newsum) < eps)
      break;
    oldsum = newsum;
  }
  
  for (j=0; j < cols; j++){
    results[j] =  t + c[j]; 
  }
  
  Free(rdelta);
  Free(cdelta);
  Free(r);
  Free(c);
  Free(z); 
}

/************************************************************************************
 **
 **  void do_RMA(double *PM, char **ProbeNames, int *rows, int * cols)
 **
 ** double *PM - matrix of dimension rows by cols (probes by chips) should already be 
 **              normalized and background corrected.
 ** char **ProbeNames - Probeset names, one for each probe.
 ** int *rows, *cols - dimensions of matrix
 **
 ** perform the multichip averaging. PM should be background corrected and normalized
 **
 ** assumed that Probes are sorted, by ProbeNames, so that we can just look at 
 ** consecutive rows in PM matrix when doing the median polish
 **
 ** each item is then used to create a matrix that is median polished to give
 ** expression estimates.
 **
 ************************************************************************************/

void do_RMA(double *PM, char **ProbeNames, int *rows, int *cols, double *results, char **outNames, int nps){
  int j = 0;
  int i = 0;
  int k = 0;
  int size;
  char *first;
  int first_ind;
  int max_nrows = 1000;


  /* buffers of size 200 should be enough. */

  int *cur_rows=Calloc(max_nrows,int);
  int nprobes=0;

  double *cur_exprs = Calloc(*cols,double);

  /* double *OLDPM = NULL; */

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
      median_polish(PM, *rows, *cols, cur_rows, cur_exprs, nprobes);
      for (k =0; k < *cols; k++){
	results[k*nps + i] = cur_exprs[k];
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
  median_polish(PM, *rows, *cols, cur_rows, cur_exprs, nprobes);
  for (k =0; k < *cols; k++){
    results[k*nps + i] = cur_exprs[k];
  } 
  size = strlen(first);
  outNames[i] = Calloc(size+1,char);
  strcpy(outNames[i],first);
  

  Free(cur_exprs);
  Free(cur_rows);
}

/********************************************************************************************
 **
 ** void rma_c_call(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_flag)
 **
 ** SEXP PMmat - matrix of Perfect-match values
 ** SEXP MMmat - matrix of Mismatch values
 ** SEXP ProbeNamesVec - vector containing names of probeset for each probe
 ** SEXP N_probes - number of PM/MM probes on an array
 ** SEXP norm_flag  - non zero for use quantile normalization, 0 for no normalization
 **
 ** a function to actually carry out the RMA method taking the R objects and manipulating
 ** into C data structures.
 **
 ** this function assumes any sort of background correction was carried out previously
 ** This function carries out the other two steps of the RMA algorithm:
 ** Normalization and Summarization.
 **
 ** In particular the data is quantile normalized and then it is
 ** summarized using median polish
 **
 *******************************************************************************************/

SEXP rma_c_call(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_flag){
  
  int rows, cols;
  double *outexpr;
  double *PM,*MM;
  char **outnames;
  char **ProbeNames;
  int i,nprobesets;
  


  SEXP dim1;
  SEXP outvec; /* ,outnamesvec; */
  SEXP dimnames,names;
  
  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1]; 

  PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
  MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));
  
  nprobesets=INTEGER(N_probes)[0];
  
  /*  printf("%i\n",nprobesets); */
  /* printf("%d ",INTEGER(norm_flag)[0]); */
  if (INTEGER(norm_flag)[0]){
  /* normalize PM using quantile normalization */
  /*  printf("Normalizing\n"); */
    Rprintf("Normalizing\n");
    qnorm_c(PM,&rows,&cols);
  }

  ProbeNames = Calloc(rows,char *);

  for (i =0; i < rows; i++)
    ProbeNames[i] = CHAR(VECTOR_ELT(ProbeNamesVec,i));
  
  
  outnames = Calloc(nprobesets,char *);

  /* PROTECT(outvec = NEW_NUMERIC(nprobesets*cols)); */
  
  PROTECT(outvec = allocMatrix(REALSXP, nprobesets, cols));


  outexpr = NUMERIC_POINTER(outvec);
 	    
  /* printf("Calculating Expression\n"); */
  Rprintf("Calculating Expression\n");


  do_RMA(PM, ProbeNames, &rows, &cols,outexpr,outnames,nprobesets);

  UNPROTECT(2);

  /* now lets put names on the matrix */

  PROTECT(dimnames = allocVector(VECSXP,2));
  PROTECT(names = allocVector(STRSXP,nprobesets));
  
  for ( i =0; i < nprobesets; i++)
    SET_VECTOR_ELT(names,i,mkChar(outnames[i]));
  
  SET_VECTOR_ELT(dimnames,0,names);
  setAttrib(outvec, R_DimNamesSymbol, dimnames);
  UNPROTECT(2);
  for (i =0; i < nprobesets; i++)
    Free(outnames[i]);
  
  Free(outnames);
  Free(ProbeNames);
  return outvec;
}

/*******************************************************************************************************************
 **
 ** SEXP rma_c_complete(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP densfunc, SEXP rho)
 **
 ** SEXP PMmat  - PM's
 ** SEXP MMmat - MM's
 ** SEXP ProbeNamesVec - names of probeset for each row
 ** SEXP N_probes - number of probesets
 ** SEXP densfunc - density function to use in computation of background
 ** SEXP rho - an R environment 
 ** SEXP norm_flag - TRUE/FALSE or 1/0 for normalize/not
 ** SEXP bg_flag - TRUE/FALSE  or 1/0 for background correct/not
 ** SEXP bg_type - integer indicating "RMA" background to use. 2 is equivalent to bg.correct.rma in affy 1.1.1
 **                all other values default to 1.0.2 "RMA" background
 ** 
 ** Main function to be called from R. Modifies the PM matrix from the parent environment. More dangerous than the
 ** function below, but less memory intensive. This is a function that implements the complete RMA method. ie
 ** background correction, quantile normalization, then expression summarization using median polish
 **
 *******************************************************************************************************************/

SEXP rma_c_complete(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP densfunc, SEXP rho,SEXP norm_flag, SEXP bg_flag, SEXP bg_type){
  if (INTEGER(bg_flag)[0]){
    Rprintf("Background correcting\n");
    PMmat = bg_correct_c(PMmat,MMmat,densfunc,rho,bg_type);
  }
  return rma_c_call(PMmat, MMmat, ProbeNamesVec,N_probes,norm_flag);
}

/********************************************************************************************************************
 **
 ** SEXP rma_c_complete_copy(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP densfunc, SEXP rho,SEXP norm_flag, SEXP bg_flag)
 **
 ** SEXP PMmat   - PM's
 ** SEXP MMmat   - MM's
 ** SEXP ProbeNamesVec - names of probeset for each row
 ** SEXP N_probes  - number of probesets
 ** SEXP densfunc - density function to use in computation of background
 ** SEXP rho - an r environment to work within when doing density call in background step
 ** SEXP norm_flag - TRUE/FALSE or 1/0 for normalize/not
 ** SEXP bg_flag - TRUE/FALSE  or 1/0 for background correct/not
 ** SEXP bg_type - integer indicating "RMA" background to use. 2 is equivalent to bg.correct.rma in affy 1.1.1
 **                all other values default to 1.0.2 "RMA" background
 *
 ** Main function to be called from R. Makes a copy of the PM matrix and then works with that. Safer than the 
 ** other function above, but more memory intensive. This is the function that implements the complete RMA method.
 ** ie background correction, quantile normalization, then expression summarization using median polish
 **
 ********************************************************************************************************************/

SEXP rma_c_complete_copy(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP densfunc, SEXP rho,SEXP norm_flag, SEXP bg_flag, SEXP bg_type){
 SEXP dim1,PMcopy,exprs;
 int rows,cols;

 if (INTEGER(bg_flag)[0]){
   Rprintf("Background correcting\n");
   PMmat = bg_correct_c_copy(PMmat,MMmat,densfunc,rho, bg_type); 
   return rma_c_call(PMmat, MMmat, ProbeNamesVec,N_probes,norm_flag);
  } else {
    PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol));
    rows = INTEGER(dim1)[0];
    cols = INTEGER(dim1)[1];
    PROTECT(PMcopy = allocMatrix(REALSXP,rows,cols));
    copyMatrix(PMcopy,PMmat,0);
    exprs = rma_c_call(PMcopy, MMmat, ProbeNamesVec,N_probes,norm_flag);
    UNPROTECT(2);
    return exprs;
  }
}
