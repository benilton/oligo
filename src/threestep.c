/************************************************************************
 **
 ** file: threestep.c (based on earlier work in rma.c, mas5.c, avglog.c)
 **
 ** aim: provide the general framework for threestep analysis
 **      and provide interface to R code
 **
 ** Copyright (C) 2003-2005 Ben Bolstad
 **
 ** created by: B. M. Bolstad   <bolstad@stat.berkeley.edu>
 ** created on: Jan 7, 2003 (rma.c dates to June 26, 2002)
 **
 ** last modified: January 17, 2003
 **
 ** License: GPL V2 or later
 **
 ** 
 ** Specific Modification History
 **
 ** Jan 7, 2003 - Initial version
 ** Jan 9, 2003 - 0.4 built and released (still needs some background methods and ability to change
 **               expression measures)
 ** Jan 9, 2003 - Made summary function pointer array
 ** Jan 11, 2003 - Clean up summary function selection by setting up own file, fix up stack imbalance call in IdealMM.
 ** Jan 17, 2003 - moved threestep_background and threestep_normalize to preprocess.c
 **                and renamed these functions pp_background  and pp_normalize
 ** Feb 24, 2003 - remove unused variable to eliminate compiler warning.
 ** Jul 23, 2003 - added a little more documentation, standard errors from
 **                threestep. eliminate the copy step. It is not
 **                required under the current setup.
 ** Oct 5, 2003 - SEXP summary_parameters added
 ** Apr 5, 2004 - all malloc/free should now be Calloc/Free
 ** Sep 13, 2005 - fix a possible gc() situation
 ** March 1, 2006 - change comments to ansi c style
 ** Oct 10, 2006 - add verbosity arguments 
 ** Jan 6, 2009 - change SET_VECTOR_ELT to SET_STRING_ELT where relevant.
 **
 ************************************************************************/

#include "rma_common.h" 


/*#include "medianpolish.h"
#include "biweight.h"
#include "avg_log.h" */

#include "preprocess.h"
#include "threestep_summary.h"
#include "threestep_summary_methods.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/********************************************************************************************
 **
 ** void threestep_summary(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP summary_type)
 **
 ** SEXP PMmat - matrix of Perfect-match values
 ** SEXP MMmat - matrix of Mismatch values
 ** SEXP ProbeNamesVec - vector containing names of probeset for each probe
 ** SEXP N_probes - number of PM/MM probes on an array
 ** SEXP summary_type  - an integer indicating the summary statistic to be used.
 **
 ** a function to manipulate the data so that the R objects are converted
 ** into C data structures. The summary statistic is then applied to the C data structures. 
 ** The summary statistic is applied only to PM data.
 **
 ** this function assumes any sort of background correction and normalization was carried out previously
 **
 *******************************************************************************************/

SEXP threestep_summary(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP summary_type, SEXP summary_parameters,SEXP verbosity){
  
  int rows, cols;
  double *outexpr, *outSE;
  double *PM,*MM;
  char **outnames;
  const char **ProbeNames;
  int i,nprobesets;
  int Method;
  int verbosity_level;

  summary_plist *summary_param = (summary_plist *)Calloc(1,summary_plist);

  SEXP dim1;
  SEXP outvec, outSEvec, output_list;        /*outnamesvec */
  SEXP dimnames,names;
  SEXP cur_param;

  
  /*
  pt2Summary funcArr[3];
  funcArr[0] = &median_polish;
  funcArr[1] = &tukeybiweight;
  funcArr[2] = &AverageLog; */

  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1]; 

  PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
  MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));
   
  verbosity_level = asInteger(verbosity);


  nprobesets=INTEGER(N_probes)[0];
  
  ProbeNames = (const char **)Calloc(rows,const char *);
  for (i =0; i < rows; i++)
    ProbeNames[i] = CHAR(STRING_ELT(ProbeNamesVec,i));
  outnames = (char **)Calloc(nprobesets,char *);

  /* PROTECT(outvec = NEW_NUMERIC(nprobesets*cols)); */
  PROTECT(outvec = allocMatrix(REALSXP, nprobesets, cols));
  outexpr = NUMERIC_POINTER(outvec);
  
  PROTECT(outSEvec = allocMatrix(REALSXP, nprobesets, cols));
  outSE = NUMERIC_POINTER(outSEvec);


  
  Method = asInteger(summary_type)-1;
  /*printf("%d ",asInteger(summary_type));*/

  cur_param =  GetParameter(summary_parameters, "psi.k");
  summary_param->psi_k = NUMERIC_POINTER(cur_param)[0];
  
  cur_param =  GetParameter(summary_parameters, "psi.type");
  summary_param->psi_method = asInteger(cur_param);


  if (verbosity_level > 0){
    Rprintf("Calculating Expression\n");
  }
  do_3summary(PM, ProbeNames, &rows, &cols,outexpr,outnames,nprobesets,SummaryMethod(Method),outSE,summary_param);
  

  /* now lets put names on the matrix */

  PROTECT(dimnames = allocVector(VECSXP,2));
  PROTECT(names = allocVector(STRSXP,nprobesets));
  for ( i =0; i < nprobesets; i++)
    SET_STRING_ELT(names,i,mkChar(outnames[i]));
  SET_VECTOR_ELT(dimnames,0,names);
  setAttrib(outvec, R_DimNamesSymbol, dimnames);
  setAttrib(outSEvec, R_DimNamesSymbol, dimnames);
  UNPROTECT(2);
  PROTECT(output_list = allocVector(VECSXP,2));

  SET_VECTOR_ELT(output_list,0,outvec);
  SET_VECTOR_ELT(output_list,1,outSEvec); 
  UNPROTECT(3);
  UNPROTECT(1);
  
  for (i =0; i < nprobesets; i++)
    Free(outnames[i]);

  Free(outnames);
  Free(ProbeNames);
  Free(summary_param);
  return output_list;
} 





/********************************************************************************************
 **
 ** SEXP R_threestep_c(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP densfunc, 
 **                    SEXP rho,SEXP norm_flag, SEXP bg_flag, SEXP bg_type,SEXP norm_type, 
 **                    SEXP summary_type,SEXP LESN_params)
 **
 ** 
 **
 ** This is the function that is called from R and interfaces with the C code. It produces
 ** expression estimates by the threestep procedure of background/signal adjustment then
 ** between chip normalization, followed by summarization.
 **
 *******************************************************************************************/


SEXP R_threestep_c(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_flag, SEXP bg_flag, SEXP bg_type,SEXP norm_type, SEXP summary_type,SEXP background_parameters,SEXP norm_parameters, SEXP summary_parameters, SEXP verbosity){
  
  SEXP dim1,PMcopy,exprs;
  int rows,cols;
  
  /*Create a copy matrix to work on. Allows us to modify data in background and normalization steps without affecting original data */
  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  /*PROTECT(PMcopy = allocMatrix(REALSXP,rows,cols));
    copyMatrix(PMcopy,PMmat,0); */
  PMcopy = PMmat;

  /*  printf("%d ",INTEGER(norm_flag)[0]);
  printf("%d ",INTEGER(bg_flag)[0]);
  printf("%d ",asInteger(bg_type));
  printf("%d ",asInteger(norm_type));
  printf("%d ",asInteger(summary_type));
  */

  /* If Background correction do it */
  if (INTEGER(bg_flag)[0]){
    PMcopy = pp_background(PMcopy, MMmat, ProbeNamesVec,N_probes,bg_type, background_parameters,verbosity); /* densfunc, rho, LESN_params); */
  }

  /* If Normalization do it */
  if (INTEGER(norm_flag)[0]){
    PMcopy = pp_normalize(PMcopy, MMmat, ProbeNamesVec,N_probes,norm_type, norm_parameters,verbosity);
  }

  /* Do Summarization */  
  exprs = threestep_summary(PMcopy, MMmat, ProbeNamesVec,N_probes,summary_type,summary_parameters,verbosity);
  UNPROTECT(1);
  return exprs;
  
}
