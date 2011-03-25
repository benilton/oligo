/*********************************************************************
 **
 ** file: preprocess.c
 **
 ** Aim: Have implementations of the background and normalization 
 **      steps 
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Jan 17, 2003
 **
 ** Last modified: Oct 12, 2006
 **
 ** This file contains the normalization and background preprocessing
 ** steps.
 **
 ** Modification history
 **
 ** Jan 17, 2003 - Initial version. Moved background and normalization
 **                steps from threestep.c to this file. renamed from
 **                threestep_* to pp_*
 **  
 ** Feb 6, 2003 - change printf to Rprintf
 ** Feb 12, 2003 - add a missing #include "preprocess.h"
 **                clean up documentation
 **                Put in some code, that should not yet be executed to allow
 **                a call to the MAS 5.0 style background.
 ** Feb 24, 2003 - cleanup declaration of variable that is not used
 ** Mar 21, 2003 - Add ability to call LESN background methods, this required
 **                adding additional parameters to pp_background function for 
 **                LESN parameters.
 ** Mar 23, 2003 - comment out printf in MAS background, this will be handled by R code.
 ** Jun 23, 2003 - Add quantiles_probeset in normalization options.
 ** Jul 24, 2003 - scaling normalization added as normalization option
 ** Jul 26, 2003 - more general framework for providing parameters to 
 **                normalization. The function GetParameter was introduced.
 **                A similar change was made to the background code
 **                so that SEXP densfunc, SEXP rho, SEXP LESN_param
 **                were all subsumed into bg_parameters
 ** Apr 5, 2004 - all malloc/free are now Calloc/Free
 ** Aug 4, 2004 - Change functions to deal with new structure and pmonly/mmonly/separate/together methodology.
 ** Apr 27, 2006 - add "quantile.robust" to normalization methods
 ** Jul 10, 2006 - add log.scalefactors support into scaling normalization.
 ** Oct 10, 2006 - add verbosity argument to functions. Higher levels of verbosity give more text output to the screen
 ** Oct 12, 2006 - add a function that does both background and normalization
 ** Mar 31, 2008 - make RMA background use preprocessCore implementation
 ** Jan 6, 2009 - change SET_VECTOR_ELT to SET_STRING_ELT where relevant.
 **
 *********************************************************************/

#include "rma_common.h" 
#include "rma_background4.h" 
//#include "qnorm.h"
#include "preprocessCore_normalization_stubs.c"
#include "idealmismatch.h"
#include "LESN.h"
#include "preprocess.h"
#include "qnorm_probeset.h"
#include "scaling.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "preprocessCore_background_stubs.c"

/********************************************************************************************
 **
 ** SEXP GetParameter(SEXP alist, char *param_name)
 **
 **
 **
 ** given a list find the parameter with the given name or halt because parameter not
 ** found.
 **
 ********************************************************************************************/

SEXP GetParameter(SEXP alist, char *param_name){

  int length,i;
  SEXP Names;

  if (isVectorList(alist) ==FALSE){
    error("Parameter list was not list.");
  }
  
  length = length(alist);
  Names = GET_NAMES(alist);

  if (length(Names) != length){
    error("Need names for all items in parameter list.");
  }

  for( i=0; i < length; i++){
    if (strcmp(CHAR(STRING_ELT(Names,i)),param_name) == 0){
      break;
    }
  }

  if (i == length){
    error("Did not find %s in parameter list.", param_name);
  }
  
  return VECTOR_ELT(alist,i);

  
}



/********************************************************************************************
 **
 ** void pp_background(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_type, SEXP verbosity)
 **
 ** SEXP PMmat - matrix of Perfect-match values
 ** SEXP MMmat - matrix of Mismatch values
 ** SEXP ProbeNamesVec - vector containing names of probeset for each probe
 ** SEXP N_probes - number of PM/MM probes on an array
 ** SEXP bg_type  - an integer indicating the background method to be used.
 ** SEXP densfunc - an R function for computing non parameteric density in RMA background
 ** SEXP rho - an R environment used in computation of RMA background
 ** SEXP LESN_param - a vector of two elements containing parameters for LESN method
 **                   first parameter is p0, the second is theta  - ADDED Mar 21, 2003
 **
 ** Carry out the background/PM adjustment preprocessing step.
 **
 ** At some point this will be generalized so that it is even easier to add new
 ** and/or different background methods.
 **
 ** bg_type definitions
 **
 ** 1 - RMA version 1 background (this is what is used in the normalization paper)
 ** 2 - RMA version 2 background (used in the NAR paper)
 ** 3 - use Ideal Mismatch Idea of MAS 5.0 to subtract IM from PM
 ** 4 - use an MAS 5.0 style within grid background correction method
 ** 5 - use both MAS 5 style followed by Ideal Mismatch correction
 ** 6 - use LESN proposal 2, half gaussian correction
 ** 7 - use LESN proposal 1, exponential correction
 ** 8 - use LESN proposal 0, just shift intensities lower so that lowest value is small
 **
 ********************************************************************************************/

SEXP pp_background(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP bg_type,SEXP background_param, SEXP verbosity){
  int i,j;
  double *PM,*MM,*PMMM;
  
  double theta, baseline;
  int rows, cols;
  const char **ProbeNames;
  
  int allrows;
  int which_lesn;
  int verbosity_level;
  
  SEXP dim1;
  /*
    SEXP densfunc;
    SEXP rho;
  */
  SEXP LESN_param;
  SEXP param;
  
  SEXP rma_bg_type;
  SEXP allPMMM;
  

  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1]; 
  UNPROTECT(1);

  PROTECT(rma_bg_type = allocVector(REALSXP,10));
  NUMERIC_POINTER(rma_bg_type)[0] = 2.0;
  

  verbosity_level = asInteger(verbosity);


  if (strcmp(CHAR(STRING_ELT(bg_type,0)),"RMA.2") == 0){
    param = GetParameter(background_param,"type");
    PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
    MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));
    

    if ((strcmp(CHAR(STRING_ELT(param,0)),"pmonly") == 0) || (strcmp(CHAR(STRING_ELT(param,0)),"separate") == 0)){
      if (verbosity_level > 0){
	Rprintf("Background correcting PM\n");
      }
      rma_bg_correct(PM, rows, cols);

    }
    if ((strcmp(CHAR(STRING_ELT(param,0)),"mmonly") == 0) || (strcmp(CHAR(STRING_ELT(param,0)),"separate") == 0)){
      if (verbosity_level > 0){
	Rprintf("Background correcting MM\n");
      }
      rma_bg_correct(MM, rows, cols);
    }
    if (strcmp(CHAR(STRING_ELT(param,0)),"together") == 0){ 
      if (verbosity_level > 0){
	Rprintf("Background correcting PM and MM together\n");
      }
      PROTECT(allPMMM = allocMatrix(REALSXP,2*rows,cols));
      allrows = 2*rows;
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  NUMERIC_POINTER(allPMMM)[j*2*rows + i] = NUMERIC_POINTER(PMmat)[j*rows + i];
	}
      }
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  NUMERIC_POINTER(allPMMM)[j*2*rows + i + rows] = NUMERIC_POINTER(MMmat)[j*rows + i];
	}
      }
      PMMM = NUMERIC_POINTER(AS_NUMERIC(allPMMM));
      rma_bg_correct(PMMM, 2*rows, cols);
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  NUMERIC_POINTER(PMmat)[j*rows + i] = NUMERIC_POINTER(allPMMM)[j*2*rows + i];
	}
      }
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  NUMERIC_POINTER(MMmat)[j*rows + i] = NUMERIC_POINTER(allPMMM)[j*2*rows + i + rows];
	}
      }
      UNPROTECT(1);
    }
  } else if ((strcmp(CHAR(STRING_ELT(bg_type,0)),"IdealMM") == 0 )|| (strcmp(CHAR(STRING_ELT(bg_type,0)),"MASIM") == 0)){
    if (verbosity_level > 0){
      Rprintf("Background correcting");
    }
    PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
    rows = INTEGER(dim1)[0];
    cols = INTEGER(dim1)[1]; 
    
    PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
    MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));

    ProbeNames =(const char **)Calloc(rows,const char *);
    for (i =0; i < rows; i++)
      ProbeNames[i] = CHAR(STRING_ELT(ProbeNamesVec,i));
    param = GetParameter(background_param,"ideal");
    if (strcmp(CHAR(STRING_ELT(param,0)),"MM") == 0){
      if (verbosity_level > 0){
	Rprintf(" PM using MM\n");
      }
      IdealMM_correct(PM,MM, &rows, &cols,ProbeNames);
    } else {
      if (verbosity_level > 0){
	Rprintf(" MM using PM\n");
      }
      IdealMM_correct(MM,PM, &rows, &cols,ProbeNames);
    }



    Free(ProbeNames);
    UNPROTECT(1);
  } else if (strncmp(CHAR(STRING_ELT(bg_type,0)),"LESN",4) == 0){
    LESN_param = GetParameter(background_param,"lesnparam");

    if (strcmp(CHAR(STRING_ELT(bg_type,0)),"LESN2") == 0){
      which_lesn =2;
    } else  if (strcmp(CHAR(STRING_ELT(bg_type,0)),"LESN1") == 0){
      which_lesn=1;
    } else if (strcmp(CHAR(STRING_ELT(bg_type,0)),"LESN0") == 0){
      which_lesn=0;
    }

   
    PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
    rows = INTEGER(dim1)[0];
    cols = INTEGER(dim1)[1]; 
    PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
    MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));
    baseline = NUMERIC_POINTER(LESN_param)[0];
    theta = NUMERIC_POINTER(LESN_param)[1];
    param = GetParameter(background_param,"type");
    if ((strcmp(CHAR(STRING_ELT(param,0)),"pmonly") == 0) || (strcmp(CHAR(STRING_ELT(param,0)),"separate") == 0)){ 
      if (verbosity_level > 0){
	Rprintf("Background correcting PM\n");
      }
      LESN_correct(PM, rows, cols, 2, baseline,theta);
    }
    if ((strcmp(CHAR(STRING_ELT(param,0)),"mmonly") == 0) || (strcmp(CHAR(STRING_ELT(param,0)),"separate") == 0)){  
      if (verbosity_level > 0){
	Rprintf("Background correcting MM\n");
      }
      LESN_correct(MM, rows, cols, 2, baseline,theta);
    }
    if (strcmp(CHAR(STRING_ELT(param,0)),"together") == 0){ 
      if (verbosity_level > 0){
	Rprintf("Background correcting PM and MM together\n");
      }
      PROTECT(allPMMM = allocMatrix(REALSXP,2*rows,cols));
      allrows = 2*rows;
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  NUMERIC_POINTER(allPMMM)[j*2*rows + i] = NUMERIC_POINTER(PMmat)[j*rows + i];
	}
      }
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  NUMERIC_POINTER(allPMMM)[j*2*rows + i + rows] = NUMERIC_POINTER(MMmat)[j*rows + i];
	}
      }
      LESN_correct(NUMERIC_POINTER(allPMMM), allrows, cols, 2, baseline,theta);
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  NUMERIC_POINTER(PMmat)[j*rows + i] = NUMERIC_POINTER(allPMMM)[j*2*rows + i];
	}
      }
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  NUMERIC_POINTER(MMmat)[j*rows + i] = NUMERIC_POINTER(allPMMM)[j*2*rows + i + rows];
	}
      }
      UNPROTECT(1);
    }
    UNPROTECT(1);
  }
  
  UNPROTECT(1);
  return PMmat;
}


/********************************************************************************************
 **
 ** void pp_normalize(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_type, SEXP verbosity)
 **
 ** SEXP PMmat - matrix of Perfect-match values
 ** SEXP MMmat - matrix of Mismatch values
 ** SEXP ProbeNamesVec - vector containing names of probeset for each probe
 ** SEXP N_probes - number of PM/MM probes on an array
 ** SEXP norm_type  - an integer indicating the normalization method to be used.
 ** SEXP verbosity - how verbose to be on the output.
 **
 ** a function to manipulate the data so that the R objects are converted
 ** into C data structures. The Normalization is then applied to the C data structures. 
 ** The summary statistic is applied only to PM data.
 **
 ** this function assumes any sort of background correction has been applied and only the
 ** PM probes need be normalized.
 **
 **
 ** norm_type  1 - Quantile Normalization
 ** norm_type  2 - Probeset Quantile Normalization
 ** norm_type  3 - Scaling normalization
 ** 
 *******************************************************************************************/

SEXP pp_normalize(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_type, SEXP norm_parameters, SEXP verbosity){


  SEXP param;
  
  int rows, cols;
  /* double *outexpr; */
  double *PM,*MM;
  double *allPMMM;
  int allrows;
  /*  char **outnames; */
  const char **ProbeNames; 
  int i,j;
  int nprobesets;
  
  double trim;
  int baseline;
  int logscale;

  int usemedian;
  int uselog2;
  int weightscheme;

  int generate_weights = 0;

  int verbosity_level;



  double *weights=0;

  SEXP dim1;
  /* SEXP outvec,outnamesvec;
     SEXP dimnames,names;*/
  
  SEXP remove_extreme, n_remove, R_weights;





  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1]; 

  PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
  MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));
  
  nprobesets=INTEGER(N_probes)[0];
  
  /*  printf("%i\n",nprobesets); */
  /* printf("%d ",INTEGER(norm_flag)[0]); */
  /* normalize PM using quantile normalization */
  /* printf("Normalizing\n"); */
  /* Rprintf("Normalizing\n"); */

  verbosity_level = asInteger(verbosity);



  if (strcmp(CHAR(STRING_ELT(norm_type,0)),"quantile") == 0){
    param = GetParameter(norm_parameters,"type");
    if ((strcmp(CHAR(STRING_ELT(param,0)),"pmonly") == 0) || (strcmp(CHAR(STRING_ELT(param,0)),"separate") == 0)){
      if (verbosity_level > 0){
	Rprintf("Normalizing PM\n");
      }
      qnorm_c(PM,&rows,&cols);
    } 
    if ((strcmp(CHAR(STRING_ELT(param,0)),"mmonly") == 0)|| (strcmp(CHAR(STRING_ELT(param,0)),"separate") == 0)){
       if (verbosity_level > 0){
	 Rprintf("Normalizing MM\n");
       }
      qnorm_c(MM,&rows,&cols);
    }
    if (strcmp(CHAR(STRING_ELT(param,0)),"together") == 0){  
      if (verbosity_level > 0){
	Rprintf("Normalizing PM and MM together\n");
      }
      allPMMM = (double *)Calloc(2*rows*cols,double);
      allrows = 2*rows;
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  allPMMM[j*2*rows + i] = PM[j*rows + i];
	}
      }
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  allPMMM[j*2*rows + i + rows] = MM[j*rows + i];
	}
      }
      qnorm_c(allPMMM,&allrows,&cols);
      
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  PM[j*rows + i] = allPMMM[j*2*rows + i];
	}
      }
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  MM[j*rows + i] = allPMMM[j*2*rows + i + rows];
	}
      }
      Free(allPMMM);

    }



  } else if (strcmp(CHAR(STRING_ELT(norm_type,0)),"quantile.probeset") == 0){
    ProbeNames = (const char **)Calloc(rows, const char *);
    for (i =0; i < rows; i++)
      ProbeNames[i] = CHAR(STRING_ELT(ProbeNamesVec,i));
    
    param = GetParameter(norm_parameters,"use.median");
    usemedian = asInteger(param);
    param = GetParameter(norm_parameters,"use.log2");
    uselog2 = asInteger(param);

    param = GetParameter(norm_parameters,"type");
    if ((strcmp(CHAR(STRING_ELT(param,0)),"pmonly") == 0) || (strcmp(CHAR(STRING_ELT(param,0)),"separate") == 0)){   
      if (verbosity_level > 0){
	Rprintf("Normalizing PM\n");
      }
      qnorm_probeset_c(PM, rows, cols, nprobesets, ProbeNames, usemedian, uselog2);
    } 
    if ((strcmp(CHAR(STRING_ELT(param,0)),"mmonly") == 0) || (strcmp(CHAR(STRING_ELT(param,0)),"separate") == 0)){ 
      if (verbosity_level > 0){
	Rprintf("Normalizing MM\n");
      }
      qnorm_probeset_c(MM, rows, cols, nprobesets, ProbeNames, usemedian, uselog2);
    }
    if (strcmp(CHAR(STRING_ELT(param,0)),"together") == 0){ 
      if (verbosity_level > 0){
	Rprintf("Normalizing PM and MM together\n");
      }
      allPMMM = (double *)Calloc(2*rows*cols,double);
      allrows = 2*rows;
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  allPMMM[j*2*rows + i] = PM[j*rows + i];
	}
      }
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  allPMMM[j*2*rows + i + rows] = MM[j*rows + i];
	}
      }
      ProbeNames = (const char **)Realloc(ProbeNames,2*rows,const char *);
      for (i =0; i < rows; i++)
	ProbeNames[rows + i] = CHAR(STRING_ELT(ProbeNamesVec,i));
      qnorm_probeset_c(allPMMM, allrows, cols, 2*nprobesets, ProbeNames, usemedian, uselog2);
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  PM[j*rows + i] = allPMMM[j*2*rows + i];
	}
      }
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  MM[j*rows + i] = allPMMM[j*2*rows + i + rows];
	}
      }
      Free(allPMMM);
    }




    Free(ProbeNames);

  } else if (strcmp(CHAR(STRING_ELT(norm_type,0)),"scaling") == 0){
    param = GetParameter(norm_parameters,"scaling.trim");
    trim = asReal(param);
    param = GetParameter(norm_parameters, "scaling.baseline");
    baseline = asInteger(param); 
    
    param = GetParameter(norm_parameters, "log.scalefactors");
    logscale = asInteger(param); 

    param = GetParameter(norm_parameters,"type");
    if ((strcmp(CHAR(STRING_ELT(param,0)),"pmonly") == 0) || (strcmp(CHAR(STRING_ELT(param,0)),"separate") == 0)){
      if (verbosity_level > 0){
	Rprintf("Normalizing PM\n");
      }
      scaling_norm(PM, rows, cols,trim, baseline,logscale);
    }
    if ((strcmp(CHAR(STRING_ELT(param,0)),"mmonly") == 0) || (strcmp(CHAR(STRING_ELT(param,0)),"separate") == 0)){
      if (verbosity_level > 0){
	Rprintf("Normalizing MM\n");
      }
      scaling_norm(MM, rows, cols,trim, baseline,logscale);
    }
    if (strcmp(CHAR(STRING_ELT(param,0)),"together") == 0){
      if (verbosity_level > 0){
	Rprintf("Normalizing PM and MM together\n");
      }
      allPMMM = (double *)Calloc(2*rows*cols,double);
      allrows = 2*rows;
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  allPMMM[j*2*rows + i] = PM[j*rows + i];
	}
      }
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  allPMMM[j*2*rows + i + rows] = MM[j*rows + i];
	}
      }
      scaling_norm(allPMMM, allrows, cols,trim, baseline,logscale);
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  PM[j*rows + i] = allPMMM[j*2*rows + i];
	}
      }
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  MM[j*rows + i] = allPMMM[j*2*rows + i + rows];
	}
      }
      Free(allPMMM);
    }
  } else if (strcmp(CHAR(STRING_ELT(norm_type,0)),"quantile.robust") == 0){      

    param = GetParameter(norm_parameters,"use.median");
    usemedian = asInteger(param);
    param = GetParameter(norm_parameters,"use.log2");
    uselog2 = asInteger(param);

    remove_extreme = GetParameter(norm_parameters,"remove.extreme");
    n_remove = GetParameter(norm_parameters,"n.remove");


    param = GetParameter(norm_parameters,"weights");

    weightscheme = 0;
    if (isString(param)){
      if (strcmp(CHAR(STRING_ELT(param,0)),"huber") == 0){
	weightscheme = 1;
      } else {
	weightscheme = 0;
      }
    } else if (isNumeric(param)){
      if (length(param) != cols){
	error("Problem with weights for quantile.robust");
      }
      weights = REAL(param);
    } else if (isNull(param)){
      /* Need to generate weights or more correctly exclude some arrays */
      generate_weights = 1;
   
    }

    param = GetParameter(norm_parameters,"type");
    if ((strcmp(CHAR(STRING_ELT(param,0)),"pmonly") == 0) || (strcmp(CHAR(STRING_ELT(param,0)),"separate") == 0)){
     if (verbosity_level > 0){
       Rprintf("Normalizing PM\n");
     }
      
      if (generate_weights){
	R_weights = R_qnorm_robust_weights(PMmat, remove_extreme, n_remove);
	weights = REAL(R_weights);
      }
      qnorm_robust_c(PM,weights, &rows, &cols, &usemedian, &uselog2, &weightscheme);
    }
    if ((strcmp(CHAR(STRING_ELT(param,0)),"mmonly") == 0) || (strcmp(CHAR(STRING_ELT(param,0)),"separate") == 0)){ 
      if (verbosity_level > 0){
	Rprintf("Normalizing MM\n");
      }
      if (generate_weights){
	R_weights = R_qnorm_robust_weights(MMmat, remove_extreme, n_remove);
	weights = REAL(R_weights);
      }
      qnorm_robust_c(MM,weights, &rows, &cols, &usemedian, &uselog2, &weightscheme);
    } 
    if (strcmp(CHAR(STRING_ELT(param,0)),"together") == 0){
      if (verbosity_level > 0){
	Rprintf("Normalizing PM and MM together\n");    
      }
      if (generate_weights){
	R_weights = R_qnorm_robust_weights(PMmat, remove_extreme, n_remove);
	weights = REAL(R_weights);
      }
      allPMMM = (double *)Calloc(2*rows*cols,double);
      allrows = 2*rows;
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  allPMMM[j*2*rows + i] = PM[j*rows + i];
	}
      }
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  allPMMM[j*2*rows + i + rows] = MM[j*rows + i];
	}
      } 
      weights = Calloc(cols,double);
      qnorm_robust_c(allPMMM,weights, &allrows, &cols, &usemedian, &uselog2, &weightscheme);
      Free(weights);
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  PM[j*rows + i] = allPMMM[j*2*rows + i];
	}
      }
      for (i=0; i < rows; i++){
	for (j=0; j < cols; j++){
	  MM[j*rows + i] = allPMMM[j*2*rows + i + rows];
	}
      }
      Free(allPMMM);
    }
  } 
  UNPROTECT(1);
  return PMmat;
}



SEXP pp_bothstages(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_flag, SEXP bg_flag, SEXP bg_type,SEXP norm_type, SEXP background_parameters,SEXP norm_parameters, SEXP verbosity){
  
  SEXP dim1,PMcopy;
  int rows,cols;
  
  /*Create a copy matrix to work on. Allows us to modify data in background and normalization steps without affecting original data */
  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  /*PROTECT(PMcopy = allocMatrix(REALSXP,rows,cols));
    copyMatrix(PMcopy,PMmat,0); */
  UNPROTECT(1);
  PMcopy = PMmat;

  
  /* If Background correction do it */
  if (INTEGER(bg_flag)[0]){
    PMcopy = pp_background(PMcopy, MMmat, ProbeNamesVec,N_probes,bg_type, background_parameters,verbosity); /* densfunc, rho, LESN_params); */
  }
  
  /* If Normalization do it */
  if (INTEGER(norm_flag)[0]){
    PMcopy = pp_normalize(PMcopy, MMmat, ProbeNamesVec,N_probes,norm_type, norm_parameters,verbosity);
  }
  
  return PMcopy;
}
