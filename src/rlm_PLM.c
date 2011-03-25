/*********************************************************************
 **
 ** file: rlm_PLM.c
 **
 ** Aim: fit robust linear models for the PLMset object.
 **
 ** Copyright (C) 2003-2008 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bmb@bmbolstad.com>
 ** 
 ** created on: Jan 17, 2003
 **
 ** Last modified: Nov 1, 2008
 **
 ** the aim here will be to fit specified robust linear models
 ** to affy data. 
 **
 ** The code in this particular file takes the objects passed from R
 ** sends what is required to the preprocessing steps (Background and Normalization)
 ** and then takes that data, forms "C" data structures and calls routines to do the 
 ** actual robust model fitting
 **
 ** Modification history
 **
 ** Jan 17, 2003 - Initial version.
 ** Jan 18, 2003 - Better setup rlm procedure to take covariates and
 **                specify different models.
 ** Jan 19, 2003 - continued implementation
 ** Jan 20, 2003 - more implementation, clean up parameter passing methodology.
 ** Jan 24, 2003 - expand the range of models that can be fit by actually
 **                making use of the chipcovariates parameter.
 ** Jan 27, 2003 - Standard error calculation
 ** Jan 28, 2003 - Ability to select different types of standard error estimation
 ** Feb 1, 2003 - remove the row naming aspect on weights and probes, this will be handled by R 
 **               more documentation
 ** Feb 6, 2003 - Change printf("Fitting models ....") to an Rprintf
 ** Feb 15, 2003 - Add a mechanism for returning intercept parameters.
 ** Feb 17, 2003 - add in a free(outnames); free(ProbeNames) to rlmPLMset;
 ** Feb 24, 2003 - get rid of unused, but declared variables.
 ** Mar 21, 2003 - modify background for LESN methods
 ** Jun 4,  2003 - Add mechanism for different psi functions
 ** Jul 26, 2003 - add normalization options parameter
 **                add background options parameter
 ** Sep 02, 2003 - we now store residuals
 ** Sep 04, 2003 - a parameter which specifies what should be outputted is now
 **                passed. This item is an R list similar to normalization
 **                and background parameter lists.
 **                Considerable clean up of how parameters are passed to 
 **                do_PLMrlm
 ** Sep 06, 2003 - Make Storage allocation routine separate.
 ** Sep 07, 2003 - output varcov
 ** Sep 12, 2003 - remove psi, psi_k etc from arguments of functions
 **                they are now in model_params
 ** Sept 14, 2003 - can intialize M estimatation starting with a fully iterated
 **                 Huber regression
 ** Apr 5, 2004   - All malloc/free are now Calloc/Free
 ** May 3, 2004   - Fixed a subtle and small memory leak.
 ** May 27, 2004  - add a way to detect that default model is being fitted
 ** July 10, 2004 - Start integrating new structure
 ** Mar 1, 2006 - change all comments to ansi style
 ** Jun 6, 2006 - fix problem with space allocation when only one probe in probeset and trying to estimate
 **               probe-effects. (Example is Soybean chips).
 ** Oct 11, 2006 - make verbosity argument get passed to summarization function
 ** Nov 30, 2007 - remove code that was commented out for being defunct. This will help maintainability.
 ** Nov 1, 2008 - comment out some defunct code (should be removed at a later date)
 ** Jan 6, 2009 - change SET_VECTOR_ELT to SET_STRING_ELT where relevant.
 **
 **
 *********************************************************************/

#include "preprocess.h"
#include "do_PLMrlm.h"
#include "common_types.h"

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>





/*********************************************************************
 **
 ** static void rlmPLM_alloc_space(PLMRoutput *Routput, PLMoutput *output,
 **                                outputsettings *store,Datagroup *data, 
 **                                PLMmodelparam *model)
 **
 ** 
 ** This function allocates all the space that is needed for storing
 ** user desired output from the PLM
 ** 
 **
 ********************************************************************/

/* static void rlmPLM_alloc_space(PLMRoutput *Routput, PLMoutput *output,outputsettings *store,Datagroup *data, PLMmodelparam *model){
 **  SEXP tmp;
 **
 **  int i;
  ** **  
  ** Routput->nprotected = 0;
 **
 **  
 **  output->outnames = (char **)Calloc(data->nprobesets,char *);
 **
 **  if (store->weights){
 **    PROTECT(Routput->weights = allocMatrix(REALSXP, data->rows, data->cols));
 **  } else {
 **    PROTECT(Routput->weights = allocMatrix(REALSXP, 0, 0));
 **  }
 **  Routput->nprotected++;
 **  output->out_weights = NUMERIC_POINTER(Routput->weights);
 **
 **
  ** PROTECT(Routput->probe_coef = allocMatrix(REALSXP,data->rows,1));
  ** Routput->nprotected++;
 **  output->out_probeparams = NUMERIC_POINTER(Routput->probe_coef);
 **
 **
 **  PROTECT(Routput->chip_coef = allocMatrix(REALSXP, data->nprobesets, model->nchipparams));
 **  Routput->nprotected++;
 **  output->out_chipparams = NUMERIC_POINTER(Routput->chip_coef);
 **  
 **  PROTECT(Routput->const_coef = allocMatrix(REALSXP, data->nprobesets, 1));
  ** Routput->nprotected++;
 **  output->out_constparams = NUMERIC_POINTER(Routput->const_coef);
 **
  ** PROTECT(Routput->chip_SE = allocMatrix(REALSXP, data->nprobesets, model->nchipparams));
 **  Routput->nprotected++;
 **  output->out_chip_SE = NUMERIC_POINTER(Routput->chip_SE);
 **
 **
 **  PROTECT(Routput->probe_SE = allocMatrix(REALSXP,data->rows,1));
  ** Routput->nprotected++;
 **  output->out_probe_SE = NUMERIC_POINTER(Routput->probe_SE);
 **
 **
  ** PROTECT(Routput->const_SE = allocMatrix(REALSXP, data->nprobesets, 1));
  ** Routput->nprotected++;
 **  output->out_const_SE = NUMERIC_POINTER(Routput->const_SE);
 **
 **
  ** if (store->residuals){
  **   PROTECT(Routput->residuals = allocMatrix(REALSXP, data->rows, data->cols));
 **  } else {
 **    PROTECT(Routput->residuals = allocMatrix(REALSXP, 0, 0));
 **  }
 **  Routput->nprotected++;
 **  output->out_resids = NUMERIC_POINTER(Routput->residuals); 
  ** 
 **
  ** if (store->residSE){
 **    PROTECT(Routput->residSE = allocMatrix(REALSXP,data->nprobesets, 2));
 **  } else {
 **    PROTECT(Routput->residSE = allocMatrix(REALSXP,0,0));
 **  }
 **  Routput->nprotected++;
 **  output->out_residSE = NUMERIC_POINTER(Routput->residSE);
 **
 **  
 **  if (store->varcov == 0){
 **    PROTECT(Routput->varcov = allocVector(VECSXP,0));
 **    output->out_varcov= NULL;
 **  } else if (store->varcov == 1){
 **    PROTECT(Routput->varcov = allocVector(VECSXP,data->nprobesets));
 **    output->out_varcov = Calloc(data->nprobesets,double*);
 **    for (i =0; i < data->nprobesets; i++){
 **      PROTECT(tmp = allocMatrix(REALSXP,model->nchipparams,model->nchipparams));
 **      SET_VECTOR_ELT(Routput->varcov,i,tmp);
 **      UNPROTECT(1);
 **      output->out_varcov[i] = NUMERIC_POINTER(VECTOR_ELT(Routput->varcov,i));
 **    }
 **  }
 **  Routput->nprotected++;
 **  
 **  
 **  
 **  
 **
 **
 **  }  
*/



static int  checkDefaultModel(const PLM_model_parameters *model){

  int i;
  int howmany=0;

  for (i=0; i < 5; i++){
    howmany+=model->which_parameter_types[i];
  }
  
  if (howmany > 2 || howmany < 2){
    return 0;
  }

  if (!(model->which_parameter_types[2] & model->which_parameter_types[4])){
    return 0;
  }

  if (model->constraints[4] == 1){
    return 0;
  }

  if (model->strata[4] !=0){
    return 0;
  }

  if (model->mmorpm_covariate != 0){
    return 0;
  }


  if (model->response_variable == 0){
    return 0;
  }
  
  if (model->psi_code > 3){
    return 0;
  }


  return 1;
  
}













static void rlm_PLM_alloc_space(PLMRoutput *Routput, PLM_output *output,PLM_outputsettings *store,PLM_Datagroup *data, PLM_model_parameters *model){
  SEXP tmp;

  int i,j;
  int n_const_col;
  int probe_multiplier =0;
  const char *first;
  int *n_probes_probeset;
  int new_nprobes=0;
  int modifier=0;
  int p;


  Routput->nprotected = 0;

  
  output->outnames = (char **)Calloc(data->n_probesets,char *);
  

  
  /**** Everything below is always stored ****/
  /**** where sample and chip-level covariate estimates are stored ****/
  if (model->n_chiplevelcovariates == 0){
    if (model->which_parameter_types[2] != 0){
      if (model->constraints[2]==0){
	PROTECT(Routput->chip_coef = allocMatrix(REALSXP, data->n_probesets, data->n_arrays));
	Routput->nprotected++;
	output->out_chipparams = NUMERIC_POINTER(Routput->chip_coef);
	PROTECT(Routput->chip_SE = allocMatrix(REALSXP, data->n_probesets, model->n_arrays));
	Routput->nprotected++;
	output->out_chip_SE = NUMERIC_POINTER(Routput->chip_SE);
      } else {
	PROTECT(Routput->chip_coef = allocMatrix(REALSXP, data->n_probesets, data->n_arrays-1));
	Routput->nprotected++;
	output->out_chipparams = NUMERIC_POINTER(Routput->chip_coef);
	PROTECT(Routput->chip_SE = allocMatrix(REALSXP, data->n_probesets, model->n_arrays-1));
	Routput->nprotected++;
	output->out_chip_SE = NUMERIC_POINTER(Routput->chip_SE);
      }
    } else {
      PROTECT(Routput->chip_coef = allocMatrix(REALSXP, 0, 0));
      Routput->nprotected++;
      output->out_chipparams = NUMERIC_POINTER(Routput->chip_coef);
      PROTECT(Routput->chip_SE = allocMatrix(REALSXP, 0, 0));
      Routput->nprotected++;
      output->out_chip_SE = NUMERIC_POINTER(Routput->chip_SE);

    }
  } else {
    PROTECT(Routput->chip_coef = allocMatrix(REALSXP, data->n_probesets, model->n_chiplevelcovariates));
    Routput->nprotected++;
    output->out_chipparams = NUMERIC_POINTER(Routput->chip_coef);
    PROTECT(Routput->chip_SE = allocMatrix(REALSXP, data->n_probesets, model->n_chiplevelcovariates));
    Routput->nprotected++;
    output->out_chip_SE = NUMERIC_POINTER(Routput->chip_SE);
  }

  /*** this is where intercepts and probe.type parameters are stored ***/

  n_const_col=0;
  if (model->which_parameter_types[0]){
    n_const_col++;
  }

  /* Check for PM or MM type covariate */
  if (model->mmorpm_covariate != 0){
    n_const_col++;
  }

  if (model->which_parameter_types[3]){
    if (model->constraints[3] == 0){
      if (model->strata[3]==0){
	n_const_col +=2;
      } else if (model->strata[3]==1){
	n_const_col +=2*model->n_arrays;
      } else if (model->strata[3]==2){
	n_const_col +=2*(model->max_probe_type_treatment_factor+1);
      }
    } else {
      if (model->strata[3]==0){	
	n_const_col +=1;
      } else if (model->strata[3]==1){
	n_const_col +=model->n_arrays;
      } else if (model->strata[3]==2){
	n_const_col +=(model->max_probe_type_treatment_factor+1);
      }
    }   
  }


  if (n_const_col != 0){
    PROTECT(Routput->const_coef = allocMatrix(REALSXP, data->n_probesets, n_const_col));
    Routput->nprotected++;
    output->out_constparams = NUMERIC_POINTER(Routput->const_coef);
    
  
    PROTECT(Routput->const_SE = allocMatrix(REALSXP, data->n_probesets, n_const_col));
    Routput->nprotected++;
    output->out_const_SE = NUMERIC_POINTER(Routput->const_SE);
  } else {
    PROTECT(Routput->const_coef = allocMatrix(REALSXP, 0, 0));
    Routput->nprotected++;
    output->out_constparams = NUMERIC_POINTER(Routput->const_coef);
    
  
    PROTECT(Routput->const_SE = allocMatrix(REALSXP, 0, 0));
    Routput->nprotected++;
    output->out_const_SE = NUMERIC_POINTER(Routput->const_SE);
  }





  i=0;
  n_probes_probeset= Calloc(data->n_probesets, int);
  new_nprobes=1;
  first = data->ProbeNames[0];
  for (j = 1; j < data->n_probes; j++){
    if ((strcmp(first,data->ProbeNames[j]) != 0) | (j == (data->n_probes -1))){
      if (j == (data->n_probes -1)){
	new_nprobes++;
      }
      n_probes_probeset[i] = new_nprobes;
      i++;
      first = data->ProbeNames[j];
      new_nprobes = 0;
    }
    new_nprobes++;
  }
  
  if (model->which_parameter_types[4]){
    PROTECT(Routput->probe_coef = allocVector(VECSXP,data->n_probesets));
    Routput->nprotected++;
    PROTECT(Routput->probe_SE = allocVector(VECSXP,data->n_probesets));
    Routput->nprotected++;
    
  

    if (model->which_parameter_types[4]){
      /* probe effect parameter figure figure out which strata and constraints */ 
      if (model->strata[4] == 0){
	/* Overall probe effect */
	probe_multiplier=1;
      } else if (model->strata[4] == 2){
	/* Probe effect within the levels of a treatment/covariate factor */
	probe_multiplier=((model->max_probe_treatment_factor+1));
      } else if (model->strata[4] == 3){
	/* Probe effect within probe type */
	probe_multiplier=2;
      } else if (model->strata[4] == 4){
	/* Probe effect within probe type within the levels of a treatment/covariate factor*/
       probe_multiplier=(2*(model->max_probe_treatment_factor+1));
      }
    }
  
    
  
    if (checkDefaultModel(model)){
      modifier =0;
    } else {    
      if (model->which_parameter_types[4]){
	if (model->constraints[4] !=0){
	  modifier = -1;
	} else {
	  modifier =0;
	}
      }
    }
    for (i=0; i < data->n_probesets; i++){
      if (n_probes_probeset[i]+modifier >= 0){
	PROTECT(tmp = allocMatrix(REALSXP,  (n_probes_probeset[i]+modifier),probe_multiplier));
	SET_VECTOR_ELT(Routput->probe_coef,i,tmp);
	UNPROTECT(1);
	PROTECT(tmp = allocMatrix(REALSXP, (n_probes_probeset[i]+modifier),probe_multiplier));
	SET_VECTOR_ELT(Routput->probe_SE,i,tmp);
	UNPROTECT(1);
      } else {
	PROTECT(tmp = allocMatrix(REALSXP,  0,probe_multiplier));
	SET_VECTOR_ELT(Routput->probe_coef,i,tmp);
	UNPROTECT(1);
	PROTECT(tmp = allocMatrix(REALSXP,  0,probe_multiplier));
	SET_VECTOR_ELT(Routput->probe_SE,i,tmp);
	UNPROTECT(1);


      }
    }


    output->out_probeparams = (double **)Calloc(data->n_probesets,double *);
    output->out_probe_SE = (double **)Calloc(data->n_probesets,double *);
    
  
    for (i=0; i < data->n_probesets; i++){
      output->out_probeparams[i] = NUMERIC_POINTER(VECTOR_ELT(Routput->probe_coef,i));
      output->out_probe_SE[i] = NUMERIC_POINTER(VECTOR_ELT(Routput->probe_SE,i));
    }
  } else {
    PROTECT(Routput->probe_coef = allocVector(VECSXP,0));
    Routput->nprotected++;
    PROTECT(Routput->probe_SE = allocVector(VECSXP,0));
    Routput->nprotected++;
  }



  /**** The Optional components ****/

  /**** Weights ****/
  if (store->weights){
    PROTECT(Routput->weights = allocVector(VECSXP,2));
    if (model->response_variable >= 0){
      PROTECT(tmp = allocMatrix(REALSXP, data->n_probes, data->n_arrays));
      SET_VECTOR_ELT(Routput->weights,0,tmp);
      UNPROTECT(1);
    } else {
      PROTECT(tmp=allocMatrix(REALSXP, 0, 0));
      SET_VECTOR_ELT(Routput->weights,0,tmp);
      UNPROTECT(1);
    }
    if (model->response_variable <=0){
      PROTECT(tmp = allocMatrix(REALSXP, data->n_probes, data->n_arrays));
      SET_VECTOR_ELT(Routput->weights,1,tmp);
      UNPROTECT(1);
    } else {
      PROTECT(tmp=allocMatrix(REALSXP, 0, 0));
      SET_VECTOR_ELT(Routput->weights,1,tmp);
      UNPROTECT(1);
    }
  } else {
    PROTECT(Routput->weights = allocVector(VECSXP,2));   
    PROTECT(tmp=allocMatrix(REALSXP, 0, 0));
    SET_VECTOR_ELT(Routput->weights,0,tmp);
    UNPROTECT(1);
    PROTECT(tmp=allocMatrix(REALSXP, 0, 0));
    SET_VECTOR_ELT(Routput->weights,1,tmp);
    UNPROTECT(1);
  }
  
  Routput->nprotected++;

  output->out_weights = (double **)Calloc(2,double *);
  for (i =0; i < 2; i++){
    output->out_weights[i] = NUMERIC_POINTER(VECTOR_ELT(Routput->weights,i));
  }


  /**** Residuals ****/
  if (store->residuals){
    PROTECT(Routput->residuals = allocVector(VECSXP,2));
    if (model->response_variable >= 0){
      PROTECT(tmp = allocMatrix(REALSXP, data->n_probes, data->n_arrays));
      SET_VECTOR_ELT(Routput->residuals,0,tmp);
      UNPROTECT(1);
    } else {
      PROTECT(tmp=allocMatrix(REALSXP, 0, 0));
      SET_VECTOR_ELT(Routput->residuals,0,tmp);
      UNPROTECT(1);
    }
    if (model->response_variable <=0){
      PROTECT(tmp = allocMatrix(REALSXP, data->n_probes, data->n_arrays));
      SET_VECTOR_ELT(Routput->residuals,1,tmp);
      UNPROTECT(1);
    } else {
      PROTECT(tmp=allocMatrix(REALSXP, 0, 0));
      SET_VECTOR_ELT(Routput->residuals,1,tmp);
      UNPROTECT(1);
    }
  } else {
    PROTECT(Routput->residuals = allocVector(VECSXP,2));   
    PROTECT(tmp=allocMatrix(REALSXP, 0, 0));
    SET_VECTOR_ELT(Routput->residuals,0,tmp);
    UNPROTECT(1);
    PROTECT(tmp=allocMatrix(REALSXP, 0, 0));
    SET_VECTOR_ELT(Routput->residuals,1,tmp);
    UNPROTECT(1);
  }
  
  Routput->nprotected++;

  output->out_resids = (double **)Calloc(2,double *);

  for (i =0; i < 2; i++){
    output->out_resids[i] = NUMERIC_POINTER(VECTOR_ELT(Routput->residuals,i));
  }
  
  /*** Residual SE ***/
  if (store->residSE){
    PROTECT(Routput->residSE = allocMatrix(REALSXP,data->n_probesets, 2));
  } else {
    PROTECT(Routput->residSE = allocMatrix(REALSXP,0,0));
  }
  Routput->nprotected++;
  output->out_residSE = NUMERIC_POINTER(Routput->residSE);

  /*** Variance Covariance matrix ***/

  if (store->varcov == 0){
    PROTECT(Routput->varcov = allocVector(VECSXP,0));
    output->out_varcov= NULL;
  } else if (store->varcov == 1){
    /* variance-covariance matrix for chip-level covariates */
    if (model->which_parameter_types[0] | model->which_parameter_types[1] | model->which_parameter_types[2]){
      PROTECT(Routput->varcov = allocVector(VECSXP,data->n_probesets));
      output->out_varcov = Calloc(data->n_probesets,double*);
      if (model->n_chiplevelcovariates == 0){
	if (model->which_parameter_types[2] != 0){
	  for (i =0; i < data->n_probesets; i++){
	    PROTECT(tmp = allocMatrix(REALSXP,data->n_arrays,data->n_arrays));
	    SET_VECTOR_ELT(Routput->varcov,i,tmp);
	    UNPROTECT(1);
	    output->out_varcov[i] = NUMERIC_POINTER(VECTOR_ELT(Routput->varcov,i));
	  }
	} else if (model->which_parameter_types[0]){
	  for (i =0; i < data->n_probesets; i++){
	    PROTECT(tmp = allocMatrix(REALSXP,1,1));
	    SET_VECTOR_ELT(Routput->varcov,i,tmp);
	    UNPROTECT(1);
	    output->out_varcov[i] = NUMERIC_POINTER(VECTOR_ELT(Routput->varcov,i));
	  }
	}
      } else {
	if(model->which_parameter_types[0]){
	  /* need to save the intercept */
	  for (i =0; i < data->n_probesets; i++){
	    PROTECT(tmp = allocMatrix(REALSXP,model->n_chiplevelcovariates+1,model->n_chiplevelcovariates+1));
	    SET_VECTOR_ELT(Routput->varcov,i,tmp);
	    UNPROTECT(1);
	    output->out_varcov[i] = NUMERIC_POINTER(VECTOR_ELT(Routput->varcov,i));
	  }
	} else {
	  for (i =0; i < data->n_probesets; i++){
	    PROTECT(tmp = allocMatrix(REALSXP,model->n_chiplevelcovariates,model->n_chiplevelcovariates));
	    SET_VECTOR_ELT(Routput->varcov,i,tmp);
	    UNPROTECT(1);
	    output->out_varcov[i] = NUMERIC_POINTER(VECTOR_ELT(Routput->varcov,i));
	  }
	}
      }
    } else {
      PROTECT(Routput->varcov = allocVector(VECSXP,0));
      output->out_varcov= NULL;
    }
  } else if (store->varcov == 2){  
    PROTECT(Routput->varcov = allocVector(VECSXP,data->n_probesets));
    output->out_varcov = (double **)Calloc(data->n_probesets,double*);
    for (i =0; i < data->n_probesets; i++){
      p = 0;
      if (model->n_chiplevelcovariates > 0){
	p+=model->n_chiplevelcovariates;
      } else if (model->which_parameter_types[2] != 0){
	if (model->which_parameter_types[0]){
	  p+=model->n_arrays-1;
	} else {
	  p+=model->n_arrays;
	}
      }
      p+= n_const_col+ (n_probes_probeset[i]+modifier)*probe_multiplier;
      PROTECT(tmp = allocMatrix(REALSXP,p,p));
      SET_VECTOR_ELT(Routput->varcov,i,tmp);
      UNPROTECT(1);
      output->out_varcov[i] = NUMERIC_POINTER(VECTOR_ELT(Routput->varcov,i));
    } 
  }
  
  Routput->nprotected++;


  Free(n_probes_probeset);


}


void rlm_PLMset_nameoutput(PLMRoutput *Routput,PLM_output *output,PLM_outputsettings *store, PLM_Datagroup *data,PLM_model_parameters *model,SEXP R_model, SEXP sampleNames, SEXP ProbeNamesVec, SEXP chipcovariates){

  SEXP chip_coef_dimnames;
  SEXP const_coef_dimnames;
  SEXP weights_dimnames;
  SEXP residuals_dimnames;
  SEXP residSE_dimnames;

  SEXP chip_coef_rownames;
  SEXP chip_coef_colnames;

  SEXP const_coef_rownames;  
  SEXP const_coef_colnames;

  SEXP residSE_rownames;
  SEXP residSE_colnames;
  
  SEXP weights_names;
  SEXP weights_rownames;
  SEXP weights_colnames;
  
  SEXP residuals_names;
  SEXP residuals_rownames;
  SEXP residuals_colnames; 

  SEXP probe_coef_names;

  SEXP probe_dimnames;
  SEXP probe_colnames;
  SEXP probe_rownames;

  SEXP varcov_names;
  SEXP varcov_dimnames;
  SEXP varcov_colnames;
  SEXP varcov_rownames;

  SEXP names;
  
  SEXP chip_covariate_names;

  SEXP dim;   /* A place to temporarily store dimensions      */


  int i,j,k;
  
  int n_const_col;
  int curcol;

  char *tmp_str;
  char *tmp_str2;
  
  SEXP probe_trt_levels,probe_type_levels;


  probe_trt_levels = GetParameter(R_model,"probe.trt.levels");
  probe_type_levels = GetParameter(R_model,"probe.type.levels");

  
  PROTECT(names = allocVector(STRSXP,data->n_probesets));
  


 


  for ( i =0; i < data->n_probesets; i++)
    SET_STRING_ELT(names,i,mkChar(output->outnames[i]));


  /** figure out how many columns in the const_coef matrix **/
  n_const_col=0;
  if (model->which_parameter_types[0]){
    n_const_col++;
  }

  /* Check for PM or MM type covariate */
  if (model->mmorpm_covariate != 0){
    n_const_col++;
  }

  if (model->which_parameter_types[3]){
    if (model->constraints[3] == 0){
      if (model->strata[3]==0){
	n_const_col +=2;
      } else if (model->strata[3]==1){
	n_const_col +=2*model->n_arrays;
      } else if (model->strata[3]==2){
	n_const_col +=2*(model->max_probe_type_treatment_factor+1);
      }
    } else {
      if (model->strata[3]==0){	
	n_const_col +=1;
      } else if (model->strata[3]==1){
	n_const_col +=model->n_arrays;
      } else if (model->strata[3]==2){
	n_const_col +=(model->max_probe_type_treatment_factor+1);
      }
    }   
  }

  
  if (n_const_col !=0){
    PROTECT(const_coef_dimnames = allocVector(VECSXP,2));
    PROTECT(const_coef_rownames = allocVector(STRSXP,data->n_probesets));
    PROTECT(const_coef_colnames = allocVector(STRSXP,n_const_col));
    
    n_const_col=0;
    if (model->which_parameter_types[0]){
      SET_STRING_ELT(const_coef_colnames,n_const_col,mkChar("Intercept"));
      n_const_col++;
    }
    if (model->mmorpm_covariate != 0){
      if (model->mmorpm_covariate < 0){
	SET_STRING_ELT(const_coef_colnames,n_const_col,mkChar("PM"));
      } else {
	SET_STRING_ELT(const_coef_colnames,n_const_col,mkChar("MM"));
      }
      n_const_col++;
    }
      
    if (model->which_parameter_types[3]){
      if (model->constraints[3] == 0){
	if (model->strata[3]==0){
	  SET_STRING_ELT(const_coef_colnames,n_const_col,mkChar("probe.type_PM"));
	  SET_STRING_ELT(const_coef_colnames,n_const_col+1,mkChar("probe.type_MM"));
	  n_const_col +=2;
	} else if (model->strata[3]==1){
	  for (i =0; i < model->n_arrays;i++){
	    tmp_str = (char *)Calloc(strlen(CHAR(STRING_ELT(sampleNames,i)))+15,char);
	    tmp_str = strcpy(tmp_str,CHAR(STRING_ELT(sampleNames,i)));
	    tmp_str = strcat(tmp_str,":probe.type_PM");
	    SET_STRING_ELT(const_coef_colnames,n_const_col+i*2,mkChar(tmp_str));
	    Free(tmp_str);
	    tmp_str = (char *)Calloc(strlen(CHAR(STRING_ELT(sampleNames,i)))+15,char);
	    tmp_str = strcpy(tmp_str,CHAR(STRING_ELT(sampleNames,i)));
	    tmp_str = strcat(tmp_str,":probe.type_MM");
	    SET_STRING_ELT(const_coef_colnames,n_const_col+i*2+1,mkChar(tmp_str));
	    Free(tmp_str);
	  }
	  n_const_col +=2*model->n_arrays;
	} else if (model->strata[3]==2){
	  for (i =0; i < model->max_probe_type_treatment_factor+1;i++){
	    tmp_str = (char *)Calloc(strlen(CHAR(STRING_ELT(getAttrib(probe_type_levels,R_NamesSymbol),0)))+ strlen(CHAR(STRING_ELT(VECTOR_ELT(probe_type_levels,0),i)))+16,char);
	    tmp_str = strcpy(tmp_str,CHAR(STRING_ELT(getAttrib(probe_type_levels,R_NamesSymbol),0)));
	    tmp_str = strcat(tmp_str,"_");
	    tmp_str = strcat(tmp_str,CHAR(STRING_ELT(VECTOR_ELT(probe_type_levels,0),i)));
	    tmp_str = strcat(tmp_str,":probe.type_PM");
	    SET_STRING_ELT(const_coef_colnames,n_const_col+2*i,mkChar(tmp_str));
	    Free(tmp_str);
	    tmp_str = (char *)Calloc(strlen(CHAR(STRING_ELT(getAttrib(probe_type_levels,R_NamesSymbol),0)))+ strlen(CHAR(STRING_ELT(VECTOR_ELT(probe_type_levels,0),i)))+16,char);
	    tmp_str = strcpy(tmp_str,CHAR(STRING_ELT(getAttrib(probe_type_levels,R_NamesSymbol),0)));
	    tmp_str = strcat(tmp_str,"_");
	    tmp_str = strcat(tmp_str,CHAR(STRING_ELT(VECTOR_ELT(probe_type_levels,0),i)));
	    tmp_str = strcat(tmp_str,":probe.type_MM");
	    SET_STRING_ELT(const_coef_colnames,n_const_col+2*i+1,mkChar(tmp_str));
	    Free(tmp_str);
	  }
	  n_const_col +=2*(model->max_probe_type_treatment_factor+1);
	}
      } else {
	if (model->strata[3]==0){
	  if (model->constraints[3] > 0){
	    SET_STRING_ELT(const_coef_colnames,n_const_col,mkChar("probe.type_MM"));
	  } else {
	    SET_STRING_ELT(const_coef_colnames,n_const_col,mkChar("probe.type_PM"));
	  }
	  n_const_col +=1;
	} else if (model->strata[3]==1){
	  if (model->constraints[3] > 0){
	    for (i =0; i < model->n_arrays;i++){
	      tmp_str = (char *)Calloc(strlen(CHAR(STRING_ELT(sampleNames,i)))+15,char);
	      tmp_str = strcpy(tmp_str,CHAR(STRING_ELT(sampleNames,i)));
	      tmp_str = strcat(tmp_str,":probe.type_MM");
	      SET_STRING_ELT(const_coef_colnames,n_const_col+i,mkChar(tmp_str));
	      Free(tmp_str);
	    }
	  } else {
	    for (i =0; i < model->n_arrays;i++){
	      tmp_str = (char *)Calloc(strlen(CHAR(STRING_ELT(sampleNames,i)))+15,char);
	      tmp_str = strcpy(tmp_str,CHAR(STRING_ELT(sampleNames,i)));
	      tmp_str = strcat(tmp_str,":probe.type_PM");
	      SET_STRING_ELT(const_coef_colnames,n_const_col+i,mkChar(tmp_str));
	      Free(tmp_str);
	    }
	  }
	  n_const_col +=model->n_arrays;
	} else if (model->strata[3]==2){ 
	  if (model->constraints[3] > 0){
	    for (i =0; i < model->max_probe_type_treatment_factor+1;i++){
	      tmp_str = (char *)Calloc(strlen(CHAR(STRING_ELT(getAttrib(probe_type_levels,R_NamesSymbol),0)))+ strlen(CHAR(STRING_ELT(VECTOR_ELT(probe_type_levels,0),i)))+16,char);
	      tmp_str = strcpy(tmp_str,CHAR(STRING_ELT(getAttrib(probe_type_levels,R_NamesSymbol),0)));
	      tmp_str = strcat(tmp_str,"_");
	      tmp_str = strcat(tmp_str,CHAR(STRING_ELT(VECTOR_ELT(probe_type_levels,0),i)));
	      tmp_str = strcat(tmp_str,":probe.type_MM");
	      SET_STRING_ELT(const_coef_colnames,n_const_col+i,mkChar(tmp_str));
	      Free(tmp_str);
	    }
	  } else {
	    for (i =0; i < model->max_probe_type_treatment_factor+1;i++){
	      tmp_str = (char *)Calloc(strlen(CHAR(STRING_ELT(getAttrib(probe_type_levels,R_NamesSymbol),0)))+ strlen(CHAR(STRING_ELT(VECTOR_ELT(probe_type_levels,0),i)))+16,char);
	      tmp_str = strcpy(tmp_str,CHAR(STRING_ELT(getAttrib(probe_type_levels,R_NamesSymbol),0)));
	      tmp_str = strcat(tmp_str,"_");
	      tmp_str = strcat(tmp_str,CHAR(STRING_ELT(VECTOR_ELT(probe_type_levels,0),i)));
	      tmp_str = strcat(tmp_str,":probe.type_PM");
	      SET_STRING_ELT(const_coef_colnames,n_const_col+i,mkChar(tmp_str));
	      Free(tmp_str);
	    }
	  }
	  n_const_col +=(model->max_probe_type_treatment_factor+1);
	}
      }   
    }

    copyVector(const_coef_rownames,names);  
    SET_VECTOR_ELT(const_coef_dimnames,0,const_coef_rownames);
    SET_VECTOR_ELT(const_coef_dimnames,1,const_coef_colnames);
  
    setAttrib(Routput->const_coef,R_DimNamesSymbol, const_coef_dimnames);
    setAttrib(Routput->const_SE,R_DimNamesSymbol, const_coef_dimnames);
    UNPROTECT(3);
  }








  /* Now the probe coef matrix */
  
  if (model->which_parameter_types[4]){
    PROTECT(probe_coef_names= allocVector(STRSXP,data->n_probesets));
    copyVector(probe_coef_names,names);
    setAttrib(Routput->probe_coef,R_NamesSymbol,probe_coef_names);
    setAttrib(Routput->probe_SE,R_NamesSymbol,probe_coef_names);
    UNPROTECT(1);
    tmp_str2 = (char *)Calloc(10,char);
    for (i =0; i < data->n_probesets;i++){
      PROTECT(dim = getAttrib(VECTOR_ELT(Routput->probe_coef,i),R_DimSymbol));	
      PROTECT(probe_dimnames = allocVector(VECSXP,2));
      PROTECT(probe_rownames = allocVector(STRSXP,INTEGER(dim)[0]));
      for (j =0; j < INTEGER(dim)[0]; j++){
	tmp_str = (char *)Calloc(15,char);
	tmp_str = strcpy(tmp_str,"probe_");
	if (model->constraints[4] == 1){
	  snprintf(tmp_str2,9,"%d",j+2);
	} else {
	  snprintf(tmp_str2,9,"%d",j+1);
	}
	tmp_str = strcat(tmp_str,tmp_str2);
	SET_STRING_ELT(probe_rownames,j,mkChar(tmp_str));
	Free(tmp_str);
      }

      PROTECT(probe_colnames = allocVector(STRSXP,INTEGER(dim)[1]));
      if (INTEGER(dim)[1]==1){
	SET_STRING_ELT(probe_colnames,0,mkChar("Overall"));
      } else if (model->strata[4] == 3){
	SET_STRING_ELT(probe_colnames,0,mkChar("probe.type_PM:"));
	SET_STRING_ELT(probe_colnames,1,mkChar("probe.type_MM:"));
      } else if (model->strata[4] == 2){
	for (j =0; j < model->max_probe_treatment_factor+1;j++){
	  tmp_str = (char *)Calloc(strlen(CHAR(STRING_ELT(getAttrib(probe_trt_levels,R_NamesSymbol),0)))+ strlen(CHAR(STRING_ELT(VECTOR_ELT(probe_trt_levels,0),j)))+16,char);
	  tmp_str = strcpy(tmp_str,CHAR(STRING_ELT(getAttrib(probe_trt_levels,R_NamesSymbol),0)));
	  tmp_str = strcat(tmp_str,"_");
	  tmp_str = strcat(tmp_str,CHAR(STRING_ELT(VECTOR_ELT(probe_trt_levels,0),j)));
	  tmp_str = strcat(tmp_str,":");
	  SET_STRING_ELT(probe_colnames,j,mkChar(tmp_str));
	  Free(tmp_str);
	}
      } else if (model->strata[4] == 4){
	for (j =0; j < model->max_probe_treatment_factor+1;j++){
	  tmp_str = (char *)Calloc(strlen(CHAR(STRING_ELT(getAttrib(probe_trt_levels,R_NamesSymbol),0)))+ strlen(CHAR(STRING_ELT(VECTOR_ELT(probe_trt_levels,0),j)))+16,char);
	  tmp_str = strcpy(tmp_str,CHAR(STRING_ELT(getAttrib(probe_trt_levels,R_NamesSymbol),0)));
	  tmp_str = strcat(tmp_str,"_");
	  tmp_str = strcat(tmp_str,CHAR(STRING_ELT(VECTOR_ELT(probe_trt_levels,0),j)));
	  tmp_str = strcat(tmp_str,":probe.type_PM:");
	  SET_STRING_ELT(probe_colnames,2*j,mkChar(tmp_str));
	  Free(tmp_str);
	  tmp_str = (char *)Calloc(strlen(CHAR(STRING_ELT(getAttrib(probe_trt_levels,R_NamesSymbol),0)))+ strlen(CHAR(STRING_ELT(VECTOR_ELT(probe_trt_levels,0),j)))+16,char);
	  tmp_str = strcpy(tmp_str,CHAR(STRING_ELT(getAttrib(probe_trt_levels,R_NamesSymbol),0)));
	  tmp_str = strcat(tmp_str,"_");
	  tmp_str = strcat(tmp_str,CHAR(STRING_ELT(VECTOR_ELT(probe_trt_levels,0),j)));
	  tmp_str = strcat(tmp_str,":probe.type_MM:");
	  SET_STRING_ELT(probe_colnames,2*j+1,mkChar(tmp_str));
	  Free(tmp_str);
	}
      }
      SET_VECTOR_ELT(probe_dimnames,0,probe_rownames);
      SET_VECTOR_ELT(probe_dimnames,1,probe_colnames);
      setAttrib(VECTOR_ELT(Routput->probe_coef,i),R_DimNamesSymbol,probe_dimnames);
      setAttrib(VECTOR_ELT(Routput->probe_SE,i),R_DimNamesSymbol,probe_dimnames);
      UNPROTECT(4);
    }

    Free(tmp_str2);
  }
  


  

  /*** Now put names on the optional components ****/  
  if (model->which_parameter_types[1] | model->which_parameter_types[2]){

    PROTECT(chip_coef_dimnames = allocVector(VECSXP,2));  
    PROTECT(chip_coef_rownames= allocVector(STRSXP,data->n_probesets));
    copyVector(chip_coef_rownames,names);
    if (model->which_parameter_types[2]){
      /* sample effects */ 
      
      if (model->constraints[2] == -1){
	PROTECT(chip_coef_colnames = allocVector(STRSXP,data->n_arrays-1));
	for (i =0; i <data->n_arrays-1;i++){
	  tmp_str=(char *)Calloc(strlen(CHAR(STRING_ELT(sampleNames,i)))+1,char);
	  tmp_str=strcpy(tmp_str,CHAR(STRING_ELT(sampleNames,i)));
	  SET_STRING_ELT(chip_coef_colnames,i,mkChar(tmp_str));
	  Free(tmp_str);
	}
      } else if (model->constraints[2] == 1){
	PROTECT(chip_coef_colnames = allocVector(STRSXP,data->n_arrays-1));
	for (i =0; i <data->n_arrays-1;i++){
	  tmp_str=(char *)Calloc(strlen(CHAR(STRING_ELT(sampleNames,i+1)))+1,char);
	  tmp_str=strcpy(tmp_str,CHAR(STRING_ELT(sampleNames,i+1)));
	  SET_STRING_ELT(chip_coef_colnames,i,mkChar(tmp_str));
	  Free(tmp_str);
	}
      } else {
	PROTECT(chip_coef_colnames = allocVector(STRSXP,data->n_arrays));
	for (i =0; i < data->n_arrays;i++){
	  tmp_str=(char *)Calloc(strlen(CHAR(STRING_ELT(sampleNames,i)))+1,char);
	  tmp_str=strcpy(tmp_str,CHAR(STRING_ELT(sampleNames,i)));
	  SET_STRING_ELT(chip_coef_colnames,i,mkChar(tmp_str));	
	  Free(tmp_str);
	}
      }
    } else if (model->which_parameter_types[1]){
      /* chip-covariates */
      PROTECT(chip_coef_colnames = allocVector(STRSXP,model->n_chiplevelcovariates));
      chip_covariate_names = VECTOR_ELT(getAttrib(chipcovariates,R_DimNamesSymbol),1);
      for (i =0; i <model->n_chiplevelcovariates ;i++){
	tmp_str=(char *)Calloc(strlen(CHAR(STRING_ELT(chip_covariate_names,i)))+1,char);
	tmp_str=strcpy(tmp_str,CHAR(STRING_ELT(chip_covariate_names,i)));
	SET_STRING_ELT(chip_coef_colnames,i,mkChar(tmp_str));
	Free(tmp_str);
      }
    } else {
      PROTECT(chip_coef_colnames = allocVector(STRSXP,0));
    }
    
    SET_VECTOR_ELT(chip_coef_dimnames,1,chip_coef_colnames);
    SET_VECTOR_ELT(chip_coef_dimnames,0,chip_coef_rownames);
    setAttrib(Routput->chip_coef, R_DimNamesSymbol, chip_coef_dimnames); 
    setAttrib(Routput->chip_SE,R_DimNamesSymbol, chip_coef_dimnames);
    UNPROTECT(3);
  }

  if (store->residSE){
    PROTECT(residSE_dimnames = allocVector(VECSXP,2));
    PROTECT(residSE_rownames= allocVector(STRSXP,data->n_probesets));
    PROTECT(residSE_colnames = allocVector(STRSXP,2));
    copyVector(residSE_rownames,names);
    SET_STRING_ELT(residSE_colnames,0,mkChar("Resid SE"));
    SET_STRING_ELT(residSE_colnames,1,mkChar("df"));
    SET_VECTOR_ELT(residSE_dimnames,1,residSE_colnames);
    SET_VECTOR_ELT(residSE_dimnames,0,residSE_rownames);
    setAttrib(Routput->residSE,R_DimNamesSymbol, residSE_dimnames);
    UNPROTECT(3);
  }
  
  if (store->varcov){ 
  
    
    /* name the elements of the varcov output */
    if (store->varcov ==1){  
      /* name the list */
      if (model->which_parameter_types[0] | model->which_parameter_types[1] | model->which_parameter_types[2]){
	PROTECT(varcov_names= allocVector(STRSXP,data->n_probesets));
	copyVector(varcov_names,names);
	setAttrib(Routput->varcov,R_NamesSymbol,varcov_names);  
	UNPROTECT(1);
      }
      /* Chip-level only - Note that this would include the intercept */
      if (model->which_parameter_types[0]){
	/* Intercept */
	if (model->which_parameter_types[2]){
	  /* with samples effects */
	  for (i=0; i < data->n_probesets;i++){
	    PROTECT(varcov_dimnames = allocVector(VECSXP,2));
	    PROTECT(varcov_colnames = allocVector(STRSXP,data->n_arrays));
	    PROTECT(varcov_rownames = allocVector(STRSXP,data->n_arrays));
	    SET_STRING_ELT(varcov_colnames,0,mkChar("Intercept"));
	    SET_STRING_ELT(varcov_rownames,0,mkChar("Intercept"));
	    if (model->constraints[2] == 1){
	      for (j =0; j <data->n_arrays-1;j++){
		tmp_str=(char *)Calloc(strlen(CHAR(STRING_ELT(sampleNames,j+1)))+1,char);
		tmp_str=strcpy(tmp_str,CHAR(STRING_ELT(sampleNames,j+1)));
		SET_STRING_ELT(varcov_rownames,j+1,mkChar(tmp_str));
		SET_STRING_ELT(varcov_colnames,j+1,mkChar(tmp_str));
		Free(tmp_str);
	      }
	    } else {
	       for (j =0; j <data->n_arrays-1;j++){
		tmp_str=(char *)Calloc(strlen(CHAR(STRING_ELT(sampleNames,j)))+1,char);
		tmp_str=strcpy(tmp_str,CHAR(STRING_ELT(sampleNames,j)));
		SET_STRING_ELT(varcov_rownames,j+1,mkChar(tmp_str));
		SET_STRING_ELT(varcov_colnames,j+1,mkChar(tmp_str));
		Free(tmp_str);
	       }
	    }
	    SET_VECTOR_ELT(varcov_dimnames,0,varcov_rownames);
	    SET_VECTOR_ELT(varcov_dimnames,1,varcov_colnames);
	    setAttrib(VECTOR_ELT(Routput->varcov,i),R_DimNamesSymbol,varcov_dimnames);
	    UNPROTECT(3);
	  }
	} else if (model->which_parameter_types[1]){
	  /* which chip-level variables */
	  for (i=0; i < data->n_probesets;i++){
	    PROTECT(varcov_dimnames = allocVector(VECSXP,2));
	    PROTECT(varcov_colnames = allocVector(STRSXP,model->n_chiplevelcovariates+1));
	    PROTECT(varcov_rownames = allocVector(STRSXP,model->n_chiplevelcovariates+1));
	    SET_STRING_ELT(varcov_colnames,0,mkChar("Intercept"));
	    SET_STRING_ELT(varcov_rownames,0,mkChar("Intercept"));
	    chip_covariate_names = STRING_ELT(getAttrib(chipcovariates,R_DimNamesSymbol),1);
	    for (j =0; j <model->n_chiplevelcovariates ;j++){
	      tmp_str=(char *)Calloc(strlen(CHAR(STRING_ELT(chip_covariate_names,j)))+1,char);
	      tmp_str=strcpy(tmp_str,CHAR(STRING_ELT(chip_covariate_names,j)));
	      SET_STRING_ELT(varcov_rownames,j+1,mkChar(tmp_str));
	      SET_STRING_ELT(varcov_colnames,j+1,mkChar(tmp_str));
	      Free(tmp_str);
	    }
	    SET_VECTOR_ELT(varcov_dimnames,0,varcov_rownames);
	    SET_VECTOR_ELT(varcov_dimnames,1,varcov_colnames);
	    setAttrib(VECTOR_ELT(Routput->varcov,i),R_DimNamesSymbol,varcov_dimnames);
	    UNPROTECT(3);
	  }
	} else {
	  /* Only the intercept */
	  for (i=0; i < data->n_probesets;i++){
	    PROTECT(varcov_dimnames = allocVector(VECSXP,2));
	    PROTECT(varcov_colnames = allocVector(STRSXP,1));
	    PROTECT(varcov_rownames = allocVector(STRSXP,1));
	    SET_STRING_ELT(varcov_colnames,0,mkChar("Intercept"));
	    SET_STRING_ELT(varcov_rownames,0,mkChar("Intercept"));
	    SET_VECTOR_ELT(varcov_dimnames,0,varcov_rownames);
	    SET_VECTOR_ELT(varcov_dimnames,1,varcov_colnames);
	    setAttrib(VECTOR_ELT(Routput->varcov,i),R_DimNamesSymbol,varcov_dimnames);
	    UNPROTECT(3);
	  }
	}
      } else {
	/* no intercept */
	if (model->which_parameter_types[2]){
	  /* with samples effects */
	  for (i=0; i < data->n_probesets;i++){
	    PROTECT(varcov_dimnames = allocVector(VECSXP,2));
	    PROTECT(varcov_colnames = allocVector(STRSXP,data->n_arrays));
	    PROTECT(varcov_rownames = allocVector(STRSXP,data->n_arrays));
	    for (j =0; j <data->n_arrays;j++){
	      tmp_str=(char *)Calloc(strlen(CHAR(STRING_ELT(sampleNames,j)))+1,char);
	      tmp_str=strcpy(tmp_str,CHAR(STRING_ELT(sampleNames,j)));
	      SET_STRING_ELT(varcov_rownames,j,mkChar(tmp_str));
	      SET_STRING_ELT(varcov_colnames,j,mkChar(tmp_str));
	      Free(tmp_str);
	    }
	    
	    SET_VECTOR_ELT(varcov_dimnames,0,varcov_rownames);
	    SET_VECTOR_ELT(varcov_dimnames,1,varcov_colnames);
	    setAttrib(VECTOR_ELT(Routput->varcov,i),R_DimNamesSymbol,varcov_dimnames);
	    UNPROTECT(3);
	  }
	} else if (model->which_parameter_types[1]){
	  for (i=0; i < data->n_probesets;i++){
	    PROTECT(varcov_dimnames = allocVector(VECSXP,2));
	    PROTECT(varcov_colnames = allocVector(STRSXP,model->n_chiplevelcovariates));
	    PROTECT(varcov_rownames = allocVector(STRSXP,model->n_chiplevelcovariates));
	    SET_STRING_ELT(varcov_colnames,0,mkChar("Intercept"));
	    SET_STRING_ELT(varcov_rownames,0,mkChar("Intercept"));
	    chip_covariate_names = VECTOR_ELT(getAttrib(chipcovariates,R_DimNamesSymbol),1);
	    for (j =0; j <model->n_chiplevelcovariates ;j++){
	      tmp_str=(char *)Calloc(strlen(CHAR(STRING_ELT(chip_covariate_names,j)))+1,char);
	      tmp_str=strcpy(tmp_str,CHAR(STRING_ELT(chip_covariate_names,j)));
	      SET_STRING_ELT(varcov_rownames,j,mkChar(tmp_str));
	      SET_STRING_ELT(varcov_colnames,j,mkChar(tmp_str));
	      Free(tmp_str);
	    }
	    SET_VECTOR_ELT(varcov_dimnames,0,varcov_rownames);
	    SET_VECTOR_ELT(varcov_dimnames,1,varcov_colnames);
	    setAttrib(VECTOR_ELT(Routput->varcov,i),R_DimNamesSymbol,varcov_dimnames);
	    UNPROTECT(3);
	  }	   
	}
      }
    } else if (store->varcov == 2){
      /* The complete varcov matrix */
      /* name the items of the list */
      PROTECT(varcov_names= allocVector(STRSXP,data->n_probesets));
      copyVector(varcov_names,names);
      setAttrib(Routput->varcov,R_NamesSymbol,varcov_names);  
      UNPROTECT(1);

      /** remember that the order is intercept/MM/chip-covariates/samples/probe.types/probes **/
      for (i=0; i < data->n_probesets;i++){
	curcol = 0;
	PROTECT(varcov_dimnames = allocVector(VECSXP,2)); 
	dim = getAttrib(VECTOR_ELT(Routput->varcov,i),R_DimSymbol);
	PROTECT(varcov_colnames = allocVector(STRSXP,INTEGER(dim)[0]));
	PROTECT(varcov_rownames = allocVector(STRSXP,INTEGER(dim)[0]));
	if (model->which_parameter_types[0]){
	  SET_STRING_ELT(varcov_colnames,curcol,mkChar("Intercept"));
	  SET_STRING_ELT(varcov_rownames,curcol,mkChar("Intercept"));
	  curcol++;
	}
	if (model->mmorpm_covariate != 0){
	  if (model->mmorpm_covariate < 0){
	    SET_STRING_ELT(varcov_colnames,curcol,mkChar("PM"));
	  } else {
	    SET_STRING_ELT(varcov_rownames,curcol,mkChar("MM"));
	  }
	  curcol++;
	}
	if (model->which_parameter_types[1]){
	  PROTECT(chip_coef_colnames = allocVector(STRSXP,model->n_chiplevelcovariates));
	  copyVector(chip_coef_colnames,VECTOR_ELT(getAttrib(chipcovariates,R_DimNamesSymbol),1));
	  for (j =0; j < model->n_chiplevelcovariates; j++){
	    SET_STRING_ELT(varcov_colnames,curcol+j,STRING_ELT(chip_coef_colnames,j));
	    SET_STRING_ELT(varcov_rownames,curcol+j,STRING_ELT(chip_coef_colnames,j));
	  }
	  curcol+=model->n_chiplevelcovariates;
	  UNPROTECT(1);
	}
	if (model->which_parameter_types[2]){
	  PROTECT(chip_coef_colnames = allocVector(STRSXP,INTEGER(getAttrib(Routput->chip_coef,R_DimSymbol))[1]));
	  copyVector(chip_coef_colnames,VECTOR_ELT(getAttrib(Routput->chip_coef,R_DimNamesSymbol),1));
	  for (j =0; j < INTEGER(getAttrib(Routput->chip_coef,R_DimSymbol))[1]; j++){
	    SET_STRING_ELT(varcov_colnames,curcol+j,STRING_ELT(chip_coef_colnames,j));
	    SET_STRING_ELT(varcov_rownames,curcol+j,STRING_ELT(chip_coef_colnames,j));
	  }
	  curcol+=INTEGER(getAttrib(Routput->chip_coef,R_DimSymbol))[1];
	  UNPROTECT(1);
	}
	if (model->which_parameter_types[3]){
	  n_const_col=0;
	  if (model->which_parameter_types[0]){
	    n_const_col++;
	  }
	  if (model->mmorpm_covariate !=0){
	    n_const_col++;
	  }
	  for (j =0; j < INTEGER(getAttrib(Routput->const_coef,R_DimSymbol))[1]-n_const_col; j++){
	    SET_STRING_ELT(varcov_colnames,curcol+j,mkChar(CHAR(STRING_ELT(VECTOR_ELT(getAttrib(Routput->const_coef,R_DimNamesSymbol),1),j+n_const_col))));
	    SET_STRING_ELT(varcov_rownames,curcol+j,mkChar(CHAR(STRING_ELT(VECTOR_ELT(getAttrib(Routput->const_coef,R_DimNamesSymbol),1),j+n_const_col))));
	  }
	  curcol+=INTEGER(getAttrib(Routput->const_coef,R_DimSymbol))[1] - n_const_col;
	}
	if (model->which_parameter_types[4]){
	  dim = getAttrib(VECTOR_ELT(Routput->probe_coef,i),R_DimSymbol);
	  if (INTEGER(dim)[1] == 1){
	    for (j=0; j < INTEGER(dim)[0]; j++){
	       SET_STRING_ELT(varcov_colnames,curcol+j,mkChar(CHAR(STRING_ELT(VECTOR_ELT(getAttrib(VECTOR_ELT(Routput->probe_coef,i),R_DimNamesSymbol),0),j))));
	       SET_STRING_ELT(varcov_rownames,curcol+j,mkChar(CHAR(STRING_ELT(VECTOR_ELT(getAttrib(VECTOR_ELT(Routput->probe_coef,i),R_DimNamesSymbol),0),j))));
	    }
	  } else {
	    for (k=0; k < INTEGER(dim)[1]; k++){
	      for (j=0; j < INTEGER(dim)[0]; j++){
		tmp_str = (char *)Calloc(strlen(CHAR(STRING_ELT(VECTOR_ELT(getAttrib(VECTOR_ELT(Routput->probe_coef,i),R_DimNamesSymbol),0),j))) + strlen(CHAR(STRING_ELT(VECTOR_ELT(getAttrib(VECTOR_ELT(Routput->probe_coef,i),R_DimNamesSymbol),1),k))) + 2,char);
		tmp_str = strcpy(tmp_str,CHAR(STRING_ELT(VECTOR_ELT(getAttrib(VECTOR_ELT(Routput->probe_coef,i),R_DimNamesSymbol),1),k)));
		tmp_str = strcat(tmp_str,CHAR(STRING_ELT(VECTOR_ELT(getAttrib(VECTOR_ELT(Routput->probe_coef,i),R_DimNamesSymbol),0),j)));
		SET_STRING_ELT(varcov_colnames,curcol+k*INTEGER(dim)[0]+j,mkChar(tmp_str));
		SET_STRING_ELT(varcov_rownames,curcol+k*INTEGER(dim)[0]+j,mkChar(tmp_str));
		Free(tmp_str);
	      }
				       
	    }
	  }
	}
	SET_VECTOR_ELT(varcov_dimnames,0,varcov_rownames);
	SET_VECTOR_ELT(varcov_dimnames,1,varcov_colnames);
	setAttrib(VECTOR_ELT(Routput->varcov,i),R_DimNamesSymbol,varcov_dimnames);
	UNPROTECT(3);
      }
    }
  }



  

  PROTECT(weights_names = allocVector(STRSXP,2));
  SET_STRING_ELT(weights_names,0,mkChar("PM.weights"));
  SET_STRING_ELT(weights_names,1,mkChar("MM.weights"));
  setAttrib(Routput->weights,R_NamesSymbol,weights_names);
  UNPROTECT(1);

  if (store->weights){
     PROTECT(weights_dimnames = allocVector(VECSXP,2));
     PROTECT(weights_rownames = allocVector(STRSXP,data->n_probes));
     PROTECT(weights_colnames = allocVector(STRSXP,data->n_arrays));
     copyVector(weights_colnames,sampleNames);
     copyVector(weights_rownames,ProbeNamesVec);
     SET_VECTOR_ELT(weights_dimnames,0,weights_rownames);
     SET_VECTOR_ELT(weights_dimnames,1,weights_colnames);
     if (model->response_variable >=0){
       setAttrib(VECTOR_ELT(Routput->weights,0),R_DimNamesSymbol, weights_dimnames);
     }
     if (model->response_variable <=0){
       setAttrib(VECTOR_ELT(Routput->weights,1),R_DimNamesSymbol, weights_dimnames);
     }
     
     UNPROTECT(3);
  }

  PROTECT(residuals_names = allocVector(STRSXP,2));
  SET_STRING_ELT(residuals_names,0,mkChar("PM.resid"));
  SET_STRING_ELT(residuals_names,1,mkChar("MM.resid"));
  setAttrib(Routput->residuals,R_NamesSymbol,residuals_names);
  UNPROTECT(1);

  if (store->residuals){ 
    PROTECT(residuals_dimnames = allocVector(VECSXP,2));
    PROTECT(residuals_rownames = allocVector(STRSXP,data->n_probes));
    PROTECT(residuals_colnames = allocVector(STRSXP,data->n_arrays));
    copyVector(residuals_colnames,sampleNames);
    copyVector(residuals_rownames,ProbeNamesVec);
    SET_VECTOR_ELT(residuals_dimnames,0,residuals_rownames);
    SET_VECTOR_ELT(residuals_dimnames,1,residuals_colnames);
    if (model->response_variable >=0){
      setAttrib(VECTOR_ELT(Routput->residuals,0),R_DimNamesSymbol, residuals_dimnames);
    }
    if (model->response_variable <=0){
      setAttrib(VECTOR_ELT(Routput->residuals,1),R_DimNamesSymbol, residuals_dimnames);
    }
    UNPROTECT(3);
  }
  
  UNPROTECT(1);
}




SEXP rlm_PLMset(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes, SEXP R_model, SEXP outputparam, SEXP modelparam, SEXP verbosity){
  
  int i;
  int verbosity_level;


  PLM_outputsettings *store = (PLM_outputsettings *)Calloc(1,PLM_outputsettings);
  PLM_Datagroup *data = (PLM_Datagroup *)Calloc(1,PLM_Datagroup);
  PLM_output *output = (PLM_output *)Calloc(1,PLM_output);
  PLM_model_parameters *model = (PLM_model_parameters *)Calloc(1,PLM_model_parameters);
  PLMRoutput *Routput = (PLMRoutput *)Calloc(1,PLMRoutput);
  
  SEXP dim1,dim2;

  SEXP output_list;
  SEXP param;
  SEXP chipcovariates;

  SEXP sampleNames;
  /* this fills up "data" */
  /* organise data to be passed to model fitting routines */
  
  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol));
  
  data->n_probes = INTEGER(dim1)[0];
  data->n_arrays = INTEGER(dim1)[1];
  
  data->PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
  data->MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));
  data->n_probesets = INTEGER(N_probes)[0];
  
  /* Get the names corresponding to each row */
    
  data->ProbeNames = (const char **)Calloc(data->n_probes,const char *);
  for (i =0; i < data->n_probes; i++){
    data->ProbeNames[i] = CHAR(STRING_ELT(ProbeNamesVec,i));
  }
  
  verbosity_level = asInteger(verbosity);


  /********************* This part of the code completely fills up "model" *************************************/
  /* figure out what the covariate matrix should be based on the rlm_model_type information and chipcovariates */
 
  /*  these are the parameters related to the rlm procedure ie the model fitting procedure */
  param = GetParameter(modelparam,"psi.type");
  model->psi_code = asInteger(param);
  param = GetParameter(modelparam,"se.type");
  model->se_method = asInteger(param);
  param = GetParameter(modelparam,"psi.k");
  model->psi_k = asReal(param);
  param = GetParameter(modelparam,"max.its");
  model->n_rlm_iterations = asInteger(param);
  param = GetParameter(modelparam,"init.method");
  if (strcmp(CHAR(STRING_ELT(param,0)),"ls") == 0){
    model->init_method = 0;
  } else if (strcmp(CHAR(STRING_ELT(param,0)),"median.polish") == 0){
    model->init_method = 1;
  } else if (strcmp(CHAR(STRING_ELT(param,0)),"Huber") == 0){
    model->init_method = 2;
  }
  
  param = GetParameter(modelparam,"weights.chip");
  if (isNull(param)){
    model->input_chip_weights = NULL;
  } else {
    model->input_chip_weights = NUMERIC_POINTER(param);
  }
  

  param = GetParameter(modelparam,"weights.probe");
  if (isNull(param)){
    model->input_probe_weights = NULL;
  } else {
    model->input_probe_weights = NUMERIC_POINTER(param);
  }

  param = GetParameter(modelparam,"trans.fn");
  if (strcmp(CHAR(STRING_ELT(param,0)),"log2") == 0){
    model->trans_fn = 0;
  } else if (strcmp(CHAR(STRING_ELT(param,0)),"ln") == 0){
    model->trans_fn = 1;
  } else if (strcmp(CHAR(STRING_ELT(param,0)),"loge") == 0){
    model->trans_fn = 1;
  } else if (strcmp(CHAR(STRING_ELT(param,0)),"log10") == 0){
    model->trans_fn = 2;
  } else if (strcmp(CHAR(STRING_ELT(param,0)),"sqrt") == 0){
    model->trans_fn = 3;
  } else if (strcmp(CHAR(STRING_ELT(param,0)),"cuberoot") == 0){
    model->trans_fn = 4;
  } else {
    error("%s is unknown transformation\n",CHAR(STRING_ELT(param,0)));
  }

  /* this code deals with the model itself */
  
  param = GetParameter(R_model,"mmorpm.covariate");
  model->mmorpm_covariate = asInteger(param);
  
  param = GetParameter(R_model,"response.variable");
  model->response_variable = asInteger(param);

  param = GetParameter(R_model,"which.parameter.types");
  model->which_parameter_types = INTEGER(param);
  param = GetParameter(R_model,"strata");
  model->strata = INTEGER(param);
  param = GetParameter(R_model,"constraints");
  model->constraints = INTEGER(param);

  param = GetParameter(R_model,"probe.type.trt.factor");
  model->probe_type_treatment_factor = INTEGER_POINTER(param);
  param = GetParameter(R_model,"max.probe.type.trt.factor");
  model->max_probe_type_treatment_factor = asInteger(param);
  
  param = GetParameter(R_model,"probe.trt.factor");
  model->probe_treatment_factor = INTEGER_POINTER(param);
  param = GetParameter(R_model,"max.probe.trt.factor");
  model->max_probe_treatment_factor = asInteger(param);

  chipcovariates = GetParameter(R_model,"chipcovariates");
  model->chiplevelcovariates = NUMERIC_POINTER(chipcovariates);
  
  PROTECT(dim2 = getAttrib(chipcovariates,R_DimSymbol));
  model->n_chiplevelcovariates = INTEGER(dim2)[1];
  model->n_arrays = data->n_arrays;


  /***************** Now load up "output" ******************************************/
  /* figure out what optional features we want to output */
  param = GetParameter(outputparam,"weights");
  store->weights = asInteger(param);
  param = GetParameter(outputparam,"residuals");
  store->residuals = asInteger(param);
  param = GetParameter(outputparam,"resid.SE");
  store->residSE = asInteger(param);
  param = GetParameter(outputparam,"varcov");
  if (strcmp(CHAR(STRING_ELT(param,0)),"none") == 0){
    store->varcov = 0;
  } else if (strcmp(CHAR(STRING_ELT(param,0)),"chiplevel") == 0){
    store->varcov =1;
  } else if (strcmp(CHAR(STRING_ELT(param,0)),"all") == 0){
    store->varcov =2;
  }


  /* Make space for output */
  rlm_PLM_alloc_space(Routput, output, store, data, model);
  
  

  /* now go actually fit the model */

  if (verbosity_level > 0){
    Rprintf("Fitting models\n");
  }
  
  do_PLM_rlm(data, model, output,store);
  
  /* now lets put names on the output matrices */
  /* First the chip coef matrix */


  sampleNames = VECTOR_ELT(getAttrib(PMmat,R_DimNamesSymbol),1);


  rlm_PLMset_nameoutput(Routput,output,store, data,model,R_model, sampleNames, ProbeNamesVec,chipcovariates);
  
  
  /*Now lets create the output_list */

  PROTECT(output_list = allocVector(VECSXP,10));

  SET_VECTOR_ELT(output_list,0,Routput->chip_coef);
  SET_VECTOR_ELT(output_list,1,Routput->probe_coef);
  SET_VECTOR_ELT(output_list,2,Routput->weights);
  SET_VECTOR_ELT(output_list,3,Routput->chip_SE);
  SET_VECTOR_ELT(output_list,4,Routput->probe_SE);
  SET_VECTOR_ELT(output_list,5,Routput->const_coef);
  SET_VECTOR_ELT(output_list,6,Routput->const_SE);
  SET_VECTOR_ELT(output_list,7,Routput->residuals);
  SET_VECTOR_ELT(output_list,8,Routput->residSE);
  SET_VECTOR_ELT(output_list,9,Routput->varcov);
  UNPROTECT(Routput->nprotected + 3);


  Free(output->out_weights);
  Free(output->out_resids);
  for ( i =0; i < data->n_probesets; i++){
    Free(output->outnames[i]);
  }
  
  Free(output->outnames);
  Free(output->out_probeparams);
  Free(output->out_probe_SE);
  Free(output->out_varcov);
  Free(data->ProbeNames);
  Free(data);
  Free(output);
  Free(Routput);
  Free(store);
  Free(model);
  
  return output_list;
  
}








SEXP R_rlm_PLMset_c(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes, SEXP R_model, SEXP outputparam, SEXP modelparam,  SEXP bg_flag, SEXP bg_type, SEXP background_parameters, SEXP norm_flag, SEXP norm_type,SEXP norm_parameters, SEXP verbosity){


  SEXP dim1,rlmPLMresults; 
  int rows,cols;

  /*Create a copy matrix to work on. Allows us to modify data in background and normalization steps without affecting original data */
  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
 
  /* If Background correction do it */
  if (INTEGER(bg_flag)[0]){
    PMmat = pp_background(PMmat, MMmat, ProbeNamesVec,N_probes,bg_type,background_parameters,verbosity);
  }

  /* If Normalization do it */
  if (INTEGER(norm_flag)[0]){
    PMmat = pp_normalize(PMmat, MMmat, ProbeNamesVec,N_probes,norm_type, norm_parameters,verbosity);
  }
  
  /* Now do RLM fit */

  rlmPLMresults = rlm_PLMset(PMmat, MMmat, ProbeNamesVec, N_probes, R_model, outputparam, modelparam,verbosity);
  
  UNPROTECT(1);

  return rlmPLMresults;
}









