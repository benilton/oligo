/*********************************************************************
 **
 ** file: do_PLMrlm.c
 **
 ** Aim: fit robust linear models for the PLMset object.
 **
 ** Copyright (C) 2003-2007 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bmb@bmbolstad.com>
 ** 
 ** created on: Jan 21, 2003
 **
 ** Last modified: Nov 30, 2007
 **
 ** the aim here will be to fit specified robust linear models
 ** to affy data. this file actually contains the code that 
 ** takes preprocessed data and actually carries out the fitting of the models
 ** 
 ** in particular it breaks the data into "blocks" where a block is
 ** a probes by chips matrix of preprocessed PM intensites. to each block it fits 
 ** the specified model.
 **
 ** the probe effects part of the design matrix is handled here.
 **
 ** Note that the function calls here suffer a little bit from "fortran-itis", ie
 ** lots of parameters passed in and out, should probably to rewritten in a 
 ** cleaner form at somepoint. (THIS HAS BEEN DONE)
 **
 ** Modification history
 **
 ** Jan 21, 2003 - Initial version.
 ** Jan 24, 2003 - Better handling of more general models.
 ** Jan 27, 2003 - Computation of standard errors
 ** Jan 28, 2003 - ability to select between se methods
 ** Jan 31, 2003 - add in a missing #include "rlm.h"
 **                note that we should be able to optimize routines, by testing 
 **                whether we really need to reallocate and initialize the X matrix 
 **                each iteration, typically almost all the probesets on a
 **                chip have same number of probesets, and X is not changed by the RLM
 **                fitting routine so we could actually keep it between blocks only changing it when
 **                the number of probes had changed.  <- TO DO WHEN WE WORK ON OPTIMIZATIONS
 **                Documentation Updates.
 **                Changed return type of rlm_PLM_block to void (At some point perhaps this should be an integer error flag)
 ** Feb 1, 2003 -  rework the design matrix code so that it is only reallocated/reinitialised if required.
 **                this results in a new function that actually allocates/assigns data values for the X matrix
 **                and changes in the parameter list passed to rlm_PLM_block
 ** Feb 11, 2003 - remove a loop that did nothing. Clean up some unused code.
 **                Optimizations in rlm_PLM_block to remove unnecessary duplication/copying of data.
 **                Clean-ups in the commenting.
 ** Feb 15, 2003 - generalize the class of models that can be fit by altering the design matrix
 **                realloc function. We also break the results copying parts of do_PLMrlm()
 **                out into a new function copy_PLM_results (this is for the model parameters,
 **                weights are still handled in the main function
 **                n,p definitions need cleaning up in this case.
 ** Feb 17, 2003 - Continue clean-up of n,p. add in the ability to fit models with no probe parameters
 ** Feb 18. 2003 - squash seg-faulting bug in rlm_designamtrix_realloc..... for method = 21
 ** Feb 23, 2003 - add implementations for method=10,11 into rlm_design_matrix_realloc(), copy_PLM_results()
 ** Feb 24, 2003 - comment out some unused variables. Trying to reduce compiler warnings.
 ** Apr 04, 2003  - make the number of rows in the PM matrix be dynamic.
 ** Jun 04, 2003  - add mechanism for more general psi functions.
 ** Jun 11, 2003  - make sure that the standard error call allows other psi functions.
 ** Jul 23, 2003 - remove one last compiler warning.
 ** Sep 02, 2003 - residuals are now stored
 ** Sep 05, 2003 - clean up in how parameters are passed to  do_PLMrlm
 ** Sep 06, 2003 - introduced the struct modelfit. It is used to group together
 **                items related to the current model (for each probeset) being
 **                fitted. More clean up in how parameters are passed between
 **                functions. residSE now outputted.
 ** Sep 07, 2003 - chip level part of variance covariance matrix returned.
 ** Sep 08, 2003 - more work on returning variance covariance
 ** Sep 13, 2003 - Modify rlm_PLM_block so it handles number of iterations
 **                and initialization method
 ** Oct 12, 2003 - fixed declaration order error             
 ** Apr 5, 2004  - Changed a malloc to a Calloc        
 ** May 11, 2004 - fixed a minor memory leak.
 ** May 27, 2004 - if the -1 + samples + probes model (with probes constraint sum to zero) use
 **                a different and faster algorithm
 ** July 9, 2004 - start integrating new structure
 ** Mar 12, 2005 - changed the loop in do_PLMrlm
 ** Aug 23, 2006 - fix a bug in checkDefaultModel. It incorrectly said that PM ~ -1 + samples + treatment:probes was a valid default model.
 ** Nov 30, 2007 - comment out old unused code. Since it is defunct it will likely be removed in the future.
 ** Nov 1, 2008 - remove old defunct code
 **               return constrained probe coef and se for default model
 **
 *********************************************************************/


#include "do_PLMrlm.h"
#include "rlm_se.h"
#include "rlm.h"
#include "psi_fns.h"
#include "common_types.h"
#include "PLM_medianpolish.h"
#include "PLM_modelmatrix.h"
#include "transfns.h"

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**********************************************************************
 **
 ** the modelfit struct is used for storing information about the 
 ** current model (the one being fitted individually to each probeset.
 **
 **********************************************************************/

typedef struct{
  double *cur_params;            /* storage */
  double *cur_se_estimates;
  double *cur_weights;
  double *cur_resids;
  double *cur_varcov;
  double *cur_residSE;
  int *cur_rows;  /* indices in the data matrix to use for current model */
  double *X;      /* design matrix */
  int n;          /* number of observations */
  int p;          /* number of parameters */
  int nprobes;    /* number of probes in current probeset */
  
} modelfit;




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

static void  rlm_PLM_probeset(const PLM_model_parameters *model, const PLM_Datagroup *data, PLM_modelfit *current, int *current_rows){
  
  int i, j;
  int isDefaultModel = 0;

  double *Y = Calloc(current->n,double);
  /*  double lg2 = log(2.0); */ /* Cache hopefully for speed :) */
  double *input_weights=NULL;
  pt2trans transfn = transFunc(model->trans_fn);
  
  /* log2 transform and create Y vector */
  
  if (model->response_variable == 1){
    /* PM response variable */
    for (j = 0; j < data->n_arrays; j++){
      for (i =0; i < current->nprobes; i++){
	Y[j*current->nprobes + i] = transfn(data->PM[j*data->n_probes + current_rows[i]]);                                /* log(data->PM[j*data->n_probes + current_rows[i]])/lg2; */
      }
    }
  } else if (model->response_variable == -1){
    /* MM response variable */
    for (j = 0; j < data->n_arrays; j++){
      for (i =0; i < current->nprobes; i++){
	Y[j*current->nprobes + i] = transfn(data->MM[j*data->n_probes + current_rows[i]]);                             /* log(data->MM[j*data->n_probes + current_rows[i]])/lg2; */
      }
    }
  } else {
    /* Both PM and MM response variables - we stack them in Y */
    for (j = 0; j < data->n_arrays; j++){
      for (i =0; i < current->nprobes; i++){
	Y[j*current->nprobes + i] = transfn(data->PM[j*data->n_probes + current_rows[i]]);       /* log(data->PM[j*data->n_probes + current_rows[i]])/lg2; */
      }
    }
    for (j = 0; j < data->n_arrays; j++){
      for (i =0; i < current->nprobes; i++){
	Y[current->nprobes*data->n_arrays + j*current->nprobes + i] = transfn(data->MM[j*data->n_probes + current_rows[i]]);  /* log(data->MM[j*data->n_probes + current_rows[i]])/lg2; */
      }
    }
  }

  if ((model->input_chip_weights == NULL) && (model->input_probe_weights == NULL)){
    if (model->init_method == 1){
      /* median polish */
      /* No longer supported */
    } else if (model->init_method ==2){
      /* fully iterated huber regression */
      rlm_fit(current->X,Y, current->n, current->p, current->cur_params, current->cur_resids, current->cur_weights, PsiFunc(0),1.345,20,0);
      
    }



    isDefaultModel = checkDefaultModel(model);


    if (isDefaultModel){
      
      current->cur_params = Realloc(current->cur_params,current->nprobes+data->n_arrays,double);
      current->cur_se_estimates = Realloc(current->cur_se_estimates,current->nprobes+data->n_arrays,double);

      /* Optimized case for the RMA style model PM ~ -1 + samples + probes  with sum to zero constraint on probes */
      rlm_fit_anova(Y, current->nprobes, data->n_arrays, current->cur_params, current->cur_resids, current->cur_weights, PsiFunc(model->psi_code),model->psi_k,model->n_rlm_iterations,model->init_method);
      rlm_compute_se_anova(Y, current->nprobes, data->n_arrays, current->cur_params, current->cur_resids, current->cur_weights, current->cur_se_estimates,current->cur_varcov, current->cur_residSE, model->se_method, PsiFunc(model->psi_code),model->psi_k);
      current->cur_params[current->nprobes+data->n_arrays - 1] = 0.0;
      for (i = data->n_arrays; i <current->nprobes+data->n_arrays -1 ; i++)
	current->cur_params[current->nprobes+data->n_arrays -1]-=current->cur_params[i];

    } else {
      /* Least Squares for general models */
      rlm_fit(current->X,Y, current->n, current->p, current->cur_params, current->cur_resids, current->cur_weights, PsiFunc(model->psi_code),model->psi_k,model->n_rlm_iterations,model->init_method);  
      rlm_compute_se(current->X,Y, current->n, current->p, current->cur_params, current->cur_resids, current->cur_weights, current->cur_se_estimates,current->cur_varcov, current->cur_residSE, model->se_method, PsiFunc(model->psi_code),model->psi_k);
    }
    Free(Y);
  } else {
    /* weighted rlm's */
    input_weights= Calloc(current->n,double);
    for (i =0; i < current->n; i++){
      input_weights[i] = 1.0;
    }
    if (model->input_chip_weights != NULL){
      for (j = 0; j < data->n_arrays; j++){
	for (i =0; i < current->nprobes; i++){
	  input_weights[j*current->nprobes + i] = input_weights[j*current->nprobes +i]*model->input_chip_weights[j];
	}
      }
    }
    
    if (model->input_probe_weights != NULL){ 
      if (model->response_variable != 0){
	for (j = 0; j < data->n_arrays; j++){
	  for (i =0; i < current->nprobes; i++){
	    input_weights[j*current->nprobes + i] = input_weights[j*current->nprobes +i]*model->input_probe_weights[current_rows[i]];
	  }
	}
      } else {
	for (j = 0; j < data->n_arrays; j++){
	  for (i =0; i < current->nprobes; i++){
	    input_weights[j*current->nprobes + i] = input_weights[j*current->nprobes +i]*model->input_probe_weights[current_rows[i]];
	  }
	}
	for (j = 0; j < data->n_arrays; j++){
	  for (i =0; i < current->nprobes; i++){
	    input_weights[data->n_arrays*current->nprobes + j*current->nprobes + i] = input_weights[data->n_arrays*current->nprobes+j*current->nprobes +i]*model->input_probe_weights[data->n_probes+current_rows[i]];
	  }
	} 
	
      }
    }

    rlm_wfit(current->X,Y,input_weights, current->n, current->p, current->cur_params, current->cur_resids, current->cur_weights, PsiFunc(model->psi_code),model->psi_k,model->n_rlm_iterations,model->init_method);  
    rlm_compute_se(current->X,Y, current->n, current->p, current->cur_params, current->cur_resids, current->cur_weights, current->cur_se_estimates,current->cur_varcov, current->cur_residSE, model->se_method, PsiFunc(model->psi_code),model->psi_k);



    Free(input_weights);
  }
}


static void copy_PLM_estimates(PLM_modelfit *current, PLM_output *output,PLM_Datagroup *data, PLM_model_parameters *model, PLM_outputsettings *store, int first_probe,int which_probeset){

  int which_const_col=0;
  int which_param=0;
  int i,k,l;
  int offset1,offset2;
  int offset=0;

  /* first copy across all the stuff that is always there */
  
  /*Parameters and SE */

  /* Check for intercept */

  if (model->which_parameter_types[0]){
    output->out_constparams[which_const_col*data->n_probesets + which_probeset] = current->cur_params[which_param];
    output->out_const_SE[which_const_col*data->n_probesets + which_probeset] = current->cur_se_estimates[which_param];
    which_const_col++;
    which_param++;
  }

  /* Check for PM or MM type covariate */
  if (model->mmorpm_covariate != 0){
    output->out_constparams[which_const_col*data->n_probesets + which_probeset] = current->cur_params[which_param]; 
    output->out_const_SE[which_const_col*data->n_probesets + which_probeset] = current->cur_se_estimates[which_param];
    which_const_col++;
    which_param++;
  }

  /* Chip level factor/covariates */
  if (model->which_parameter_types[1]){
    offset = which_param;
    for (i = 0; i < model->n_chiplevelcovariates; i++){
      output->out_chipparams[i*data->n_probesets + which_probeset] = current->cur_params[which_param+i];
      output->out_chip_SE[i*data->n_probesets + which_probeset] = current->cur_se_estimates[which_param+i];
    }

    which_param+=model->n_chiplevelcovariates;
  }

  /* Sample effects */
  if (model->which_parameter_types[2]){
    offset = which_param;
    if (model->constraints[2] == 0){
      for (i = 0; i < data->n_arrays; i++){
	output->out_chipparams[i*data->n_probesets + which_probeset] = current->cur_params[which_param+i];
	output->out_chip_SE[i*data->n_probesets + which_probeset] = current->cur_se_estimates[which_param+i];
      }
      which_param+=data->n_arrays;
    } else {
      for (i = 0; i < data->n_arrays-1; i++){
	output->out_chipparams[i*data->n_probesets + which_probeset] = current->cur_params[which_param+i];
	output->out_chip_SE[i*data->n_probesets + which_probeset] = current->cur_se_estimates[which_param+i];
      }
      which_param+=data->n_arrays-1;

    }
  }

  /* Probe type effects */


  if (model->which_parameter_types[3]){
    if (model->constraints[3] == 0){
      if (model->strata[3]==0){
	for (i =0; i < 2; i++){
	  output->out_constparams[(which_const_col+i)*data->n_probesets + which_probeset] = current->cur_params[which_param+i]; 
	  output->out_const_SE[(which_const_col+i)*data->n_probesets + which_probeset] = current->cur_se_estimates[which_param+i]; 
	}
	which_param+=2;
	which_const_col +=2;
      } else if (model->strata[3]==1){
	for ( i=0; i < 2*model->n_arrays; i++){
	  output->out_constparams[(which_const_col + i)*data->n_probesets + which_probeset] = current->cur_params[which_param+i];  
	  output->out_const_SE[(which_const_col+i)*data->n_probesets + which_probeset] = current->cur_se_estimates[which_param+i]; 
	}
	which_param+=2*model->n_arrays;
	which_const_col +=2*model->n_arrays;
      } else if (model->strata[3]==2){
	for ( i=0; i < 2*(model->max_probe_type_treatment_factor+1); i++){
	  output->out_constparams[(which_const_col + i)*data->n_probesets + which_probeset] = current->cur_params[which_param+i];  
	  output->out_const_SE[(which_const_col+i)*data->n_probesets + which_probeset] = current->cur_se_estimates[which_param+i]; 
	}
	which_param+=2*(model->max_probe_type_treatment_factor+1);
	which_const_col +=2*(model->max_probe_type_treatment_factor+1);
      }
    } else {
      if (model->strata[3]==0){	
	for (i =0; i < 1; i++){
	  output->out_constparams[(which_const_col+i)*data->n_probesets + which_probeset] = current->cur_params[which_param+i]; 
	  output->out_const_SE[(which_const_col+i)*data->n_probesets + which_probeset] = current->cur_se_estimates[which_param+i]; 
	}
	which_param+=1;
	which_const_col +=1;
      } else if (model->strata[3]==1){
	for ( i=0; i < model->n_arrays; i++){
	  output->out_constparams[(which_const_col + i)*data->n_probesets + which_probeset] = current->cur_params[which_param+i];  
	  output->out_const_SE[(which_const_col+i)*data->n_probesets + which_probeset] = current->cur_se_estimates[which_param+i]; 
	}
	which_param+=model->n_arrays;
	which_const_col +=model->n_arrays;
      } else if (model->strata[3]==2){
	for ( i=0; i < (model->max_probe_type_treatment_factor+1); i++){
	  output->out_constparams[(which_const_col + i)*data->n_probesets + which_probeset] = current->cur_params[which_param+i];  
	  output->out_const_SE[(which_const_col+i)*data->n_probesets + which_probeset] = current->cur_se_estimates[which_param+i]; 
	}
	which_param+=(model->max_probe_type_treatment_factor+1);
	which_const_col +=(model->max_probe_type_treatment_factor+1);
      }
    }   
  }
  
  /* Probe effects */

  if (checkDefaultModel(model)){
    for (i=0; i < current->nprobes;i++){
      output->out_probeparams[which_probeset][i] = current->cur_params[which_param+i];
      output->out_probe_SE[which_probeset][i] = current->cur_se_estimates[which_param+i];
    }
  } else {
    if (model->which_parameter_types[4]){
      /* probe effect parameter figure figure out which strata and constraints */ 
      if (model->constraints[4] ==0){
	/* Unconstrainted */
	if (model->strata[4] == 0){
	  /* Overall probe effect */
	  for (i=0; i < current->nprobes;i++){
	    output->out_probeparams[which_probeset][i] = current->cur_params[which_param+i];
	    output->out_probe_SE[which_probeset][i] = current->cur_se_estimates[which_param+i];
	  }
	  which_param+=current->nprobes;
	} else if (model->strata[4] == 2){
	  /* Probe effect within the levels of a treatment/covariate factor */
	  for (i=0; i < current->nprobes*(model->max_probe_treatment_factor+1);i++){
	    output->out_probeparams[which_probeset][i] = current->cur_params[which_param+i];
	    output->out_probe_SE[which_probeset][i] = current->cur_se_estimates[which_param+i];
	  }
	  which_param+=(current->nprobes*(model->max_probe_treatment_factor+1));
	} else if (model->strata[4] == 3){
	  /* Probe effect within probe type */
	  for (i=0; i < 2*current->nprobes;i++){
	    output->out_probeparams[which_probeset][i] = current->cur_params[which_param+i];
	    output->out_probe_SE[which_probeset][i] = current->cur_se_estimates[which_param+i];
	  }
	  which_param+=(2*current->nprobes);
	} else if (model->strata[4] == 4){
	  /* Probe effect within probe type within the levels of a treatment/covariate factor*/
	  for (i=0; i < 2*current->nprobes*(model->max_probe_treatment_factor+1);i++){
	    output->out_probeparams[which_probeset][i] = current->cur_params[which_param+i];
	    output->out_probe_SE[which_probeset][i] = current->cur_se_estimates[which_param+i];
	  }
	  which_param+=(2*current->nprobes*(model->max_probe_treatment_factor+1));
	}
      } else {
	/* Constrained */
	if (model->strata[4] == 0){
	  /* Overall probe effect */
	  for (i=0; i < current->nprobes-1;i++){
	    output->out_probeparams[which_probeset][i] = current->cur_params[which_param+i];
	    output->out_probe_SE[which_probeset][i] = current->cur_se_estimates[which_param+i];
	  }
	  which_param+=current->nprobes-1;
	} else if (model->strata[4] == 2){
	  /* Probe effect within the levels of a treatment/covariate factor */
	  for (i=0; i < ((current->nprobes-1)*(model->max_probe_treatment_factor+1));i++){
	    output->out_probeparams[which_probeset][i] = current->cur_params[which_param+i];
	    output->out_probe_SE[which_probeset][i] = current->cur_se_estimates[which_param+i];
	  }
	  which_param+=((current->nprobes-1)*(model->max_probe_treatment_factor+1));
	} else if (model->strata[4] == 3){
	  /* Probe effect within probe type */
	  for (i=0; i < (2*(current->nprobes-1));i++){
	    output->out_probeparams[which_probeset][i] = current->cur_params[which_param+i];
	    output->out_probe_SE[which_probeset][i] = current->cur_se_estimates[which_param+i];
	  }
	  which_param+=(2*(current->nprobes-1));
	} else if (model->strata[4] == 4){
	  /* Probe effect within probe type within the levels of a treatment/covariate factor*/
	  for (i=0; i <(2*(current->nprobes-1)*(model->max_probe_treatment_factor+1));i++){
	    output->out_probeparams[which_probeset][i] = current->cur_params[which_param+i];
	    output->out_probe_SE[which_probeset][i] = current->cur_se_estimates[which_param+i];
	  }
	  which_param+=(2*(current->nprobes-1)*(model->max_probe_treatment_factor+1));
	}
      }
    }
  }

  /* now check for the optional stuff */
  /* Weights */
  if (store->weights){
    offset1 =0;
    if (model->response_variable ==0){
      offset1 =  current->nprobes; /*  2*current->nprobes; */
      offset2 = data->n_arrays*current->nprobes;
    } else {
      offset1 = current->nprobes;
      offset2 = 0;
    }
    if (model->response_variable >=0){
      for(k=0; k < data->n_arrays; k++){
	for (l=0; l < current->nprobes; l++){
	  output->out_weights[0][k*(data->n_probes) + (first_probe + l)] = current->cur_weights[k*(offset1) + l];
	}
      }
    }
    if (model->response_variable <=0){
      for(k=0; k < data->n_arrays; k++){
	for (l=0; l < current->nprobes; l++){
	  output->out_weights[1][k*(data->n_probes) + (first_probe + l)] = current->cur_weights[k*(offset1) + l+ offset2];
	}
      }
    }
  }

  /* Residuals */
  if (store->residuals){
    offset1 =0; 
    if (model->response_variable ==0){
      offset1 =   current->nprobes;  /* 2*current->nprobes; */
      offset2 = data->n_arrays*current->nprobes;
    } else {
      offset1 = current->nprobes;
      offset2 = 0;
    }
    if (model->response_variable >=0){
      for(k=0; k < data->n_arrays; k++){
	for (l=0; l < current->nprobes; l++){
	  output->out_resids[0][k*(data->n_probes) + (first_probe + l)] = current->cur_resids[k*(offset1) + l];
	}
      }
    }
    if (model->response_variable <=0){
      for(k=0; k < data->n_arrays; k++){
	for (l=0; l < current->nprobes; l++){
	  output->out_resids[1][k*(data->n_probes) + (first_probe + l)] = current->cur_resids[k*(offset1) + l+offset2];
	}
      }
    }
  }

  /* Residual SE */
  if (store->residSE){
    output->out_residSE[which_probeset] = current->cur_residSE[0];
    output->out_residSE[data->n_probesets+which_probeset] = current->n - current->p;
  }

  /* Covariance Matrix */
  if (store->varcov){
    /*  error("varcov option all not currently supported"); */
    if (store->varcov == 1){
      /* chip-level (note that this includes the intercept) */
      if (model->which_parameter_types[2]){
	/* sample effect model */
	if (model->which_parameter_types[0]){
	  /* need to include the intercept parameter */
	  output->out_varcov[which_probeset][0] = current->cur_varcov[0];
	  for (l=0; l < model->n_arrays-1; l++){
	    output->out_varcov[which_probeset][0*model->n_arrays + (l+1)] = current->cur_varcov[(l+offset)*current->p + 0];
	    output->out_varcov[which_probeset][(l+1)*model->n_arrays + 0] = current->cur_varcov[(l+offset)*current->p + 0];
	  }
	  for (k = 0; k < model->n_arrays-1; k++){
	    for (l = 0; l <= k ; l++){
	      output->out_varcov[which_probeset][(k+1)*model->n_arrays + (l+1)] = current->cur_varcov[(k+offset)*current->p + (l+offset)];
	      output->out_varcov[which_probeset][(l+1)*model->n_arrays + (k+1)] = output->out_varcov[which_probeset][(k+1)*model->n_arrays + (l+1)];
	    }
	  }
	} else {
	  /* no intercept */
	  for (k = 0; k < model->n_arrays; k++){
	    for (l = 0; l <= k ; l++){
	      output->out_varcov[which_probeset][k*model->n_arrays + l] = current->cur_varcov[(k+offset)*current->p + (l+offset)];
	      output->out_varcov[which_probeset][l*model->n_arrays + k] = output->out_varcov[which_probeset][k*model->n_arrays + l];
	    }
	  }
	}
      } else {
	/* treatment covariates model */
	if (model->which_parameter_types[0]){
	  /* need to include the intercept parameter */ 
	  output->out_varcov[which_probeset][0] = current->cur_varcov[0];
	  for (l=0; l < model->n_chiplevelcovariates; l++){
	    output->out_varcov[which_probeset][0*(model->n_chiplevelcovariates+1) + (l+1)] =current->cur_varcov[(l+offset)*current->p + 0];
	    output->out_varcov[which_probeset][(l+1)*(model->n_chiplevelcovariates+1) + 0] =current->cur_varcov[(l+offset)*current->p + 0];
	  }
	  for (k = 0; k < model->n_chiplevelcovariates; k++){
	    for (l = 0; l <= k ; l++){
	      output->out_varcov[which_probeset][(k+1)*(model->n_chiplevelcovariates+1) + (l+1)] = current->cur_varcov[(k+offset)*current->p + (l+offset)];
	      output->out_varcov[which_probeset][(l+1)*(model->n_chiplevelcovariates+1) + (k+1)] = output->out_varcov[which_probeset][(k+1)*(model->n_chiplevelcovariates+1) + (l+1)];
	    }
	  }
	} else {
	  for (k = 0; k < model->n_chiplevelcovariates; k++){
	    for (l = 0; l <= k ; l++){
	      output->out_varcov[which_probeset][k*model->n_chiplevelcovariates + l] = current->cur_varcov[(k+offset)*current->p + (l+offset)];
	      output->out_varcov[which_probeset][l*model->n_chiplevelcovariates + k] = output->out_varcov[which_probeset][k*model->n_chiplevelcovariates + l];
	    }
	  }
	}
      }
    } else if (store->varcov ==2){
      for (k = 0; k < current->p; k++){
	for (l = 0; l <= k ; l++){
	  output->out_varcov[which_probeset][k*current->p + l] = current->cur_varcov[k*current->p + l];
	  output->out_varcov[which_probeset][l*current->p + k] = output->out_varcov[which_probeset][k*current->p + l];
	}
      }
    }


  }


}







void do_PLM_rlm(PLM_Datagroup *data,  PLM_model_parameters *model, PLM_output *output, PLM_outputsettings *store){

  int i,j,k;
  int start;
  int new_nprobes=0;
  int size;
  int first_ind;
  int max_nrows = 1000;
  int *cur_rows= (int *)Calloc(max_nrows,int);
  const char *first;
  PLM_modelfit *current= new_PLMmodelfit();
  
  first = data->ProbeNames[0];
  first_ind = 0;
  i = 0;     /* indexes current probeset */
  j = 0;    /* indexes current row in PM matrix */
  k = 0;    /* indexes current probe in probeset */
  
  /*  new_nprobes=1;
      for (j = 1; j < data->n_probes; j++){
      if ((strcmp(first,data->ProbeNames[j]) != 0) | (j == (data->n_probes -1))){
      if (j == (data->n_probes -1)){
      new_nprobes++;
      for (k = 0; k < new_nprobes; k++){
      if (k >= max_nrows){
      max_nrows = 2*max_nrows;
      cur_rows = Realloc(cur_rows, max_nrows, int);
	  }
          cur_rows[k] = (j+1 - new_nprobes)+k;
	  }
	  start = j+1 - new_nprobes;
	  } else {
	  for (k = 0; k < new_nprobes; k++){
	  if (k >= max_nrows){
	  max_nrows = 2*max_nrows;
	  cur_rows = Realloc(cur_rows, max_nrows, int);
	  }
          cur_rows[k] = (j - new_nprobes)+k;
	}
	start = j - new_nprobes;
	}
	
	PLM_build_model_matrix(model, data, current, cur_rows, new_nprobes);
	rlm_PLM_probeset(model,data,current,cur_rows);
	//printf("%d %d %d\n",j,j-(new_nprobes),new_nprobes);
	copy_PLM_estimates(current, output, data, model, store, start, i);  // j-(new_nprobes),i);
	
	size = strlen(first);
	output->outnames[i] = Calloc(size+1,char);
	strcpy(output->outnames[i],first);  
	i++;
	first = data->ProbeNames[j];
	first_ind = j;
	new_nprobes = 0;
	}
	new_nprobes++;
	} */

  while ( j < data->n_probes){
    if (strcmp(first,data->ProbeNames[j]) == 0){
      if (k >= max_nrows){
	max_nrows = 2*max_nrows;
	cur_rows = Realloc(cur_rows, max_nrows, int);
      }
      cur_rows[k] = j;
      k++;
      j++;
      
    } else {
      new_nprobes = k;
      start = j - new_nprobes;
      PLM_build_model_matrix(model, data, current, cur_rows, new_nprobes);
      rlm_PLM_probeset(model,data,current,cur_rows);
      copy_PLM_estimates(current, output, data, model, store, start, i);  /* j-(new_nprobes),i); */
	
      size = strlen(first);
      output->outnames[i] = Calloc(size+1,char);
      strcpy(output->outnames[i],first);  
      i++;
      first = data->ProbeNames[j];


      k = 0;
    }
  }
  new_nprobes = k;
  start = j - new_nprobes;
  PLM_build_model_matrix(model, data, current, cur_rows, new_nprobes);
  rlm_PLM_probeset(model,data,current,cur_rows);
  copy_PLM_estimates(current, output, data, model, store, start, i);  /* j-(new_nprobes),i); */
  
  size = strlen(first);
  output->outnames[i] = Calloc(size+1,char);
  strcpy(output->outnames[i],first);  
  






  Free(cur_rows);
  free_PLMmodelfit(current);

}






