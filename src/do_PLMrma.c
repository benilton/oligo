/*********************************************************************
 **
 ** file: do_PLMrlm.c
 **
 ** Aim: fit rma model as a PLMset object.
 **
 ** Copyright (C) 2003-2005 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Sept 14, 2003
 **
 ** History
 ** Sept 14, 2003 - Initial version
 ** Oct 5, 2003 - add missing #include
 ** Oct 12, 2003 - fix declaration order error
 ** May 11, 2004 - fix a memory leak
 ** Mar 13, 2005 - change the main loop
 ** Mar 1, 2006 - change comments to ansi style comments
 **
 *********************************************************************/

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

#include "PLM_medianpolish.h"
#include "common_types.h"

#include "rmaPLM_pseudo.h"



#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



void rma_PLM_block(const Datagroup *data, const PLMmodelparam *model, modelfit *current){
  
  int i;

  double *probeparam;
  double *chipparam;
  double constparam;


  /* median polish */
  probeparam = Calloc(current->nprobes,double);
  chipparam = Calloc(data->cols,double);
  median_polishPLM(data->PM,data->rows, data->cols, current->cur_rows, probeparam, chipparam, &constparam, current->nprobes, current->cur_resids); 
  for (i =0; i < (current->nprobes); i++){
    current->cur_params[i] = probeparam[i];
  }
  for (i = 0; i < data->cols; i++){
    current->cur_params[i+(current->nprobes)] = chipparam[i];
  }
  current->cur_params[(current->nprobes)+(data->cols)] = constparam;
  
  Free(probeparam);
  Free(chipparam);
}

/*********************************************************************
 **
 ** void copy_PLM_results(modelfit *current, PLMoutput *output, Datagroup *data,const PLMmodelparam *model, int j, int i)
 **
 ** This function should copy the results of the current fit into 
 ** the appropiate output areas.
 **
 ********************************************************************/

void copy_rmaPLM_results(modelfit *current, PLMoutput *output, Datagroup *data,const PLMmodelparam *model, const outputsettings *store, int j, int i){

  int k,l;
  /*  int offset; */

  /* depending on what model was fit copy the parameter estimates and standard errors to the appropriate places */
  /* Parameter estimates */

  if (j == (data->rows -1)){
    for (k = 0; k < current->nprobes; k++){
      output->out_probeparams[j+1  - (current->nprobes) + k] =  current->cur_params[k];
    }
  } else {
    for (k = 0; k < current->nprobes; k++){
      output->out_probeparams[j  - (current->nprobes) + k] =  current->cur_params[k];
    }
  }

  for (k = 0; k < data->cols; k++){
    output->out_chipparams[k*data->nprobesets +i] = current->cur_params[k+ current->nprobes] +  current->cur_params[data->cols +current->nprobes ];
  }


  /* Standard errors */

  if (!store->pseudoSE){
    if (j == (data->rows -1)){
      for (k = 0; k < current->nprobes; k++){
	output->out_probe_SE[j+1  - (current->nprobes) + k] = R_NaN;
      }
    } else {
      for (k = 0; k < current->nprobes; k++){
	output->out_probe_SE[j  - (current->nprobes) + k] = R_NaN;
      }
    }
    
    for (k = 0; k < data->cols; k++){
      output->out_chip_SE[k*data->nprobesets +i] =  R_NaN;
    }
  } else {
    /* we will compute pseudoSE based upon the weighting system specified */
    
    compute_pseudoSE(current->cur_resids, current->cur_se_estimates,current->nprobes, data->cols,model->psi_code,model->psi_k);
    if (j == (data->rows -1)){
      for (k = 0; k < current->nprobes; k++){
	output->out_probe_SE[j+1  - (current->nprobes) + k] = current->cur_se_estimates[k];
      }
    } else {
      for (k = 0; k < current->nprobes; k++){
	output->out_probe_SE[j  - (current->nprobes) + k] = current->cur_se_estimates[k];
      }
    }
    
    for (k = 0; k < data->cols; k++){
      output->out_chip_SE[k*data->nprobesets +i] =  current->cur_se_estimates[k + (current->nprobes)];
    }
  }
  
  /* copy the weights and residuals into output */
  /* note that we use the values in "store"
     to determine whether to save what has been returned
     for everything that follows                       */
    
  if (store->weights){

    compute_pseudoweights(current->cur_resids,current->cur_weights,current->nprobes, data->cols,model->psi_code,model->psi_k);



    if (j == (data->rows -1)){
      for(k=0; k < data->cols; k++){
	for (l=0; l < current->nprobes; l++){
	  output->out_weights[k*(data->rows) + (j+1 - (current->nprobes) + l)] = current->cur_weights[k*(current->nprobes) + l];
	}
      }
    } else {
      for(k=0; k < data->cols; k++){
	for (l=0; l < current->nprobes; l++){
	  output->out_weights[k*(data->rows) + (j - (current->nprobes) + l)] = current->cur_weights[k*(current->nprobes) + l];
	  /* printf("%d ",(j - (current->nprobes) + l));
	  **  printf("% f",current->cur_weights[k*(current->nprobes) + l]); */
	}
	/* printf("\n"); */
      }
    }
    
  }

  if (store->residuals){
    if (j == (data->rows -1)){
      for(k=0; k < data->cols; k++){
	for (l=0; l < current->nprobes; l++){
	  output->out_resids[k*(data->rows) + (j+1 - (current->nprobes) + l)] = current->cur_resids[k*(current->nprobes) + l];
	}
      }
    } else {
      for(k=0; k < data->cols; k++){
	for (l=0; l < current->nprobes; l++){
	  output->out_resids[k*(data->rows) + (j - (current->nprobes) + l)] = current->cur_resids[k*(current->nprobes) + l];
	}
      }
    }
    
  }

  if (store->residSE){
    output->out_residSE[i] = current->cur_residSE[0];
    output->out_residSE[data->nprobesets+i] = current->n - current->p;
  }
}





void do_PLMrma(Datagroup *data,  PLMmodelparam *model, PLMoutput *output, outputsettings *store){
  int i = 0,j=0,k=0;
  int size;
  const char *first;
  int first_ind;
  int max_nrows = 1000;
  int old_nprobes =0;
  
  /* buffers of size 200 should be enough. */

  modelfit *current = (modelfit *)Calloc(1,modelfit);

  current->cur_rows=Calloc(max_nrows,int);
  current->cur_weights = Calloc(data->cols,double);
  current->cur_params = Calloc(data->cols+100,double);
  current->cur_se_estimates = Calloc(data->cols+100,double);
  current->cur_resids = Calloc(data->cols,double);
  current->p = 0;
  current->nprobes = 0;
  current->n = 0;
  current->cur_residSE = 0; /* Calloc(2,double); */
  current->cur_varcov = 0; /* Calloc(4,double); */
  current->X = 0;
  

  first = data->ProbeNames[0];
  first_ind = 0;
  i =0;
  current->nprobes = 1;
  /*  for (j = 1; j < data->rows; j++){
    if ((strcmp(first,data->ProbeNames[j]) != 0) | (j == (data->rows -1))){
      if (j == (data->rows -1)){
        current->nprobes++;
        for (k = 0; k < current->nprobes; k++){
	  if (k >= max_nrows){
	    max_nrows = 2*max_nrows;
	    current->cur_rows = Realloc(current->cur_rows, max_nrows, int);
	  }
          current->cur_rows[k] = (j+1 - current->nprobes)+k;
        }
      } else {
        for (k = 0; k < current->nprobes; k++){
	  if (k >= max_nrows){
	    max_nrows = 2*max_nrows;
	    current->cur_rows = Realloc(current->cur_rows, max_nrows, int);
	  }
          current->cur_rows[k] = (j - current->nprobes)+k;
	}
      }
      if (old_nprobes != current->nprobes){
	current->n = current->nprobes*(data->cols);
	current->p = current->nprobes + (data->cols) + 1;
	current->cur_weights = Realloc(current->cur_weights,current->n,double);
	current->cur_resids = Realloc(current->cur_resids,current->n,double);
	current->cur_params  = Realloc(current->cur_params,current->p,double);
	current->cur_se_estimates  = Realloc(current->cur_se_estimates,current->p,double);
	// current->cur_varcov = Realloc(current->cur_varcov,current->p*current->p, double);
	// current->X = Realloc(current->X,current->n*current->p,double);
	// rlm_design_matrix_realloc(current->X, current->nprobes, data->cols, current->p, model->input_chipcovariates, model->method);
	//
	old_nprobes = current->nprobes;
      }


      rma_PLM_block(data, model, current);

      copy_rmaPLM_results(current, output, data, model, store, j,i);
      
    
      size = strlen(first);
      output->outnames[i] = Calloc(size+1,char);
      strcpy(output->outnames[i],first);  
      i++;
      first = data->ProbeNames[j];
      first_ind = j;
      current->nprobes = 0;
    }
    current->nprobes++;
  } */

 current->nprobes = 0;
  i = 0;     /* indexes current probeset */
  j = 0;    /* indexes current row in PM matrix */
  k = 0;    /* indexes current probe in probeset */
   while ( j < data->rows){
    if (strcmp(first,data->ProbeNames[j]) == 0){
      if (k >= max_nrows){
	max_nrows = 2*max_nrows;
	current->cur_rows = Realloc(current->cur_rows, max_nrows, int);
      }
      current->cur_rows[k] = j;
      k++;
      j++;
      current->nprobes++;
      
    } else {
      if (old_nprobes != current->nprobes){
	current->n = current->nprobes*(data->cols);
	current->p = current->nprobes + (data->cols) + 1;
	current->cur_weights = Realloc(current->cur_weights,current->n,double);
	current->cur_resids = Realloc(current->cur_resids,current->n,double);
	current->cur_params  = Realloc(current->cur_params,current->p,double);
	current->cur_se_estimates  = Realloc(current->cur_se_estimates,current->p,double);
	/* current->cur_varcov = Realloc(current->cur_varcov,current->p*current->p, double);
	**current->X = Realloc(current->X,current->n*current->p,double);
	**rlm_design_matrix_realloc(current->X, current->nprobes, data->cols, current->p, model->input_chipcovariates, model->method);
	*/
	old_nprobes = current->nprobes;
      }
      rma_PLM_block(data, model, current);
      
      copy_rmaPLM_results(current, output, data, model, store, j,i);
      
      
      size = strlen(first);
      output->outnames[i] = Calloc(size+1,char);
      strcpy(output->outnames[i],first);  
      i++;
      first = data->ProbeNames[j];
      first_ind = j;
      current->nprobes = 0;
      k = 0;
    }
   }
   j--;
   
   if (old_nprobes != current->nprobes){
     current->n = current->nprobes*(data->cols);
     current->p = current->nprobes + (data->cols) + 1;
     current->cur_weights = Realloc(current->cur_weights,current->n,double);
     current->cur_resids = Realloc(current->cur_resids,current->n,double);
     current->cur_params  = Realloc(current->cur_params,current->p,double);
     current->cur_se_estimates  = Realloc(current->cur_se_estimates,current->p,double);
     /*current->cur_varcov = Realloc(current->cur_varcov,current->p*current->p, double);
     **current->X = Realloc(current->X,current->n*current->p,double);
     **rlm_design_matrix_realloc(current->X, current->nprobes, data->cols, current->p, model->input_chipcovariates, model->method);
     */
     old_nprobes = current->nprobes;
   } 
   rma_PLM_block(data, model, current);
   
   copy_rmaPLM_results(current, output, data, model, store, j,i); 
   
   size = strlen(first);
   output->outnames[i] = Calloc(size+1,char);
   strcpy(output->outnames[i],first);  
   i++;
   /*first = data->ProbeNames[j]; */
   






/*  Free(current->X); */
/* Free(current->cur_varcov); */
  Free(current->cur_resids);
  Free(current->cur_se_estimates);
  Free(current->cur_params);
  Free(current->cur_weights);
  Free(current->cur_rows);
  Free(current);
}
