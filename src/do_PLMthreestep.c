/*********************************************************************
 **
 ** file: do_PLMthreestep.c
 **
 ** Aim: do a threestep summary as a PLMset object.
 **
 ** Copyright (C) 2003-2005 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Oct 9, 2003
 **
 ** History
 ** Oct 9, 2003 - Initial version
 ** Oct 10, 2003 - Introduce a general mechanism for threestep PLM
 ** Oct 12, 2003 - Fix declaration order problem
 ** Mar 13, 2005 - change loop
 ** Sep 18, 2005 - fix a malloc/Free pair that was causing seg faults on windows
 ** Mar 1, 2006 - change all comments to ansi style
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
#include "threestep_summary_methods_param.h"



#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void threestep_PLM_block(Datagroup *data, PLMmodelparam *model, modelfit *current){

  summary_plist summary_param;
  
  summary_param.psi_method = model->psi_code;
  summary_param.psi_k = model->psi_k;

  model->PLM3stepSummary(data->PM, data->rows, data->cols, current->cur_rows, current->cur_params, current->nprobes,current->cur_se_estimates,current->cur_resids,&summary_param);

}



void copy_threestepPLM_results(modelfit *current, PLMoutput *output, Datagroup *data,const PLMmodelparam *model, const outputsettings *store, int j, int i){

  int k,l;
  
  for (k = 0; k < data->cols; k++){
    output->out_chipparams[k*data->nprobesets +i] = current->cur_params[k];
  }
    
  for (k = 0; k < data->cols; k++){
    output->out_chip_SE[k*data->nprobesets +i] =  current->cur_se_estimates[k];
  }
  
  
  /* copy the weights and residuals into output */
  /* note that we use the values in "store"
     to determine whether to save what has been returned
     for everything that follows                       */
    
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
}












void do_PLMthreestep(Datagroup *data,  PLMmodelparam *model, PLMoutput *output, outputsettings *store){
  int i = 0,j=0,k=0;
  int size;
  const char *first;
  int first_ind;
  int max_nrows = 1000;
  int old_nprobes =0;
  
  /* buffers of size 200 should be enough. */

  modelfit *current = Calloc(1,modelfit);

  current->cur_rows=Calloc(max_nrows,int);
  current->cur_weights = 0; /* weights are not returned by threestep routines */
  current->cur_params = Calloc(data->cols,double);
  current->cur_se_estimates = Calloc(data->cols,double);
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
  /*current->nprobes = 1;
  for (j = 1; j < data->rows; j++){
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

      // Check last number of probes and only Realloc when needed  
      if (old_nprobes != current->nprobes){
	current->n = current->nprobes*(data->cols);
	current->cur_resids = Realloc(current->cur_resids,data->cols*current->nprobes,double);
	old_nprobes = current->nprobes;
      }
      
      current->cur_resids = Realloc(current->cur_resids,data->cols*current->nprobes,double);


      

      threestep_PLM_block(data, model, current);

      copy_threestepPLM_results(current, output, data, model, store, j,i);
      
    
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
      
    } else{
      if (old_nprobes != current->nprobes){
	current->n = current->nprobes*(data->cols);
	current->cur_resids = Realloc(current->cur_resids,data->cols*current->nprobes,double);
	old_nprobes = current->nprobes;
      }
      
      /* current->cur_resids = Realloc(current->cur_resids,data->cols*current->nprobes,double); */


      

      threestep_PLM_block(data, model, current);

      copy_threestepPLM_results(current, output, data, model, store, j,i);
      
    
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
    current->cur_resids = Realloc(current->cur_resids,data->cols*current->nprobes,double);
    old_nprobes = current->nprobes;
  }
  
/* current->cur_resids = Realloc(current->cur_resids,data->cols*current->nprobes,double); */
  
  threestep_PLM_block(data, model, current);
  copy_threestepPLM_results(current, output, data, model, store, j,i);
  size = strlen(first);
  output->outnames[i] = Calloc(size+1,char);
  strcpy(output->outnames[i],first);  
  



/*  Free(current->X);
**  Free(current->cur_varcov); */
  Free(current->cur_resids);
  Free(current->cur_se_estimates);
  Free(current->cur_params);
/* Free(current->cur_weights); */
  Free(current->cur_rows);
  Free(current);
}
