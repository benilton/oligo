/*********************************************************************
 **
 ** file: PLM_modelmatrix.c
 **
 ** Aim: build a model matrix
 **
 ** Copyright (C) 2004 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: July 6, 2004
 **
 ** Some assumptions:
 **      1. Observations are stacked in probes by array by probetype order
 **      2. the trt_cov variables take on values between 0 and max_trt_cov -1
 **
 ** History
 ** March 1, 2006 - change all comments to ansi style
 **
 *********************************************************************/

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include "PLM_modelmatrix.h"
#include "transfns.h"

/*********************************************************************
 **
 ** int PLM_matrix_intercept(double *X,
 **                          int n_arrays, 
 **                          int n_probes, 
 **                          int n_probetypes, 
 **                          int start_column)
 **
 ** creates a vector of 1's in an assigned column
 **
 *********************************************************************/


static int PLM_matrix_intercept(double *X,int n_arrays, int n_probes, int n_probetypes, int start_column){

  int i;
  int n_row=n_arrays*n_probes*n_probetypes;


  for (i=0;i < n_row; i++)
    X[start_column*n_row+ i] = 1.0;

  
  return 1;
}

/*********************************************************************
 **
 ** A wrapper to test the function from R using .C()
 **
 *********************************************************************/


void R_PLM_matrix_intercept(double *X,int *n_arrays, int *n_probes, int *n_probetypes, int *start_column){
  
  PLM_matrix_intercept(X,*n_arrays, *n_probes, *n_probetypes,*start_column);
}


static int PLM_matrix_MM(double *X,int n_arrays, int n_probes, int n_probetypes, int start_column, double *MM){

  int i;
  int n_row=n_arrays*n_probes*n_probetypes;

  for (i=0;i < n_row; i++)
    X[start_column*n_row+ i] = MM[i];

  
  return 1;

}



/*********************************************************************
 **
 ** A wrapper to test the function from R using .C()
 **
 *********************************************************************/

void R_PLM_matrix_MM(double *X,int *n_arrays, int *n_probes, int *n_probetypes, int *start_column, double *MM){
  
  PLM_matrix_MM(X,*n_arrays, *n_probes, *n_probetypes,*start_column, MM);
}

static int PLM_matrix_sample_effect(double *X,int n_arrays, int n_probes, int n_probetypes, int start_column, int constraint_type){

  int i,j,k, currow=0, curcol;
  int n_row=n_arrays*n_probes*n_probetypes;

  if (constraint_type == 0){
    for (k=0; k < n_probetypes; k++){
      for (i=0;i < n_arrays; i++){
	for (j=0; j < n_probes; j++){
	  X[(start_column+i)*n_row+ currow] = 1.0;
	  currow++;
	}
      }
    }
    return n_arrays;
  } else if (constraint_type == 1){
    for (k=0; k < n_probetypes; k++){
      for (i=0;i < n_arrays; i++){
	for (j=0; j < n_probes; j++){
	  if (i !=0){
	    X[(start_column+i - 1)*n_row+ currow] = 1.0;
	  }
	  currow++;
	}
      }
    }
    return n_arrays-1;
  } else if (constraint_type == -1){
     for (k=0; k < n_probetypes; k++){
      for (i=0;i < n_arrays; i++){
	for (j=0; j < n_probes; j++){
	  if (i != n_arrays -1){
	    X[(start_column+i)*n_row+ currow] = 1.0;
	  } else {
	    for (curcol=0; curcol < n_arrays-1; curcol++)
	      X[(start_column+curcol)*n_row+ currow] = -1.0;
	  }
	  currow++;
	}
      }
     }
     return n_arrays-1;
  }
  return 1;





}


/*********************************************************************
 **
 ** A wrapper to test the function from R using .C()
 **
 *********************************************************************/


void R_PLM_matrix_sample_effect(double *X,int *n_arrays, int *n_probes, int *n_probetypes, int *start_column, int *constraint_type){
  
  PLM_matrix_sample_effect(X,*n_arrays, *n_probes, *n_probetypes,*start_column, *constraint_type);
}


/**********************************************************************
 **
 ** static int PLM_matrix_probe_type_effect(double *X,
 **                                         int n_arrays, 
 **                                         int n_probes, 
 **                                         int n_probetypes, 
 **                                         int start_column, 
 **                                         int constraint_type, 
 **                                         int strata, 
 **                                         int *trt_cov, 
 **                                         int max_trt_cov)
 **
 **
 ** int strata - 0 - overall
 **            - 1 - within arrays
 **	       - 2 - which levels of a factor variable
 **
 **
 *********************************************************************/


static int PLM_matrix_probe_type_effect(double *X,int n_arrays, int n_probes, int n_probetypes, int start_column, int constraint_type, int strata, int *trt_cov, int max_trt_cov){

  int i,j,k, currow=0,curcol=0;
  int n_row=n_arrays*n_probes*n_probetypes;

  if (n_probetypes !=2){
    return 0;
  }
  if (strata == 0){
    if (constraint_type == 0){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    X[(start_column+k)*n_row+ currow] = 1.0;
	    currow++;
	  }
	}
      }
      return n_probetypes;
    } else if (constraint_type == 1){
      for (i=0;i < n_arrays; i++){
	for (j=0; j < n_probes; j++){
	  X[(start_column)*n_row+ currow + (n_arrays*n_probes)] = 1.0;
	  currow++;
	}
      }
      return 1;
    } else if (constraint_type == -1){
      for (i=0;i < n_arrays; i++){
	for (j=0; j < n_probes; j++){
	  X[(start_column)*n_row+ currow] = 1.0;
	  currow++;
	}
      }
      for (i=0;i < n_arrays; i++){
	for (j=0; j < n_probes; j++){
	  X[(start_column)*n_row+ currow] = -1.0;
	  currow++;
	}
      }
      

      return 1;
    }
  } else if (strata == 1){
      if (constraint_type == 0){
	for (k=0; k < n_probetypes; k++){
	  curcol = 0;
	  for (i=0;i < n_arrays; i++){
	    for (j=0; j < n_probes; j++){
	      X[(start_column+ k+2*curcol)*n_row+ currow] = 1.0;
	      currow++;
	    }
	    curcol++;
	   
	  }
	}
	return n_probetypes*n_arrays;
      } else if (constraint_type == 1){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    X[(start_column+curcol)*n_row+ currow + (n_arrays*n_probes)] = 1.0;
	    currow++;
	  }
	  curcol++;
	}
	return n_arrays;
      } else if (constraint_type == -1){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    X[(start_column+curcol)*n_row+ currow] = 1.0;
	    currow++;
	  }
	  curcol++;
	}
	curcol=0;
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    X[(start_column+curcol)*n_row+ currow] = -1.0;
	    currow++;
	  }
	  curcol++;
	}
	return n_arrays;
      }
  } else if (strata == 2){
    if (constraint_type == 0){
      for (i=0;i < n_arrays; i++){
	for (j=0; j < n_probes; j++){
	  X[(start_column+trt_cov[i]*2)*n_row+ currow] = 1.0;
	  currow++;
	}
      }
      for (i=0;i < n_arrays; i++){
	for (j=0; j < n_probes; j++){
	  X[(start_column+ trt_cov[i]*2 + 1)*n_row+ currow] = 1.0;
	  currow++;
	}
      }
      return n_probetypes*(max_trt_cov + 1);
    } else if (constraint_type == 1){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    X[(start_column+trt_cov[i])*n_row+ currow + (n_arrays*n_probes)] = 1.0;
	    currow++;
	  }
	}
	return (max_trt_cov + 1);
    } else if (constraint_type == -1){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    X[(start_column+trt_cov[i])*n_row+ currow] = 1.0;
	    currow++;
	  }
	}
	curcol=0;
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    X[(start_column+trt_cov[i])*n_row+ currow] = -1.0;
	    currow++;
	  }
	}
	return (max_trt_cov + 1);
    }
  }
  return 0;
}




void R_PLM_matrix_probe_type_effect(double *X,int *n_arrays, int *n_probes, int *n_probetypes, int *start_column, int *constraint_type, int *strata,int *trt_cov, int *max_trt_cov){
  
  PLM_matrix_probe_type_effect(X,*n_arrays, *n_probes, *n_probetypes,*start_column, *constraint_type, *strata,trt_cov, *max_trt_cov);
}

/*
 *
 *
 * int strata - 0  - overall
 *            - 2  - within the different levels of a treatment 
 *	      - 3  - within the two different  probe.types
 *	      - 4  - within probe-types  within treatment factor
 *
 *
 */


static int PLM_matrix_probe_effect(double *X,int n_arrays, int n_probes, int n_probetypes, int start_column,int constraint_type, int strata, int *trt_cov, int max_trt_cov){
  

  int i,j,k, currow=0,curcol=0;
  int n_row=n_arrays*n_probes*n_probetypes;
  
  if (strata == 0){
    if (constraint_type == 0){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    X[(start_column+j)*n_row+ currow] = 1.0;
	    currow++;
	  }
	}
      }
      return n_probes;
    } else if (constraint_type == 1){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    if (j !=0){
	      X[(start_column+(j-1))*n_row+ currow] = 1.0;
	    }
	    currow++;
	  }
	}
      }
      return n_probes -1;
    } else if (constraint_type ==-1){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    if (j != n_probes -1){
	      X[(start_column+ j)*n_row+ currow] = 1.0;
	    } else {
	      for (curcol=0; curcol < n_probes -1; curcol++){
		X[(start_column+ curcol)*n_row+ currow] = -1.0;
	      }
	    }
	    currow++;
	  }
	}
      }
      return n_probes -1;
    }
  } else if (strata == 2){
    if (constraint_type == 0){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    X[(start_column + j + trt_cov[i]*(n_probes))*n_row+ currow] = 1.0;
	    currow++;
	  }
	}
      }
      return (max_trt_cov+1)*n_probes;
    } else if (constraint_type== 1){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    if (j !=0){
	      X[(start_column+(j-1) + trt_cov[i]*(n_probes-1))*n_row+ currow] = 1.0;
	    }
	    currow++;
	  }
	}
      }
      return (max_trt_cov+1)*(n_probes -1);
    } else if (constraint_type== -1){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    if (j != n_probes -1){
	      X[(start_column + j + trt_cov[i]*(n_probes-1))*n_row+ currow] = 1.0;
	    } else {
	      for (curcol=0; curcol < n_probes -1; curcol++){
		X[(start_column+ curcol+ trt_cov[i]*(n_probes-1))*n_row+ currow] = -1.0;
	      }
	    }
	    currow++;
	  }
	}
      }
      return (max_trt_cov+1)*(n_probes -1);

    }

  } else if (strata == 3){
    /* probes within probe-types */
    if (constraint_type == 0){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    X[(start_column+j + k*n_probes)*n_row+ currow] = 1.0;
	    currow++;
	  }
	}
      }
      return n_probetypes*n_probes;
    } else if (constraint_type == 1){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    if (j !=0){
	      X[(start_column+(j-1) + k*(n_probes-1))*n_row+ currow] = 1.0;
	    }
	    currow++;
	  }
	}
      }
      return n_probetypes*(n_probes -1);


    } else if (constraint_type == -1){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    if (j != n_probes -1){
	      X[(start_column+ j + k*(n_probes-1))*n_row+ currow] = 1.0;
	    } else {
	      for (curcol=0; curcol < n_probes -1; curcol++){
		X[(start_column+ curcol+ k*(n_probes-1))*n_row+ currow] = -1.0;
	      }
	    }
	    currow++;
	  }
	}
      }
      return n_probetypes*(n_probes -1);
    }
  } else if (strata == 4){
    /* probes within probe-types and within treatment factor*/
    if (constraint_type == 0){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    X[(start_column+j +  trt_cov[i]*n_probes + k*(max_trt_cov+1)*n_probes)*n_row+ currow] = 1.0;
	    currow++;
	  }
	}
      }
      return (max_trt_cov+1)*n_probetypes*n_probes;

    } else if (constraint_type == 1){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    if (j !=0){
	      X[(start_column+(j-1) +trt_cov[i]*(n_probes-1) + k*(max_trt_cov+1)*(n_probes-1))*n_row+ currow] = 1.0;
	    }
	    currow++;
	  }
	}
      }
      return  (max_trt_cov+1)*n_probetypes*(n_probes -1);
    }  else if (constraint_type == -1){
      for (k=0; k < n_probetypes; k++){
	for (i=0;i < n_arrays; i++){
	  for (j=0; j < n_probes; j++){
	    if (j != n_probes -1){
	      X[(start_column+ j+trt_cov[i]*(n_probes-1) + k*(max_trt_cov+1)*(n_probes-1))*n_row+ currow] = 1.0;
	    } else {
	      for (curcol=0; curcol < n_probes -1; curcol++){
		X[(start_column+ curcol + trt_cov[i]*(n_probes-1) + k*(max_trt_cov+1)*(n_probes-1))*n_row+ currow] = -1.0;
	      }
	    }
	    currow++;
	  }
	}
      }
      return (max_trt_cov+1)*n_probetypes*(n_probes -1);
    }
  }

  return 0;
}



/*********************************************************************
 **
 ** A wrapper to test the function from R using .C()
 **
 *********************************************************************/

void R_PLM_matrix_probe_effect(double *X,int *n_arrays, int *n_probes, int *n_probetypes, int *start_column, int *constraint_type, int *strata,int *trt_cov, int *max_trt_cov){
  
  PLM_matrix_probe_effect(X,*n_arrays, *n_probes, *n_probetypes,*start_column, *constraint_type, *strata,trt_cov, *max_trt_cov);
}




static int PLM_matrix_chiplevel(double *X,int n_arrays, int n_probes, int n_probetypes, int start_column,double *covariates, int n_covariates){


  int i,j,k, currow=0,curcol=0;
  int n_row=n_arrays*n_probes*n_probetypes;
  
  for (k=0; k < n_probetypes; k++){
    for (i=0;i < n_arrays; i++){
      for (j=0; j < n_probes; j++){
	for (curcol=0; curcol < n_covariates; curcol++)
	  X[(start_column+curcol)*n_row+currow] = covariates[curcol*n_arrays + i];
	currow++;
      }
    }
  }
  return n_covariates;

}



/*********************************************************************
 **
 ** A wrapper to test the function from R using .C()
 **
 *********************************************************************/

void R_PLM_matrix_chiplevel(double *X,int *n_arrays, int *n_probes, int *n_probetypes, int *start_column,double *covariates,int *n_covariates){
  
  PLM_matrix_chiplevel(X,*n_arrays, *n_probes, *n_probetypes,*start_column, covariates,*n_covariates);
}


/*********************************************************************
 **
 ** A wrapper to test the function from R using .C()
 **
 *********************************************************************/


void R_PLM_Matrix_constructtest(double *X, int *n_arrays, int *n_probes, int *n_probetypes, int *hasintercept,int *hassampleeffects, int *hasprobetypes, int *hasprobeeffects,int *constrainttype){

  int curcol=0;
  int sumtozero=*constrainttype;
  int noconstraint=0;

  if (*hasintercept){
    curcol+=PLM_matrix_intercept(X,*n_arrays, *n_probes, *n_probetypes, curcol);
  }

  if (*hassampleeffects){
    if (*hasintercept){
      curcol+= PLM_matrix_sample_effect(X,*n_arrays,*n_probes,*n_probetypes,curcol,sumtozero);
    } else {
      curcol+= PLM_matrix_sample_effect(X,*n_arrays,*n_probes,*n_probetypes,curcol,noconstraint);
    }
  }
  
  if (*hasprobetypes){
    if (*hasintercept|| *hassampleeffects){
      curcol+= PLM_matrix_probe_type_effect(X,*n_arrays, *n_probes, *n_probetypes, curcol,sumtozero, 0, 0, 0);
    } else {
      curcol+= PLM_matrix_probe_type_effect(X,*n_arrays, *n_probes, *n_probetypes, curcol,noconstraint, 0, 0, 0);
    }
  }

  
  if (*hasprobeeffects){
    if (*hasintercept || *hassampleeffects || *hasprobetypes){
      curcol+=PLM_matrix_probe_effect(X,*n_arrays, *n_probes, *n_probetypes, curcol,sumtozero, 0, 0, 0);
    } else {
      curcol+=PLM_matrix_probe_effect(X,*n_arrays, *n_probes, *n_probetypes, curcol,noconstraint, 0, 0, 0);
    }
  }
}











void PLM_current_model_update_space(PLM_modelfit *current, int new_nprobes,int n,int p){
  int i;
  current->X = Realloc(current->X,n*p, double);
  for (i=0; i<n*p; i++){
    current->X[i] = 0.0;
  }
  current->cur_params = Realloc(current->cur_params,p,double);
  current->cur_se_estimates = Realloc(current->cur_se_estimates,p,double);
  current->cur_weights = Realloc(current->cur_weights,n,double);
  current->cur_resids = Realloc(current->cur_resids,n,double);
  current->cur_varcov = Realloc(current->cur_varcov,p*p,double);
  current->cur_residSE= Realloc(current->cur_residSE,2,double);
  current->n=n;
  current->p=p;
  current->nprobes = new_nprobes;
}


static int PLM_compute_n_params(const PLM_model_parameters *model, int new_nprobes){
  
  int p=0;

  /* check for intercept */
  if (model->which_parameter_types[0]){
    p++;
  }

  if (model->mmorpm_covariate !=0){
    p++;

  }
  if (model->which_parameter_types[1]){
    p+=model->n_chiplevelcovariates;
  }
  if (model->which_parameter_types[2]){
    if (model->constraints[2] == 0){
      p+=model->n_arrays;
    } else {
      p+=model->n_arrays-1;
    }
  }
  if (model->which_parameter_types[3]){
    /* probe.types parameter figure out which strata and constraints */
    if (model->constraints[3] == 0){
      if (model->strata[3]==0){
	p+=2;
      } else if (model->strata[3]==1){
	p+=2*model->n_arrays;
      } else if (model->strata[3]==2){
	p+=2*(model->max_probe_type_treatment_factor+1);
      }
    } else {
      if (model->strata[3]==0){
	p+=1;
      } else if (model->strata[3]==1){
	p+=model->n_arrays;
      } else if (model->strata[3]==2){
	p+=(model->max_probe_type_treatment_factor+1);
      }
    }
  }
  if (model->which_parameter_types[4]){
    /* probe effect parameter figure figure out which strata and constraints */ 
    if (model->constraints[4] ==0){
      if (model->strata[4] == 0){
	p+=new_nprobes;
      } else if (model->strata[4] == 2){
	p+=(new_nprobes*(model->max_probe_treatment_factor+1));
      } else if (model->strata[4] == 3){
	p+=(2*new_nprobes);
      } else if (model->strata[4] == 4){
	p+=(2*new_nprobes*(model->max_probe_treatment_factor+1));
      }
    } else {
      if (model->strata[4] == 0){
	p+=new_nprobes-1;
      } else if (model->strata[4] == 2){
	p+=((new_nprobes-1)*(model->max_probe_treatment_factor+1));
      } else if (model->strata[4] == 3){
	p+=(2*(new_nprobes-1));
      } else if (model->strata[4] == 4){
	p+=(2*(new_nprobes-1)*(model->max_probe_treatment_factor+1));
      }
    }
  }
  return p;
  
}


/**********************************************************************
 ** 
 ** This function will build the model matrix allocating whatever memory
 ** is needed and the space needed for storage of results
 **
 **
 **********************************************************************/

void PLM_build_model_matrix(const PLM_model_parameters *model, const PLM_Datagroup *data, PLM_modelfit *current, int *current_rows, int new_nprobes){
  
  int n,p;
  int i,j;
  int n_probetypes;
  int curcol=0;
  double *MMcovariates;
  pt2trans transfn = transFunc(model->trans_fn);

  /* check to see whether we need to do any reallocation */
  if ((model->mmorpm_covariate == 0) && (new_nprobes == current->nprobes)){
    /* Not a model with a PM/MM covariate and the number of probes in the probeset
    ** did not change so no need to reallocate anything 
    */
    return;
  }
  

  /* figure out how many observations */

  if (model->response_variable==0){
    /* both PM and MM are being used as response variables */
    n = 2*new_nprobes*data->n_arrays;
    n_probetypes = 2;
  } else {
    n = new_nprobes*data->n_arrays;
    n_probetypes = 1;
  }


  /* check to see if all we need to do is copy across the current PM or MM's */
  if((model->mmorpm_covariate !=0)  && (new_nprobes == current->nprobes)){
    /* there is a covariate MM or PM but the number of probes did not change */
    if (model->which_parameter_types[0]){
      /* intercept so second column contains the covariate */
      if (model->response_variable > 0){
	/* PM response MM covariate */
	for (i=0; i< data->n_arrays; i++){
	  for (j=0; j < new_nprobes; j++)
	    current->X[n + i*new_nprobes + j] = transfn(data->MM[i*data->n_probes + current_rows[j]]);         /* log(data->MM[i*data->n_probes + current_rows[j]])/log(2.0); */
	}
      } else {

	/* MM response PM covariate */
	for (i=0; i< data->n_arrays; i++){
	  for (j=0; j < new_nprobes; j++)
	    current->X[n + i*new_nprobes + j] = transfn(data->PM[i*data->n_probes + current_rows[j]]);            /* log(data->PM[i*data->n_probes + current_rows[j]])/log(2.0); */
	}
      }
    } else {
      /* no intercept so covariate goes in first column */
     if (model->response_variable > 0){
	for (i=0; i< data->n_arrays; i++){
	  for (j=0; j < new_nprobes; j++)
	    current->X[i*new_nprobes + j] = transfn(data->MM[i*data->n_probes + current_rows[j]]); /* log(data->MM[i*data->n_probes + current_rows[j]])/log(2.0);; */
	}
      } else {
	for (i=0; i< data->n_arrays; i++){
	  for (j=0; j < new_nprobes; j++)
	    current->X[i*new_nprobes + j] = transfn(data->PM[i*data->n_probes + current_rows[j]]);  /*  log(data->PM[i*data->n_probes + current_rows[j]])/log(2.0);; */
	}
      }


    }
    return;
  }
  /* from here on it is assumed that everything needs to be reallocated */
  

  
  /* figure out how many parameters are in current model */
  p = PLM_compute_n_params(model, new_nprobes);

  /* allocate the storage space */
  PLM_current_model_update_space(current, new_nprobes, n, p);
  
  /* start building the model matrix */
  

  /* check for intercept parameter */
  if (model->which_parameter_types[0]){
    curcol+=PLM_matrix_intercept(current->X,data->n_arrays, current->nprobes, n_probetypes, curcol);
  }  
  if (model->mmorpm_covariate !=0){
    MMcovariates = Calloc(n,double);
    if (model->response_variable < 0){
      /* MM response PM covariate */
      for (i=0; i< data->n_arrays; i++){
	for (j=0; j < new_nprobes; j++)
	  MMcovariates[i*new_nprobes + j] = log(data->PM[i*data->n_probes + current_rows[j]])/log(2.0);
      }
    } else {
      /* PM response MMcovariate */
      for (i=0; i< data->n_arrays; i++){
	for (j=0; j < new_nprobes; j++)
	  MMcovariates[i*new_nprobes + j] = log(data->MM[i*data->n_probes + current_rows[j]])/log(2.0);
      }
    }
    curcol+=PLM_matrix_MM(current->X,data->n_arrays, current->nprobes,n_probetypes, curcol,MMcovariates);
    Free(MMcovariates);
  }
  if (model->which_parameter_types[1]){
    curcol+=PLM_matrix_chiplevel(current->X,data->n_arrays, current->nprobes,n_probetypes, curcol,model->chiplevelcovariates,model->n_chiplevelcovariates);
  }
  if (model->which_parameter_types[2]){
    curcol+= PLM_matrix_sample_effect(current->X,data->n_arrays, current->nprobes,n_probetypes, curcol, model->constraints[2]);
  }
  if (model->which_parameter_types[3]){
    curcol+=PLM_matrix_probe_type_effect(current->X,data->n_arrays, current->nprobes,n_probetypes,curcol,model->constraints[3], model->strata[3], model->probe_type_treatment_factor, model->max_probe_type_treatment_factor);
  }
  if (model->which_parameter_types[4]){
    curcol+=PLM_matrix_probe_effect(current->X,data->n_arrays, current->nprobes,n_probetypes,curcol,model->constraints[4], model->strata[4], model->probe_treatment_factor, model->max_probe_treatment_factor);
  }

}
