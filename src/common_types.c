#include "common_types.h"

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

PLM_modelfit *new_PLMmodelfit(){
  
  PLM_modelfit *x = Calloc(1,PLM_modelfit);
  x->cur_params = Calloc(1,double);
  x->cur_weights= Calloc(1,double);
  x->cur_resids= Calloc(1,double);
  x->cur_varcov= Calloc(1,double);
  x->cur_residSE= Calloc(1,double);
  x->X= Calloc(1,double);   
  x->n =0;
  x->p=0;
  x->nprobes=0;
    


  return x;

}

void free_PLMmodelfit(PLM_modelfit *x){
  
  Free(x->cur_params);
  Free(x->cur_weights);
  Free(x->cur_resids);
  Free(x->cur_varcov);
  Free(x->cur_residSE);
  Free(x->cur_se_estimates);
  Free(x->X);
  Free(x);
}



