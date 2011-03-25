/******************************************************************
 **
 ** file: SCAB.c
 **
 ** Copyright (C) 2003 B M Bolstad
 **
 ** aim: implement parts of the SCAB (Standardized Curve Adjusts
 **      Background) method for adjusting signal intensities
 **
 ** Created on: Apr 11, 2003
 ** 
 ** Created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 **
 ** Last modified by: B. M. Bolstad
 **
 ** History
 ** Apr 11, 2003 - Initial version
 ** Jun 4, 2003 - add new parameter to rlm_fit call
 ** Jul 22, 2003 - can select various M estimators, parameters
 **                added to function calls.
 ** Jul 23, 2003 - removed a variable that was declared but not
 **                used. This removed another compiler warning.
 ** Jul 27, 2003 - add ability to fit model without probe pair effects
 ** Jul 28, 2003 - robust average distance between PM and MM pairs
 **                see finct fit_Difference_model
 ** Apr 5, 2004 - All calloc/free are now Calloc/Free
 **
 ******************************************************************/

#include <stdlib.h>
#include <math.h>
#include "rlm.h"
#include "psi_fns.h"
#include "threestep_common.h"

#include <lm.h>
#include <R.h> 
#include <Rdefines.h>

/*******************************************************************
 **
 ** double fit_Probeset_model(double *PM, double *MM, int length, int robust)
 **
 ** double *PM - vector of Perfect Match Intensities 
 ** double *MM - vector of Mismatch Intensities
 ** int length - number of probe pairs in probeset
 ** int robust - 0 if fit using standard linear regression, 1 if robust fit is 
 **              to be used.
 ** 
 **
 ** RETURNS value of parameter fitted in model for difference in effect between
 **         PM and MM probes
 **
 ** This function fits a joint model to PM and MM probes in particular to
 **
 ** y = (log2(PM_1),...,log2(PM_n),log2(MM_1),...,log2(MM_n))
 **
 ** fit the model
 **
 ** y = probepair effect + probetype effect + epsilon
 ** 
 ** where there are n probepair effects and 1 probetype effect (either PM or not/MM or not)
 **
 ** The probetype effect is returned by this function
 **
 ******************************************************************/

double fit_Probeset_model(double *PM, double *MM, int length, int robust, int psicode, double psi_k, int probepair_effects){

  int i; /*,j; */

  int nparams;

  double parameter;
  double tol = 1e-7;
  
  double *y = Calloc(2*length,double);
  double *weights= Calloc(2*length,double);
  double *x; /* = Calloc(2*length*(length +1),double); */
  double *out_beta = Calloc(length +1,double);
  double *out_resids = Calloc(2*length,double);

  

  /* create Y vector */

  for (i=0; i < length; i++){
    y[i] = log(PM[i])/log(2.0);
    weights[i] = 1.0;
  }

  for (i=0; i < length; i++){
    y[i+length] = log(MM[i])/log(2.0);
    weights[i+length] = 1.0;
  }
  
  /* make design matrix */
  

  if (probepair_effects){
    nparams = length +1;
    x = Calloc(2*length*(length +1),double);
    for (i=0; i < length; i++){
      x[2*length*(i+1) + i%length] = 1.0;
      x[2*length*(i+1) + (length) + i%length] = 1.0;  
    }
    for(i=0; i < length; i++){
      x[i] = 1.0;
    }

  
  } else {
    nparams =1;
    x = Calloc(2*length*nparams,double);  
    for(i=0; i < length; i++){
      x[i] = 0.5;
    }
    for(i=length; i < 2*length; i++){
      x[i] = -0.5;
    }
  }
  

  /*  for (i=0; i< 2*length; i++){
    for (j =0; j < length+1; j++){
      printf("%1.1f ",x[j*2*length + i]);
    }
    printf("\n");
    } */

  if (robust){
    rlm_fit(x, y, 2*length,nparams, out_beta, out_resids, weights,PsiFunc(psicode),psi_k,20,0);
  } else {    
    lm_wfit(x, y, weights, 2*length,nparams,tol,out_beta,out_resids);
  }

  
  parameter = out_beta[0];

  /*  printf("%f\n",parameter);*/
  Free(y);
  Free(weights);
  Free(x);
  Free(out_beta);
  Free(out_resids);

  return parameter;

}


double fit_Difference_model(double *PM, double *MM, int length, int robust, int psicode, double psi_k, int probepair_effects){


  int i; /*,j; */

  int nparams;

  double parameter;
  double tol = 1e-7;
  
  double *y = Calloc(length,double);
  double *weights= Calloc(length,double);
  double *x; /* = Calloc(2*length*(length +1),double); */
  double *out_beta = Calloc(1,double);
  double *out_resids = Calloc(length,double);

  

  /* create Y vector */

  for (i=0; i < length; i++){
    y[i] = log(PM[i])/log(2.0) -log(MM[i])/log(2.0) ;
    weights[i] = 1.0;
  }

  /* make design matrix */
  
  nparams =1;
  x = Calloc(length*nparams,double);  
  for(i=0; i < length; i++){
    x[i] = 1;
  }
  
  

  /*  for (i=0; i< 2*length; i++){
    for (j =0; j < length+1; j++){
      printf("%1.1f ",x[j*2*length + i]);
    }
    printf("\n");
    } */

  if (robust){
    rlm_fit(x, y, length,nparams, out_beta, out_resids, weights,PsiFunc(psicode),psi_k,20,0);
  } else {    
    lm_wfit(x, y, weights, length,nparams,tol,out_beta,out_resids);
  }

  
  parameter = out_beta[0];

  /*  printf("%f\n",parameter);*/
  Free(y);
  Free(weights);
  Free(x);
  Free(out_beta);
  Free(out_resids);

  return parameter;


}



double median_Difference(double *PM, double *MM, int length){
  
  int i;

  double *y = (double *)Calloc(length,double);
  double med_diff;
  /* create Y vector */

  for (i=0; i < length; i++){
    y[i] = log(PM[i])/log(2.0) -log(MM[i])/log(2.0) ;
  }

  med_diff = median(y,length);

  Free(y);

  return med_diff;
}



/*******************************************************************
 **
 ** void R_fit_Probeset_model(double *PM, double *MM, int *length, 
 **                           int *robust, int *psicode, double *psi_k)
 **
 **
 ** A function that can be called via the .C() interface from R
 ** interface to fit_Probeset_model. and fit_Difference_model
 **
 ******************************************************************/


void R_fit_Probeset_model(double *PM, double *MM, int *length, int *robust, int *psicode, double *psi_k, int *probepair_effects, int *difference, int *median){
  double result;
  if (!(*difference)){
    result = fit_Probeset_model(PM, MM, *length, *robust, *psicode, *psi_k,*probepair_effects);
  } else if(!(*median)) {
    result = fit_Difference_model(PM,MM,*length, *robust, *psicode, *psi_k, *difference);
  } else {
    result = median_Difference(PM,MM,*length);
  }

  PM[0] = result;
}
