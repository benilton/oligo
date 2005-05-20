/******************************************************************************
 **
 ** file: rma_background.c
 **
 ** Copyright (C) 2002 - 2003 B. M. Bolstad
 ** 
 ** Written by: B. M. Bolstad  <bolstad@stat.berkeley.edu>
 ** Implementation dates: 2002-2003
 **
 ** implement background correction/PM adjustment for rma calculation
 ** note this is meant to be called from within R.
 ** 
 ** Note that the density estimation is carried out by a user supplied R function
 **
 ** this implements the same background as that included within affy (1.0.2) 
 ** and also the background in affy 1.1 and later releases.
 **
 ** Last changed: Apr 22, 2003 
 **
 ** History
 **
 ** Pre Nov 2002 - original version
 ** Nov 4, 2002 - make changes for Affy2. Add in Alternate bg computation.
 ** Nov 8, 2002 - testing
 ** Dec 31, 2002 - integrate changes into affy package
 ** Jan 2, 2003 - Documentation clean up
 ** Dec 26, 2002 - fixed non-ANSI C way of commenting out lines (Laurent)  
 ** Jan 6, 2003 actually // is correct according to standards. Some compilers are 
 **             of course not fully standards compliant :)
 **
 ** Jan 9, 2003 - check that background version switching happens
 ** Feb 6, 2003 - change two printf to Rprintf (so Windows users actually see some verbage)
 ** Feb 17,2003 - change a free to Free in find_max()
 ** Feb 25, 2003 - Fixes to remove some compiler warnings (show up when using -Wall with gcc)
 ** Apr 22, 2003 - fix error in bg_parameters2 so that can more closely duplicate the results in R
 ** Mar 06, 2004 - Changed last remaining malloc to a Calloc
 ** Aug 04, 2004 - Remove "Background correcting" message
 **
 *****************************************************************************/

#include "rma_background2.h"
#include "rma_common.h"
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h> 
#include <R_ext/Applic.h>

/***********************************************************
 **
 ** double find_max(double *x, int length)
 **
 ** this function returns the max of x
 **
 ************************************************************/

double find_max(double *x,int length){
  int i;
  double *buffer = Calloc(length,double);
  double max;
  
  for (i=0; i < length; i++){
    buffer[i] = x[i];
  }

  qsort(buffer,length,sizeof(double),(int(*)(const void*, const void*))sort_double);

  max = buffer[length-1];
  
  /* printf("max is %f \n", max); */
  Free(buffer);
  return max;
}

/***************************************************************
 **
 ** double get_sd(double *MM, double MMmax, int rows, int cols, int column)
 ** 
 ** double *PM - pm matrix
 ** double PMmax - value of mode of PMs for column
 ** int rows,cols - dimensions of matrix
 ** int column - column (chip) of interest
 **
 ** estimate the sigma parameter given vector MM value of maximum of density
 ** of MM, dimensions of MM matrix and column of interest
 **
 **
 ***************************************************************/

double get_sd(double *MM, double MMmax, int rows, int cols, int column){
  double sigma;
  double tmpsum = 0.0;
  int numtop=0;
  int i;


  for (i=0; i < rows; i++){
    if (MM[column*rows + i] < MMmax){
      tmpsum = tmpsum + (MM[column*rows + i] - MMmax)*(MM[column*rows + i] - MMmax);
      numtop++;
    }
  }
  sigma = sqrt(tmpsum/(numtop -1))*sqrt(2.0)/0.85;
  return sigma;
  
}

/*********************************************************************************
 **
 ** double  get_alpha(double *PM,double PMmax, int rows,int cols,int column)
 **
 ** double *PM - pm matrix
 ** double PMmax - value of mode of PMs for column
 ** int rows,cols - dimensions of matrix
 ** int column - column (chip) of interest
 **
 ** estimate the alpha parameter given vector PM value of maximum of density
 ** of PM, dimensions of MM matrix and column of interest
 **
 **
 ***********************************************************************/

double get_alpha(double *PM,double PMmax, int rows,int cols,int column){
  double alpha;
  double tmpsum = 0.0;
  int numtop=0;
  int i;

  for (i=0; i < rows; i++){
    if (PM[column*rows + i] > PMmax){
      tmpsum = tmpsum + (PM[column*rows + i] - PMmax);
      numtop++;
    }
  }
  alpha = numtop/tmpsum;
  return alpha ;
  
}

/**************************************************************************************
 **
 ** double max_density(double *z,int rows,int cols,int column, SEXP fn,SEXP rho)
 **
 ** double *z - matrix of dimension rows*cols
 ** int cols - matrix dimension
 ** int rows - matrix dimension
 ** int column - column of interest
 ** SEXP fn - R function for estimation of density
 ** SEXP rho - an R environment to work within
 **
 *************************************************************************************/

double max_density(double *z,int rows,int cols,int column, SEXP fn,SEXP rho){

  int i;
  int N;
  SEXP x;
  SEXP results;

  double *dens_x;
  double *dens_y;
  double max_y,max_x;
  
  PROTECT(x = allocVector(REALSXP,rows));

  for (i=0; i< rows; i++){
    REAL(x)[i] = z[column*rows +i];
  }

  /* for(i=0; i < rows; i++)
     printf("%f\n",REAL(x)[i]);*/
  
  defineVar(install("x"),x,rho);
  PROTECT(results = eval(fn,rho));
  
  N = INTEGER(VECTOR_ELT(results,3))[0];
  
  dens_x = NUMERIC_POINTER(VECTOR_ELT(results,0));
  dens_y = NUMERIC_POINTER(VECTOR_ELT(results,1));

  max_y = find_max(dens_y,16384);
  
  i = 0;
  do {
    if (dens_y[i] == max_y)
      break;
    i++;

  } while(1);
  
  max_x = dens_x[i];
  

  /*  PROTECT(names = getAttrib(results,R_NamesSymbol)); */
  

  /* for (i = 0; i < length(results); i++)
  //  printf("%S \n",CHAR(STRING_ELT(names,i)));

  //printf("max_x: %f\n",max_x);
  */

  UNPROTECT(2);


  return max_x;

}

/********************************************************************************
 **
 ** void bg_parameters(double *PM,double *MM, double *param, int rows, int cols, int column,SEXP fn,SEXP rho)
 **
 ** estimate the parameters for the background, they will be returned in *param
 ** param[0] is alpha, param[1] is mu, param[2] is sigma.
 **
 ** parameter estimates are same as those given by affy_1.1.1
 **
 **
 *******************************************************************************/

void bg_parameters(double *PM,double *MM, double *param, int rows, int cols, int column,SEXP fn,SEXP rho){
  
  double PMmax;
  double MMmax;
  double sd,alpha;
  
  
  PMmax = max_density(PM,rows, cols, column,fn,rho);
  MMmax = max_density(MM,rows, cols, column,fn,rho);

  sd = get_sd(MM,MMmax,rows,cols,column);
  alpha = get_alpha(PM,PMmax,rows,cols,column);

  param[0] = alpha;
  param[1] = MMmax;
  param[2] = sd;

}



/**********************************************************************************
 **
 ** double Phi(double x)
 **
 ** Compute the standard normal distribution function
 **
 *********************************************************************************/

double Phi(double x){
   return pnorm5(x,0.0,1.0,1,0);
}

/***********************************************************************************
 **
 ** double phi(double x)
 **
 ** compute the standard normal density.
 **
 **
 **********************************************************************************/

double phi(double x){
  return dnorm4(x,0.0,1.0,0);
}

/************************************************************************************
 **
 ** void bg_adjust(double *PM,double *MM, double *param, int rows, int cols, int column)
 **
 ** double *PM - PM matrix of dimension rows by cols
 ** double *MM - MM matrix of dimension rows by cols
 ** double *param - background model parameters
 ** int rows, cols - dimension of matrix
 ** int column - which column to adjust
 **
 ** note we will assume that param[0] is alpha, param[1] is mu, param[2] is sigma
 **
 ***********************************************************************************/

void bg_adjust(double *PM,double *MM, double *param, int rows, int cols, int column){
  int i;
  double a;
  
  for (i=0; i < rows; i++){
    a = PM[column*rows + i] - param[1] - param[0]*param[2]*param[2]; 
    PM[column*rows + i] = a + param[2] * phi(a/param[2])/Phi(a/param[2]);
  }
  
}

/***************************************************************
 **
 ** double  get_alpha2(double *PM,double PMmax, int rows,int cols,int column)
 **
 ** estimate the alpha parameter given vector PM value of maximum of density
 ** of PM, dimensions of MM matrix and column of interest using method proposed
 ** in affy2
 **
 **
 ***************************************************************/

double get_alpha2(double *PM, double PMmax, int length,SEXP fn,SEXP rho){
  double alpha;
  /* double tmpsum = 0.0;
     int numtop=0; */
  int i;

  for (i=0; i < length; i++){
    PM[i] = PM[i] - PMmax;
  }

  alpha = max_density(PM,length, 1,0,fn,rho);

  alpha = 1.0/alpha;
  return alpha ;  
}

/********************************************************************************
 **
 ** void bg_parameters2(double *PM,double *MM, double *param, int rows, int cols, int column,SEXP fn,SEXP rho)
 **
 ** estimate the parameters for the background, they will be returned in *param
 ** param[0] is alpha, param[1] is mu, param[2] is sigma.
 **
 ** parameter estimates are same as those given by affy in bg.correct.rma (Version 1.1 release of affy)
 **
 *******************************************************************************/


void bg_parameters2(double *PM,double *MM, double *param, int rows, int cols, int column,SEXP fn,SEXP rho){
  int i = 0;
  double PMmax;
  /* double MMmax; */
  double sd,alpha;
  int n_less=0,n_more=0;
  double *tmp_less = (double *)Calloc(rows,double);
  double *tmp_more = (double *)Calloc(rows,double);
  
  
  PMmax = max_density(PM,rows, cols, column,fn,rho);
  
  for (i=0; i < rows; i++){
    if (PM[column*rows +i] < PMmax){
      tmp_less[n_less] = PM[column*rows +i];
      n_less++;
    }

  }  

  PMmax = max_density(tmp_less,n_less,1,0,fn,rho);
  sd = get_sd(PM,PMmax,rows,cols,column)*0.85; 

  for (i=0; i < rows; i++){
    if (PM[column*rows +i] > PMmax) {
      tmp_more[n_more] = PM[column*rows +i];
      n_more++;
    }
  }

  /* the 0.85 is to fix up constant in above */
  alpha = get_alpha2(tmp_more,PMmax,n_more,fn,rho);

  param[0] = alpha;
  param[1] = PMmax;
  param[2] = sd;

  /* printf("%f %f %f\n",param[0],param[1],param[2]); */


  Free(tmp_less);
  Free(tmp_more);
}

/************************************************************************************
 **
 ** SEXP bg_correct_c(SEXP PMmat, SEXP MMmat, SEXP densfunc, SEXP rho)
 ** 
 ** given R matricies PMmat and MMmat background correct the columns of PMmat
 **
 ** SEXP PMmat - matrix of PM's
 ** SEXP MMmat - matrix of MM's
 ** SEXP densfunc - function for computing density function
 ** SEXP rho - r enviroment to work within.
 **
 ** this is the function to be called using .Call() from within R.
 **
 ** this function can be dangerous since it changes PM.
 **
 ***********************************************************************************/

SEXP bg_correct_c(SEXP PMmat, SEXP MMmat, SEXP densfunc, SEXP rho, SEXP bgtype){
  
  SEXP dim1;
  int j;
  int rows;
  int cols;
  double *PM,*MM;
  double param[3];

  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol));

  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  
  PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
  MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));
  /* printf("Background correcting\n"); */
  /* Rprintf("Background correcting\n"); */
  for (j=0; j < cols; j++){
    if (asInteger(bgtype) == 2){
      bg_parameters2(PM,MM,param,rows,cols,j,densfunc,rho);
    } else {
      bg_parameters(PM,MM,param,rows,cols,j,densfunc,rho);
    } 
    bg_adjust(PM,MM,param,rows,cols,j);
  }
  
  UNPROTECT(1);

  return PMmat;
}

/************************************************************************************
 **
 ** SEXP bg_correct_c_copy(SEXP PMmat, SEXP MMmat, SEXP densfunc, SEXP rho)
 ** 
 ** given R matricies PMmat and MMmat background correct the columns of PMmat
 **
 ** SEXP PMmat - matrix of PM's
 ** SEXP MMmat - matrix of MM's
 ** SEXP densfunc - function for computing density function
 ** SEXP rho - r enviroment to work within.
 **
 ** this is the function to be called using .Call() from within R.
 **
 ** this function copies the PM matrix and works on the copy.
 ** 
 ***********************************************************************************/

SEXP bg_correct_c_copy(SEXP PMmat, SEXP MMmat, SEXP densfunc, SEXP rho, SEXP bgtype){
  
  SEXP dim1,PMcopy;
  int j;
  int rows;
  int cols;
  double *PM,*MM;
  double param[3];

  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol));

  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  
  PROTECT(PMcopy = allocMatrix(REALSXP,rows,cols));
    
  copyMatrix(PMcopy,PMmat,0);

  PM = NUMERIC_POINTER(AS_NUMERIC(PMcopy));
  MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));
  /* printf("Background correcting\n"); */
  /* Rprintf("Background correcting\n"); */
  /* printf("%d \n", INTEGER(bgtype)[0]); */
  for (j=0; j < cols; j++){
    if (asInteger(bgtype)==2){
      /* printf("using type 2\n");*/ 
      bg_parameters2(PM,MM,param,rows,cols,j,densfunc,rho);
    } else {
      bg_parameters(PM,MM,param,rows,cols,j,densfunc,rho);
    } 
    bg_adjust(PM,MM,param,rows,cols,j);
  }
  
  UNPROTECT(2);

  return PMcopy;
}





