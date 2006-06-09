/**********************************************************
 **
 ** file: qnorm.c
 **
 ** aim: A c implementation of the quantile normalization method 
 **
 ** Copyright (C) 2002-2006    Ben Bolstad
 **
 ** written by: B. M. Bolstad  <bmb@bmbolstad.com>
 **
 ** written: Feb 2, 2002
 ** last modified: Jun 4, 2006
 ** 
 ** This c code implements the quantile normalization method
 ** for normalizing high density oligonucleotide data as discussed
 ** in
 **
 ** Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003)(2003) 
 ** A Comparison of Normalization Methods for High 
 ** Density Oligonucleotide Array Data Based on Bias and Variance.
 ** Bioinformatics 19,2,pp 185-193
 **
 ** History
 ** Feb 2, 2002 - Intial c code version from original R code
 ** Apr 19, 2002 - Update to deal more correctly with ties (equal rank)
 ** Jan 2, 2003 - Documentation/Commenting updates reformating
 ** Feb 17, 2003 - add in a free(datvec) to qnorm(). clean up freeing of dimat
 ** Feb 25, 2003 - try to reduce or eliminate compiler warnings (with gcc -Wall)
 ** Feb 28, 2003 - update reference to normalization paper in comments
 ** Mar 25, 2003 - ability to use median, rather than mean in so called "robust" method
 ** Aug 23, 2003 - add ability to do normalization on log scale in "robust" method.
 **                also have added .Call() interface c functions which may be called
 **                now from R as alterative to traditonal means.
 **                Fixed a bug where use_median was not being dereferenced in "robust method"
 ** Oct 7, 2003 - fix a bug with length is qnorm_robust
 ** Mar 6, 2004 - change malloc/free pairs to Calloc/Free
 ** Mar 3, 2005 - port across the low memory quantile normalization from RMAExpress (and make it the new qnorm_c (previous version made qnorm_c_old)
 ** Mar 12, 2006 - make some internal functions static
 ** Mar 13, 2006 - re-working of the "robust" quantile normalizer. The old function is
 **                still here with a _old added to the name. Also now
 **                have a .Call() interface for the robust method
 ** Apr 27-28, 2006 - Add C level functionality for determining which outliers
 **                to exclude for the "robust" quantile normalizer.
 ** Jun 2, 2006  - Add a quantile normalization function that accepts a target
 **                distribution. Improve/add a few comments
 ** Jun 4, 2006 - Add a .Call interface for target based quantile normalization.
 **               Add a function for determing target distribution.
 ** Jun 5, 2006 - Re-organize code blocks
 **               Add normalization within blocks functions
 ** Jun 9, 2006 - change nearbyint to floor(x +0.5) (to fix problems on Sparc Solaris builds)
 **
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rma_common.h"
#include "qnorm.h"


#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
 
/*****************************************************************************************************
 *****************************************************************************************************
 **
 ** This section defines utility functions and data types
 **
 **
 *****************************************************************************************************
 *****************************************************************************************************/


/*************************************************************
 **
 ** the dataitem record is used to keep track of data indicies 
 ** along with data value when sorting and unsorting in the 
 ** quantile algorithm.
 **
 ************************************************************/

typedef struct{
  double data;
  int rank;
} dataitem;
  


/*************************************************************
 **
 ** the dataitem_block record is used to keep track of data indicies 
 ** along with data value when sorting and unsorting in the 
 ** quantile algorithm in blocks
 **
 ************************************************************/

typedef struct{
  double data;
  int rank;
  int block;
} dataitem_block;


/***********************************************************
 **  
 ** int min(int x1, int x2)							    
 **
 ** returns the minimum of x1 and x2
 **		    
 **********************************************************/

static int min(int x1,int x2){
  if (x1 > x2)
    return x2;
  else
    return x1;
}

/**********************************************************
 **
 ** int sort_fn(const void *a1,const void *a2)
 **
 ** a comparison function for sorting objects of the dataitem type.
 **
 **
 **********************************************************/

static int sort_fn(const void *a1,const void *a2){
  dataitem *s1, *s2;
  s1 = (dataitem *)a1;
  s2 = (dataitem *)a2;
  
  if (s1->data < s2->data)
    return (-1);
  if (s1 ->data > s2->data)
    return (1);
  return 0;
}


/**********************************************************
 **
 ** int sort_fn_blocks(const void *a1,const void *a2)
 **
 ** a comparison function for sorting objects of the dataitem_blocks type.
 **
 **
 **********************************************************/

static int sort_fn_blocks(const void *a1,const void *a2){
  dataitem_block *s1, *s2;
  s1 = (dataitem_block *)a1;
  s2 = (dataitem_block *)a2;
  
  if (s1->block < s2->block){
    return (-1);
  } else if (s1->block > s2->block){
    return (1);
  } else {
    if (s1->data < s2->data)
      return (-1);
    if (s1 ->data > s2->data)
      return (1);
    return 0;
  }
}





/************************************************************
 **
 ** dataitem **get_di_matrix(double *data, int rows, int cols)
 **
 ** given data  form a matrix of dataitems, each element of
 ** matrix holds datavalue and original index so that 
 ** normalized data values can be resorted to the original order
 **
 ***********************************************************/

static dataitem **get_di_matrix(double *data, int rows, int cols){
  int i,j;
  dataitem **dimat;
  /* dataitem *xtmp; */
  
  dimat = (dataitem **)Calloc((cols),dataitem *);
  
  if (dimat == NULL){
    printf("\nERROR - Sorry the normalization routine could not allocate adequate memory\n       You probably need more memory to work with a dataset this large\n");
  }

  /* xtmp = malloc(cols*rows*sizeof(dataitem));
     for (j=0; j < cols; j++, xtmp +=rows) dimat[j] = xtmp; */
  
  for (j=0; j < cols; j++)
    dimat[j] = Calloc(rows,dataitem);



  for (j =0; j < cols; j++)
    for (i =0; i < rows; i++){
      dimat[j][i].data = data[j*rows + i];
      dimat[j][i].rank = i;
    }

  return(dimat); 
}

/************************************************************
 **
 ** double *get_ranks(dataitem *x,int n)
 **
 ** get ranks in the same manner as R does. Assume that *x is
 ** already sorted
 **
 *************************************************************/

static void get_ranks(double *rank, dataitem *x,int n){
  int i,j,k;
   
  i = 0;

  while (i < n) {
    j = i;
    while ((j < n - 1) && (x[j].data  == x[j + 1].data))
      j++;
    if (i != j) {
      for (k = i; k <= j; k++)
	rank[k] = (i + j + 2) / 2.0;
    }
    else
      rank[i] = i + 1;
    i = j + 1;
  }
  /*return rank;*/
}


/************************************************************
 **
 ** double *get_ranks_blocks(dataitem *x,int n)
 **
 ** get ranks in the same manner as R does. Assume that *x is
 ** already sorted
 **
 *************************************************************/

static void get_ranks_blocks(double *rank, dataitem_block *x,int n){
  int i,j,k;
   
  i = 0;

  while (i < n) {
    j = i;
    while ((j < n - 1) && (x[j].data  == x[j + 1].data) && (x[j].block  == x[j + 1].block))
      j++;
    if (i != j) {
      for (k = i; k <= j; k++)
	rank[k] = (i + j + 2) / 2.0;
    }
    else
      rank[i] = i + 1;
    i = j + 1;
  }
  /*return rank;*/
}









/*************************************************************************
 **
 ** static double weights_huber(double u, double k)
 **
 ** double u - standardized residuals
 ** doubke k - tuning parameter
 **
 ** Used to get weights for M-estimation.
 **
 *************************************************************************/

static double weights_huber(double u, double k){
  
  if ( 1 < k/fabs(u)){
    return 1.0;
  } else {
    return  k/fabs(u);
  }
}


/**************************************************************************
 **
 ** double median(double *x, int length)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the median of *x
 **
 *************************************************************************/

static double median(double *x, int length){
  int i;
  int half;
  double med;
  double *buffer = Calloc(length,double);
  
  for (i = 0; i < length; i++)
    buffer[i] = x[i];
  
  qsort(buffer,length,sizeof(double), (int(*)(const void*, const void*))sort_double);
  half = (length + 1)/2;
  if (length % 2 == 1){
    med = buffer[half - 1];
  } else {
    med = (buffer[half] + buffer[half-1])/2.0;
  }
  
  Free(buffer);
  return med;
}



/**************************************************************************
 **
 ** static double med_abs(double *x, int length)
 **
 ** double *x - a data vector
 ** int length - length of x
 **
 ** Compute the median absolute value of a data vector
 **
 *************************************************************************/

static double med_abs(double *x, int length){
  int i;
  double med_abs;
  double *buffer = Calloc(length,double);

  for (i = 0; i < length; i++)
    buffer[i] = fabs(x[i]);

  med_abs = median(buffer,length);

  Free(buffer);
  return(med_abs);
}







/*****************************************************************************************************
 *****************************************************************************************************
 **
 ** THE FOLLOWING ARE DEFUNCT. HISTORICAL REMINDERS OF EARLIER IMPLEMENTATIONS
 **
 **
 *****************************************************************************************************
 *****************************************************************************************************/
 

/*********************************************************
 **
 ** void qnorm_c_old(double *data, int *rows, int *cols)
 **
 **  this is the function that actually implements the 
 ** quantile normalization algorithm. It is called from R
 **
 ** Previous implementation, replaced with lower memory overhead version below
 ** Remains here for historical interest only.
 ** 
 ********************************************************/

static void qnorm_c_old(double *data, int *rows, int *cols){
  int i,j,ind;
  dataitem **dimat;
  double sum;
  double *row_mean = (double *)Calloc(*rows,double);
  double *datvec = (double *)Calloc(*cols,double);
  double *ranks = (double *)Calloc(*rows,double);

  /*# sort original columns */
  
  dimat = get_di_matrix(data, *rows, *cols);
  
  for (j=0; j < *cols; j++){
    qsort(dimat[j],*rows,sizeof(dataitem),sort_fn);
  }

  /*# calculate means */
  
  for (i =0; i < *rows; i++){
    sum = 0.0;
    for (j=0; j < *cols; j++)
      datvec[j] = dimat[j][i].data;
    /*qsort(datvec,*cols,sizeof(double),(int(*)(const void*, const void*))sort_double); */
    for (j=0; j < *cols; j++){
      sum +=datvec[j]/(double)*cols;;
    }
    row_mean[i] = sum; /*/(double)*cols;*/
  }
  
  /*# unsort mean columns */
  for (j =0; j < *cols; j++){
    get_ranks(ranks,dimat[j],*rows);
    for (i =0; i < *rows; i++){
      ind = dimat[j][i].rank;
      data[j*(*rows) + ind] = row_mean[(int)floor(ranks[i])-1];
    }
  }
  Free(ranks);
  Free(datvec);   

  for (j=0; j < *cols; j++){
    Free(dimat[j]);
  }

  Free(dimat);
  Free(row_mean); 
}


/*********************************************************
 **
 ** void qnorm_robust_c_old(double *data,double *weights, int *rows, int *cols, int *use_median)
 **
 ** double *data - datamatrix
 ** double *weights - weights to give each chip when computing normalization chip
 ** int *rows, *cols - matrix dimensions
 ** int *use_median - 0 if using weighted mean, otherwise use median.
 ** int *use_log2  - 0 if natural scale, 1 if we should log2 the data and normalize
 **                  on that scale
 **
 ** this is the function that actually implements the 
 ** quantile normalization algorithm. It is called from R. 
 ** this function allows the user to downweight particular
 ** chips, in the calculation of the mean or to use the median
 ** rather than the mean. If median is used weights are ignored.
 **
 ** note that log scale with mean is equivalent to geometric mean.
 **
 ** SPECIAL NOTE. This is the old qnorm_robust_c function. 
 ** Remains here for historical interest only.
 **
 **
 ********************************************************/

static void qnorm_robust_c_old(double *data,double *weights, int *rows, int *cols, int *use_median, int *use_log2){
  int i,j,ind;
  int half,length;
  dataitem **dimat;
  double sum,sumweights;
  double *row_mean = Calloc(*rows,double);
  double *datvec = Calloc(*cols,double);
  double *ranks = Calloc(*rows,double);

  /* Log transform the data if needed */

  if ((*use_log2)){
    for (j =0; j < *cols; j++)
      for (i =0; i < *rows; i++){
	data[j *(*rows) + i] = log(data[j*(*rows) + i])/log(2.0);
      }
  }
  

  dimat = get_di_matrix(data, *rows, *cols);

  
  for (j=0; j < *cols; j++){
    qsort(dimat[j],*rows,sizeof(dataitem),sort_fn);
  }

  for (i =0; i < *rows; i++){
    sum = 0.0;
    for (j=0; j < *cols; j++)
      datvec[j] = dimat[j][i].data;
    /* qsort(datvec,*cols,sizeof(double),(int(*)(const void*, const void*))sort_double); */
    if (!(*use_median)){
      for (j=0; j < (*cols); j++){
	sum +=weights[j]*datvec[j];
      }
      sumweights = 0.0;
      for (j=0; j < (*cols); j++){
	sumweights = sumweights + weights[j];
      }
      row_mean[i] = sum/sumweights;
    } else {
       qsort(datvec,*cols,sizeof(double),(int(*)(const void*, const void*))sort_double);
       half = (*cols + 1)/2;
       length = *cols;
       if (length % 2 == 1){
	 row_mean[i] = datvec[half - 1];
       } else {
	 row_mean[i] = (datvec[half] + datvec[half-1])/2.0;
       }
    }

  }
  
  /*# unsort mean columns */
  for (j =0; j < *cols; j++){
    get_ranks(ranks,dimat[j],*rows);

    for (i =0; i < *rows; i++){
      ind = dimat[j][i].rank;
      data[j*(*rows) + ind] = row_mean[(int)floor(ranks[i])-1];
    }
  }

  /* antilog  return to the natural scale*/

   if ((*use_log2)){
    for (j =0; j < (*cols); j++)
      for (i =0; i < (*rows); i++){
	data[j *(*rows) + i] = pow(2.0,data[j*(*rows) + i]);
      }
   }

  Free(datvec);
  Free(ranks); 

  for (j=0; j < *cols; j++){
    Free(dimat[j]);
  }

  Free(dimat);
  Free(row_mean); 
}


/*****************************************************************************************************
 *****************************************************************************************************
 **
 ** The following block implements the standard quantile normalization function
 **
 **
 *****************************************************************************************************
 *****************************************************************************************************/

/*********************************************************
 **
 ** void qnorm_c(double *data, int *rows, int *cols)
 **
 **  this is the function that actually implements the
 ** quantile normalization algorithm. It is called from R.
 **
 ** returns 1 if there is a problem, 0 otherwise
 **
 ********************************************************/

int qnorm_c(double *data, int *rows, int *cols){
  int i,j,ind;
  dataitem **dimat;
  /*  double sum; */
  double *row_mean = (double *)Calloc((*rows),double);
  double *datvec; /* = (double *)Calloc(*cols,double); */
  double *ranks = (double *)Calloc((*rows),double);
  
  datvec = (double *)Calloc(*rows,double);
  
  for (i =0; i < *rows; i++){
    row_mean[i] = 0.0;
  }
  
  /* first find the normalizing distribution */
  for (j = 0; j < *cols; j++){
    for (i =0; i < *rows; i++){
      datvec[i] = data[j*(*rows) + i];
    }
    qsort(datvec,*rows,sizeof(double),(int(*)(const void*, const void*))sort_double);
    for (i =0; i < *rows; i++){
      row_mean[i] += datvec[i]/((double)*cols);
    }
  }
  
  /* now assign back distribution */
  dimat = (dataitem **)Calloc(1,dataitem *);
  dimat[0] = (dataitem *)Calloc(*rows,dataitem);
  
  for (j = 0; j < *cols; j++){
    for (i =0; i < *rows; i++){
      dimat[0][i].data = data[j*(*rows) + i];
      dimat[0][i].rank = i;
    }
    qsort(dimat[0],*rows,sizeof(dataitem),sort_fn);
    get_ranks(ranks,dimat[0],*rows);
    for (i =0; i < *rows; i++){
      ind = dimat[0][i].rank;
      data[j*(*rows) +ind] = row_mean[(int)floor(ranks[i])-1];
      }
  }
  
  Free(ranks);
  Free(datvec);
  Free(dimat[0]);
  
  Free(dimat);
  Free(row_mean);
  return 0;
}





/*********************************************************
 **
 ** SEXP R_qnorm_c(SEXP X)
 **
 ** SEXP X      - a matrix
 ** SEXP copy   - a flag if TRUE then make copy
 **               before normalizing, if FALSE work in place
 **               note that this can be dangerous since
 **               it will change the original matrix.
 **
 ** returns a quantile normalized matrix.
 **
 ** This is a .Call() interface for quantile normalization
 **
 *********************************************************/

SEXP R_qnorm_c(SEXP X, SEXP copy){

  SEXP Xcopy,dim1;
  double *Xptr;
  int rows,cols;
  
  PROTECT(dim1 = getAttrib(X,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  if (asInteger(copy)){
    PROTECT(Xcopy = allocMatrix(REALSXP,rows,cols));
    copyMatrix(Xcopy,X,0);
  } else {
    Xcopy = X;
  }
  Xptr = NUMERIC_POINTER(AS_NUMERIC(Xcopy));
  
  qnorm_c(Xptr, &rows, &cols);
  if (asInteger(copy)){
    UNPROTECT(2);
  } else {
    UNPROTECT(1);
  }
  return Xcopy;
}




/*****************************************************************************************************
 *****************************************************************************************************
 **
 ** The following block of code provides the "robust" quantile normalization. In addition it tries to 
 ** give the equivalent to the R code functionality for selecting arrays to remove before determining
 **
 **
 *****************************************************************************************************
 *****************************************************************************************************/

/*********************************************************
 **
 ** void qnorm_robust_c(double *data,double *weights, int *rows, int *cols, int *use_median,int *use_log2,int *weight_scheme)
 ** 
 ** double *data
 ** double *weights
 ** int *rows
 ** int *cols
 ** int *use_median
 ** int *use_log2
 ** int *weight_scheme
 **
 ** This function implements the "robust" quantile normalizer
 **
 ********************************************************/

int qnorm_robust_c(double *data,double *weights, int *rows, int *cols, int *use_median, int *use_log2, int *weight_scheme){
  
  int i,j,ind,rep;  
  int half,length;
  dataitem **dimat;
  double *row_mean = (double *)Calloc((*rows),double);
  double *datvec=0; /* = (double *)Calloc(*cols,double); */
  double *ranks = (double *)Calloc((*rows),double);
  
  double sum_weights = 0.0;
  double mean, scale; /* used in M-estimation routine */


  
  for (i =0; i < *rows; i++){
    row_mean[i] = 0.0;
  }


  if ((*weight_scheme == 0) && !(*use_median)){
    datvec = (double *)Calloc(*rows,double);
    
    if (!(*use_log2)){
      for (j = 0; j < *cols; j++){
	sum_weights+=weights[j];
      }


      for (j = 0; j < *cols; j++){
	for (i =0; i < *rows; i++){
	  datvec[i] = data[j*(*rows) + i];
	}
	qsort(datvec,*rows,sizeof(double),(int(*)(const void*, const void*))sort_double);
	if (weights[j] > 0.0){
	  for (i =0; i < *rows; i++){
	    row_mean[i] += weights[j]*datvec[i]/sum_weights;
	  }
	}
      } 
    } else {
      for (j = 0; j < *cols; j++){
	sum_weights+=weights[j];
      }


      for (j = 0; j < *cols; j++){
	for (i =0; i < *rows; i++){
	  datvec[i] = data[j*(*rows) + i];
	}
	qsort(datvec,*rows,sizeof(double),(int(*)(const void*, const void*))sort_double);
	if (weights[j] > 0.0){
	  for (i =0; i < *rows; i++){
	    row_mean[i] += weights[j]*(log(datvec[i])/log(2.0))/sum_weights;
	  }
	}
      } 
      for (i =0; i < *rows; i++){
	row_mean[i] = pow(2.0,row_mean[i]);
      }
    } 
  } else if ((*weight_scheme == 1) && !(*use_median)){
    /** row-wise huber weights **/
    dimat = get_di_matrix(data, *rows, *cols);
   
    datvec = Calloc(*cols,double);
    
    for (j=0; j < *cols; j++){
      qsort(dimat[j],*rows,sizeof(dataitem),sort_fn);
    }
    

    if (!(*use_log2)){
      for (i=0; i < *rows; i++){
	for (j=0; j < *cols; j++)
	  datvec[j] = dimat[j][i].data;
	
	/* five step huber estimate of location */
	mean = 0.0;
	for (j=0; j < *cols; j++){
	  mean += datvec[j]/(double)(*cols);
	}
	
	for (rep = 0; rep < 5; rep++){
	  for (j=0; j < *cols; j++){
	    datvec[j] = datvec[j] - mean;
	  }
	  scale = med_abs(datvec,*cols)/0.6745;
	  if (scale == 0.0){
	    break;
	  }
	  
	  for (j=0; j < *cols; j++){
	    datvec[j] = (datvec[j] - mean)/scale;
	  }
	  
	  mean = 0.0;
	  sum_weights=0.0;
	  for (j=0; j < *cols; j++){
	    mean+= weights_huber(datvec[j],1.345) * dimat[j][i].data;
	    sum_weights+=weights_huber(datvec[j],1.345);
	  }
	  mean/=sum_weights;
	  for (j=0; j < *cols; j++)
	    datvec[j] = dimat[j][i].data;
	  /* Rprintf("rep %d %f %f\n",rep,mean,scale); */
	}
	row_mean[i] = mean;
      }
    } else {
      for (i=0; i < *rows; i++){
	for (j=0; j < *cols; j++)
	  datvec[j] = log(dimat[j][i].data)/log(2.0);
	
	/* five step huber estimate of location */
	mean = 0.0;
	for (j=0; j < *cols; j++){
	  mean += datvec[j]/(double)(*cols);
	}
	
	for (rep = 0; rep < 5; rep++){
	  for (j=0; j < *cols; j++){
	    datvec[j] = datvec[j] - mean;
	  }
	  scale = med_abs(datvec,*cols)/0.6745;
	  if (scale == 0.0){
	    break;
	  }
	  
	  for (j=0; j < *cols; j++){
	    datvec[j] = (datvec[j] - mean)/scale;
	  }
	  
	  mean = 0.0;
	  sum_weights=0.0;
	  for (j=0; j < *cols; j++){
	    mean+= weights_huber(datvec[j],1.345) * log(dimat[j][i].data)/log(2.0);
	    sum_weights+=weights_huber(datvec[j],1.345);
	  }
	  mean/=sum_weights;
	  for (j=0; j < *cols; j++)
	    datvec[j] = log(dimat[j][i].data)/log(2.0);
	  /* Rprintf("rep %d %f %f\n",rep,mean,scale); */
	}
	row_mean[i] = pow(2.0,mean);
      }
    }
    for (j=0; j < *cols; j++){
      Free(dimat[j]);
    }
    
    Free(dimat);




  } else if ((*use_median)){
    dimat = get_di_matrix(data, *rows, *cols);
   
    datvec = Calloc(*cols,double);
    
    for (j=0; j < *cols; j++){
      qsort(dimat[j],*rows,sizeof(dataitem),sort_fn);
    }
    
    for (i=0; i < *rows; i++){
      for (j=0; j < *cols; j++)
	datvec[j] = dimat[j][i].data;
      
      
      qsort(datvec,*cols,sizeof(double),(int(*)(const void*, const void*))sort_double);
      half = (*cols + 1)/2;
      length = *cols;
      if (length % 2 == 1){
	row_mean[i] = datvec[half - 1];
      } else {
	row_mean[i] = (datvec[half] + datvec[half-1])/2.0;
      }
    }
    for (j=0; j < *cols; j++){
      Free(dimat[j]);
    }
    
    Free(dimat);
  } else {
    error("Not sure that these inputs are recognised for the robust quantile normalization routine.\n");


  }
	     




  /* now assign back distribution */
  dimat = (dataitem **)Calloc(1,dataitem *);
  dimat[0] = (dataitem *)Calloc(*rows,dataitem);
  
  for (j = 0; j < *cols; j++){
    for (i =0; i < *rows; i++){
      dimat[0][i].data = data[j*(*rows) + i];
      dimat[0][i].rank = i;
    }
    qsort(dimat[0],*rows,sizeof(dataitem),sort_fn);
    get_ranks(ranks,dimat[0],*rows);
    for (i =0; i < *rows; i++){
      ind = dimat[0][i].rank;
      data[j*(*rows) +ind] = row_mean[(int)floor(ranks[i])-1];
      }
  }
  
  Free(ranks);
  Free(datvec);
  Free(dimat[0]);
  
  Free(dimat);
  Free(row_mean);
  return 0;
  
}




/*********************************************************
 **
 ** SEXP R_qnorm_robust_c(SEXP X)
 **
 ** SEXP X      - a matrix
 ** SEXP copy   - a flag if TRUE then make copy
 **               before normalizing, if FALSE work in place
 **               note that this can be dangerous since
 **               it will change the original matrix.
 **
 ** returns a quantile normalized matrix.
 **
 ** This is a .Call() interface for quantile normalization (of the robust variety)
 **
 *********************************************************/


SEXP R_qnorm_robust_c(SEXP X, SEXP copy, SEXP R_weights, SEXP R_use_median, SEXP R_use_log2, SEXP R_weight_scheme){

  SEXP Xcopy,dim1;
  double *Xptr;
  int rows,cols;
  
  double *weights;
  int use_median;
  int use_log2;
  int weight_scheme;


  
  PROTECT(dim1 = getAttrib(X,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  if (asInteger(copy)){
    PROTECT(Xcopy = allocMatrix(REALSXP,rows,cols));
    copyMatrix(Xcopy,X,0);
  } else {
    Xcopy = X;
  }
  Xptr = NUMERIC_POINTER(AS_NUMERIC(Xcopy));
  
  weights =  NUMERIC_POINTER(AS_NUMERIC(R_weights));
  
  use_median = INTEGER(R_use_median)[0];
  use_log2 = INTEGER(R_use_log2)[0];
  weight_scheme = INTEGER(R_weight_scheme)[0];


  qnorm_robust_c(Xptr,weights, &rows, &cols, &use_median, &use_log2, &weight_scheme);
  if (asInteger(copy)){
    UNPROTECT(2);
  } else {
    UNPROTECT(1);
  }
  return Xcopy;
}



/*****************************************************************
 **
 ** static double compute_var(double *x, int length)
 **
 ** double *x - data vector
 ** int length - length of x
 **
 ** compute the sample variance of a data vector
 **
 *****************************************************************/

static double compute_var(double *x, int length){

  int i;
  double sum=0.0,sum2=0.0;

  for (i = 0; i < length; i++){
    sum+=x[i];
  }

  sum = sum/(double)length;
  for (i=0; i < length; i++){
    sum2+=(x[i]-sum)*(x[i] - sum);
  }

  return(sum2/(double)(length-1));
}

/*****************************************************************
 **
 ** static double compute_means(double *x, int length)
 **
 ** double *x - data vector
 ** int length - length of x
 **
 ** compute the sample mean of a data vector
 **
 **
 *****************************************************************/

static double compute_means(double *x, int length){

  int i;
  double sum=0.0; 

  for (i = 0; i < length; i++){
    sum+=x[i];
  }

  sum = sum/(double)length;
  
  return(sum);
}

/*****************************************************************
 **
 ** static void remove_order_variance(double *x, int rows, int cols, int n_remove, double *weights)
 **
 ** double *x 
 ** int rows
 ** int cols
 ** int n_remove
 ** double *weights
 **
 *****************************************************************/

static void remove_order_variance(double *x, int rows, int cols, int n_remove, double *weights){

  double *vars = Calloc(cols,double);
  double *vars_row = Calloc(cols,double);
  double *vars_col = Calloc(cols,double);

  double *results = Calloc(cols*cols,double);
  
  int i,j;


  for (j=0; j < cols; j++){
    vars[j] = compute_var(&x[j*rows],rows);
  }
  
  for (i = 0; i < cols -1; i++){
    for (j = i+1; j < cols; j++){
      results[j*cols + i] = vars[i]/vars[j];
      results[i*cols + j] = vars[j]/vars[i];
    }
  }

  for (i = 0; i < cols; i++){
    vars_row[i] = 0.0;
    for (j=0; j < cols; j++){
      vars_row[i]+=results[j*cols + i];
    }

  }


  for (j = 0; j < cols; j++){
    vars_col[j] = 0.0;
    for (i=0; i < cols; i++){
      vars_col[j]+=results[j*cols + i];
    }
  }
  

  for (j=0; j < cols; j++){
    vars_row[j] = vars[j] = vars_row[j] + vars_col[j];
  }

  qsort(vars_row,cols,sizeof(double),(int(*)(const void*, const void*))sort_double);
  
  for (i=cols-1; i >= cols - n_remove; i--){
    for (j=0; j < cols; j++){

      if (vars[j] == vars_row[i]){
	weights[j] =0.0;
	break;
      }
    }
  }
  
  Free(results);
  Free(vars);
  Free(vars_row);
  Free(vars_col);


}


/*****************************************************************
 **
 ** static void remove_order_mean(double *x, int rows, int cols, int n_remove, double *weights)
 **
 ** double *x 
 ** int rows
 ** int cols
 ** int n_remove
 ** double *weights
 **
 *****************************************************************/



static void remove_order_mean(double *x, int rows, int cols, int n_remove, double *weights){

  
  double *means = Calloc(cols,double);
  double *means_row = Calloc(cols,double);
  double *means_col = Calloc(cols,double);

  double *results = Calloc(cols*cols,double);
  
  int i,j;

  for (j=0; j < cols; j++){
    means[j] = compute_means(&x[j*rows],rows);
  }
  
  for (i = 0; i < cols -1; i++){
    for (j = i+1; j < cols; j++){
      results[j*cols + i] = means[i] - means[j];
      results[i*cols + j] = means[j]- means[i];
    }
  }


  for (j = 0; j < cols; j++){
    means_col[j] = 0.0;
    for (i=0; i < cols; i++){
      means_col[j]+=results[j*cols + i];
    }
  }
  

  for (j=0; j < cols; j++){
    means_row[j] = means[j] = fabs(means_col[j]);
  }

  qsort(means_row,cols,sizeof(double),(int(*)(const void*, const void*))sort_double);
  
  for (i=cols-1; i >= cols - n_remove; i--){
    for (j=0; j < cols; j++){
      if (means[j] == means_row[i]){
	weights[j] =0.0;
	break;
      }
    }
  }
  
  Free(results);
  Free(means);
  Free(means_row);
  Free(means_col);





}

/*****************************************************************
 **
 ** static void remove_order_both(double *x, int rows, int cols, int n_remove, double *weights)
 **
 ** double *x 
 ** int rows
 ** int cols
 ** int n_remove
 ** double *weights
 **
 *****************************************************************/


static void remove_order_both(double *x, int rows, int cols, int n_remove, double *weights){

  double *means = Calloc(cols,double);
  double *means_row = Calloc(cols,double);
  double *means_col = Calloc(cols,double);

  double *vars = Calloc(cols,double);
  double *vars_row = Calloc(cols,double);
  double *vars_col = Calloc(cols,double);



  double *results = Calloc(cols*cols,double);
  
  int i,j;

  int n_remove_mean;
  int n_remove_var;


  if (n_remove % 2 ==0){
    n_remove_var = n_remove/2;
    n_remove_mean = n_remove/2;
  } else {
    n_remove_var = n_remove/2 + 1;
    n_remove_mean = n_remove/2;
  }


  /* Work out all the stuff for excluding means */


  for (j=0; j < cols; j++){
    means[j] = compute_means(&x[j*rows],rows);
  }
  
  for (i = 0; i < cols -1; i++){
    for (j = i+1; j < cols; j++){
      results[j*cols + i] = means[i] - means[j];
      results[i*cols + j] = means[j]- means[i];
    }
  }


  for (j = 0; j < cols; j++){
    means_col[j] = 0.0;
    for (i=0; i < cols; i++){
      means_col[j]+=results[j*cols + i];
    }
  }
  

  for (j=0; j < cols; j++){
    means_row[j] = means[j] = fabs(means_col[j]);
  }

  qsort(means_row,cols,sizeof(double),(int(*)(const void*, const void*))sort_double);



  /* Work out all the stuff for excluding variances */


 for (j=0; j < cols; j++){
    vars[j] = compute_var(&x[j*rows],rows);
  }
  
  for (i = 0; i < cols -1; i++){
    for (j = i+1; j < cols; j++){
      results[j*cols + i] = vars[i]/vars[j];
      results[i*cols + j] = vars[j]/vars[i];
    }
  }

  for (i = 0; i < cols; i++){
    vars_row[i] = 0.0;
    for (j=0; j < cols; j++){
      vars_row[i]+=results[j*cols + i];
    }

  }


  for (j = 0; j < cols; j++){
    vars_col[j] = 0.0;
    for (i=0; i < cols; i++){
      vars_col[j]+=results[j*cols + i];
    }
  }
  

  for (j=0; j < cols; j++){
    vars_row[j] = vars[j] = vars_row[j] + vars_col[j];
  }

  qsort(vars_row,cols,sizeof(double),(int(*)(const void*, const void*))sort_double);
  
  for (i=cols-1; i >= cols - n_remove_var; i--){
    for (j=0; j < cols; j++){
      if (vars[j] == vars_row[i]){
	weights[j] =0.0;
	break;
      }
    }
  }

  for (i=cols-1; i >= cols - n_remove_mean; i--){
    for (j=0; j < cols; j++){
      if (means[j] == means_row[i]){
	if (weights[j] ==0.0){
	  /* means it has already been excluded by variance rule. So need to look one more along */
	  n_remove_mean+=1;
	} else {
	  weights[j] =0.0;
	  break;
	}
      }
    }
  }


}











SEXP R_qnorm_robust_weights(SEXP X, SEXP remove_extreme, SEXP n_remove){


  SEXP weights,dim1;


  int rows, cols;
  int j;

  PROTECT(dim1 = getAttrib(X,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];

  PROTECT(weights = allocVector(REALSXP,cols));

  for (j=0; j < cols; j++){
    REAL(weights)[j] = 1.0;
  }

  if (strcmp(CHAR(VECTOR_ELT(remove_extreme,0)),"variance") == 0){
    remove_order_variance(REAL(X), rows, cols, INTEGER(n_remove)[0], REAL(weights));
  }

  if (strcmp(CHAR(VECTOR_ELT(remove_extreme,0)),"mean") == 0){
    remove_order_mean(REAL(X), rows, cols, INTEGER(n_remove)[0], REAL(weights));
  }

  if (strcmp(CHAR(VECTOR_ELT(remove_extreme,0)),"both") == 0){
    remove_order_both(REAL(X), rows, cols, INTEGER(n_remove)[0], REAL(weights));
  }



  
  UNPROTECT(2);
  return weights;

}


/*****************************************************************************************************
 *****************************************************************************************************
 **
 ** The following block of code provides quantile normalization where a specified target vector is given.
 ** In addition it deals with cases of un equal length by estimating the appropriate quantiles
 **
 *****************************************************************************************************
 *****************************************************************************************************/



/*****************************************************************
 **
 ** int qnorm_c_using_target(double *data, int *rows, int *cols, double *target, int *targetrows)
 **
 ** double *data - a matrix of data to be normalized
 ** int *rows - dimensions of data
 ** int *cols - dimensions of data
 ** double *target - vector containing target distribution (ie distribution to be
 **                  normalized to)
 ** int *targetrows - length of target distribution vector
 **
 ** Note that it is assumed that there is no NA or Inf type values in the vectors. (ie no missing data)
 **
 ** if targetrows == rows then the standard methodology is used.
 ** 
 ** in other cases the appropriate quantiles to be normalized to are determined in a method
 ** equivalent to what you get using "type 8" with the quantile function
 **
 ** Note sample percentiles are calculated using i/(n+1)  (ie if there is 
 ** only 2 observations, the first sample percentile is 1/3 = 0.333,
 ** the second sample percentile will be 2/3 = 0.6666
 **
 ** 
 **
 *****************************************************************/



int qnorm_c_using_target(double *data, int *rows, int *cols, double *target, int *targetrows){


  
  int i,j,ind,target_ind;
  dataitem **dimat;

  double *row_mean; 
  double *datvec;
  double *ranks = (double *)Calloc((*rows),double);

  double samplepercentile;
  double target_ind_double,target_ind_double_floor;


  
  row_mean = (double *)Calloc(*targetrows,double);
  
  datvec = (double *)Calloc(*rows,double);
  


  /* first find the normalizing distribution */
  for (i =0; i < *targetrows; i++){
    row_mean[i] = target[i];
  }

  qsort(row_mean,*targetrows,sizeof(double),(int(*)(const void*, const void*))sort_double);

  
  if (*rows == *targetrows){
    /* now assign back distribution */
    /* this is basically the standard story */

    dimat = (dataitem **)Calloc(1,dataitem *);
    dimat[0] = (dataitem *)Calloc(*rows,dataitem);
    
    for (j = 0; j < *cols; j++){
      for (i =0; i < *rows; i++){
	dimat[0][i].data = data[j*(*rows) + i];
	dimat[0][i].rank = i;
      }
      qsort(dimat[0],*rows,sizeof(dataitem),sort_fn);
      get_ranks(ranks,dimat[0],*rows);
      for (i =0; i < *rows; i++){
	ind = dimat[0][i].rank;
	data[j*(*rows) +ind] = row_mean[(int)floor(ranks[i])-1];
      }
    }
  } else {
    /** need to estimate quantiles **/
    dimat = (dataitem **)Calloc(1,dataitem *);
    dimat[0] = (dataitem *)Calloc(*rows,dataitem);
    
    for (j = 0; j < *cols; j++){
      for (i =0; i < *rows; i++){
	dimat[0][i].data = data[j*(*rows) + i];
	dimat[0][i].rank = i;
      }
      qsort(dimat[0],*rows,sizeof(dataitem),sort_fn);
      get_ranks(ranks,dimat[0],*rows);
      for (i =0; i < *rows; i++){

	samplepercentile = (double)ranks[i]/(double)(*rows +1);
	target_ind_double = 1.0/3.0 + ((double)(*targetrows) + 1.0/3.0) * samplepercentile;
	target_ind_double_floor = floor(target_ind_double + 4*DOUBLE_EPS);
	
	target_ind_double = target_ind_double - target_ind_double_floor;

	if (fabs(target_ind_double) <=  4*DOUBLE_EPS){
	  target_ind_double = 0.0;
	}

	
	if (target_ind_double  == 0.0){
	  target_ind = (int)floor(target_ind_double_floor + 0.5); /* nearbyint(target_ind_double_floor); */	
	  ind = dimat[0][i].rank;
	  data[j*(*rows) +ind] = row_mean[target_ind-1];
	} else if (target_ind_double == 1.0){
	  target_ind = (int)floor(target_ind_double_floor + 1.5); /* (int)nearbyint(target_ind_double_floor + 1.0); */ 
	  ind = dimat[0][i].rank;
	  data[j*(*rows) +ind] = row_mean[target_ind-1];
	} else {
	  target_ind = (int)floor(target_ind_double_floor + 0.5); /* nearbyint(target_ind_double_floor); */	
	  ind = dimat[0][i].rank;
	  if ((target_ind < *targetrows) && (target_ind > 0)){
	    data[j*(*rows) +ind] = (1.0- target_ind_double)*row_mean[target_ind-1] + target_ind_double*row_mean[target_ind];
	  } else if (target_ind >= *targetrows){
	    data[j*(*rows) +ind] = row_mean[*rows-1];
	  } else {
	    data[j*(*rows) +ind] = row_mean[0];
	  }
	}
	
      }
    }
  }



  Free(ranks);
  Free(datvec);
  Free(dimat[0]);
  
  Free(dimat);
  Free(row_mean);
  return 0;



}





int qnorm_c_determine_target(double *data, int *rows, int *cols, double *target, int *targetrows){


  int i,j,ind,row_mean_ind;
  dataitem **dimat;
  /*  double sum; */
  double *row_mean = (double *)Calloc((*rows),double);
  double *datvec;
  double *ranks = (double *)Calloc((*rows),double);  
  
  double row_mean_ind_double,row_mean_ind_double_floor;
  double samplepercentile;


  
  datvec = (double *)Calloc(*rows,double);
  
  for (i =0; i < *rows; i++){
    row_mean[i] = 0.0;
  }
  
  /* first find the normalizing distribution */
  for (j = 0; j < *cols; j++){
    for (i =0; i < *rows; i++){
      datvec[i] = data[j*(*rows) + i];
    }
    qsort(datvec,*rows,sizeof(double),(int(*)(const void*, const void*))sort_double);
    for (i =0; i < *rows; i++){
      row_mean[i] += datvec[i]/((double)*cols);
    }
  }
  
  if (*rows == *targetrows){
    for (i =0; i < *rows; i++){
      target[i] = row_mean[i];
    }
  } else {
    /* need to estimate quantiles */
    for (i =0; i < *targetrows; i++){
      samplepercentile = (double)(i+1)/(double)(*targetrows +1);

      row_mean_ind_double = 1.0/3.0 + ((double)(*rows) + 1.0/3.0) * samplepercentile;
      row_mean_ind_double_floor = floor(row_mean_ind_double + 4*DOUBLE_EPS);
	
      row_mean_ind_double = row_mean_ind_double - row_mean_ind_double_floor;

      if (fabs(row_mean_ind_double) <=  4*DOUBLE_EPS){
	row_mean_ind_double = 0.0;
      }


      if (row_mean_ind_double  == 0.0){
	row_mean_ind = (int)floor(row_mean_ind_double_floor + 0.5);  /* (int)nearbyint(row_mean_ind_double_floor); */	
	target[i] = row_mean[row_mean_ind-1];
      } else if (row_mean_ind_double == 1.0){
	row_mean_ind = (int)floor(row_mean_ind_double_floor + 1.5);  /* (int)nearbyint(row_mean_ind_double_floor + 1.0); */ 
	target[i] = row_mean[row_mean_ind-1];
      } else {
	row_mean_ind =  (int)floor(row_mean_ind_double_floor + 0.5); /* (int)nearbyint(row_mean_ind_double_floor); */

	if ((row_mean_ind < *rows) && (row_mean_ind > 0)){
	  target[i] = (1.0- row_mean_ind_double)*row_mean[row_mean_ind-1] + row_mean_ind_double*row_mean[row_mean_ind];
	} else if (row_mean_ind >= *rows){
	  target[i] = row_mean[row_mean_ind-1];
	} else {
	  target[i] = row_mean[0];
	}
      }
    } 






  }
}




SEXP R_qnorm_using_target(SEXP X, SEXP target,SEXP copy){


  SEXP Xcopy,dim1,dim2;
  int rows, cols;
  int target_rows, target_cols;
  double *Xptr;
  double *targetptr;


  PROTECT(dim1 = getAttrib(X,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);
  if (asInteger(copy)){
    PROTECT(Xcopy = allocMatrix(REALSXP,rows,cols));
    copyMatrix(Xcopy,X,0);
  } else {
    Xcopy = X;
  }
  Xptr = NUMERIC_POINTER(AS_NUMERIC(Xcopy));

  if (isVector(target)){
    target_rows = length(target);
  } else if (isMatrix(target)){
    PROTECT(dim1 = getAttrib(X,R_DimSymbol));
    target_rows = INTEGER(dim1)[0];
    target_cols = INTEGER(dim1)[1];
    UNPROTECT(1);

    target_rows = target_rows*target_cols;
  } 

  
  targetptr = NUMERIC_POINTER(AS_NUMERIC(target));


  qnorm_c_using_target(Xptr, &rows, &cols,targetptr,&target_rows);


  UNPROTECT(1);
  return Xcopy;
}





SEXP R_qnorm_determine_target(SEXP X, SEXP targetlength){


  SEXP dim1,target;
  int rows, cols;
  int length;
  double *Xptr;
  double *targetptr;


  PROTECT(dim1 = getAttrib(X,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  length = asInteger(targetlength);

  Rprintf("%d\n",length);

  PROTECT(target=allocVector(REALSXP,length));


  Xptr = NUMERIC_POINTER(AS_NUMERIC(X));
  targetptr = NUMERIC_POINTER(target);
    
  qnorm_c_determine_target(Xptr,&rows,&cols,targetptr, &length);


  UNPROTECT(1);
  return target;

}




/*****************************************************************************************************
 *****************************************************************************************************
 **
 ** The following block of code implements quantile normalization within blocks.
 ** What this means is that the normalization is still carried out across arrrays (or columns)
 ** but separate subsets of rows (these are blocks) each get there own normalization
 **
 *****************************************************************************************************
 *****************************************************************************************************/


/*****************************************************************
 **
 ** int qnorm_c_within_blocks(double *x, int *rows, int *cols, int *blocks)
 ** 
 ** double *x - matrix to be normalized
 ** int *rows - dimensions of the matrix
 ** int *cols -
 ** int *blocks - labeling telling which block each row belongs to.
 **
 *****************************************************************/


int qnorm_c_within_blocks(double *x, int *rows, int *cols, int *blocks){


  int i,j,ind;
  dataitem_block **dimat_block;
  /*  double sum; */
  double *row_mean = (double *)Calloc((*rows),double);
  double *ranks = (double *)Calloc((*rows),double);
  

  dimat_block = (dataitem_block **)Calloc(1,dataitem_block *);
  dimat_block[0] = (dataitem_block *)Calloc(*rows,dataitem_block);
  
  for (i =0; i < *rows; i++){
    row_mean[i] = 0.0;
  }
  
  /* first find the normalizing distribution */
  for (j = 0; j < *cols; j++){
    for (i =0; i < *rows; i++){
      dimat_block[0][i].data = x[j*(*rows) + i];
      dimat_block[0][i].block = blocks[i];
    }
    qsort(dimat_block[0],*rows,sizeof(dataitem_block),sort_fn_blocks);
    /*   for (i=0; i < *rows; i++){
      Rprintf("%f %d\n",dimat_block[0][i].data,dimat_block[0][i].block);
      } */
    

    for (i =0; i < *rows; i++){
      row_mean[i] += dimat_block[0][i].data/((double)*cols);
    }
  }
  
  /* now assign back distribution */

  
  for (j = 0; j < *cols; j++){
    for (i =0; i < *rows; i++){
      dimat_block[0][i].data = x[j*(*rows) + i];
      dimat_block[0][i].block = blocks[i];
      dimat_block[0][i].rank = i;
    }
    qsort(dimat_block[0],*rows,sizeof(dataitem_block),sort_fn_blocks);
    get_ranks_blocks(ranks,dimat_block[0],*rows);
    for (i =0; i < *rows; i++){
      ind = dimat_block[0][i].rank;
      x[j*(*rows) +ind] = row_mean[(int)floor(ranks[i])-1];
    }
  }
  
  Free(ranks);
  
  Free(dimat_block[0]);
  
  Free(dimat_block);
  Free(row_mean);
  return 0;
  


}



SEXP R_qnorm_within_blocks(SEXP X,SEXP blocks,SEXP copy){

  SEXP Xcopy,dim1,blocksint;
  double *Xptr;
  int *blocksptr;
  int rows,cols;
  
  PROTECT(dim1 = getAttrib(X,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  UNPROTECT(1);

  if (asInteger(copy)){
    PROTECT(Xcopy = allocMatrix(REALSXP,rows,cols));
    copyMatrix(Xcopy,X,0);
  } else {
    Xcopy = X;
  }

  PROTECT(blocksint = coerceVector(blocks,INTSXP));
    

  Xptr = NUMERIC_POINTER(AS_NUMERIC(Xcopy));
  blocksptr  = INTEGER_POINTER(blocksint);


  
  qnorm_c_within_blocks(Xptr, &rows, &cols,blocksptr);
  if (asInteger(copy)){
    UNPROTECT(2);
  } else {
    UNPROTECT(1);
  }
  return Xcopy;



}
