#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

/**************************************************************************
 ** ORIGINAL CODE BY BEN BOLSTAD
 ** THIS IS A MODIFICATION BY BENILTON, PROBLEMS HERE ARE BENILTON'S FAULT :-)
 ** (except for any problems that are caused by Ben)
 **************************************************************************/
 
/**************************************************************************
 **
 ** double median_nocopy(double *x, int length)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the median of *x. note x is not order preserved when this function
 ** is called.
 **
 *************************************************************************/

double  median_nocopy(double *x, int length){
  int i;
  int half;
  double med;
  double *buffer = x;  //Calloc(length,double);
  


  half = (length + 1)/2;
  /*  
      qsort(buffer,length,sizeof(double), (int(*)(const void*, const void*))sort_double);  
      
      if (length % 2 == 1){
      med = buffer[half - 1];
      } else {
      med = (buffer[half] + buffer[half-1])/2.0;
      }
  */

  rPsort(buffer, length, half-1);
  med = buffer[half-1];
  if (length % 2 == 0){
    rPsort(buffer, length, half);
    med = (med + buffer[half])/2.0;
  }
  
 
  return med;
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
  med_abs = median_nocopy(buffer,length);

  Free(buffer);

  return(med_abs);
}


static double med_abs_indices(double *x, int length, int *indices, int which_index, int length_index){
  int i,j;
  double med_abs;
  double *buffer = Calloc(length_index,double);
  /*  Rprintf("%d %d\n",which_index,length_index); */
  j=0;
  for (i = 0; i < length; i++){
    if (indices[i] == which_index){
      buffer[j] = fabs(x[i]);
      j++;
    }
  }
  med_abs = median_nocopy(buffer,length_index);

  Free(buffer);

  return(med_abs);
}






static void huber_rows2(double *data, double *m1, double *m2, int *m3, int *class, int rows, int cols, double k){
  int i, j, rep, n1, n2, n3, l;
  double mean1, mean2, mean3, scale1, scale2, scale3, sum_weights1, sum_weights2, sum_weights3;
  int done1, done2, done3;

  double *datvec=Calloc(cols,double);
  int *classvec=Calloc(cols,int);
 
  for (i=0; i < rows; i++){
    n1=0;
    n2=0;
    n3=0;
    mean1=0.0;
    mean2=0.0;
    mean3=0.0;
    done1 = done2 = done3 = 0;
    
    for (j=0; j < cols; j++){
      if (class[j*rows + i] == 1){
	datvec[j]=data[j*rows + i];
	mean1+=datvec[j];
	++n1;
	classvec[j] = 1;
      } else if (class[j*rows + i] == 2){
	datvec[j]=data[j*rows + i];
	mean2+=datvec[j];
	++n2;
	classvec[j] = 2;
      } else if (class[j*rows + i] == 3){
	datvec[j]=data[j*rows + i];
	mean3+=datvec[j];
	++n3;
	classvec[j] = 3;
      } else {
	/* Should be the NA's */
	classvec[j] = class[j*rows + i];
      }
    }

    if (n1 > 0){
      mean1/=n1;
    } else {
      mean1 = R_NaReal;
      scale1 = R_NaReal;
      done1 = 1;
    }

    if (n2 > 0){
      mean2/=n2;
    } else {
      mean2 = R_NaReal;
      scale2 = R_NaReal;
      done2 = 1;
    }

    if (n3 > 0){
      mean3/=n3;
    } else {
      mean3 = R_NaReal;
      scale3 = R_NaReal;
      done3 = 1;
    }

    /* ten step huber estimate of location */

    for (rep=0; rep < 10; rep++){
      for (j=0; j < cols; j++){
	if (classvec[j] == 1){
	  datvec[j]=data[j*rows+i]-mean1;
	}
	if (classvec[j] == 2){
	  datvec[j]=data[j*rows +i]-mean2;
	}
	if (classvec[j]== 3){
	  datvec[j]=data[j*rows + i]-mean3;
	}
      }
      
      if (n1 > 0)
	scale1=med_abs_indices(datvec, cols, classvec, 1, n1)/0.6744908;
      if (n2 > 0)
	scale2=med_abs_indices(datvec, cols, classvec, 2, n2)/0.6744908;
      if (n3 >  0)
	scale3=med_abs_indices(datvec, cols, classvec, 3, n3)/0.6744908;
      
      
      for (j=0; j < cols; j++){
	if (classvec[j] == 1){
	  if (n1 > 1){
	    datvec[j]/=scale1;
	  } else {
	    mean1 = datvec[j];
	    scale1 = 0.0;
	    done1 = 1;
	  }

	} else if (classvec[j] == 2){
	  if (n2 > 1){
	    datvec[j]/=scale2;
	  } else {
	    mean2 = datvec[j];
	    scale2 = 0.0;
	    done2 = 1;
	  }
	} else if (classvec[j] == 3){
	  if (n3 > 1){
	    datvec[j]/=scale3;
	  } else {
	    mean3 = datvec[j];
	    scale3 = 0.0;
	    done3 = 1;
	  }

	  if (!done1)
	    mean1 = 0.0;
	  if (!done2)
	    mean2 = 0.0;
	  if (!done3)
	    mean3 = 0.0;

	  sum_weights1 = sum_weights2 = sum_weights3 = 0.0;
	  for (j=0; j < cols; j++){
	    if (classvec[j] == 1 && !done1){
	      mean1+=weights_huber(datvec[j],k) * data[j*rows + i];
	      sum_weights1+=weights_huber(datvec[j],k);
	    } else if (classvec[j] == 2 && !done2){
	      mean2+=weights_huber(datvec[j],k) * data[j*rows + i];
	      sum_weights2+=weights_huber(datvec[j],k);
	    } else if (classvec[j] == 3 && !done3){
	      mean3+=weights_huber(datvec[j],k) * data[j*rows + i];
	      sum_weights3+=weights_huber(datvec[j],k);
	    }
	  }
	  if (!done1)
	    mean1/=sum_weights1;
	  if (!done2)
	    mean2/=sum_weights2;
	  if (!done3)
	    mean3/=sum_weights3;

	}
      }

      /* Rprintf("%f %f %f %f %f %f %d %d %d %f %f %f \n",mean1,mean2,mean3,scale1,scale2,scale3,n1,n2,n3,sum_weights1,sum_weights2,sum_weights3); */

    }
    m1[i]=mean1;
    m1[rows + i]=mean2;
    m1[2*rows + i]=mean3;
    m2[i]=scale1;
    m2[rows + i]=scale2;
    m2[2*rows + i]=scale3;
    m3[i]=n1;
    m3[rows + i]=n2;
    m3[2*rows + i]=n3;
  }
  Free(datvec);
}

SEXP R_HuberMatrixRows2(SEXP X, SEXP Y, SEXP K){
  SEXP dim1;
  SEXP center, scale, output, sizes;
  double *Xptr, *Mptr1, *Mptr2;
  int *Yptr, *Mptr3;
  int rows, cols;
  // SEXP dimnames;

  double k;
  
  PROTECT(dim1 = getAttrib(X,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  
  Xptr = NUMERIC_POINTER(AS_NUMERIC(X));
  Yptr = INTEGER_POINTER(AS_INTEGER(Y));

  PROTECT(center = allocMatrix(REALSXP, rows, 3));
  PROTECT(scale = allocMatrix(REALSXP, rows, 3));
  PROTECT(sizes = allocMatrix(INTSXP, rows, 3));

  Mptr1 = NUMERIC_POINTER(center);
  Mptr2 = NUMERIC_POINTER(scale);
  Mptr3 = INTEGER_POINTER(sizes);

  k = NUMERIC_POINTER(K)[0];  
  
  huber_rows2(Xptr, Mptr1, Mptr2, Mptr3, Yptr, rows, cols, k);

  /*  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, getAttrib(X, R_RowNamesSymbol));
  SET_VECTOR_ELT(dimnames, 1, R_NilValue);
  setAttrib(center, R_DimNamesSymbol, dimnames);
  setAttrib(scale, R_DimNamesSymbol, dimnames);
  setAttrib(sizes, R_DimNamesSymbol, dimnames); */
     
  PROTECT(output = allocVector(VECSXP,3));
  SET_VECTOR_ELT(output, 0, center);
  SET_VECTOR_ELT(output, 1, scale);
  SET_VECTOR_ELT(output, 2, sizes);

  UNPROTECT(5);
  
  return output;
}
