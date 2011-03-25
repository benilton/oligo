/*********************************************************************
 **
 ** file: threestep_common.c
 **
 ** Aim: Commonly used routines for threestep methods
 **
 ** Copyright (C) 2003 Ben Bolstad 
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Jan 13, 2003
 **
 ** Last modified: Jan 13, 2003
 **
 ** History
 **
 ** Jan 13, 2003 - Initial version
 ** Jul 25, 2003 - Added median_low
 ** Apr 5, 2004 - All malloc/free are now Calloc/Free
 ** Nov 13, 2006 - remove median which is now in rma_common.c/.h
 **
 ********************************************************************/


#include "rma_common.h"
#include "threestep_common.h"

#include <stdlib.h>
#include <R.h> 
#include <Rdefines.h>


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

double median_low(double *x, int length){
  int half;
  double med;
  double *buffer = Calloc(length,double);

  /* for (i = 0; i < length; i++)
     buffer[i] = x[i];

     qsort(buffer,length,sizeof(double), (int(*)(const void*, const void*))sort_double);
     half = (length + 1)/2;
     if (length % 2 == 1){
     med = buffer[half - 1];
     } else {
     med = buffer[half-1];
     }
  */
  memcpy(buffer,x,length*sizeof(double));

  half = (length + 1)/2;
  
  rPsort(buffer, length, half-1); 
  med = buffer[half - 1];
  Free(buffer);
  return med;
}
