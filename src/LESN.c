/***************************************************************************
 **
 ** file: LESN.c  (Originally this file was called AdHoc.c, changed to LESN.c,
 **                to keep consistant naming across R/C code. Offical renaming was on Mar 21, 2003)
 **
 ** aim: implement LESN - Low End Signal is Noise background/signal adjustment methods
 **                       implement stretching, shifting algorithms background/signal adjustments
 **
 ** Copyright (C) 2002-2004   Ben Bolstad
 **
 **
 ** created by: B. M. Bolstad  on  Dec 12, 2002
 **
 ** last modified: Feb 27, 2003
 **
 ** The functions here implement some AdHoc background methods
 ** one I call "shifting" which is to shift the distribution of probe intensities
 ** down to some level, the second I call "stretching" where intensities are
 ** strecthed so that the lowest intensities get moved the most and the higher intensities are
 ** hardly altered at all (down to some baseline level).
 **
 ** When I say baseline I mean a minimmum level (the lowest intensity probe on a chip will be 
 ** changed to have that value)
 **
 ** History:
 ** 
 ** Dec 12, 2002 - Initial version
 ** Feb 12, 2003 - Resume working with this method. Prepare for integration
 **                into AffyExtensions. Clean up documentation. rename bw_exponential2 to bw_gaussian
 ** Feb 13, 2003 - make a more general stretchdown function that takes more arbitrary
 **                weighting function. Integrate into AffyExtensions.
 **                R interfaces will be changed to .Call in future.
 ** Feb 14, 2003 - modify to deal with cases where minimum is less than baseline value
 **                the rule will be that all values less than the baseline will be floored at the baseline and
 **                no futher adjustment made.
 **                if at some later date it appears like there is a number of cases where there is one low outlier
 **                 and the rest of the intensities are higher and need to be adjusted then this should be changed.
 ** Feb 27, 2003 - Propose calling methods LESR for Low End Signal Reduction (or "Lesser"). Another 
 **                possibility would be LESRM  for  Low End Signal Reduce Magnitude or even
 **                LESRMA  for  Low End Signal Reduce Magnitude and Attenutation.
 **                Could easily replace E with I for Intensity. Yes another possibility is
 **                LESN (ie Lessen) for Low End Signal is Noise.
 ** Mar 21, 2003 - Rename this file to LESN.c. Add function LESN_correct to be called from outside
 ** Mar 22, 2003 - a few documentation touch ups
 ** Apr 5,  2004 - All malloc/free are now Calloc/Free
 ** May 6,  2005 - changed log2 to affyPLMlog2
 **
 *************************************************************************/

#include "LESN.h"
#include "rma_common.h"


#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/***************************************************************************
 **
 ** void shift_down(double *data, double P0, int rows, int cols)
 **
 ** double *data - the data to be adjusted
 ** double P0 - minimum to be adjusted too 
 ** int rows,cols - dimensions of Data.
 **
 ** Shift all intensities down so that minimum on chip is P0
 **
 ** P0 should really be the less than or equal to the minimum already on the chip
 **
 **************************************************************************/

static void shift_down(double *data, double P0, int rows, int cols){

  double *buffer = (double *)Calloc(rows,double);
  double Pmin;
  int i,j;

  for (j=0; j < cols; j++){
    
    for (i=0; i < rows; i++){
      buffer[i] = data[j*rows + i];
    }
    
    qsort(buffer,rows,sizeof(double), (int(*)(const void*, const void*))sort_double);
    
    Pmin = buffer[0];

    for (i=0; i < rows; i++){
      data[j*rows + i] = data[j*rows + i] - (Pmin - P0);
    }
  }

  Free(buffer);

}

/***************************************************************************
 **
 ** void R_shift_down(double *data, double P0, int rows, int cols)
 **
 ** double *data - the data to be adjusted
 ** double *P0 - minimum to be adjusted too
 ** int *rows,*cols - dimensions
 **
 ** shift all intensities down so that minimum on chip is P0
 **
 ** wrapper function that can be called from R
 **
 **************************************************************************/

void R_shift_down(double *data, double *P0, int *rows, int *cols){

  shift_down(data, *P0, *rows, *cols);
}



/******************************************************************************
 **
 ** double bw_linear(double x, double pmin, double pmax)
 **
 ** double x - current probe intensity
 ** double pmin - minimum probe intensity
 ** double pmax - maximum probe intensity
 **
 ** this function is the linear background weighting function
 **
 **
 ******************************************************************************/

static double bw_linear(double x, double pmin, double pmax, double theta){

  return(x - pmax)/(pmin - pmax);

}

/******************************************************************************
 **
 ** double bw_exponential(double x, double pmin, double pmax,double theta)
 **
 ** double x - current probe intensity
 ** double pmin - minimum probe intensity
 ** double pmax - maximum probe intensity
 ** double theta - a parameter in the weighting process
 **
 ** this function is the exponential background weighting function
 **
 **
 ******************************************************************************/

static double bw_exponential(double x, double pmin, double pmax, double theta){

  return exp(-(x - pmin)/theta);
}

/******************************************************************************
 **
 ** double bw_gaussian(double x, double pmin, double pmax)
 **
 ** double x - current probe intensity
 ** double pmin - minimum probe intensity
 ** double pmax - maximum probe intensity
 ** double theta - a parameter in the weighting process
 **
 ** this function is the half gaussian background weighting function
 **
 **
 ******************************************************************************/

static double bw_gaussian(double x, double pmin, double pmax, double theta){

  return exp(-(x - pmin)*(x-pmin)/theta);
}


/*******************************************************************************
 **
 **  static double affyPLMlog2(double x)
 **
 ** internal function for computing log base 2
 **
 **
 *******************************************************************************/

static double affyPLMlog2(double x){

  return log(x)/log(2.0);

}


/****************************************************************************
 **
 ** void stretch_down(double *data, double baseline, int rows, int cols, double theta, 
 **                  int use_logs, double (*BackWeighting)(double,double, double, double))
 **
 ** double *data
 ** double baseline
 ** int rows, int cols 
 ** double theta
 ** int use_logs
 ** double (*BackWeighting)(double,double, double, double)
 **
 *****************************************************************************/



void stretch_down(double *data, double baseline, int rows, int cols, double theta, int use_logs, double (*BackWeighting)(double,double, double, double)){

  double *buffer = (double *)Calloc(rows,double);
  double Pmin,Pmax;
  int i,j;

  for (j=0; j < cols; j++){
    
    for (i=0; i < rows; i++){
      buffer[i] = data[j*rows + i];
    }
    
    qsort(buffer,rows,sizeof(double), (int(*)(const void*, const void*))sort_double);
    
    Pmin = buffer[0];
    Pmax = buffer[rows-1];
    
    if (use_logs){
      if (Pmin < baseline){
	for (i=0; i < rows; i++){
	  if (data[j*rows + i] < baseline){
	    data[j*rows + i] = baseline;
	  }
	}
      } else {
	for (i=0; i < rows; i++){
	  data[j*rows + i] = pow(2.0,affyPLMlog2(data[j*rows + i]) -  BackWeighting(affyPLMlog2(data[j*rows + i]),affyPLMlog2(Pmin),affyPLMlog2(Pmax),theta)*(affyPLMlog2(Pmin) - affyPLMlog2(baseline)));
	}
      }
    } else {
      for (i=0; i < rows; i++){
	data[j*rows + i] = data[j*rows + i] -  BackWeighting(data[j*rows + i],Pmin,Pmax,theta)*(Pmin - baseline);
      }
    }
  }


  Free(buffer);
 
}

/***************************************************************************
 **
 ** void LESN_correct(double *data, int rows, int cols, int method, double baseline, double theta)
 **
 ** double *data - PM intensities to be corrected
 ** int rows - number of rows in matrix
 ** int cols - number of cols in matrix
 ** int method - integer indicating which of the LESN methods to use
 ** double baseline - baseline value. Minimum value which we will allow probe intensites to take
 ** double theta - tuning parameter used in LESN methods 
 **
 ** This function is c code interface that should be called from the outside when
 ** background correcting/signal adjusting using LESN methods.
 **
 ****************************************************************************/

void LESN_correct(double *data, int rows, int cols, int method, double baseline, double theta){

  /*  printf("Baseline:theta %f %f\n",baseline,theta); */

  if (method == 2){
    /* half gaussian */
    stretch_down(data,baseline,rows,cols,2*theta*theta,1,bw_gaussian);
  } else if (method == 1){
    /* exponential decay */
    stretch_down(data,baseline,rows,cols,theta,1,bw_exponential);
  } else {
    /* a simple shifting */
     shift_down(data, baseline, rows, cols);
  }
}


/***************************************************************************
 **
 ** void R_stretch_down(double *data, double P0, int rows, int cols)
 **
 ** double *data - the data to be adjusted
 ** double *P0 - minimum to be adjusted too
 ** int *rows,*cols - dimensions
 ** int type 0 is linear, type 1 is 
 ** shift all intensities down so that minimum on chip is P0
 **
 ** wrapper function that can be called from R
 **
 **************************************************************************/

void R_stretch_down(double *data, double *baseline, int *rows, int *cols,int *type,double *theta){


  if (*type ==1){
    stretch_down(data,*baseline,*rows,*cols,*theta,0,bw_linear);
  } else if (*type ==2) {
    stretch_down(data,*baseline,*rows,*cols,*theta,0,bw_exponential);
  } else if (*type ==3) {
    stretch_down(data,*baseline,*rows,*cols,*theta,1,bw_linear);
  }  else if (*type ==4) {
    stretch_down(data,*baseline,*rows,*cols,*theta,1,bw_exponential);
  } else if (*type == 5){
    stretch_down(data,*baseline,*rows,*cols,*theta,1,bw_gaussian);
  }
}




/***************************************************************************
 **
 ** void shift_down(double *data, double P0, int rows, int cols)
 **
 ** double *data - the data to be adjusted
 ** double P0 - minimum to be adjusted too
 ** int rows,cols - dimensions
 **
 ** shift all intensities down so that minimum on chip is P0
 **
 **
 **
 **************************************************************************/

void shift_down_log(double *data, double baseline, int rows, int cols){

  double *buffer = (double *)Calloc(rows,double);
  double Pmin;
  int i,j;

  for (j=0; j < cols; j++){
    
    for (i=0; i < rows; i++){
      buffer[i] = data[j*rows + i];
    }
    
    qsort(buffer,rows,sizeof(double), (int(*)(const void*, const void*))sort_double);
    
    Pmin = buffer[0];

    if (Pmin < baseline){
      for (i=0; i < rows; i++){
	if (data[j*rows + i] < baseline){
	  data[j*rows + i] = baseline;
	}
      }
    } else {
      for (i=0; i < rows; i++){
	data[j*rows + i] = pow(2.0,(affyPLMlog2(data[j*rows + i]) - (affyPLMlog2(Pmin) - affyPLMlog2(baseline))));
      }
    }
  }
  Free(buffer);
}



void R_shift_down_log(double *data, double *P0, int *rows, int *cols){

  shift_down_log(data, *P0, *rows, *cols);
}


