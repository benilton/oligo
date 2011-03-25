/***************************************************************************
 **
 ** file: idealmismatch.c
 **
 ** aim: Implement Affymetrix ideal mismatch as documented in the
 **      Statistical Algorithms Description Document
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** written by: B. M. Bolstad   <bolstad@stat.berkeley.edu>
 ** 
 ** created: Oct 3, 2002
 **
 ** last Modified: Jan 11, 2003
 **
 ** History:
 ** Oct 3, 2002 - Initial version
 ** Oct 4, 2002 - R wrapper, bug test
 ** Oct 14, 2002 - A little c function to test IM on a single probeset
 **                with a R wrapper. Found and fixed an indexing error  
 ** Jan 2, 2003 - Clean up commenting, prepare for integration into 
 **               AffyExtensions 0.4           
 ** Jan 7, 2003 - More significant reworking, function renaming etc.
 **               remove R testing stuff
 ** Jan 11, 2003 - check up on some stack imbalance problems.(turns out to be in threestep.c)
 ** Feb 24, 2003 - add in missing #include <string.h>, remove or comment out unused variables
 ** Aug  7, 2003 - add a specific biweight only correction
 ** Apr  5, 2004 - all calloc/free are now Calloc/Free. (Note that dynamic array sizing still needs implementation)
 ** June 28, 2004 - dynamic array sizing implemented
 **
 ***************************************************************************/

#include "rma_common.h"
#include "biweight.h"
#include "idealmismatch.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h> 
#include <Rdefines.h>


/**************************************************************************************
 **
 ** void IdealMM_correct(double *PM, double *MM, int nprobes, const char** ProbeNames)
 **
 ** Correct PM by using the idea of an Ideal MM as documented in Affymetrix Statistical
 ** Algorithms Reference guide.
 **
 ** double *PM - PM probe intensities
 ** double *MM - MM probe intensities
 ** int nprobes - total number of PM probes
 ** char** ProbeNames - name for each probes indicating probeset membership.
 **
 *************************************************************************************/

static void IdealMM_correct_single(double *PM, double *MM, int rows, const char** ProbeNames){
  int i,j;
  double IM,SB;
  int current_buf_size=200;
  double contrast_tau = 0.03;
  double scale_tau = 10.0;
  char *curname =Calloc(200,char);
  int *cur_rows =Calloc(200,int); 
  double *buffer = Calloc(200,double);
  int k = 0;
  /*  int size; */
  const char *first;
  int first_ind;
  int nprobes;

  /*  int nprobes_set=0; */
  
  /**** build some structures so we can work with probes in the context of
	probesets ***/

  
  first = ProbeNames[0];
  first_ind = 0;
  i =0;
  nprobes = 1;
  for (j = 1; j < rows; j++){
    if ((strcmp(first,ProbeNames[j]) != 0) | (j == (rows -1))){
      if (nprobes > current_buf_size){
	current_buf_size = nprobes;
	cur_rows = Realloc(cur_rows,current_buf_size,int);
	buffer = Realloc(buffer,current_buf_size,double);
	  
      }
      if (j == (rows -1)){
        nprobes++;
        for (k = 0; k < nprobes; k++){
          cur_rows[k] = (j+1 - nprobes)+k;
          /* printf("%d ", (j+1 - nprobes)+k); */
        }
      } else {
        for (k = 0; k < nprobes; k++){
          cur_rows[k] = (j - nprobes)+k;
          /* printf("%d ", (j - nprobes)+k); */
        }
      }

      for (i = 0; i < nprobes; i++){
	buffer[i] = log(PM[cur_rows[i]])/log(2.0) - log(MM[cur_rows[i]])/log(2.0);
      }
      SB = Tukey_Biweight(buffer,nprobes);
      /* printf("%f \n ", SB); */
      for (i=0; i < nprobes; i++){
	if (PM[cur_rows[i]] > MM[cur_rows[i]]){
	  IM = MM[cur_rows[i]];
	} else if (SB > contrast_tau){
	  IM = PM[cur_rows[i]]/pow(2.0,SB);
	} else {
	  IM = PM[cur_rows[i]]/pow(2.0,contrast_tau/(1.0 + (contrast_tau - SB)/scale_tau));
	}
	/* printf("%f ",PM[cur_rows[i]]); */
	PM[cur_rows[i]] = PM[cur_rows[i]] - IM;
	/* printf("%f %f\n",PM[cur_rows[i]],IM); */
      }
      first = ProbeNames[j];
      first_ind = j;
      nprobes = 0;
    }
    nprobes++;
  }
    
  Free(buffer);
  Free(curname);
  Free(cur_rows);

}



/**************************************************************************************
 **
 ** void IdealMM_correct_R(double *PM, double *MM, int nprobes, const char** ProbeNames)
 **
 ** Wrapper to be called from R .Correct PM by using the idea of an Ideal MM as 
 ** documented in Affymetrix Statistical Algorithms Reference guide.
 **
 ** double *PM - PM probe intensities
 ** double *MM - MM probe intensities
 ** int nprobes - total number of PM probes
 ** const char** ProbeNames - name for each probes indicating probeset membership.
 **
 *************************************************************************************/

void IdealMM_correct(double *PM, double *MM, int *nprobes, int *nchips, const char** ProbeNames){
  int i =0;
  for (i =0; i < *nchips; i++){
    /* printf("%d \n",i); */
    IdealMM_correct_single(&PM[(*nprobes)*i],&MM[(*nprobes)*i], *nprobes, ProbeNames);
  }
}


/**************************************************************************************
 **************************************************************************************
 ** 
 ** Note that what follows is not a standard affymetrix method.
 **
 **************************************************************************************
 *************************************************************************************/




static void SpecificBiweightCorrect_single(double *PM, double *MM, int rows, const char** ProbeNames){
  int i,j;
  double IM,SB;
  /*  double contrast_tau = 0.03; */
  /*double scale_tau = 10.0; */
  int current_buf_size=200;
  /* char *curname =Calloc(200,char); */
  int *cur_rows =Calloc(200,int); 
  double *buffer = Calloc(200,double);
  int k = 0;
  /*  int size; */
  const char *first;
  int first_ind;
  int nprobes;

  /*  int nprobes_set=0; */
  
  /**** build some structures so we can work with probes in the context of
	probesets ***/

  

  first = ProbeNames[0];
  first_ind = 0;
  i =0;
  nprobes = 1;
  for (j = 1; j < rows; j++){
    if ((strcmp(first,ProbeNames[j]) != 0) | (j == (rows -1))){   
      if (nprobes > current_buf_size){
	current_buf_size = nprobes;
	cur_rows = Realloc(cur_rows,current_buf_size,int);
	buffer = Realloc(buffer,current_buf_size,double);  
      }
      if (j == (rows -1)){
        nprobes++;
        for (k = 0; k < nprobes; k++){
          cur_rows[k] = (j+1 - nprobes)+k;
          /* printf("%d ", (j+1 - nprobes)+k); */
        }
      } else {
        for (k = 0; k < nprobes; k++){
          cur_rows[k] = (j - nprobes)+k;
          /* printf("%d ", (j - nprobes)+k); */
        }
      }

      for (i = 0; i < nprobes; i++){
	buffer[i] = log(PM[cur_rows[i]])/log(2.0) - log(MM[cur_rows[i]])/log(2.0);
      }
      SB = Tukey_Biweight(buffer,nprobes);
      /* printf("%f \n ", SB); */
      for (i=0; i < nprobes; i++){
	/* if (SB > contrast_tau){ */
	  IM = PM[cur_rows[i]]/pow(2.0,SB);
	  /* } else {
	  ** IM = PM[cur_rows[i]]/pow(2.0,contrast_tau/(1.0 + (contrast_tau - SB)/scale_tau));
	  ** } 
	  */
      /* printf("%f ",PM[cur_rows[i]]); */
	PM[cur_rows[i]] = PM[cur_rows[i]] - IM;
	/* printf("%f %f\n",PM[cur_rows[i]],IM); */
      }
      first = ProbeNames[j];
      first_ind = j;
      nprobes = 0;
    }
    nprobes++;
  }
    
  Free(buffer);
/*  Free(curname); */
  Free(cur_rows);

}





void SpecificBiweightCorrect(double *PM, double *MM, int *nprobes, int *nchips, const char** ProbeNames){

  int i =0;
  for (i =0; i < *nchips; i++){
    /* printf("%d \n",i); */
    SpecificBiweightCorrect_single(&PM[(*nprobes)*i],&MM[(*nprobes)*i], *nprobes, ProbeNames);
  }
}
