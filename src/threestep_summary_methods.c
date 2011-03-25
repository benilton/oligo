/********************************************************************
 **
 ** file: threestep_summary_methods.c
 **
 ** written by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** Created on: Jan 11, 2003
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** aim: general interface to threestep summary methods
 **      to add another summary method only the header file should need 
 **      be modified.
 **
 ** last modified: Feb 6, 2003
 **
 ** History
 ** 
 ** Jan 11, 2003 - Initial version
 ** Jan 13, 2003 - added rlm_threestep method.
 ** Feb 6, 2003 - added four new methods: LogAverage, LogMedianPM, MedianLogPM, LogNthLargestPM
 ** Jul 23, 2003 - a three step method should return a SE estimate
 ** Oct 5, 2003 - some of the function names have changed.
 ** Oct 10, 2003 - added equivalent functionality for threestepPLM summaries
 **
 ********************************************************************/

#include "threestep_summary_methods.h"
#include "stdlib.h"

int number_summary_methods=9;

pt2Summary funcArr[9];

pt2Summary SummaryMethod(int code){
  
  funcArr[0] = &median_polish_threestep;
  funcArr[1] = &TukeyBiweight_threestep;
  funcArr[2] = &AverageLog_threestep;
  funcArr[3] = &rlm_threestep;
  funcArr[4] = &LogAverage_threestep;
  funcArr[5] = &LogMedianPM_threestep;
  funcArr[6] = &MedianLogPM_threestep;
  funcArr[7] = &LogNthLargestPM;
  funcArr[8] = &lm_threestep;

  return funcArr[code];
}


pt2PLMSummary funcArr2[9];

pt2PLMSummary PLMSummaryMethod(int code){
  
  funcArr2[0] = &median_polish_threestep_PLM;
  funcArr2[1] = &TukeyBiweight_PLM;
  funcArr2[2] = &AverageLog_PLM;
  funcArr2[3] = &rlm_threestep_PLM;
  funcArr2[4] = &LogAverage_threestep_PLM;
  funcArr2[5] = &LogMedianPM_threestep_PLM;
  funcArr2[6] = &MedianLogPM_PLM;
  funcArr2[7] = &LogNthLargestPM_PLM;
  funcArr2[8] = &lm_threestep_PLM;

  return funcArr2[code];
}
