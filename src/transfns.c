/*********************************************************************
 **
 ** file: trans_fns.c
 **
 ** Aim: transformation functions to be applied to the data
 **      
 **
 ** Copyright (C) 2004 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 **
 ** Created on July 27, 2004
 **
 ** History
 **
 ********************************************************************/

#include "transfns.h"
#include "math.h"

#define SGN(A) ((A) < 0.0 ? -1.0 : ((A) > 0.0 ? 1.0 : 0.0))




double trans_log2(double x){
  return log(x)/log(2.0);
}



double trans_loge(double x){
  return log(x);
}

double trans_log10(double x){
  return log10(x);
}


double trans_sqrt(double x){
  return sqrt(x);
}


double trans_cuberoot(double x){
  return SGN(x)*pow(fabs(x),0.3333333333333333333);

}



static pt2trans transArr[5];

pt2trans transFunc(int code){
  
  transArr[0] = &trans_log2;
  transArr[1] = &trans_loge;
  transArr[2] = &trans_log10;
  transArr[3] = &trans_sqrt;
  transArr[4] = &trans_cuberoot;
  
  return transArr[code];
}
