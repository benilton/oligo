#ifndef NTHLARGESTPM_H
#define NTHLARGESTPM_H 1


#include "threestep_summary_methods_param.h"

void LogNthLargestPM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes,double *resultsSE, summary_plist *summary_param);
void LogNthLargestPM_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes,double *resultsSE, double *residuals, summary_plist *summary_param);


#endif
