#ifndef PLM_MEDIAN_LOGPM_H
#define PLM_MEDIAN_LOGPM_H 1


#include "threestep_summary_methods_param.h"

void MedianLogPM_threestep(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param);
void MedianLogPM_noSE(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes);


void MedianLogPM_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param);


#endif
