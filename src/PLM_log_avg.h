#ifndef LOG_AVG_H
#define LOG_AVG_H 1

#include "threestep_summary_methods_param.h"

void LogAverage(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE);
void LogAverage_threestep(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param);
void LogAverage_threestep_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals,summary_plist *summary_param);


#endif
