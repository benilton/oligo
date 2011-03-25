#ifndef LM_THREESTEP_H
#define LM_THREESTEP_H 1

#include  "threestep_summary_methods_param.h"

void lm_threestep(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param);
void lm_threestep_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param);



#endif
