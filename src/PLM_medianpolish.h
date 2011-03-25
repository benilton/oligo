#ifndef PLM_MEDIANPOLISH_H
#define PLM_MEDIANPOLISH_H 1


#include "threestep_summary_methods_param.h"

void median_polish_threestep(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE,summary_plist *summary_param);
void median_polishPLM(double *data, int rows, int cols, int *cur_rows, double *probe_param, double *chip_param, double *intercept_param, int nprobes, double *residuals);
void median_polish_threestep_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param);

#endif
