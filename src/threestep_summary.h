#ifndef THREESTEP_SUMMARY_H 
#define THREESTEP_SUMMARY_H 1

#include "threestep_summary_methods.h"

void do_3summary(double *PM, const char **ProbeNames, int *rows, int *cols, double *results, char **outNames, int nps,void (* SummaryMeth)(double*, int, int, int *,double *, int, double *, summary_plist *),double *resultsSE, summary_plist *summary_param);

#endif
