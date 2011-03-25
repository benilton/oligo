#ifndef THREESTEP_SUMMARY_METHODS_H 
#define THREESTEP_SUMMARY_METHODS_H 1

#include "PLM_medianpolish.h"
#include "PLM_biweight.h"
#include "PLM_avg_log.h"
#include "rlm_threestep.h"
#include "PLM_log_avg.h"
#include "PLM_medianPM.h"
#include "PLM_median_logPM.h"
#include "nthLargestPM.h"
#include "lm_threestep.h"

#include "threestep_summary_methods_param.h"

typedef void (*pt2Summary)(double*, int, int, int *,double *, int, double *, summary_plist*);

pt2Summary SummaryMethod(int code);
int threestep_summary_code(char *Name);

typedef void (*pt2PLMSummary)(double*, int, int, int *,double *, int, double *, double*, summary_plist*);

pt2PLMSummary PLMSummaryMethod(int code);

#endif
