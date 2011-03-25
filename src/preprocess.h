#ifndef PREPROCESS_H
#define PREPROCESS_H 1

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


SEXP GetParameter(SEXP alist, char *param_name);
SEXP pp_background(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP bg_type,SEXP background_parameters,SEXP verbosity);
SEXP pp_normalize(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_type, SEXP norm_parameters,SEXP verbosity);

#endif
