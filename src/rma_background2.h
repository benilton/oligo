#ifndef RMA_BACKGROUND2
#define RMA_BACKGROUND2 1

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>


SEXP bg_correct_c(SEXP PMmat, SEXP MMmat, SEXP densfunc, SEXP rho, SEXP bgtype);
SEXP bg_correct_c_copy(SEXP PMmat, SEXP MMmat, SEXP densfunc, SEXP rho, SEXP bytype);

#endif
