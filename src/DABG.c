#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP R_DABG_P(SEXP X, SEXP Y, SEXP G){
  /************************************
   ** X is a N x S intensity matrix
   ** Y is a list of length L (=0:length(max(G))) of matrices with M_l rows x S columns
   ** G is an integer vector of length N
   ** Objective: mean( X[i, j] < Y[[G]][, j])
   ************************************/
  int N, S, M, i, j, k, count;
  double *ref;
  SEXP dim1, dim2, R;
  int *Gp;
  double *Xp, *Rp;

  Xp = NUMERIC_POINTER(AS_NUMERIC(X));
  Gp = INTEGER_POINTER(AS_INTEGER(G));


  PROTECT(dim1 = getAttrib(X, R_DimSymbol));
  N = INTEGER(dim1)[0];
  S = INTEGER(dim1)[1];
  PROTECT(R = allocMatrix(REALSXP, N, S));
  Rp = NUMERIC_POINTER(AS_NUMERIC(R));

  for (i=0; i < N; i++){
    ref = REAL(VECTOR_ELT(Y, Gp[i]));
    PROTECT(dim2 = getAttrib(VECTOR_ELT(Y, Gp[i]), R_DimSymbol));
    M = INTEGER(dim2)[0];
    for (j=0; j < S; j++){
      count = 0;
      for (k = 0; k < M; k++){
	if (Xp[i + j*N] < ref[k + j*M])
	  count++;
      }
      Rp[i + j*N] = (double) count/M;
    }
    UNPROTECT(1);
  }
  UNPROTECT(2);
  return R;
}
