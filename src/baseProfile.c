#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/RConverters.h>
#include <R_ext/Rdynload.h>
#include <string.h>
#define STR(SE) CHAR(STRING_ELT(SE,0))

SEXP gcrma_getSeq2(SEXP,SEXP,SEXP);
SEXP gcrma_getSeq2(SEXP psequence, SEXP x, SEXP length) {
    SEXP outMatrix;
    char *pseq;
    int k,i;
    R_len_t K,seql;
    K=INTEGER(x)[0];
    pseq=STR(psequence);
    seql=INTEGER(length)[0];

    PROTECT(outMatrix = allocMatrix(INTSXP, K, 3*seql));
        for(k=0;k<K;k++){
      for (i = 0; i < seql; i++) {
	if (pseq[k*seql+i] == 'A')
	  INTEGER(outMatrix)[i*K+k]=1;
	else
	  INTEGER(outMatrix)[i*K+k]=0;
	
	if (pseq[k*seql+i] == 'C')
	  INTEGER(outMatrix)[(i+seql)*K+k] = 1;
	else
	  INTEGER(outMatrix)[(i+seql)*K+k] = 0;
	
	if (pseq[k*seql+i] == 'G')
	  INTEGER(outMatrix)[(i+2*seql)*K+k] = 1;
	else
	  INTEGER(outMatrix)[(i+2*seql)*K+k] = 0;
      }
      }
    /*  for (i=0;i<strlen(pseq);i++){
	  INTEGER(outMatrix)[i]=100;
	  }*/
    UNPROTECT(1);
    return(outMatrix);
}
