#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/RConverters.h>
#include <R_ext/Rdynload.h>
#include <string.h>
#define STR(SE) CHAR(STRING_ELT(SE,0))

SEXP gcrma_getSeq2(SEXP,SEXP);
SEXP gcrma_getSeq2(SEXP psequence, SEXP x) {
    SEXP outMatrix;
    char *pseq;
    int k,i;
    R_len_t nx,K;
    K=INTEGER(x)[0];
    pseq=STR(psequence);

    PROTECT(outMatrix = allocMatrix(INTSXP, K, 75));
        for(k=0;k<K;k++){
      for (i = 0; i < 25; i++) {
	if (pseq[k*25+i] == 'A')
	  INTEGER(outMatrix)[i*K+k]=1;
	else
	  INTEGER(outMatrix)[i*K+k]=0;
	
	if (pseq[k*25+i] == 'C')
	  INTEGER(outMatrix)[(i+25)*K+k] = 1;
	else
	  INTEGER(outMatrix)[(i+25)*K+k] = 0;
	
	if (pseq[k*25+i] == 'G')
	  INTEGER(outMatrix)[(i+50)*K+k] = 1;
	else
	  INTEGER(outMatrix)[(i+50)*K+k] = 0;
      }
      }
    /*  for (i=0;i<strlen(pseq);i++){
	  INTEGER(outMatrix)[i]=100;
	  }*/
    UNPROTECT(1);
    return(outMatrix);
}
