/*------------------------------------------------*/
/* get the CGAT content of a sequence             */
/*------------------------------------------------*/

#include <R.h>
#include <Rinternals.h>
#include "R_ext/Arith.h"
#include "R_ext/Error.h"
#include "R_ext/Applic.h" /* machar */

#include <string.h>
#include <stdlib.h>

//#include "strutils.h"

char errmess[256];

SEXP basecontent(SEXP x)
{
  SEXP rv, rownames, colnames, dimnames, dim;
  const char *seq;
  int i, j, n, ia, ic, ig, it;

  if( !isString(x) )
    error("argument must be a string");

  n = length(x);
  PROTECT(rv = allocVector(INTSXP, n*4));

  for(i=0; i<n; i++) {
    seq = CHAR(STRING_ELT(x, i));
    ia = ic = ig = it = 0;

    for(j=0; j<strlen(seq); j++) {
      switch(seq[j]) {
      case 'a': 
      case 'A':
	ia++;
	break;
      case 't':
      case 'T':
	it++;
	break;
      case 'c':
      case 'C':
	ic++;
	break;
      case 'g':
      case 'G':
	ig++;
	break;
      default:
	sprintf(errmess, "Unknown base %c in row %d, column %d.", seq[j], i+1, j+1);
	error(errmess);
      }
    }
    INTEGER(rv)[i    ] = ia; 
    INTEGER(rv)[i+n  ] = it; 
    INTEGER(rv)[i+n*2] = ic; 
    INTEGER(rv)[i+n*3] = ig; 
  }

  /* dim */
  PROTECT(dim = allocVector(INTSXP, 2));
  INTEGER(dim)[0] = n;
  INTEGER(dim)[1] = 4;
  setAttrib(rv, R_DimSymbol, dim);

  /* dim names */
  PROTECT(colnames = allocVector(STRSXP, 4));
  SET_STRING_ELT(colnames, 0, mkChar("A"));
  SET_STRING_ELT(colnames, 1, mkChar("T"));
  SET_STRING_ELT(colnames, 2, mkChar("C"));
  SET_STRING_ELT(colnames, 3, mkChar("G"));

  /* dim names */
  PROTECT(rownames = allocVector(STRSXP, n));
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, rownames);
  SET_VECTOR_ELT(dimnames, 1, colnames);
  setAttrib(rv, R_DimNamesSymbol, dimnames);

  UNPROTECT(5);
  return(rv);
}
