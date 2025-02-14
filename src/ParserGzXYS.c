#define R_NO_REMAP
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <zlib.h>

#define LINEMAX 20001
#define WORD 32

/***************************************************************
 ** countLines: counts lines, returns integer
***************************************************************/
static int gzcountLines(gzFile file){
  int lines = 0;
  char buffer[1000];
  char *token;
  int fields = 3;

  while (!gzeof(file)){
    gzgets(file, buffer, 1000);
    token = strtok(buffer, " \t");
    if (token == NULL) fields--;
    token = strtok(NULL, " \t");
    if (token == NULL) fields--;
    token = strtok(NULL, " \t");
    if (token == NULL) fields--;
    if (fields == 3) lines++;
  }

  return lines;
}


/***************************************************************
 ** xys_header_field: gets a field in the header of an XYS file.
 **
 ** NOTE: user must deallocate memory
 **       after using this function!!!!
***************************************************************/

static char *gzxys_header_field(const char *currentFile, const char *field){
  gzFile fp;
  int j;
  char *result, *final;
  char buffer[LINEMAX];

  fp = gzopen(currentFile, "rb");
  if (fp == NULL)
    Rf_error("Can't open %s.\n", currentFile);

  gzgets(fp, buffer, LINEMAX);
  gzclose(fp);

  j = strlen(buffer)-1;
  if (buffer[j] == '\n')
    buffer[j] = '\0';

  result = strstr(buffer, field);
  if (result == NULL)
    Rf_error("Can't find \'%s\' field. %s corrupted?", field, currentFile);
  result = strtok(result, "=");
  result = strtok(NULL, "\t");
  final = R_Calloc(strlen(result)+1, char);
  strcpy(final, result);
  return final;
}

/***************************************************************
 ** read_one_xys: reads XYS file. Stores XY coordinates on 'xy'
 **               and intensities on 'signal'
***************************************************************/
static void gzread_one_xys(const char *filename, double *signal,
			   int *xy, int i, int nrow, int verbose){
  int count, n, x, y, j;
  char buffer[LINEMAX], sc[WORD], ss[WORD], *endps;
  gzFile fp=NULL;

  if (verbose) Rprintf("Reading %s.\n", filename);
  fp = gzopen(filename, "rb");
  if (fp == NULL)
    Rf_error("Can't open %s.\n", filename);

  // Header - 2 lines - skip
  while (gzgetc(fp) != '\n');
  while (gzgetc(fp) != '\n');
  count = 0;

  while (count < nrow){
    gzgets(fp, buffer, LINEMAX);

    j = strlen(buffer)-1;
    if (buffer[j] == '\n')
      buffer[j] = '\0';

    n = sscanf(buffer, "%d\t%d\t%s\t%s", &x, &y, ss, sc);

    // If it's the end of file, we're done.
    if (n == EOF) break;

    // Last line is a CR, so must match n == 4
    if (n == 4){
      // Read XY coordinates only for first file
      // and trust the rest...
      if (i == 0){
	xy[count] = x;
	xy[count + nrow] = y;
      }
      if (ss[0] != 'N'){
	signal[count + i*nrow] = strtod(ss, &endps);
      } else { // NA_REAL is R-specific
	signal[count + i*nrow] = NA_REAL;
      }
    } else {
      gzclose(fp);
      Rf_error("Line %d of %s has an unexpected format.\n", count, filename);
    }
    count++;
  }
  gzclose(fp);
  if (count != nrow)
    Rf_error("%s: Expected %d lines. Found %d lines. Corrupted file?",
	  filename, nrow, count);
}

SEXP R_read_gzxys_files(SEXP filenames, SEXP verbosity){
  int nfiles, nrows, i, verbose, *ptr2xy;
  double *ptr2signal;
  gzFile fp;
  SEXP signal, xy, output;
  SEXP dimnames, dimnamesxy, fnames, colnamesxy, namesout, dates;
  char *d0, *d1;

  verbose = Rf_asLogical(verbosity);
  nfiles = Rf_length(filenames);
  fp = gzopen(CHAR(STRING_ELT(filenames, 0)), "rb");
  if (fp == NULL)
    Rf_error("Can't open %s.\n", CHAR(STRING_ELT(filenames, 0)));
  nrows = gzcountLines(fp)-2;
  gzclose(fp);

  // Test files are of the same type here
  if (verbose) Rprintf("Checking designs for each XYS file... ");
  d0 = gzxys_header_field(CHAR(STRING_ELT(filenames, 0)), "designname=");
  if (nfiles > 1)
    for (i = 1; i < nfiles; i++){
      d1 = gzxys_header_field(CHAR(STRING_ELT(filenames, i)), "designname=");
      if(strcasecmp(d1, d0) != 0){
	R_Free(d0);
	R_Free(d1);
	Rf_error("\'%s\' and \'%s\' use different designs.\n",
	      CHAR(STRING_ELT(filenames, 0)),
	      CHAR(STRING_ELT(filenames, i)));
      }
      R_Free(d1); // Missed: 12/02/09
    }
  R_Free(d0);
  if (verbose) Rprintf("Done.\n");

  // Allocating memory in R
  if (verbose) Rprintf("Allocating memory... ");
  PROTECT(signal = Rf_allocMatrix(REALSXP, nrows, nfiles));
  PROTECT(xy = Rf_allocMatrix(INTSXP, nrows, 2));
  PROTECT(dates = Rf_allocVector(STRSXP, nfiles));
  if (verbose) Rprintf("Done.\n");
  ptr2xy = INTEGER_POINTER(xy);
  ptr2signal = NUMERIC_POINTER(signal);

  // Parsing files here
  for (i=0; i < nfiles; i++){
    gzread_one_xys(CHAR(STRING_ELT(filenames, i)), ptr2signal,
		 ptr2xy, i, nrows, verbose);
    d0 = gzxys_header_field(CHAR(STRING_ELT(filenames, i)), "date=");
    SET_STRING_ELT(dates, i, Rf_mkChar(d0));
    R_Free(d0);
  }

  PROTECT(output = Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(output, 0, xy);
  SET_VECTOR_ELT(output, 1, signal);
  SET_VECTOR_ELT(output, 2, dates);

  // Dimnames +5 PROTECTs
  PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
  PROTECT(fnames = Rf_allocVector(STRSXP, nfiles));
  for (i=0; i < nfiles; i++)
    SET_STRING_ELT(fnames, i, Rf_mkChar(CHAR(STRING_ELT(filenames, i))));
  SET_VECTOR_ELT(dimnames, 1, fnames);
  SET_VECTOR_ELT(dimnames, 0, R_NilValue);
  Rf_setAttrib(signal, R_DimNamesSymbol, dimnames);
  Rf_setAttrib(dates, R_NamesSymbol, fnames);

  PROTECT(dimnamesxy = Rf_allocVector(VECSXP, 2));
  PROTECT(colnamesxy = Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(colnamesxy, 0, Rf_mkChar("X"));
  SET_STRING_ELT(colnamesxy, 1, Rf_mkChar("Y"));
  SET_VECTOR_ELT(dimnamesxy, 0, R_NilValue);
  SET_VECTOR_ELT(dimnamesxy, 1, colnamesxy);
  Rf_setAttrib(xy, R_DimNamesSymbol, dimnamesxy);

  PROTECT(namesout = Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(namesout, 0, Rf_mkChar("coordinates"));
  SET_STRING_ELT(namesout, 1, Rf_mkChar("intensities"));
  SET_STRING_ELT(namesout, 2, Rf_mkChar("date"));
  Rf_setAttrib(output, R_NamesSymbol, namesout);

  UNPROTECT(9);
  return(output);
}
