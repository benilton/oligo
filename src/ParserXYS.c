/***************************************************************
**  --------------------------------------
**    PARSER FOR XYS FILES - B Carvalho - 2009
**  --------------------------------------
**  
**  I) HISTORY
**  
**  April 28, 2009 - Initial working version
**                   (earlier versions didn't
**                   quite match R's results)
**  April 30, 2009 - Added dimnames to matrices
**                 - Added output list
**  May    7, 2009 - Added parser for header
**                   Added check for designname
**                   Calling this version 1.0
**  
**  II) TODO
**  None, for the moment
**  
**  III) CREDITS
**  Great help from Harris A. Jaffee <hj@jhu.edu>
**  countLines from ShortRead (Martin Morgan)
**  
***************************************************************/

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#define LINEMAX 20001
#define WORD 32

/***************************************************************
 ** countLines: counts lines, returns integer
***************************************************************/
static int countLines(FILE *file){
  size_t bytes_read;
  char buf[LINEMAX + 1];
  int lines = 0;

  while ((bytes_read = fread(buf, sizeof(char), LINEMAX, file)) != 0){
    char *p = buf;
    while ((p = memchr(p, '\n', (buf + bytes_read) - p))) {
      ++p;
      ++lines;
    }
    bytes_read += bytes_read;
  }
  return lines;
}

/***************************************************************
 ** Tokenset: stores pairs <key=value> and
 **           the number of such pairs.
***************************************************************/

typedef struct{
  char **key;
  char **value;
  int n;
} tokenset;


/***************************************************************
 ** tokenizer: stores <key=value> on a
 **            tokenset
***************************************************************/

static void tokenizer(tokenset *header, const char *key, const char *value){
  int i;
  i = header->n;
  header->n++;
  header->key = Realloc(header->key, header->n, char*);
  header->value = Realloc(header->value, header->n, char*);
  header->key[i] = Calloc(strlen(key)+1, char);
  header->value[i] = Calloc(strlen(value)+1, char);
  strcpy(header->key[i], key);
  strcpy(header->value[i], value);
}

/***************************************************************
 ** untokenizer: deallocates memory used
 **              by a tokenset
***************************************************************/

static void untokenizer(tokenset *header){
  int i;
  for (i=0; i < header->n; i++){
    Free(header->key[i]);
    Free(header->value[i]);
  }
  Free(header->key);
  Free(header->value);
  Free(header);
}

/***************************************************************
 ** buffer2tokenset: parses tags <key=value>
 **                  separated by tabs
 **                  as tokensets
***************************************************************/

static tokenset *buffer2tokenset(char *buffer){
  tokenset *header = Calloc(1, tokenset);
  char *eq, *key, *value;
  header->n = 0;
  header->key = NULL;
  header->value = NULL;

  // Header starts with "# "
  key = strtok(buffer, " ");
  // Keys are tab-delimited
  while (1){
    key = strtok(NULL, "\t");
    if (key == NULL) break;
    eq = strchr(key, '=');
    *eq = '\0';
    value = ++eq;
    tokenizer(header, key, value);
  }
  return header;
}

/***************************************************************
 ** xys_header_field: gets a field in the header of an XYS file.
 **
 ** NOTE: user must deallocate memory
 **       after using this function!!!!
***************************************************************/

static char *xys_header_field(const char *currentFile, const char *field){
  FILE *fp;
  int j;
  char *result, *final;
  char buffer[LINEMAX];

  fp = fopen(currentFile, "r");
  if (fp == NULL)
    error("Can't open %s.\n", currentFile);

  fgets(buffer, LINEMAX, fp);
  fclose(fp);

  j = strlen(buffer)-1;
  if (buffer[j] == '\n')
    buffer[j] = '\0';
  
  result = strstr(buffer, field);
  if (result == NULL)
    error("Can't find \'%s\' field. %s corrupted?", field, currentFile);
  result = strtok(result, "=");
  result = strtok(NULL, "\t");
  final = Calloc(strlen(result)+1, char);
  strcpy(final, result);
  return final;
}

/***************************************************************
 ** read_one_xys: reads XYS file. Stores XY coordinates on 'xy'
 **               and intensities on 'signal'
***************************************************************/

static void read_one_xys(const char *filename, double *signal,
			 int *xy, int i, int nrow, int verbose){
  int n, count;
  char sx[WORD], sy[WORD], ss[WORD], sc[WORD], *endpx, *endpy, *endps;
  FILE *fp;
  if (verbose) Rprintf("Reading %s.\n", filename);
  fp = fopen(filename, "r");
  if (fp == NULL)
    error("Can't open %s.\n", filename);

  // Header - 2 lines - skip
  while (fgetc(fp) != '\n');
  while (fgetc(fp) != '\n');
  count = 0;
  while (!feof(fp)){
    n = fscanf(fp, "%s %s %s %s", sx, sy, ss, sc);
    
    // Last line is a CR, so must match n == 4
    if (n == 4){
      // Read XY coordinates only for first file
      // and trust the rest...
      if (i == 0){
	xy[count] = strtol(sx, &endpx, 0);
	xy[count + nrow] = strtol(sy, &endpy, 0);
      }
      if (ss[0] != 'N'){
	signal[count + i*nrow] = strtod(ss, &endps);
      }else{ // NA_REAL is R-specific
	signal[count + i*nrow] = NA_REAL;
      }
    }

    count++;
  }
  fclose(fp);
}

/***************************************************************
 ** R_read_xys_files: reads XYS files.
 **    filenames: character vector with filenames
 **    verbosity: TRUE/FALSE flag
 **
 ** Function returns a list with two elements:
 **    coordinates: XY coordinates (N x 2 matrix)
 **    intensities: intensities (N x N_FILES matrix)
 **    date: date for each XYS file (character vector N_FILES)
***************************************************************/

SEXP R_read_xys_files(SEXP filenames, SEXP verbosity){
  int nfiles, nrows, i, verbose, *ptr2xy;
  double *ptr2signal;
  FILE *fp;
  SEXP signal, xy, output;
  SEXP dimnames, dimnamesxy, dimnamesout, fnames, colnamesxy, namesout, dates;
  char *d0, *d1;

  verbose = asLogical(verbosity);
  nfiles = length(filenames);
  fp = fopen(CHAR(STRING_ELT(filenames, 0)), "r");
  if (fp == NULL)
    error("Can't open %s.\n", CHAR(STRING_ELT(filenames, 0)));
  nrows = countLines(fp)-2;
  fclose(fp);

  // Test files are of the same type here
  if (verbose) Rprintf("Checking designs for each XYS file... ");
  d0 = xys_header_field(CHAR(STRING_ELT(filenames, 0)), "designname=");
  if (nfiles > 1)
    for (i = 1; i < nfiles; i++){
      d1 = xys_header_field(CHAR(STRING_ELT(filenames, i)), "designname=");
      if(strcasecmp(d1, d0) != 0){
	Free(d0);
	Free(d1);
	error("\'%s\' and \'%s\' use different designs.\n",
	      CHAR(STRING_ELT(filenames, 0)),
	      CHAR(STRING_ELT(filenames, i)));
      }
    }
  Free(d0);
  if (verbose) Rprintf("Done.\n");

  // Allocating memory in R
  if (verbose) Rprintf("Allocating memory... ");
  PROTECT(signal = allocMatrix(REALSXP, nrows, nfiles));
  PROTECT(xy = allocMatrix(INTSXP, nrows, 2));
  PROTECT(dates = allocVector(STRSXP, nfiles));
  if (verbose) Rprintf("Done.\n");
  ptr2xy = INTEGER_POINTER(xy);
  ptr2signal = NUMERIC_POINTER(signal);

  // Parsing files here
  for (i=0; i < nfiles; i++){
    read_one_xys(CHAR(STRING_ELT(filenames, i)), ptr2signal,
		 ptr2xy, i, nrows, verbose);
    d0 = xys_header_field(CHAR(STRING_ELT(filenames, i)), "date=");
    SET_STRING_ELT(dates, i, mkChar(d0));
  }
  Free(d0);

  PROTECT(output = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(output, 0, xy);
  SET_VECTOR_ELT(output, 1, signal);
  SET_VECTOR_ELT(output, 2, dates);

  // Dimnames +5 PROTECTs
  PROTECT(dimnames = allocVector(VECSXP, 2));
  PROTECT(fnames = allocVector(STRSXP, nfiles));
  for (i=0; i < nfiles; i++)
    SET_STRING_ELT(fnames, i, mkChar(CHAR(STRING_ELT(filenames, i))));
  SET_VECTOR_ELT(dimnames, 1, fnames);
  SET_VECTOR_ELT(dimnames, 0, R_NilValue);
  setAttrib(signal, R_DimNamesSymbol, dimnames);
  setAttrib(dates, R_NamesSymbol, fnames);

  PROTECT(dimnamesxy = allocVector(VECSXP, 2));
  PROTECT(colnamesxy = allocVector(STRSXP, 2));
  SET_STRING_ELT(colnamesxy, 0, mkChar("X"));
  SET_STRING_ELT(colnamesxy, 1, mkChar("Y"));
  SET_VECTOR_ELT(dimnamesxy, 0, R_NilValue);
  SET_VECTOR_ELT(dimnamesxy, 1, colnamesxy);
  setAttrib(xy, R_DimNamesSymbol, dimnamesxy);

  PROTECT(namesout = allocVector(STRSXP, 3));
  SET_STRING_ELT(namesout, 0, mkChar("coordinates"));
  SET_STRING_ELT(namesout, 1, mkChar("intensities"));
  SET_STRING_ELT(namesout, 2, mkChar("date"));
  setAttrib(output, R_NamesSymbol, namesout);

  UNPROTECT(9);
  return(output);
}


/***************************************************************
 ** R_read_xys_header: reads header of 1 XYS file.
 **    filename:  filenames
 **
 ** Function returns a list with several elements.
***************************************************************/

SEXP R_read_xys_header(SEXP filename){
  int j;
  FILE *fp;
  char buffer[LINEMAX];
  const char *currentFile;
  tokenset *header;
  SEXP output, namesout;

  currentFile = CHAR(STRING_ELT(filename, 0));
  fp = fopen(currentFile, "r");
  if (fp == NULL)
    error("Can't open %s.\n", currentFile);
  fgets(buffer, LINEMAX, fp);
  fclose(fp);

  j = strlen(buffer)-1;
  if (buffer[j] == '\n')
    buffer[j] = '\0';

  header = buffer2tokenset(buffer);
  PROTECT(output = allocVector(VECSXP, header->n));
  PROTECT(namesout = allocVector(STRSXP, header->n));
  for (j=0; j < header->n; j++){
    SET_VECTOR_ELT(output, j, mkString(header->value[j]));
    SET_STRING_ELT(namesout, j, mkChar(header->key[j]));
  }
  setAttrib(output, R_NamesSymbol, namesout);

  UNPROTECT(2);
  untokenizer(header);
  return(output);
}
