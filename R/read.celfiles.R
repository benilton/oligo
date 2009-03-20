read.celfiles <- function( ..., filenames, pkgname, phenoData,
                          featureData, experimentData, notes,
                          verbose=TRUE, sampleNames, rm.mask=FALSE,
                          rm.outliers=FALSE, rm.extra=FALSE,
                          sd=FALSE, checkType=TRUE){

  if (!missing(filenames)){
    filenames <- c(filenames, unlist(list(...)))
  }else{
    filenames <- unlist(list(...))
  }
  checkValidFilenames(filenames)
  if (checkType) stopifnot(checkChipTypes(filenames, verbose, "affymetrix"))
  
  ## Read in the first Array details
  headdetails <- readCelHeader(filenames[1])
  chiptype <- headdetails[["chiptype"]]
  
  if (missing(pkgname))
    pkgname <- cleanPlatformName(chiptype)
  
  if (require(pkgname, character.only=TRUE)){
    if (verbose)
      message("Platform design info loaded.")
  }else{
    stop("Must install the ", pkgname, " package.")
  }
  
  arrayType <- kind(get(pkgname))
  tmpExprs <- readCelIntensities2(filenames, rm.outliers=rm.outliers,
                                  rm.masked=rm.mask,
                                  rm.extra=rm.extra, verbose=verbose)
  dimnames(tmpExprs) <- NULL

  metadata <- getMetadata(tmpExprs, filenames, phenoData, featureData,
                          experimentData, notes, sampleNames)
  
  if (sd) warning("Reading in Standard Errors not yet implemented.\n")
  theClass <- switch(arrayType,
                     tiling="TilingFeatureSet",
                     expression="ExpressionFeatureSet",
                     SNP="SnpFeatureSet",
                     SNPCNV="SnpCnvFeatureSet",
                     exon="ExonFeatureSet",
                     gene="GeneFeatureSet",
                     stop("Unknown array type: ", arrayType))
  
  return(new(theClass, exprs=tmpExprs, manufacturer="Affymetrix",
             annotation=pkgname, phenoData=metadata[["phenoData"]],
             experimentData=metadata[["experimentData"]],
             featureData=metadata[["featureData"]]))
  
}

list.celfiles <-   function(...){
    files <- list.files(...)
    return(files[grep("\\.[cC][eE][lL]$", files)])
}

## reimplementation of readCelIntensites from affxparser, allows for rm.mask, 
## rm.outliers and rm.extra arguments as used in the affy package

readCelIntensities2 <- function(filenames, indices=NULL,
                                rm.masked=FALSE, rm.outliers=FALSE,
                                rm.extra=FALSE, verbose=0){
	if(rm.extra) rm.outliers=TRUE;rm.masked=TRUE;
	if (length(filenames) == 0) 
		stop("Argument 'filenames' is empty.")
	filenames <- file.path(dirname(filenames), basename(filenames))
	missing <- !file.exists(filenames)
	if (any(missing)) {
		missing <- paste(filenames[missing], collapse = ", ")
		stop("Cannot read CEL files. Some files not found: ", 
				missing)
	}
	if (length(verbose) != 1) {
		stop("Argument 'verbose' must be a single integer.")
	}
	if (!is.finite(as.integer(verbose))) {
		stop("Argument 'verbose' must be an integer: ", verbose)
	}
	verbose <- as.integer(verbose)
	if (verbose > 1) {
		cat("Entering readCelIntensities2()\n ... reading headers\n")
	}
	all.headers <- lapply(as.list(filenames), readCelHeader)
	chiptype <- unique(sapply(all.headers, function(x) x$chiptype))
	if (length(chiptype) != 1) {
		warning("The CEL files do not have the same chiptype.")
	}
	nrows <- unique(sapply(all.headers, function(x) x$rows))
	ncols <- unique(sapply(all.headers, function(x) x$cols))
	if (length(nrows) != 1 || length(ncols) != 1) {
		stop("The CEL files dimension do not match.")
	}
	nfiles <- length(filenames)
	if (verbose > 0) {
		cat(" ... allocating memory for intensity matrix\n")
	}
	if (is.null(indices)) {
		intensities <- matrix(NA, nrow = nrows * ncols, ncol = nfiles)
	}
	else {
		intensities <- matrix(NA, nrow = length(indices), ncol = nfiles)
	}
	colnames(intensities) <- filenames
	for (i in 1:nfiles) {
		if (verbose > 0) 
			cat(" ... reading", filenames[i], "\n")
		tmp <- readCel(filename = filenames[i], 
				indices = indices, readIntensities = TRUE, readHeader = FALSE, 
				readStdvs = FALSE, readPixels = FALSE, readXY = FALSE, 
				readOutliers = rm.outliers, readMasked = rm.masked, verbose = (verbose - 1))
		if(rm.outliers || rm.masked) tmp$intensities[c(tmp$outliers,tmp$masked)] <- NA
		intensities[,i] <- tmp$intensities
	}
	intensities
}

