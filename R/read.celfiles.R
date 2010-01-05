read.celfiles <- function( ..., filenames, pkgname, phenoData,
                          featureData, experimentData, notes,
                          verbose=TRUE, sampleNames, rm.mask=FALSE,
                          rm.outliers=FALSE, rm.extra=FALSE,
                          sd=FALSE, checkType=TRUE, useAffyio=TRUE){
  
  filenames <- getFilenames(filenames=filenames, ...)
  checkValidFilenames(filenames)
  if (checkType)
    stopifnot(checkChipTypes(filenames, verbose, "affymetrix", useAffyio))
  
  ## Read in the first Array details
  chiptype <- getCelChipType(filenames[1], useAffyio)
  
  if (missing(pkgname))
    pkgname <- cleanPlatformName(chiptype)
  
  if (requireAnnotation(pkgname, verbose=verbose)){
    if (verbose)
      message("Platform design info loaded.")
  }else{
    stop("The annotation package, ", pkgname, ", could not be loaded.")
  }
  
  arrayType <- kind(get(pkgname))
  if (useAffyio){
    headdetails <- .Call("ReadHeader", as.character(filenames[1]),
                         PACKAGE="affyio")
    tmpExprs <- .Call("read_abatch", filenames, rm.mask, rm.outliers,
                      rm.extra, headdetails[[1]], headdetails[[2]],
                      verbose, PACKAGE="affyio")
    rm(headdetails)
  }else{
    tmpExprs <- readCelIntensities2(filenames,
                                    rm.outliers=rm.outliers,
                                    rm.masked=rm.mask,
                                    rm.extra=rm.extra,
                                    verbose=verbose)
  }
  datetime <- GetAffyTimeDateAsString(filenames, useAffyio=useAffyio)

  metadata <- getMetadata(tmpExprs, filenames, phenoData, featureData,
                          experimentData, notes, sampleNames, AffyDate2Posix(datetime))
  colnames(tmpExprs) <- Biobase::sampleNames(metadata[["phenoData"]])

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
    return(files[grep("\\.[cC][eE][lL]\\.[gG][zZ]$|\\.[cC][eE][lL]$", files)])
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

## TilingFeatureSet2

read.celfiles2 <- function(channel1, channel2, pkgname, phenoData,
                          featureData, experimentData, notes,
                          verbose=TRUE, sampleNames, rm.mask=FALSE,
                          rm.outliers=FALSE, rm.extra=FALSE,
                          sd=FALSE, checkType=TRUE, useAffyio=TRUE){

  filenames <- c(channel1, channel2)
  checkValidFilenames(filenames)
  if (checkType) stopifnot(checkChipTypes(filenames, verbose, "affymetrix", useAffyio))
  
  ## Read in the first Array details
##   headdetails <- readCelHeader(filenames[1])
##   chiptype <- headdetails[["chiptype"]]
  chiptype <- getCelChipType(filenames[1], useAffyio)
  
  if (missing(pkgname))
    pkgname <- cleanPlatformName(chiptype)
  
  if (requireAnnotation(pkgname, verbose=verbose)){
    if (verbose)
      message("Platform design info loaded.")
  }else{
    stop("Must install the ", pkgname, " package.")
  }
  
  arrayType <- kind(get(pkgname))

  if (useAffyio){
    headdetails <- .Call("ReadHeader", as.character(filenames[1]),
                         PACKAGE="affyio")
    channel1Intensities <- .Call("read_abatch", channel1, rm.mask,
                      rm.outliers, rm.extra, headdetails[[1]],
                      headdetails[[2]], verbose, PACKAGE="affyio")
    channel2Intensities <- .Call("read_abatch", channel2, rm.mask,
                      rm.outliers, rm.extra, headdetails[[1]],
                      headdetails[[2]], verbose, PACKAGE="affyio")
    rm(headdetails)
  }else{
    channel1Intensities <- readCelIntensities2(channel1,
                                  rm.outliers=rm.outliers,
                                  rm.masked=rm.mask,
                                  rm.extra=rm.extra, verbose=verbose)
  
    channel2Intensities <- readCelIntensities2(channel2,
                                  rm.outliers=rm.outliers,
                                  rm.masked=rm.mask,
                                  rm.extra=rm.extra, verbose=verbose)
  }
  dimnames(channel1Intensities) <- NULL
  dimnames(channel2Intensities) <- NULL
  date1 <- GetAffyTimeDateAsString(channel1, useAffyio=useAffyio)
  date2 <- GetAffyTimeDateAsString(channel2, useAffyio=useAffyio)
  date1 <- AffyDate2Posix(date1)
  date2 <- AffyDate2Posix(date2)

  metadata <- getMetadata2(channel1Intensities, channel2Intensities,
                           channel1, channel2,
                           phenoData, featureData, experimentData, notes, sampleNames,
                           date1, date2)
  colnames(channel1Intensities) <- colnames(channel2Intensities) <- Biobase::sampleNames(metadata[["phenoData"]])

  out <- new("TilingFeatureSet",
             channel1=channel1Intensities,
             channel2=channel2Intensities,
             manufacturer="Affymetrix",
             annotation=pkgname,
             phenoData=metadata[["phenoData"]],
             experimentData=metadata[["experimentData"]],
             featureData=metadata[["featureData"]])
  return(out)
}
