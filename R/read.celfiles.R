read.celfiles <- function( ..., filenames, pkgname, phenoData,
                          featureData, experimentData, notes,
                          verbose=TRUE, sampleNames, rm.mask=FALSE,
                          rm.outliers=FALSE, rm.extra=FALSE,
                          sd=FALSE, checkType=TRUE){

  filenames <- c(filenames, unlist(list(...)))
  checkValidFilenames(filenames)
  if (checkType) stopifnot(checkChipTypes(filenames, verbose))
  
  ## Read in the first Array details
  headdetails <- readCelHeader(filenames[1])
  chiptype <- chips[1]
  
  if (missing(pkgname))
    pkgname <- cleanPlatformName(chiptype)
  
  if (require(pkgname, character.only=TRUE)){
    if (is(get(pkgname), "platformDesign"))
      stop("Create a pdInfo package using the 'pdInfoBuilder' package")    
    if (verbose)
      message("Platform design info loaded.")
  }else{
    stop("Must install the ", pkgname, " package.")
  }
  
  arrayType <- kind(get(pkgname))
  tmpExprs <- readCelIntensities2(filenames, rm.outliers=rm.outliers,
                                  rm.masked=rm.masked,
                                  rm.extra=rm.extra, verbose=verbose)
  dimnames(tmpExprs) <- NULL
  metadata <- getMetadataForMatrix(filenames, phenoData, featureData, experimentData, notes, sampleNames)
  
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

## made read.celfiles look more similar to read.affybatch in affy --MS

read.celfiles.ms <- function( ..., filenames = character(0),
                          pkgname,
                          phenoData=new("AnnotatedDataFrame"),
                          featureData, experimentData=NULL, notes="",
                          verbose=TRUE, sampleNames=NULL,
                          rm.mask=FALSE, rm.outliers=FALSE,
                          rm.extra=FALSE, sd=FALSE){

  filenames <- c(filenames, unlist(list(...)))
  checkValidFilenames(filenames)
  n <- length(filenames)
  
  pdata <- pData(phenoData)
  ## try to read sample names form phenoData. if not there use CEL
  ## filenames
  if(dim(pdata)[1] != n) {
    ## if empty pdata filename are samplenames
    warning("Incompatible phenoData object. Created a new one.\n")
		
    samplenames <- sub("^/?([^/]*/)*", "", filenames, extended=TRUE)
    pdata <- data.frame(sample=1:n, row.names=samplenames)
    phenoData <- new("AnnotatedDataFrame",
                     data=pdata,
                     varMetadata=data.frame(
                       labelDescription="arbitrary numbering",
                       row.names="sample"))
  } else samplenames <- rownames(pdata)
  if (!is.null(samplenames)) samplenames <- sampleNames
  if (is.null(experimentData))
    {
      experimentData <- new("MIAME")
      preproc(experimentData)$filenames <- filenames
      preproc(experimentData)$oligoversion <- packageDescription("oligo")$Version
    }
  if (length(notes)!=0) notes(experimentData) <- notes
  
  chips <- sapply(filenames, function(x) readCelHeader(x)$chiptype)
  if (length(unique(chips)) > 1){
    print(table(chips))
    stop("All the CEL files must be of the same type.")
  }

  ## Read in the first Array details
  headdetails <- readCelHeader(filenames[1])
  ##dim.intensity <- c(headdetails$rows, headdetails$cols)
  chiptype <- chips[1]
    
  if (missing(pkgname))
    pkgname <- cleanPlatformName(chiptype)
  
  if (require(pkgname, character.only=TRUE)){
    if (verbose)
      message("Platform design info loaded.")
  } else {
    stop("Must install the ", pkgname, " package.")
  }
  if (is(get(pkgname), "platformDesign"))
    stop("Create a pdInfo package using the 'pdInfoBuilder' package")    
  
  arrayType <- kind(get(pkgname))
  
  tmpExprs <- readCelIntensities2(filenames,rm.outliers=rm.outliers,rm.masked=rm.masked,rm.extra=rm.extra,verbose=verbose)
  colnames(tmpExprs) <- samplenames
  rownames(tmpExprs) <- 1:nrow(tmpExprs)
  
  if (missing(featureData)) featureData <- annotatedDataFrameFrom(tmpExprs,byrow=TRUE)
  if (sd) warning("Reading in Standard Errors not yet impletmented.\n")
  theClass <- switch(arrayType,
                     tiling="TilingFeatureSet",
                     expression="ExpressionFeatureSet",
                     SNP="SnpFeatureSet",
                     SNPCNV="SnpCnvFeatureSet",
                     exon="ExonFeatureSet",
                     gene="GeneFeatureSet",
                     stop("unknown array type: ", arrayType))
  return(new(theClass,
             exprs  = tmpExprs,
             platform=pkgname, ## why do we have this why not just use annotation
             manufacturer="Affymetrix",
             annotation=pkgname,
             phenoData  = phenoData,
             experimentData = experimentData,
             featureData = featureData))
	 ##se.exprs = array(NaN, dim=dim.sd),
	 ##cdfName    = cdfname,   ##cel@cdfName,
}

## read.celfiles.old <- function(filenames,
##                           pkgname=NULL,
##                           phenoData=new("AnnotatedDataFrame"),
##                           featureData=NULL,
##                           description=NULL,
##                           notes="",
##                           verbose = FALSE,
##                           rm.mask = FALSE,
##                           rm.outliers=FALSE,
##                           rm.extra=FALSE,
##                           tmpdir=getwd(),
##                           memory.bound=FALSE){
## 
##   ## FIXME: this is not an acceptable function name
##   tmp <- stuffForXYSandCELreaders(filenames, phenoData, description, notes,
##                                   verbose)
##   filenames <- tmp$filenames
##   n <- length(filenames)
##   if (n == 0)
##     stop("no file name given")
## 
##   chips <- sapply(filenames, function(x) readCelHeader(x)$chiptype)
##   if (length(unique(chips)) > 1){
##     print(table(chips))
##     stop("All the CEL files must be of the same type.")
##   }
##   headdetails <- readCelHeader(filenames[1])
##   dim.intensity <- c(headdetails$rows, headdetails$cols)
##   ref.cdfName <- chips[1]
## 
##   
##   ## read the first file to see what we have
## ##  headdetails <- read.celfile.header(filenames[1])
##   
##   ##now we use the length
## ##  dim.intensity <- headdetails[[2]]
##   
##   ##and the cdfname as ref
## ##  ref.cdfName <- headdetails[[1]]
## 
##   if (!memory.bound){
## ##    tmpExprs <- .Call("read_abatch", filenames, rm.mask, rm.outliers,
## ##                      rm.extra, ref.cdfName, dim.intensity, verbose, PACKAGE="affyio")
## 
##     tmpExprs <- readCelIntensities(filenames)
##     
##   }else{
##     tmpExprs <- createBufferedMatrix(prod(dim.intensity), 0, directory=tmpdir)
##     set.buffer.dim(tmpExprs, 100, 1)
##     for (i in 1:length(filenames)){
##       AddColumn(tmpExprs)
## ##       tmpExprs[,i] <- .Call("read_abatch", filenames[i], FALSE, FALSE, FALSE,
## ##                             headdetails$cdfName, headdetails[["CEL dimensions"]],
## ##                             FALSE, PACKAGE="affyio")
## 
##       tmpExprs[, i] <- readCel(filenames[i], readHeader=FALSE, readOutliers=FALSE, readMasked = FALSE)$intensities
##     }
##   }
##     
##   ## RI: We should check if the pd package is available here. If not try
##   ## to install it
##   ## load pdInfo
##   if (is.null(pkgname))
##     pkgname <- cleanPlatformName(ref.cdfName)
##   if (require(pkgname, character.only=TRUE)){
##     message("Platform design info loaded.")
##   }else{
##     stop("Must install the ", pkgname, " package.")
##   }
##   arrayType <- kind(get(pkgname))
## 
##   ## BC: Nov 15-16 2005, the PDenv is ordered already in the way
##   ##     we want PM/MMs to be. So, we need to reorder the cel
##   ##     input, so the pm/mm methods are faster.
##   pdInfo <- get(pkgname)
##   if (is(pdInfo, "platformDesign")) {
##     order_index <- pdInfo$order_index
##     for (i in 1:ncol(tmpExprs))
##       tmpExprs[,i] <- tmpExprs[order_index, i]
##   }
## 
##   rownames(tmpExprs) <- 1:nrow(tmpExprs)
##   colnames(tmpExprs) <- basename(filenames)
## 
##   theClass <- switch(arrayType,
##                      tiling="TilingFeatureSet",
##                      expression="ExpressionFeatureSet",
##                      SNP="SnpFeatureSet",
##                      SNPCNV="SnpCnvFeatureSet",
##                      exon="ExonFeatureSet",
##                      stop("unknown array type: ", arrayType))
## 
##   if (theClass %in% c("SnpFeatureSet", "SnpCnvFeatureSet") & is.matrix(tmpExprs))
##     tmpExprs <- as.BufferedMatrix(tmpExprs, 5000, 1, directory=tmpdir)
## 
##   if (is.null(featureData))
##     featureData <- new("AnnotatedDataFrame",
##                        data=data.frame(idx=1:nrow(tmpExprs),
##                          row.names=1:nrow(tmpExprs)),
##                        varMetadata=data.frame(varLabels="idx",
##                          row.names="idx"))
##   out <- new(theClass,
##              exprs=tmpExprs,
##              platform=pkgname,
##              manufacturer="Affymetrix",
##              phenoData=tmp$phenoData,
##              featureData=featureData,
##              experimentData=tmp$description,
##              annotation=pkgname)
##   out
## }

list.celfiles <-   function(...){
    files <- list.files(...)
    return(files[grep("\\.[cC][eE][lL]$", files)])
}

## reimplementation of readCelIntensites from affxparser, allows for rm.mask, 
## rm.outliers and rm.extra arguments as used in the affy package
readCelIntensities2 <- 
function (filenames, indices = NULL,rm.masked=FALSE,rm.outliers=FALSE,rm.extra=FALSE, verbose = 0) 
{
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

