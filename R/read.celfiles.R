read.celfiles <- function(filenames,
                          pkgname=NULL,
                          phenoData=NULL,
                          featureData=NULL,
                          experimentData=NULL,
                          notes=NULL,
                          verbose = TRUE){

  n <- length(filenames)
  if (n == 0)
    stop("no file name given")

  chips <- sapply(filenames, function(x) readCelHeader(x)$chiptype)
  if (length(unique(chips)) > 1){
    print(table(chips))
    stop("All the CEL files must be of the same type.")
  }

  ## Array details
  headdetails <- readCelHeader(filenames[1])
  dim.intensity <- c(headdetails$rows, headdetails$cols)
  ref.cdfName <- chips[1]
    
  ## RI: We should check if the pd package is available here. If not try
  ## to install it
  ## load pdInfo
  if (is.null(pkgname))
    pkgname <- cleanPlatformName(ref.cdfName)
  if (require(pkgname, character.only=TRUE)){
    if (verbose)
      message("Platform design info loaded.")
  }else{
    stop("Must install the ", pkgname, " package.")
  }
  arrayType <- kind(get(pkgname))

  if (arrayType %in% c("SNP", "SNPCNV"))
    warning("SNP chips should be handled using justSNPRMA or justCRLMM.")

  if (verbose)
    message(sprintf("The intensity matrix will require %3.2f MB of RAM.",
                    prod(dim.intensity)*8/(1024^2)*length(filenames)))
    
  tmpExprs <- readCelIntensities(filenames)

  ## BC: Nov 15-16 2005, the PDenv is ordered already in the way
  ##     we want PM/MMs to be. So, we need to reorder the cel
  ##     input, so the pm/mm methods are faster.
  pdInfo <- get(pkgname)
  if (is(pdInfo, "platformDesign")) {
    order_index <- pdInfo$order_index
    tmpExprs <- tmpExprs[order_index, ]
  }

  rownames(tmpExprs) <- 1:nrow(tmpExprs)
  colnames(tmpExprs) <- basename(filenames)

  theClass <- switch(arrayType,
                     tiling="TilingFeatureSet",
                     expression="ExpressionFeatureSet",
                     SNP="SnpFeatureSet",
                     SNPCNV="SnpCnvFeatureSet",
                     exon="ExonFeatureSet",
                     stop("unknown array type: ", arrayType))

  out <- new(theClass,
             exprs=tmpExprs,
             platform=pkgname,
             manufacturer="Affymetrix",
             annotation=pkgname)

  if (is.null(featureData)){
    featureData(out) <- annotatedDataFrameFrom(assayData(out), byrow=TRUE)
  }else{
    featureData(out) <- featureData
  }
  if (is.null(phenoData)){
    phenoData(out) <- annotatedDataFrameFrom(assayData(out), byrow=FALSE)
  }else{
    phenoData(out) <- phenoData
  }
  if (is.null(experimentData)){
    ed <- new("MIAME")
    preproc(ed)$filenames <- filenames
    preproc(ed)$oligoversion <- packageDescription("oligo")$Version
    if (!is.null(notes)) notes(ed) <- notes
    experimentData(out) <- ed
  }else{
    experimentData(out) <- experimentData
  }
  out
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
    return(files[grep("\\.[cC][eE][lL]\\.gz$|\\.[cC][eE][lL]$", files)])
}

