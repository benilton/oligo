oligoReadXys <- function(cols, headdetails, filenames, out){
  if (length(cols) > 0){
    grpCols <- splitIndicesByLength(cols, ocSamples())
    open(out)
    for (theCols in grpCols)
      out[, theCols] <- .Call("R_read_xys_files", filenames[theCols], FALSE)[["intensities"]]
    close(out)
    rm(grpCols, out)
    gc()
  }
  TRUE
}


smartReadXYS <- function(filenames, sampleNames, verbose=TRUE){
  ## this runs on the master node
  if (isPackageLoaded("ff")){
    ## set filename here and location?

    ## wasting one read
    ## just to get nrow
    ## FIX ME
    nr <- nrow(.Call("R_read_xys_files", filenames[1],
                     FALSE)[["intensities"]])

    tmpExprs <- createFF("intensities-", dim=c(nr, length(filenames)),
                         dimnames=list(as.character(1:nr), sampleNames))
    intensityFile <- filename(tmpExprs)
    
    samplesByNode <- splitIndicesByNode(1:length(filenames))
    oLapply(samplesByNode, oligoReadXys, NULL, filenames, tmpExprs)
  }else{
    intensityFile <- NA_character_
    tmp <- .Call("R_read_xys_files", filenames, verbose)
    tmpExprs <- tmp[["intensities"]]
    datetime <- tmp[["date"]]
    rm(tmp)
    dimnames(tmpExprs) <- list(as.character(1:nrow(tmpExprs)), sampleNames)
  }
  return(list(intensityFile=intensityFile, exprMatrix=tmpExprs))
}



list.xysfiles <-   function(...){
  files <- list.files(...)
  return(files[grep("\\.[xX][yY][sS]\\.gz$|\\.[xX][yY][sS]$", files)])
}

readxysHeader <- function(filename) scan(filename,nlines=1,quiet=TRUE, what=character(0))

readonexysfile <- function(filename)
  read.delim(filename, comment.char="#")

readXysMatrix <- function(filenames){
  tmpExprs <- NULL
  for (i in seq(along=filenames)){
    ## Read XYS "as is"
    tmpE <- readonexysfile(filenames[i])
    if (length(filenames) > 1){
      tmpExprs <- cbind(tmpExprs, tmpE$SIGNAL)
    }else{
      tmpExprs <- matrix(tmpE$SIGNAL, ncol=1)
    }
  }
  list(intensities=tmpExprs, X=tmpE[["X"]], Y=tmpE[["Y"]])
}

read.xysfiles <- function(..., filenames, pkgname, phenoData,
                          featureData, experimentData, protocolData, notes,
                          verbose=TRUE, sampleNames, checkType=TRUE) {

  filenames <- getFilenames(filenames=filenames, ...)
  checkValidFilenames(filenames)
  if (checkType) stopifnot(checkChipTypes(filenames, verbose, "nimblegen"))
  if (!missing(sampleNames))
    stopifnot(length(sampleNames) == length(filenames))

  ## Get design name from the first
  firstline <- readxysHeader(filenames[1])
  designname <- unlist(strsplit(firstline[grep("designname",
                firstline, fixed=TRUE, useBytes=TRUE)], "="))[2]
    
  ## Load PDenv for the XYS files
  if (missing(pkgname))
    pkgname <- cleanPlatformName(designname)
  if (requireAnnotation(pkgname, verbose=verbose)){
    if (verbose)
      message("Platform design info loaded.")
  }else{
    stop("The annotation package, ", pkgname, ", could not be loaded.")
  }

  if (missing(sampleNames))
    sampleNames <- basename(filenames)

  results <- smartReadXYS(filenames, sampleNames)
  tmpExprs <- results[["exprMatrix"]]
  intensityFile <- results[["intensityFile"]]
  rm(results)
  
  arrayType <- kind(get(pkgname))
  theClass <- switch(arrayType,
                     tiling="TilingFeatureSet",
                     expression="ExpressionFeatureSet",
                     SNP="SnpFeatureSet",
                     SNPCNV="SnpCnvFeatureSet",
                     exon="ExonFeatureSet",
                     gene="GeneFeatureSet",
                     stop("Unknown array type: ", arrayType))


  out <- new(theClass)
  slot(out, "assayData") <- assayDataNew(exprs=tmpExprs)
  if (missing(phenoData))
    phenoData <- basicPhenoData(tmpExprs, filenames)
  slot(out, "phenoData") <- phenoData
  rm(phenoData)
  if (missing(featureData))
    featureData <- basicFeatureData(tmpExprs)
  slot(out, "featureData") <- featureData
  rm(featureData)
  if (missing(protocolData))
    protocolData <- basicProtocolData(tmpExprs)
  slot(out, "protocolData") <- protocolData
  rm(protocolData)
  slot(out, "manufacturer") <- "Nimblegen"
  slot(out, "annotation") <- pkgname
  slot(out, "intensityFile") <- intensityFile
  if (validObject(out)){
    return(out)
  }else{
    stop("Resulting object is invalid.")
  }
}


## For 2 channels - Tiling
read.xysfiles2 <- function(channel1, channel2, pkgname, phenoData,
                          featureData, experimentData, protocolData,
                          notes, verbose=TRUE, sampleNames,
                          checkType=TRUE) {
  
  filenames <- c(channel1, channel2)
  checkValidFilenames(filenames)
  if (checkType) stopifnot(checkChipTypes(filenames, verbose, "nimblegen"))
  if (!missing(sampleNames))
    stopifnot(length(sampleNames) == length(channel1), length(sampleNames) == length(channel2))

  ## Get design name from the first
  firstline <- readxysHeader(filenames[1])
  designname <- unlist(strsplit(firstline[grep("designname",
                firstline, fixed=TRUE, useBytes=TRUE)], "="))[2]
    
  ## Load PDenv for the XYS files
  if (missing(pkgname))
    pkgname <- cleanPlatformName(designname)
  if (requireAnnotation(pkgname, verbose=verbose)){
    if (verbose)
      message("Platform design info loaded.")
  }else{
    stop("Must install the ", pkgname, " package.")
  }

  arrayType <- kind(get(pkgname))

  if (missing(sampleNames))
    sampleNames <- basename(channel1)
  results <- smartReadXYS(channel1, sampleNames)
  channel1Intensities <- results[["exprMatrix"]]
  intensityFile1 <- results[["intensityFile"]]
  rm(results)
  results <- smartReadXYS(channel2, sampleNames)
  channel2Intensities <- results[["exprMatrix"]]
  intensityFile2 <- results[["intensityFile"]]
  rm(results)

  theClass <- "TilingFeatureSet"
  out <- new(theClass)
  slot(out, "assayData") <- assayDataNew(channel1=channel1Intensities,
                                         channel2=channel2Intensities)
  if (missing(phenoData))
    phenoData <- basicPhenoData2(channel1Intensities,
                                 channel2Intensities,
                                 channel1, channel2)
  slot(out, "phenoData") <- phenoData
  rm(phenoData)
  if (missing(featureData))
    featureData <- basicFeatureData(channel1Intensities)
  slot(out, "featureData") <- featureData
  rm(featureData)
  if (missing(protocolData))
    protocolData <- basicProtocolData(channel1Intensities)
  slot(out, "protocolData") <- protocolData
  rm(protocolData)
  slot(out, "manufacturer") <- "Nimblegen"
  slot(out, "annotation") <- pkgname
  slot(out, "intensityFile") <- c(intensityFile1, intensityFile2)
  if (validObject(out)){
    return(out)
  }else{
    stop("Resulting object is invalid.")
  }
}
