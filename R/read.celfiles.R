smartReadCEL <- function(filenames, sampleNames, headdetails,
                         verbose=TRUE){

  nr <- prod(headdetails[[2]])
  dns <- list(as.character(1:nr), sampleNames)
  
  if (isPackageLoaded("ff")){
    tmpExprs <- createFF("intensities-", dim=c(nr, length(filenames)))
    intensityFile <- filename(tmpExprs)
    samplesByNode <- splitIndicesByNode(1:length(filenames))
    datetime <- ocLapply(samplesByNode, oligoReadCels, headdetails,
                         filenames, tmpExprs, neededPkgs="oligo")
    datetime <- unlist(datetime)
  }else{
    intensityFile <- NA_character_
    tmpExprs <- .Call("read_abatch", filenames, FALSE, FALSE, FALSE,
                      headdetails[[1]], headdetails[[2]], verbose,
                      PACKAGE="affyio")
    datetime <- sapply(filenames, getCelDateTime)
  }
  dimnames(tmpExprs) <- dns
  rm(headdetails, dns)
  return(list(exprMatrix=tmpExprs, intensityFile=intensityFile, datetime=datetime))
}

getCelDateTime <- function(celFile)
    read.celfile.header(celFile, info="full")[["DatHeader"]]

oligoReadCels <- function(cols, headdetails, filenames, out){
  ## runs on the nodes
  if (length(cols) > 0){
    grpCols <- splitIndicesByLength(cols, ocSamples())
    dates <- vector("list", length(grpCols))
    open(out)
    i <- 1
    for (theCols in grpCols){
      out[, theCols] <- .Call("read_abatch", filenames[theCols], FALSE,
                              FALSE, FALSE, headdetails[[1]],
                              headdetails[[2]], FALSE, PACKAGE="affyio")
      dates[[i]] <- sapply(filenames(theCols), getCelDateTime)
      i <- i + 1
  }
    close(out)
    rm(grpCols, out)
    gc()
    return(unlist(dates))
  }
  TRUE
}

read.celfiles <- function( ..., filenames, pkgname, phenoData,
                          featureData, experimentData, protocolData,
                          notes, verbose=TRUE, sampleNames,
                          rm.mask=FALSE, rm.outliers=FALSE,
                          rm.extra=FALSE, checkType=TRUE){
  ##       add protocolData with scandate
  filenames <- getFilenames(filenames=filenames, ...)
  checkValidFilenames(filenames)
  if (checkType)
    stopifnot(checkChipTypes(filenames, verbose, "affymetrix",
                             TRUE))
  
  ## Read in the first Array details
  chiptype <- getCelChipType(filenames[1], TRUE)
  
  if (missing(pkgname))
    pkgname <- cleanPlatformName(chiptype)
  
  if (requireAnnotation(pkgname, verbose=verbose)){
    if (verbose)
      message("Platform design info loaded.")
  }else{
    stop("The annotation package, ", pkgname, ", could not be loaded.")
  }

  headdetails <- .Call("ReadHeader", as.character(filenames[1]),
                       PACKAGE="affyio")

  if (missing(sampleNames))
    sampleNames <- basename(filenames)

  results <- smartReadCEL(filenames, sampleNames, headdetails=headdetails)
  tmpExprs <- results[["exprMatrix"]]
  intensityFile <- results[["intensityFile"]]
  datetime <- results[["datetime"]]
  rm(results)
##  datetime <- GetAffyTimeDateAsString(filenames, useAffyio=TRUE)

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
      phenoData <- basicPhData1(tmpExprs, filenames)
  slot(out, "phenoData") <- phenoData
  rm(phenoData)
  if (missing(featureData))
      featureData <- basicAnnotatedDataFrame(tmpExprs, TRUE)
  slot(out, "featureData") <- featureData
  rm(featureData)
  if (missing(protocolData))
      protocolData <- basicPData(tmpExprs, filenames, datetime)
  slot(out, "protocolData") <- protocolData
  rm(protocolData)
  slot(out, "manufacturer") <- "Affymetrix"
  slot(out, "annotation") <- pkgname
  slot(out, "intensityFile") <- intensityFile
  if (validObject(out)){
    return(out)
  }else{
    stop("Resulting object is invalid.")
  }

}

## TilingFeatureSet2

read.celfiles2 <- function(channel1, channel2, pkgname, phenoData,
                           featureData, experimentData, protocolData, notes,
                           verbose=TRUE, sampleNames, rm.mask=FALSE,
                           rm.outliers=FALSE, rm.extra=FALSE,
                           checkType=TRUE){

  filenames <- c(channel1, channel2)
  checkValidFilenames(filenames)
  if (checkType) stopifnot(checkChipTypes(filenames, verbose, "affymetrix", TRUE))
  
  ## Read in the first Array details
  headdetails <- readCelHeader(filenames[1])
  chiptype <- headdetails[["chiptype"]]
  
  if (missing(pkgname))
    pkgname <- cleanPlatformName(chiptype)
  
  if (requireAnnotation(pkgname, verbose=verbose)){
    if (verbose)
      message("Platform design info loaded.")
  }else{
    stop("Must install the ", pkgname, " package.")
  }
  
  arrayType <- kind(get(pkgname))
  if (missing(sampleNames))
    sampleNames <- basename(channel1)

  results <- smartReadCEL(channel1, sampleNames, headdetails=headdetails)
  channel1Intensities <- results[["exprMatrix"]]
  intensityFile1 <- results[["intensityFile"]]
  datetime1 <- results[["datetime"]]
  rm(results)
  results <- smartReadCEL(channel2, sampleNames, headdetails=headdetails)
  channel2Intensities <- results[["exprMatrix"]]
  intensityFile2 <- results[["intensityFile"]]
  datetime2 <- results[["datetime"]]
  rm(results, headdetails)

  theClass <- "TilingFeatureSet"
  out <- new(theClass)
  slot(out, "assayData") <- assayDataNew(channel1=channel1Intensities,
                                         channel2=channel2Intensities)
  if (missing(phenoData))
      phenoData <- basicPhData2(channel1Intensities,
                                channel2Intensities,
                                channel1, channel2)
  slot(out, "phenoData") <- phenoData
  rm(phenoData)
  if (missing(featureData))
      featureData <- basicAnnotatedDataFrame(channel1Intensities, TRUE)
  slot(out, "featureData") <- featureData
  rm(featureData)
  if (missing(protocolData))
    protocolData <- basicPData2(channel1Intensities,
                                channel2Intensities,
                                channel1, channel2,
                                datetime1, datetime2)
  slot(out, "protocolData") <- protocolData
  rm(protocolData)
  slot(out, "manufacturer") <- "Affymetrix"
  slot(out, "annotation") <- pkgname
  slot(out, "intensityFile") <- c(intensityFile1, intensityFile2)
  if (validObject(out)){
    return(out)
  }else{
    stop("Resulting object is invalid.")
  }
}
