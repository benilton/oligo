read.celfiles <- function( ..., filenames, pkgname, phenoData,
                          featureData, experimentData, protocolData,
                          notes, verbose=TRUE, sampleNames,
                          rm.mask=FALSE, rm.outliers=FALSE,
                          rm.extra=FALSE, checkType=TRUE,
                          useAffyio=TRUE, intensityFile){
  ## TODO: remove useAffyio
  ##       add protocolData with scandate
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

  headdetails <- .Call("ReadHeader", as.character(filenames[1]),
                       PACKAGE="affyio")
  if (missing(sampleNames)){
    dns <- list(as.character(1:prod(headdetails[[2]])),
                basename(filenames))
  }else{
    dns <- list(as.character(1:prod(headdetails[[2]])),
                sampleNames)
  }
  
  if (isPackageLoaded("bigmemory")){
    stopifnot(!missing(intensityFile))
    tmpExprs <- filebacked.big.matrix(prod(headdetails[[2]]),
                                      length(filenames), type="double",
                                      backingfile=basename(intensityFile),
                                      backingpath=dirname(intensityFile),
                                      descriptorfile=paste(basename(intensityFile),
                                        "desc", sep="."), dimnames=dns)
    for (i in 1:length(filenames)){
      if (verbose) message("Reading ", filenames[i])
      tmpExprs[,i] <- read.celfile(filenames[i], TRUE)[["INTENSITY"]][["MEAN"]]
    }
    rm(i)
  }else{
    intensityFile <- NULL
    tmpExprs <- .Call("read_abatch", filenames, rm.mask, rm.outliers,
                      rm.extra, headdetails[[1]], headdetails[[2]],
                      verbose, PACKAGE="affyio")
    dimnames(tmpExprs) <- dns
  }
  rm(headdetails, dns)

  datetime <- GetAffyTimeDateAsString(filenames, useAffyio=useAffyio)

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

  ## assayData
  slot(out, "assayData") <- assayDataNew(exprs=tmpExprs)

  ## phenoData
  if (missing(phenoData)){
    pdd <- data.frame(exprs=filenames)
    vmd <- data.frame(labelDescription="Filenames",
                      channel=factor("exprs", levels=c("exprs", "_ALL_")))
    phenoData <- new("AnnotatedDataFrame", data=pdd, varMetadata=vmd)
    sampleNames(phenoData) <- colnames(tmpExprs)
    rm(pdd, vmd)
  }
  slot(out, "phenoData") <- phenoData
  rm(phenoData)
  
  ## featureData
  if (missing(featureData))
    featureData <- Biobase:::annotatedDataFrameFromMatrix(tmpExprs, byrow=TRUE)
  slot(out, "featureData") <- featureData
  rm(featureData)
  
  ## protocolData
  if (missing(protocolData))
    protocolData <- Biobase:::annotatedDataFrameFromMatrix(tmpExprs, FALSE)
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

## Now in oligoClasses
## list.celfiles <-   function(...){
##     files <- list.files(...)
##     return(files[grep("\\.[cC][eE][lL]\\.[gG][zZ]$|\\.[cC][eE][lL]$", files)])
## }

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
