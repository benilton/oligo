oligoReadXys <- function(cols, headdetails, filenames, desc, path=oligoBigObjectPath()){
  if (length(cols) > 0){
    grpCols <- splitIndicesByLength(cols, oligoSamples())
    out <- attach.big.matrix(desc, backingpath=path)
    for (theCols in grpCols)
      out[, theCols] <- .Call("R_read_xys_files", filenames[theCols], FALSE)[["intensities"]]
    rm(grpCols, out)
  }
  TRUE
}


smartReadXYS <- function(filenames, sampleNames, prefix="intensities-", path=oligoBigObjectPath(), uid=genDatasetUID(filenames), verbose=TRUE){
  ## this runs on the master node
  if (isPackageLoaded("bigmemory")){
    bn <- paste(prefix, uid, sep="")
    intensityFile <- file.path(path, bn)

    ## missing one read
    ## just to get nrow
    ## FIX ME
    tmp <- .Call("R_read_xys_files", filenames[1], FALSE)[["intensities"]]
    nr <- nrow(tmp)
    rm(tmp)
    dns <- list(as.character(1:nr), sampleNames)
    tmpExprs <- big.matrix(nr, length(filenames), type="double",
                           backingfile=basename(intensityFile),
                           backingpath=path,
                           descriptorfile=paste(basename(intensityFile),
                           "desc", sep="."), dimnames=dns,
                           separated=TRUE)
    samplesByNode <- splitIndicesByNode(1:length(filenames))
    oLapply(samplesByNode, oligoReadXys, NULL, filenames,
            describe(tmpExprs), path=path)
  }else{
    intensityFile <- NA_character_
    tmp <- .Call("R_read_xys_files", filenames, verbose)
    tmpExprs <- tmp[["intensities"]]
    datetime <- tmp[["date"]]
    rm(tmp)
    dimnames(tmpExprs) <- list(as.character(1:nrow(tmpExprs)), sampleNames)
  }
  return(list(exprMatrix=tmpExprs, intensityFile=intensityFile))
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


## stuffForXYSandCELreaders <- function(filenames,
##                                      phenoData=new("AnnotatedDataFrame"),
##                                      description=NULL,
##                                      notes="",
##                                      verbose = FALSE,
##                                      nwells = 1,
##                                      designname=NULL) {
##   
##   nfiles <- length(filenames)
##   n <- nfiles*nwells
##   
##   ## error if no file name !
## 
##   if (n == 0)
##     stop("No file name given !")
##   if(mode(filenames)!="character")
##     stop("filenames must be of type character!")
##   
##   pdata <- pData(phenoData)
##   if(dim(pdata)[1] != n) {
## 
##     ##if empty pdata filename are samplenames
## 
##     cat("Incompatible phenoData object. Created a new one.\n")
##     samplenames <- basename(unlist(filenames))
## 
##     if(nwells>1){
##       wells <- paste(".",as.character(unique(get(designname)@lookup$container)),sep="")
##       samplenames <- as.character(t(outer(samplenames,wells,paste,sep="")))
##     }
##     
##     pdata <- data.frame(sample=1:n, row.names=samplenames)
##     phenoData <- new("AnnotatedDataFrame",
##                      data=pdata,
##                      varMetadata=data.frame(labelDescription="arbitrary numbering", row.names="sample"))
##   }
##   else samplenames <- rownames(pdata)
##   
##   if (is.null(description))
##     {
##       description <- new("MIAME")
##       description@preprocessing$filenames <- filenames
##       description@preprocessing$oligoversion <- packageDescription("oligo")$Version
##     }
##   
##   return(list(filenames=filenames,samplenames=samplenames, phenoData=phenoData, description=description))
## }

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

  results <- smartReadXYS(filenames, sampleNames, path=oligoBigObjectPath())
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
  uid <- genDatasetUID(filenames)
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
  results <- smartReadXYS(channel1, sampleNames, prefix="intensities1-",
                          path=oligoBigObjectPath(), uid=uid)
  channel1Intensities <- results[["exprMatrix"]]
  intensityFile1 <- results[["intensityFile"]]
  rm(results)
  results <- smartReadXYS(channel2, sampleNames, prefix="intensities2-",
                          path=oligoBigObjectPath(), uid=uid)
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








  
##   tmp <- .Call("R_read_xys_files", channel1, verbose)
##   channel1Intensities <- tmp[["intensities"]]
##   date1 <- tmp[["date"]]
##   rm(tmp)
##   tmp <- .Call("R_read_xys_files", channel2, verbose)
##   channel2Intensities <- tmp[["intensities"]]
##   date2 <- tmp[["date"]]
##   rm(tmp)
  ## must get dates from here

##   metadata <- getMetadata2(channel1Intensities, channel2Intensities,
##                            channel1, channel2, phenoData, featureData,
##                            experimentData, notes, sampleNames,
##                            NgsDate2Posix(date1), NgsDate2Posix(date2))
  ## colnames(channel1Intensities) <- colnames(channel2Intensities) <- Biobase::sampleNames(metadata[["phenoData"]])

##   out <- new("TilingFeatureSet",
##              channel1=channel1Intensities,
##              channel2=channel2Intensities,
##              manufacturer="NimbleGen",
##              annotation=pkgname,
##              phenoData=metadata[["phenoData"]],
##              experimentData=metadata[["experimentData"]],
##              featureData=metadata[["featureData"]])
##   return(out)
}
