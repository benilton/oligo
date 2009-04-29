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


stuffForXYSandCELreaders <- function(filenames,
                                     phenoData=new("AnnotatedDataFrame"),
                                     description=NULL,
                                     notes="",
                                     verbose = FALSE,
                                     nwells = 1,
                                     designname=NULL) {
  
  nfiles <- length(filenames)
  n <- nfiles*nwells
  
  ## error if no file name !

  if (n == 0)
    stop("No file name given !")
  if(mode(filenames)!="character")
    stop("filenames must be of type character!")
  
  pdata <- pData(phenoData)
  if(dim(pdata)[1] != n) {

    ##if empty pdata filename are samplenames

    cat("Incompatible phenoData object. Created a new one.\n")
    samplenames <- sub("^/?([^/]*/)*", "", unlist(filenames), extended=TRUE)

    if(nwells>1){
      wells <- paste(".",as.character(unique(get(designname)@lookup$container)),sep="")
      samplenames <- as.character(t(outer(samplenames,wells,paste,sep="")))
    }
    
    pdata <- data.frame(sample=1:n, row.names=samplenames)
    phenoData <- new("AnnotatedDataFrame",
                     data=pdata,
                     varMetadata=data.frame(labelDescription="arbitrary numbering", row.names="sample"))
  }
  else samplenames <- rownames(pdata)
  
  if (is.null(description))
    {
      description <- new("MIAME")
      description@preprocessing$filenames <- filenames
      description@preprocessing$oligoversion <- packageDescription("oligo")$Version
    }
  
  return(list(filenames=filenames,samplenames=samplenames, phenoData=phenoData, description=description))
}

read.xysfiles <- function(..., filenames, pkgname, phenoData,
                          featureData, experimentData, notes,
                          verbose=TRUE, sampleNames, checkType=TRUE) {

  filenames <- getFilenames(filenames=filenames, ...)
  checkValidFilenames(filenames)
  if (checkType) stopifnot(checkChipTypes(filenames, verbose, "nimblegen"))

  ## Get design name from the first
  firstline <- readxysHeader(filenames[1])
  designname <- unlist(strsplit(firstline[grep("designname",
                firstline, fixed=TRUE, useBytes=TRUE)], "="))[2]
    
  ## Load PDenv for the XYS files
  if (missing(pkgname))
    pkgname <- cleanPlatformName(designname)
  if (require(pkgname, character.only=TRUE)){
    if (verbose)
      message("Platform design info loaded.")
  }else{
    stop("Must install the ", pkgname, " package.")
  }

  arrayType <- kind(get(pkgname))
  tmpExprs <- NULL
  if (verbose) cat("Reading XYS files...\n")
  for (i in seq(along=filenames)){
    ## Read XYS "as is"
    if (verbose) cat(" ...  ", filenames[i], ".\n")
    tmpE <- readonexysfile(filenames[i])
    if (length(filenames) > 1){
      tmpExprs <- cbind(tmpExprs, tmpE$SIGNAL)
    }else{
      tmpExprs <- matrix(tmpE$SIGNAL, ncol=1)
    }
  }
  tmpE$index <- tmpE$X + (tmpE$Y-1)*geometry(get(pkgname))[2]
  idx <- order(tmpE$index)
  tmpExprs <- tmpExprs[idx,, drop=FALSE]
  rm(tmpE, idx)
  dimnames(tmpExprs) <- NULL
  
  metadata <- getMetadata(tmpExprs, filenames, phenoData, featureData,
                          experimentData, notes, sampleNames)

  theClass <- switch(arrayType,
                     tiling="TilingFeatureSet",
                     expression="ExpressionFeatureSet",
                     SNP="SnpFeatureSet",
                     SNPCNV="SnpCnvFeatureSet",
                     exon="ExonFeatureSet",
                     gene="GeneFeatureSet",
                     stop("Unknown array type: ", arrayType))

  return(new(theClass, exprs=tmpExprs, manufacturer="NimbleGen",
             annotation=pkgname, phenoData=metadata[["phenoData"]],
             experimentData=metadata[["experimentData"]],
             featureData=metadata[["featureData"]]))
}


## For 2 channels - Tiling
read.xysfiles2 <- function(channel1, channel2, pkgname, phenoData,
                          featureData, experimentData, notes,
                          verbose=TRUE, sampleNames, checkType=TRUE) {
  
  filenames <- c(channel1, channel2)
  checkValidFilenames(filenames)
  if (checkType) stopifnot(checkChipTypes(filenames, verbose, "nimblegen"))

  ## Get design name from the first
  firstline <- readxysHeader(filenames[1])
  designname <- unlist(strsplit(firstline[grep("designname",
                firstline, fixed=TRUE, useBytes=TRUE)], "="))[2]
    
  ## Load PDenv for the XYS files
  if (missing(pkgname))
    pkgname <- cleanPlatformName(designname)
  if (require(pkgname, character.only=TRUE)){
    if (verbose)
      message("Platform design info loaded.")
  }else{
    stop("Must install the ", pkgname, " package.")
  }

  arrayType <- kind(get(pkgname))
  channel1Intensities <- readXysMatrix(channel1)
  channel2Intensities <- readXysMatrix(channel2)
  idxChannel1 <- channel1Intensities[["X"]]+(channel1Intensities[["Y"]]-1)*geometry(get(pkgname))[2]
  idxChannel2 <- channel2Intensities[["X"]]+(channel2Intensities[["Y"]]-1)*geometry(get(pkgname))[2]
  stopifnot(identical(sort(idxChannel1), sort(idxChannel2)))
  channel1Intensities <- channel1Intensities[["intensities"]][order(idxChannel1),, drop=FALSE]
  channel2Intensities <- channel2Intensities[["intensities"]][order(idxChannel2),, drop=FALSE]
  rm(idxChannel1, idxChannel2)

  metadata <- getMetadata(channel1Intensities, channel1, phenoData, featureData, experimentData, notes, sampleNames)

  out <- new("TilingFeatureSet2",
             channel1=channel1Intensities,
             channel2=channel2Intensities,
             manufacturer="NimbleGen",
             annotation=pkgname,
             phenoData=metadata[["phenoData"]],
             experimentData=metadata[["experimentData"]],
             featureData=metadata[["featureData"]])
  return(out)
}
