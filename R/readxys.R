# Reading NimbleGen Data
# Author: Benilton Carvalho
# Last Modification: May 28, 2005

list.xysfiles <-   function(...){
  files <- list.files(...)
  return(files[grep("\\.[xX][yY][sS]\\.gz$|\\.[xX][yY][sS]$", files)])
}

readxysHeader <- function(filename) scan(filename,nlines=1,quiet=TRUE, what=character(0))

readonexysfile <- function(filename)
  read.delim(filename, comment.char="#")

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

read.xysfiles <- function(filenames,
                          pkgname=NULL,
                          phenoData=NULL,
                          featureData=NULL,
                          experimentData=NULL,
                          notes=NULL,
                          verbose = FALSE) {
  
  ## Create space to store the design names
  designnamelist <- NULL

  ## Read the header for all XYS files, get design name for each
  for (xysfile in filenames){
    firstline <- readxysHeader(xysfile)
    designname <- unlist(strsplit(firstline[grep("designname",firstline,fixed=TRUE,useBytes=TRUE)],"="))[2]
    designnamelist <- rbind(designnamelist,designname)
  }

  ## How many different designs?
  numberdesigns <- length(unique(designnamelist))

  ## All XYS files should point to one NDF file
  if(numberdesigns != 1)
    stop("XYS Files do not refer to the same design!")
    
  ## Load PDenv for the XYS files
  if (is.null(pkgname))
    pkgname <- cleanPlatformName(designname)
  if (require(pkgname, character.only=TRUE)){
    if (verbose)
      message("Platform design info loaded.")
  }else{
    stop("Must install the ", pkgname, " package.")
  }

  tmpExprs <- NULL
  for (i in seq(along=filenames)){
    if (verbose) cat(".")

    ## Read XYS "as is"
    tmpE <- readonexysfile(filenames[i])
    if (length(filenames) > 1){
      tmpExprs <- cbind(tmpExprs, tmpE$SIGNAL)
    }else{
      tmpExprs <- matrix(tmpE$SIGNAL, ncol=1)
    }
    if(verbose) cat(" Done.\n")
  }
  tmpE$index <- tmpE$X + (tmpE$Y-1)*geometry(get(pkgname))[2]
  idx <- order(tmpE$index)
  tmpExprs <- tmpExprs[idx,, drop=FALSE]
  rm(tmpE, idx)

  arrayType <- kind(get(pkgname))
  theClass <- switch(arrayType,
                     tiling="TilingFeatureSet",
                     expression="ExpressionFeatureSet",
                     SNP="SnpFeatureSet",
                     exon="ExonFeatureSet",
                     stop("unknown array type: ", arrayType))

  out <- new(theClass,
             exprs=tmpExprs,
             platform=pkgname,
             manufacturer="NimbleGen",
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
    ed <- experimentData
    experimentData(out) <- ed
  }
  out
}
