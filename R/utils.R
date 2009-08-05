cleanPlatformName <- function(x)
  gsub("[_-]", ".", paste("pd.", tolower(x), sep=""))

HuberAllRowsByGroup <- function(X, Y, k=1.5){
  .Call("R_HuberMatrixRows2",X, Y, k, PACKAGE="oligo")
}

trigammaInverse <- function(x) {
#	Solve trigamma(y) = x for y
#	Gordon Smyth
#	8 Sept 2002.  Last revised 12 March 2004.

#	Non-numeric or zero length input
	if(!is.numeric(x)) stop("Non-numeric argument to mathematical function")
	if(length(x)==0) return(numeric(0))

#	Treat out-of-range values as special cases
	omit <- is.na(x)
	if(any(omit)) {
		y <- x
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x < 0)
	if(any(omit)) {
		y <- x
		y[omit] <- NaN
		warning("NaNs produced")
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x > 1e7)
	if(any(omit)) {
		y <- x
		y[omit] <- 1/sqrt(x[omit])
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x < 1e-6)
	if(any(omit)) {
		y <- x
		y[omit] <- 1/x[omit]
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}

#	Newton's method
#	1/trigamma(y) is convex, nearly linear and strictly > y-0.5,
#	so iteration to solve 1/x = 1/trigamma is monotonically convergent
	y <- 0.5+1/x
	iter <- 0
	repeat {
		iter <- iter+1
		tri <- trigamma(y)
		dif <- tri*(1-tri/x)/psigamma(y,deriv=2)
		y <- y+dif
		if(max(-dif/y) < 1e-8) break
		if(iter > 50) {
			warning("Iteration limit exceeded")
			break
		}
	}
	y
}


checkValidFilenames <- function(filenames) {
  ## must be character
  stopifnot(is.character(filenames))

  ## must exist
  dirs <- file.info(filenames)[["isdir"]]
  if (any(is.na(dirs))){
    msg <- paste("These do not exist:",
                 paste("\t", filenames[is.na(dirs)], collapse="\n"), sep="\n")
    stop(msg, call.=FALSE)
  }
  
  ## must be files, not dir
  if (any(dirs)){
    msg <- paste("These are directories:",
                 paste("\t", filenames[dirs], collapse="\n"), sep="\n")
    stop(msg, call.=FALSE)
  }

  ## must be readable
  readable <- file.access(filenames, 4) == 0
  if (any(!readable)){
    msg <- paste("These are not readable:",
                 paste("\t", filenames[!readable], collapse="\n"), sep="\n")
    stop(msg, call.=FALSE)
  }
  TRUE
}

checkValidPhenodataForFiles <- function(filenames, pd){
  if (missing(pd)) return(FALSE)
  if (class(pd) != "AnnotatedDataFrame") return(FALSE)
  return(length(filenames) == nrow(pd))
}

createDefaultMiame <- function(filenames, notes){
  experimentData <- new("MIAME")
  preproc(experimentData)$filenames <- filenames
  preproc(experimentData)$oligoversion <- packageDescription("oligo", field="Version")
  if (!missing(notes))
    notes(experimentData) <- notes
  experimentData
}

getCelChipType <- function(x, useAffyio){
  ifelse(useAffyio,
         read.celfile.header(x)[["cdfName"]],
         readCelHeader(x)[["chiptype"]])
}

checkChipTypes <- function(filenames, verbose=TRUE, manufacturer, useAffyio){
  if (missing(manufacturer)) stop("'checkChipTypes' needs 'manufacturer'")
  if (manufacturer == "affymetrix"){
    chips <- sapply(filenames, getCelChipType, useAffyio)
    ok <- length(unique(chips)) == 1
    if(!ok & verbose) message("All the CEL files must be of the same type.")
  }else if(manufacturer == "nimblegen"){
    designnamelist <- NULL
    for (xysfile in filenames){
      firstline <- readxysHeader(xysfile)
      designname <- unlist(strsplit(firstline[grep("designname",firstline,fixed=TRUE,useBytes=TRUE)],"="))[2]
      designnamelist <- rbind(designnamelist,designname)
    }
    ok <- length(unique(designnamelist)) == 1
    if(!ok & verbose) message("All the XYS files must be of the same type.")
  }else{
    stop("'manufacturer' ", manufacturer, " unknown")
  }
  ok
}

getMetadata <- function(theMatrix, filenames, phenoData, featureData,
                        experimentData, notes, sampleNames, datetime){
  stopifnot(!missing(theMatrix), !missing(filenames),
            is.matrix(theMatrix), is.character(filenames),
            ncol(theMatrix) == length(filenames))
  if (!checkValidPhenodataForFiles(filenames, phenoData)){
    phenoData <- new("AnnotatedDataFrame",
                     data=data.frame(filenames = filenames),
                     varMetadata=data.frame(labelDescription=c("names of files used to create object")))
  }
  if (missing(sampleNames))
    sampleNames <- basename(filenames)
  sampleNames(phenoData) <- sampleNames
  
  if (!missing(datetime)){
    pdDateTime <- new("AnnotatedDataFrame",
                      data=data.frame(DateTime = datetime),
                      varMetadata=data.frame(labelDescription="date/time from raw files"))
    sampleNames(pdDateTime) <- sampleNames(phenoData)
    phenoData <- combine(phenoData, pdDateTime)
    rm(pdDateTime)
  }
  
  if (missing(experimentData))
    experimentData <- createDefaultMiame(filenames)
  if (missing(featureData))
    featureData <- annotatedDataFrameFrom(theMatrix,byrow=TRUE)

  out <- list(filenames=filenames,
              phenoData=phenoData,
              featureData=featureData,
              experimentData=experimentData)
}


getMetadata2 <- function(theMatrix1, theMatrix2,
                         filesChannel1, filesChannel2,
                         phenoData, featureData,
                         experimentData, notes, sampleNames,
                         datetime1, datetime2){
  stopifnot(!missing(theMatrix1), !missing(theMatrix2),
            !missing(filesChannel1), !missing(filesChannel2),
            is.matrix(theMatrix1), is.matrix(theMatrix2),
            is.character(filesChannel1), is.character(filesChannel2),
            ncol(theMatrix1) == length(filesChannel1),
            ncol(theMatrix2) == length(filesChannel2),
            length(filesChannel1) == length(filesChannel2))

  if (missing(phenoData)){
    theDF <- data.frame(filenamesChannel1=filesChannel1,
                        filenamesChannel2=filesChannel2)
    vmDF <- data.frame(labelDescription=c(
                         "names of files used to create 'channel1'",
                         "names of files used to create 'channel2'"),
                       channel=factor(c("channel1", "channel2"),
                         levels=c("channel1", "channel2", "_ALL_")))
    phenoData <- new("AnnotatedDataFrame",
                     data=theDF,
                     varMetadata=vmDF)
    rm(theDF, vmDF)
  }

  if (!missing(datetime1)){
    theDF <- data.frame(channel1DateTime=datetime1)
    vmDF <- data.frame(labelDescription=c("date/time from raw files"),
                       channel=factor("channel1", levels=c("channel1", "channel2", "_ALL_")))
    tmpPD <- new("AnnotatedDataFrame",
                 data=theDF,
                 varMetadata=vmDF)
    sampleNames(tmpPD) <- sampleNames(phenoData)
    phenoData <- combine(phenoData, tmpPD)
    rm(theDF, vmDF, tmpPD)
  }
  if (!missing(datetime2)){
    theDF <- data.frame(channel2DateTime=datetime2)
    vmDF <- data.frame(labelDescription=c("date/time from raw files"),
                       channel=factor("channel2", levels=c("channel1", "channel2", "_ALL_")))
    tmpPD <- new("AnnotatedDataFrame",
                 data=theDF,
                 varMetadata=vmDF)
    sampleNames(tmpPD) <- sampleNames(phenoData)
    phenoData <- combine(phenoData, tmpPD)
    rm(theDF, vmDF, tmpPD)
  }

  if (missing(sampleNames))
    sampleNames <- basename(filesChannel1)
  sampleNames(phenoData) <- sampleNames
  
  if (missing(experimentData))
    experimentData <- createDefaultMiame(filesChannel1)
  if (missing(featureData))
    featureData <- annotatedDataFrameFrom(theMatrix1,byrow=TRUE)

  out <- list(filesChannel1=filesChannel1,
              filesChannel2=filesChannel2,
              phenoData=phenoData,
              featureData=featureData,
              experimentData=experimentData)
}


basicRMA <- function(pmMat, pnVec, normalize=TRUE, background=TRUE,
                     bgversion=2, destructive=FALSE, verbose=TRUE, ...){
  pns <- unique(pnVec)
  nPn <- length(unique(pnVec))
  pnVec <- split(0:(length(pnVec)-1), pnVec)
  
  if (destructive){
    theExprs <- .Call("rma_c_complete", pmMat, pnVec, nPn, normalize,
                      background, bgversion, verbose, PACKAGE="oligo")
  }else{
    theExprs <- .Call("rma_c_complete_copy", pmMat, pnVec, nPn,
                      normalize, background, bgversion,
                      verbose, PACKAGE="oligo")
  }
  colnames(theExprs) <- colnames(pmMat)
  rownames(theExprs) <- pns
  return(theExprs)
}

sequenceDesignMatrix <- function(seqs){
  if(length(unique(sapply(seqs,nchar)))!=1) stop("Sequences must be of same length.")
  oligolength <- nchar(seqs[1])
  mat <- .Call("gcrma_getSeq2",paste(seqs,collapse=""),length(seqs),oligolength,PACKAGE="oligo")
  colnames(mat) <- paste(rep(c("A","C","G"),rep(oligolength,3)),position=rep(1:oligolength,3),sep="_")
  return(mat)
}

basecontent <- function(seq) {
  good   = !is.na(seq)
  havena = !all(good)
  if(havena)
    seq = seq[good]
  
  rv = .Call("basecontent", seq, PACKAGE="oligo")

  if(havena) {
    z = rv
    rv = matrix(NA, nrow=length(good), ncol=ncol(z))
    colnames(rv) = colnames(z)
    rv[good, ] = z
  }
  
  return(rv)
}

## Helper to parse the dots and combine with filenames
getFilenames <- function(filenames, ...){
  if (!missing(filenames)){
    filenames <- c(filenames, unlist(list(...)))
  }else{
    filenames <- unlist(list(...))
  }
  filenames
}

getNgsColorsInfo <- function(path=".", pattern1="_532", pattern2="_635", ...){
  files1 <- list.xysfiles(path, pattern=pattern1, ...)
  files2 <- list.xysfiles(path, pattern=pattern2, ...)
  prefix1 <- sapply(strsplit(basename(files1), pattern1), "[[", 1)
  prefix2 <- sapply(strsplit(basename(files2), pattern2), "[[", 1)
  sugg <- identical(prefix1, prefix2)
  stopifnot(length(files1) == length(files2))
  out <- data.frame(color1=files1, color2=files2, stringsAsFactors=FALSE)
  if (sugg)
    out[["sampleNames"]] <- prefix1
  return(out)
}

## Date/time extractors
GetAffyTimeDateAsString <- function(filenames, useAffyio=TRUE){
  if (useAffyio){
    f <- function(x) read.celfile.header(x, "full")[["DatHeader"]]
  }else{
    f <- function(x) readCelHeader(x)[["datheader"]]
  }
  DatHeader = sapply(filenames, f)
  gsub("(.* )([0-9]{2}/[0-9]{2}/[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2})(.*)",
       "\\2", DatHeader)
}

NgsDate2Posix <- function(txt){
  original <- txt
  tmp <- gsub("(.* )(.{1,3}) ([0-9]{4})", "\\1\\3", txt)
  out <- as.POSIXct(tmp, format="%A %B %d %H:%M:%S %Y")
  if (any(is.na(out))){
    warning("Returning dates/times as strings due to incompatibility.")
    out <- original
  }
  out
}

AffyDate2Posix <- function(txt){
  original <- txt
  out <- as.POSIXct(txt, format="%m/%d/%y %H:%M:%S")
  if (any(is.na(out))){
    warning("Returning dates/times as strings - format not recognized.")
    out <- original
  }
  out
}
