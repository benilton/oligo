cleanPlatformName <- function(x)
  gsub("[_-]", ".", paste("pd.", tolower(x), sep=""))

i2xy <- function(i,obatch){
  xy <- mget(c("X","Y"),envir=featureInfo(getPD(obatch)))
  return(cbind(xy$X[i],xy$Y[i]))
}

xy2i <- function(x,y,obatch){
  xy <- mget(c("X","Y"),envir=featureInfo(getPD(obatch)))
  xy1 <- xy$X+(xy$Y-1)*max(xy$X)
  xy2 <- x+(y-1)*max(xy$X)
  match(xy2,xy1)
}

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

checkChipTypes <- function(filenames, verbose=TRUE){
  chips <- sapply(filenames, function(x) readCelHeader(x)[["chiptype"]])
  ok <- length(unique(chips)) == 1
  if(!ok & verbose) message("All the CEL files must be of the same type.")
  ok
}

getMetadata <- function(theMatrix, filenames, phenoData, featureData,
                        experimentData, notes, sampleNames){
  stopifnot(!missing(theMatrix), !missing(filenames),
            is.matrix(theMatrix), is.character(filenames))
  if (!checkValidPhenodataForFiles(filenames, phenoData))
    phenoData <- annotatedDataFrameFrom(theMatrix, byrow=FALSE)
  if (!missing(sampleNames)){
    sampleNames(phenoData) <- sampleNames
  }else{
    sampleNames(phenoData) <- sub("^/?([^/]*/)*", "", filenames, extended=TRUE)
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
