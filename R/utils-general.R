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

## checkValidPhenodataForFiles <- function(filenames, pd){
##   if (missing(pd)) return(FALSE)
##   if (class(pd) != "AnnotatedDataFrame") return(FALSE)
##   return(length(filenames) == nrow(pd))
## }

## createDefaultMiame <- function(filenames, notes){
##   experimentData <- new("MIAME")
##   preproc(experimentData)$filenames <- filenames
##   preproc(experimentData)$oligoversion <- packageDescription("oligo", field="Version")
##   if (!missing(notes))
##     notes(experimentData) <- notes
##   experimentData
## }

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

## Selectors for probes - useful for exon/gene arrays

getFidProbeset <- function(object){
  conn <- db(object)
  sql <- "SELECT fid, fsetid FROM pmfeature"
  featureInfo <- dbGetQuery(conn, sql)
  featureInfo <- featureInfo[order(featureInfo[["fsetid"]]),]
  rownames(featureInfo) <- NULL
  return(featureInfo)
}

getFidMetaProbesetCore <- function(object){
  conn <- db(object)
  sql <- "SELECT fid, meta_fsetid as fsetid FROM pmfeature INNER JOIN core_mps USING(fsetid)"
  featureInfo <- dbGetQuery(conn, sql)
  featureInfo <- featureInfo[order(featureInfo[["fsetid"]]),]
  rownames(featureInfo) <- NULL
  return(featureInfo)
}

getFidMetaProbesetFull <- function(object){
  conn <- db(object)
  sql <- "SELECT fid, meta_fsetid as fsetid FROM pmfeature INNER JOIN full_mps USING(fsetid)"
  featureInfo <- dbGetQuery(conn, sql)
  featureInfo <- featureInfo[order(featureInfo[["fsetid"]]),]
  rownames(featureInfo) <- NULL
  return(featureInfo)
}

getFidMetaProbesetExtended <- function(object){
  conn <- db(object)
  sql <- "SELECT fid, meta_fsetid as fsetid FROM pmfeature INNER JOIN extended_mps USING(fsetid)"
  featureInfo <- dbGetQuery(conn, sql)
  featureInfo <- featureInfo[order(featureInfo[["fsetid"]]),]
  rownames(featureInfo) <- NULL
  return(featureInfo)
}


## Date/time extractors
GetAffyTimeDateAsString <- function(filenames){
    f <- function(x) read.celfile.header(x, "full")[["DatHeader"]]
    DatHeader = sapply(filenames, f)
    gsub("(.* )([0-9]{2}/[0-9]{2}/[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2})(.*)",
         "\\2", DatHeader)
}

NgsDate2Posix <- function(txt){
  if (!missing(txt)){
    original <- txt
    tmp <- gsub("(.* )(.{1,3}) ([0-9]{4})", "\\1\\3", txt)
    out <- as.POSIXct(tmp, format="%A %B %d %H:%M:%S %Y")
    if (any(is.na(out))){
      warning("Returning dates/times as strings due to incompatibility.")
      out <- original
    }
  }else{
    out <- NULL
  }
  out
}

AffyDate2Posix <- function(txt){
  if (!missing(txt)){
    original <- txt
    out <- as.POSIXct(txt, format="%m/%d/%y %H:%M:%S")
    if (any(is.na(out))){
      warning("Returning dates/times as strings - format not recognized.")
      out <- original
    }
  }else{
    out <- NULL
  }
  out
}

basicPData <- function(mat, filenames, dates){
  if (missing(dates))
      dates <- rep(NA, length(filenames))
  stopifnot(length(filenames) == length(dates))
  pdd <- data.frame(exprs=filenames, dates=dates)
  vmd <- data.frame(labelDescription=c(
                    "Names of files used in 'exprs'",
                    "Run dates for files used in 'exprs'"),
                    channel=factor("_ALL_", levels=c("exprs", "_ALL_")))
  phenoData <- new("AnnotatedDataFrame", data=pdd, varMetadata=vmd)
  sampleNames(phenoData) <- colnames(mat)
  return(phenoData)
}

basicPData2 <- function(mat1, mat2, filenames1, filenames2,
                            dates1, dates2){
  if (missing(dates1)) dates1 <- rep(NA, length(filenames1))
  if (missing(dates2)) dates2 <- rep(NA, length(filenames2))
  stopifnot(identical(colnames(mat1), colnames(mat2)),
            length(filenames1) == length(filenames2),
            length(dates1) == length(dates2))
  pdd <- data.frame(filenamesChannel1=filenames1,
                    filenamesChannel2=filenames2,
                    dates1=dates1, dates2=dates2)
  vmd <- data.frame(labelDescription=c(
                      "Names of files used in 'channel1'",
                      "Names of files used in 'channel2'",
                      "Run dates for files used in 'channel1'",
                      "Run dates for files used in 'channel2'"),
                    channel=factor(rep(c("channel1", "channel2"), 2),
                      levels=c("channel1", "channel2", "_ALL_")))
  phenoData <- new("AnnotatedDataFrame", data=pdd, varMetadata=vmd)
  sampleNames(phenoData) <- colnames(mat1)
  return(phenoData)
}

basicPhData1 <- function(mat){
  pdd <- data.frame(index=1:ncol(mat))
  vmd <- data.frame(labelDescription=c("Index"),
                    channel=factor("_ALL_", levels=c("exprs", "_ALL_")))
  phenoData <- new("AnnotatedDataFrame", data=pdd, varMetadata=vmd)
  sampleNames(phenoData) <- colnames(mat)
  return(phenoData)
}

basicPhData2 <- function(mat1, mat2){
  stopifnot(identical(colnames(mat1), colnames(mat2)))
  pdd <- data.frame(index=1:ncol(mat1))
  vmd <- data.frame(labelDescription=c("Index"),
                    channel=factor("_ALL_",
                      levels=c("channel1", "channel2", "_ALL_")))
  phenoData <- new("AnnotatedDataFrame", data=pdd, varMetadata=vmd)
  sampleNames(phenoData) <- colnames(mat1)
  return(phenoData)
}

basicAnnotatedDataFrame <- function(mat, byrow=FALSE){
    Biobase:::annotatedDataFrameFromMatrix(mat, byrow=byrow)
}

## colors I like in oligo
darkColors <- function(n){
  cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
            "#E6AB02", "#A6761D", "#666666")
  fff <- colorRampPalette(cols)
  fff(n)
}

seqColors <- function(n){
  cols <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
            "#4292C6", "#2171B5", "#08519C", "#08306B")
  fff <- colorRampPalette(cols)
  fff(n)
}

## gcrma-like stuff
## from jean

getAffinitySplineCoefficients <- function(intensities, sequences){
  if (is(sequences, "DNAStringSet"))
    sequences <- as.character(sequences)
  stopifnot(is.character(sequences))
  stopifnot(is.matrix(intensities))
  design <- sequenceDesignMatrix(sequences)
  rm(sequences)
  B <- ns(1:25, df=5)
  design <- cbind(design[,  1:25] %*% B,
                  design[, 26:50] %*% B,
                  design[, 51:75] %*% B)
  fits <- lm(intensities~design)
  coefs <- coef(fits)[-1,]
  rownames(coefs) <- rep(c("A", "C", "G"), each=5)
  coefs
}

getBaseProfile <- function(coefs, probeLength=25, plot=FALSE, ...){
  stopifnot(is.vector(coefs))
  P <- as.integer(length(coefs)/3)
  B <- ns(1:probeLength, df=5)
  effects <- matrix(0, nrow=probeLength, ncol=4)
  colnames(effects) <- c("A", "C", "G", "T")
  for (i in 1:3)
    effects[, i] <- B %*% coefs[(i-1)*P + (1:P)]
  effects <- sweep(effects, 1, rowMeans(effects))
  if(plot) matplot(1:probeLength, effects, ...)
  invisible(effects)
}

## SnpSuperSet method?
plotM <- function(x, snp, ...){
  crlmmInfo <- file.path(system.file("extdata", package=annotation(x)),
                         paste(annotation(x), "CrlmmInfo.rda", sep=""))
  infoObj <- load(crlmmInfo)
  f0 <- get(infoObj)$params$f0
  rm(list=infoObj)
  obj <- x[snp,]
  theM <- getM(obj)[,,]
  theF <- cbind(antisense=assayDataElement(obj, "antisenseF")[,],
                sense=assayDataElement(obj, "senseF")[,])
  theCalls <- calls(obj)[,]
  correctedM <- theM + (theCalls - 2) * (f0-theF)
  plot(correctedM, ...)
  invisible(correctedM)
}
