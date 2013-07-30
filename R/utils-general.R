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

getFidProbeset <- function(object, sortBy='fsetid'){
    if (!is.null(sortBy))
        sortBy <- match.arg(sortBy, c('fid', 'fsetid'))
  conn <- db(object)
  sql <- "SELECT fid, fsetid FROM pmfeature"
  featureInfo <- dbGetQuery(conn, sql)
  if (!is.null(sortBy)){
      featureInfo <- featureInfo[order(featureInfo[[sortBy]]),]
      rownames(featureInfo) <- NULL
  }
  return(featureInfo)
}

getFidMetaProbesetCore <- function(object, sortBy='fsetid'){
    if (!is.null(sortBy))
        sortBy <- match.arg(sortBy, c('fid', 'fsetid'))
  conn <- db(object)
  sql <- "SELECT fid, meta_fsetid as fsetid FROM pmfeature INNER JOIN core_mps USING(fsetid)"
  featureInfo <- dbGetQuery(conn, sql)
  if (!is.null(sortBy)){
      featureInfo <- featureInfo[order(featureInfo[[sortBy]]),]
      rownames(featureInfo) <- NULL
  }
  return(featureInfo)
}

getFidMetaProbesetFull <- function(object, sortBy='fsetid'){
    if (!is.null(sortBy))
        sortBy <- match.arg(sortBy, c('fid', 'fsetid'))
  conn <- db(object)
  sql <- "SELECT fid, meta_fsetid as fsetid FROM pmfeature INNER JOIN full_mps USING(fsetid)"
  featureInfo <- dbGetQuery(conn, sql)
  if (!is.null(sortBy)){
      featureInfo <- featureInfo[order(featureInfo[[sortBy]]),]
      rownames(featureInfo) <- NULL
  }
  return(featureInfo)
}

getFidMetaProbesetExtended <- function(object, sortBy='fsetid'){
    if (!is.null(sortBy))
        sortBy <- match.arg(sortBy, c('fid', 'fsetid'))
  conn <- db(object)
  sql <- "SELECT fid, meta_fsetid as fsetid FROM pmfeature INNER JOIN extended_mps USING(fsetid)"
  featureInfo <- dbGetQuery(conn, sql)
  if (!is.null(sortBy)){
      featureInfo <- featureInfo[order(featureInfo[[sortBy]]),]
      rownames(featureInfo) <- NULL
  }
  return(featureInfo)
}

stArrayPmInfo <- function(object, target='core', sortBy='fsetid'){
    ## *PmInfo returns a data.frame with 'fid' and 'fsetid'
    target <- match.arg(target, c('probeset', 'core', 'full', 'extended'))
    theFun <- switch(target,
                     probeset=getFidProbeset,
                     core=getFidMetaProbesetCore,
                     full=getFidMetaProbesetFull,
                     extended=getFidMetaProbesetExtended)
    theFun(object, sortBy=sortBy)
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

seqColors2 <- function(n){
    cols <- c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59",
              "#EF6548", "#D7301F", "#B30000", "#7F0000")
    fff <- colorRampPalette(cols)
    fff(n)
}

divColors <- function(n){
    cols <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
              "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
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


## MA plot functions
maplot <- function (x, transfo=log2, groups, refSamples,
                    which, pch=".", summaryFun=rowMedians,
                    plotFun=smoothScatter,
                    main="vs pseudo-median reference chip",
                    pairs=FALSE, ...){
    transfo <- match.fun(transfo)
    stopifnot(is.function(transfo))
    summaryFun <- match.fun(summaryFun)
    stopifnot(is.function(summaryFun))

    ## Checking for errors
    if (!missing(refSamples)){
        if (pairs)
            stop("Cannot combine pairs=TRUE with 'refSamples'")
        if (!is.numeric(refSamples))
            stop("'refSamples' must be integer(s)")
        refSamples <- as.integer(refSamples)
        rg <- range(refSamples)
        if (min(refSamples) < 1 | max(refSamples) > ncol(x))
            stop("Ensure 'refSamples' are integers between 1 and ", ncol(x), ".")
    }

    if (!missing(which)){
        if (!is.numeric(which))
            stop("'which' must be integer(s)")
        which <- as.integer(which)
        rg <- range(which)
        if (min(which) < 1 | max(which) > ncol(x))
            stop("Ensure 'which' contains integers between 1 and ", ncol(x), ".")
    }

    if (!missing(groups)){
        stopifnot(is.factor(groups))
        stopifnot(length(groups) == ncol(x))
        ## if groups is given, 'which' refers to group
        if (!missing(which))
            if (max(which) > length(levels(groups)))
                stop("'which' must be smaller than ", length(levels(groups)))
        ## if groups is given, 'refSamples' refers to group(s) of
        ## reference
        if (!missing(refSamples))
            if (max(refSamples) > length(levels(groups)))
                stop("'refSamples' must be smaller than ", length(levels(groups)))
    }
    ## END OF ERROR CHECK

##    x <- transfo(exprs(object))
    x <- transfo(x)

    if (missing(groups)){
        if (missing(which))
            which <- 1:ncol(x)
        if (!pairs) {
            if (missing(refSamples)) {
                medianchip <- summaryFun(x)
            } else if (length(refSamples) > 1) {
                medianchip <- summaryFun(x[, refSamples])
            } else {
                medianchip <- x[, refSamples]
            }

            M <- sweep(x, 1, medianchip, FUN = "-")
            A <- 1/2 * sweep(x, 1, medianchip, FUN = "+")
            for (i in which) {
                if (missing(refSamples)){
                    themain <- paste(colnames(x)[i], main)
                }else{
                    if (length(refSamples) == 1) {
                        if (i != refSamples)
                            themain <- paste(colnames(x)[i], "vs", colnames(x)[refSamples])
                    } else {
                        themain <- paste(colnames(x)[i], main)
                    }
                }
                basicMvAplot(A[, i], M[, i], main=themain, xlab="A",
                             ylab="M", plotFun=plotFun, pch=pch, ...)
            }
        } else {
            basicMvApairsPlot(x[, which], transfo=transfo, plotFun=plotFun, ...)
        } ## end if (!pairs)
    } else { ## if groups
        groups.list <- split(1:ncol(x), groups)
        grouped.data <- matrix(0, nrow(x), length(groups.list))
        colnames(grouped.data) <- names(groups.list)
        for (i in 1:length(groups.list))
            grouped.data[, i] <- rowMeans(x[, groups.list[[i]], drop = FALSE])
        if (missing(which))
            which <- 1:length(groups.list)

        if (!pairs) {
            if (missing(refSamples)) {
                medianchip <- summaryFun(grouped.data)
            } else if (length(refSamples) == 1) {
                medianchip <- grouped.data[, refSamples]
            } else {
                medianchip <- summaryFun(grouped.data[, refSamples])
            }

            M <- sweep(grouped.data, 1, medianchip, FUN = "-")
            A <- 1/2 * sweep(grouped.data, 1, medianchip, FUN = "+")
            for (i in which) {
                if (missing(refSamples)){
                    main <- paste(levels(groups)[i], main)
                }else{
                    if (length(refSamples) == 1) {
                        if (i != refSamples)
                            main <- paste(levels(groups)[i], "vs", levels(groups)[refSamples])
                    } else {
                        main <- paste(levels(groups)[i], main)
                    }
                }
                basicMvAplot(A[, i], M[, i], main=main, xlab="A",
                             ylab="M", pch=pch, plotFun=plotFun, ...)
            }
        } else {
            basicMvApairsPlot(grouped.data[, which], transfo=transfo, plotFun=plotFun, ...)
        } ## end if(!pairs)
    }
}

basicMvAplot <- function(A, M, subset=sample(length(M), min(c(1e4, length(M)))),
                    show.statistics=TRUE, span=2/3,
                    family.loess="gaussian", cex=2,
                    plotFun=smoothScatter, addLoess=TRUE, lwd=1, lty=1,
                    loess.col="red", ...){

    fn.call <- list(...)
    nmdots <- names(fn.call)
    sigma <- IQR(M)
    mean <- median(M)
    yloc <- ifelse(is.element("ylim", nmdots), max(fn.call[['ylim']]), max(M))
    xloc <- ifelse(is.element("xlim", nmdots), max(fn.call[['xlim']]), max(A))

    plotFun(A, M, cex=cex, ...)

    if (addLoess) {
        aux <- loess(M[subset]~A[subset], degree=1, span=span,
                     family=family.loess)[['fitted']]
        o <- order(A[subset])
        A <- A[subset][o]
        M <- aux[o]
        o <- which(!duplicated(A))
        lines(approx(A[o], M[o]), col=loess.col, lwd=lwd, lty=lty)
    }
    abline(0, 0, col="blue")
    if (show.statistics){
        txt <- paste("IQR:", format(sigma, digits = 3))
        txt2 <- paste("Median:", format(mean, digits = 3))
        text(xloc, yloc, paste(txt2, txt, sep="\n"), cex=cex, adj=c(1, 1))
    }
}

basicMvApairsPlot <- function (x, labels=colnames(x), transfo=log2, span=2/3,
                       family.loess="gaussian", digits=3,
                       main="MVA plot", xlab="A", ylab="M",
                       cex=2, plotFun=smoothScatter, addLoess=TRUE,
                       parParams=list(mgp=c(0, .2, 0), mar=rep(1, 4), oma=c(1, 2, 2, 1)), ...){

    x <- transfo(x)
    J <- ncol(x)
    frame()
    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    par(mfrow=c(J, J))
    do.call(par, parParams)
    for (j in 1:(J - 1)) {
        par(mfg=c(j, j))
        plot(1, 1, type="n", xaxt="n", yaxt="n", xlab="", ylab="")
        text(1, 1, labels[j], cex=cex)
        for (k in (j + 1):J) {
            par(mfg=c(j, k))
            yy <- x[, j] - x[, k]
            xx <- (x[, j] + x[, k])/2
            sigma <- IQR(yy)
            mean <- median(yy)
            basicMvAplot(xx, yy, tck=0, plotFun=plotFun, addLoess=addLoess,
                    show.statistics=FALSE,  pch=".", xlab="", ylab="", tck=0, span=span,  ...)
            par(mfg=c(k, j))
            qs <- unique(quantile(xx, seq(0, 1, .1)))
            grps <- cut(xx, qs, include.lowest=TRUE)
            ng <- length(qs)
            boxplot(yy~grps, range=0, xaxt='n', yaxt='n', xlab='', ylab='', xlim=c(0, ng), ...)
            ## txt <- paste("IQR:", format(sigma, digits=digits))
            ## txt2 <- paste("Median:", format(mean, digits=digits))
            ## plot(c(0, 1), c(0, 1), type="n", ylab="", xlab="", xaxt="n", yaxt="n")
            ## text(0.5, 0.5, paste(txt2, txt, sep="\n"), cex=cex)
        }
    }
    par(mfg=c(J, J))
    plot(1, 1, type="n", xaxt="n", yaxt="n", xlab="",  ylab="")
    text(1, 1, labels[J], cex=cex)
    mtext(xlab, 1, outer=TRUE, cex=1.5)
    mtext(ylab, 2, outer=TRUE, cex=1.5)
    mtext(main, 3, outer=TRUE, cex=1.5)
    invisible()
}

setMethod("MAplot", "matrix",
          function(object, what=identity, transfo=identity, groups, refSamples, which,
                   pch=".", summaryFun=rowMedians, plotFun=smoothScatter,
                   main="vs pseudo-median reference chip", pairs=FALSE, ...){
              stopifnot(is.function(what))
              maplot(x=what(object), transfo=transfo, groups=groups,
                     refSamples=refSamples, which=which, pch=pch,
                     summaryFun=summaryFun, main=main, pairs=pairs, ...)
          })

readCEL <- function(cels, fid){
    g <- function(x) read.celfile(x, TRUE)[['INTENSITY']][['MEAN']]
    if (missing(fid)){
        sapply(cels, g)
    }else{
        sapply(cels, function(x) g(x)[fid])
    }
}

txtMsg <- function(..., domain=NULL, appendLF=FALSE)
    message(..., domain=domain, appendLF=appendLF)

msgOK <- function()
    message("OK")

cloneFS <- function(fs){
    nms <- sort(assayDataElementNames(fs))
    if (is(assayDataElement(fs, nms[1]), "ff_matrix")){
        if (nms[1] == "exprs" & length(nms) == 1){
            c1 <- clone(assayDataElement(fs, nms),
                           pattern=file.path(ldPath(), 'oligo-exprs-'))
            ad <- assayDataNew(exprs=c1)
            slot(fs, "assayData") <- ad
        }else if (all (nms == c("channel1", "channel2"))){
            c1 <- clone(assayDataElement(fs, 'channel1'),
                        pattern=file.path(ldPath(), 'oligo-channel1-'))
            c2 <- clone(assayDataElement(fs, 'channel2'),
                        pattern=file.path(ldPath(), 'oligo-channel2-'))
            ad <- assayDataNew(channel1=c1, channel2=c2)
            slot(fs, "assayData") <- ad
        }else{
            stop("Don't know how to clone elements: ", nms)
        }
    }
    return(fs)
}


#### Pair -> XYS
pair2xys <- function(pairFile){
  hdr <- readLines(pairFile, n=1)
  pair <- read.table(pairFile, header=TRUE, comment.char='#', sep='\t', stringsAsFactors=FALSE)
  ## For the moment, this gets only the PM probes
  out <- merge(expand.grid(X=1:1050, Y=1:4200), pair[, c('X', 'Y', 'PM')], all.x=TRUE)
  out$COUNT <- with(out, ifelse(is.na(PM), NA, 1L))
  out <- out[order(out$Y, out$X),]
  rownames(out) <- NULL
  names(out) <- c('X', 'Y', 'SIGNAL', 'COUNT')
  fout <- gsub("\\.pair$", "\\.xys", pairFile)
  writeLines(hdr, fout)
  suppressWarnings(write.table(out, file=fout, sep='\t', quote=FALSE, append=TRUE, row.names=FALSE))
  fout
}


### parallel stuff
### iExprsProbesets <- function(x, ...){
###     rns <- rownames(x)
###     grps <- split(1:length(rns), rns)
###     it <- idiv(length(grps), ...)
###     i <- 1
###     nextEl <- function(){
###         n <- nextElem(it)
###         set <- seq(i, length=n)
###         i <<- i+n
###         x[unlist(grps[set]),, drop=FALSE]
###     }
###     obj <- list(nextElem=nextEl)
###     class(obj) <- c('iExprsProbesets', 'abstractiter', 'iter')
###     obj
### }
