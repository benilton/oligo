############################################
## utilities
############################################
getSomeRowsAllCols <- function(cols, rows, inObj, outObj,
                               path=oligoBigObjectPath(), inMatLoaded,
                               outMatLoaded){
  ## cols inObj match cols outObj
  ## rows inObj do NOT match rows outObj
  envir <- .oligoPkgEnv
  if (length(cols) > 0){
    if (missing(inMatLoaded)){
      inMat <- attach.big.matrix(inObj, backingpath=path)
    }else{
      stopifnot(is.character(inMatLoaded))
      inMat <- get(inMatLoaded, envir=envir)
    }
    if (missing(outMatLoaded)){
      outMat <- attach.big.matrix(outObj, backingpath=path)
    }else{
      stopifnot(is.character(outMatLoaded))
      outMat <- get(outMatLoaded, envir=envir)
    }
    grpCols <- splitIndicesByLength(cols, oligoSamples())
    for (theCols in grpCols)
      outMat[, theCols] <- inMat[rows, theCols, drop=FALSE]
    rm(inMat, outMat, grpCols, theCols)
    gc()
  }
  TRUE
}

getAllRowsSomeCols <- function(cols, allCols, inObj, outObj,
                               path=oligoBigObjectPath(), inMatLoaded,
                               outMatLoaded){
  ## cols inObj do NOT match cols outObj
  ## rows inObj match rows outObj
  envir <- .oligoPkgEnv
  if (length(cols) > 0){
    if (missing(inMatLoaded)){
      inMat <- attach.big.matrix(inObj, backingpath=path)
    }else{
      stopifnot(is.character(inMatLoaded))
      inMat <- get(inMatLoaded, envir=envir)
    }
    if (missing(outMatLoaded)){
      outMat <- attach.big.matrix(outObj, backingpath=path)
    }else{
      stopifnot(is.character(outMatLoaded))
      outMat <- get(outMatLoaded, envir=envir)
    }
    grpCols <- splitIndicesByLength(cols, oligoSamples())
    for (theCols in grpCols){
      idx <- match(theCols, allCols)
      outMat[, idx] <- inMat[, theCols, drop=FALSE]
    }
    rm(inMat, outMat, grpCols, idx, theCols)
    gc()
  }
  TRUE
}

getSomeRowsSomeCols <- function(cols, allCols, rows, inObj, outObj,
                                path=oligoBigObjectPath(), inMatLoaded,
                                outMatLoaded){
  if (length(cols) > 0){
    envir <- .oligoPkgEnv
    if (missing(inMatLoaded)){
      inMat <- attach.big.matrix(inObj, backingpath=path)
    }else{
      stopifnot(is.character(inMatLoaded))
      inMat <- get(inMatLoaded, envir=envir)
    }
    if (missing(outMatLoaded)){
      outMat <- attach.big.matrix(outObj, backingpath=path)
    }else{
      stopifnot(is.character(outMatLoaded))
      outMat <- get(outMatLoaded, envir=envir)
    }
    grpCols <- splitIndicesByLength(cols, oligoSamples())
    for (theCols in grpCols){
      idx <- match(theCols, allCols)
      outMat[, idx] <- inMat[rows, theCols, drop=FALSE]
    }
    rm(inMat, outMat, grpCols, idx, theCols)
    gc()
  }
  TRUE
}

subsetBO <- function(rows, cols, object, path=oligoBigObjectPath(),
                      fname=basename(tempfile("oligo-", path)),
                      dfname=paste(fname, "desc", sep="."),
                      returnIfExistent=FALSE, nameInEnv="subMatrix",
                      clean=TRUE){
  ## hypothesis: nr >> nc
  objDesc <- file.path(path, dfname)
  objFile <- file.path(path, fname)
  if (file.exists(objDesc) && file.exists(objFile) && returnIfExistent){
    out <- attach.big.matrix(dfname, backingpath=path)
  }else{
    rm(objDesc, objFile)
    stopifnot(validWorkDir(path))
    stopifnot(!file.exists(file.path(path, fname)))
    stopifnot(!file.exists(file.path(path, dfname)))
    inMat <- attach.big.matrix(object, backingpath=path)
    ## if cluster is up and running
    ## send 'out' to nodes
    ## so we don't need to load it
    ## every time
    sendBO2PkgEnv(object, path, "inMat")
    
    rns <- rownames(inMat)
    cns <- colnames(inMat)
    if (missing(rows)){
      nr <- nrow(inMat)
    }else{
      nr <- length(rows)
      rns <- rns[rows]
    }
    if (missing(cols)){
      nc <- ncol(inMat)
    }else{
      nc <- length(cols)
      cns <- cns[cols]
    }
    dns <- list(rns, cns)
    rm(cns, rns)
    out <- big.matrix(nr, nc, type=bigmemory::typeof(inMat),
                      backingpath=path, backingfile=fname,
                      descriptorfile=dfname, dimnames=dns,
                      separated=TRUE)
    rm(dns)

    ## if cluster is up and running
    ## send 'out' to nodes
    ## so we don't need to load it
    ## every time
    sendBO2PkgEnv(describe(out), path, nameInEnv)
    
    if ((!missing(rows)) && missing(cols)){
      ## cols iObj match cols oOubj
      samplesByNode <- splitIndicesByNode(1:ncol(out))
      oLapply(samplesByNode, getSomeRowsAllCols, rows, object,
              describe(out), path=path, inMatLoaded="inMat",
              outMatLoaded=nameInEnv)
    }else if (missing(rows) && (!(missing(cols)))){
      samplesByNode <- splitIndicesByNode(cols)
      oLapply(samplesByNode, getAllRowsSomeCols, cols, object,
              describe(out), path=path, inMatLoaded="inMat",
              outMatLoaded=nameInEnv)
    }else if ((!missing(rows)) && (!missing(cols))){
      samplesByNode <- splitIndicesByNode(cols)
      oLapply(samplesByNode, getSomeRowsSomeCols, cols, rows, object,
              describe(out), path=path, inMatLoaded="inMat",
              outMatLoaded=nameInEnv)
    }else{
      stop("Must specify at least one of 'rows'/'cols'")
    }
  }
  rmFromPkgEnv("inMat")
  if (clean)
    rmFromPkgEnv(nameInEnv)
  return(out)
}


############################################
## BackgroundCorrection
############################################
rmaBgCorrectBO <- function(cols, object, path=oligoBigObjectPath(),
                           matInEnv){
  ## this runs on the node
  if (length(cols) > 0 ){
    if (missing(matInEnv)){
      object <- attach.big.matrix(object, backingpath=path)
    }else{
      stopifnot(is.character(matInEnv))
      object <- get(matInEnv, envir=.oligoPkgEnv)
    }
    grpCols <- splitIndicesByLength(cols, oligoSamples())
    for (theCols in grpCols)
      object[, theCols] <- rma.background.correct(object[, theCols, drop=FALSE])
    rm(object, grpCols)
  }
  TRUE
}

rma.background.correct.bigObject <-  function(object, copy=TRUE, path=oligoBigObjectPath()){
  ## This runs on the master node
  stopifnot(oligoBigObjectSupport())
  if (copy){
    out <- deepcopy(object, type=bigmemory::typeof(object))
  }else{
    out <- object
  }
  rm(object)
  samplesByNode <- splitIndicesByNode(1:ncol(out))
  oLapply(samplesByNode, rmaBgCorrectBO, describe(out), path=path)
  rm(samplesByNode)
  return(out)
}

## This is visible for user
setMethod("backgroundCorrect", "matrix",
          function(object, method="rma", copy=TRUE, verbose=TRUE){
            method <- match.arg(method, "rma")
            if (verbose) cat("Background correcting... ")
            if (method == "rma"){
              out <- rma.background.correct(object, copy=copy)
            }
            if (verbose) cat("OK\n")
            return(out)
          })

setMethod("backgroundCorrect", "big.matrix",
          function(object, method="rma", copy=TRUE, verbose=TRUE){
            method <- match.arg(method, "rma")
            if (verbose) cat("Background correcting... ")
            if (method == "rma"){
              out <- rma.background.correct.bigObject(object, copy=copy)
            }
            if (verbose) cat("OK\n")
            return(out)
          })

############################################
## Normalization
############################################
qnTargetStats <- function(cols, object, path=oligoBigObjectPath(), matInEnv){
  ## this runs on the node
  if (missing(matInEnv)){
    object <- attach.big.matrix(object, backingpath=path)
  }else{
    stopifnot(is.character(matInEnv))
    object <- get(matInEnv, envir=.oligoPkgEnv)
  }
  total <- rep(0, nrow(object))
  if (length(cols) > 0){
    for (i in cols)
      total <- total+sort(object[,i])
  }
  rm(object)
  list(total=total, n=length(cols))
}

qnToTargetBO <- function(cols, target, object, path=oligoBigObjectPath(), matInEnv){
  ## this runs on the node
  if (length(cols) > 0){
    if (missing(matInEnv)){
      object <- attach.big.matrix(object, backingpath=path)
    }else{
      stopifnot(is.character(matInEnv))
      object <- get(matInEnv, envir=.oligoPkgEnv)
    }
    grpCols <- splitIndicesByLength(cols, oligoSamples())
    for (theCols in grpCols)
      object[, theCols] <- normalize.quantiles.use.target(object[, theCols, drop=FALSE], target)
    rm(object, grpCols)
  }
  TRUE
}

quantileNormalizationBO <- function(object, path=oligoBigObjectPath()){
  ## this runs on the master node
  samplesByNode <- splitIndicesByNode(1:ncol(object))
  stats <- oLapply(samplesByNode, qnTargetStats, describe(object), path)
  totalN <- sum(sapply(stats, "[[", "n"))
  total <- rowSums(sapply(stats, "[[", "total"))
  target <- total/totalN
  rm(stats, total, totalN)
  oLapply(samplesByNode, qnToTargetBO, target, describe(object), path)
  return(object)
}

setMethod("normalize", "matrix",
          function(object, method="quantile", copy=TRUE, verbose=TRUE){
            method = match.arg(method, "quantile")
            if (verbose) cat("Normalizing... ")
            if (method == "quantile"){
              out <- normalize.quantiles(object, copy=copy)
            }
            if (verbose) cat("OK\n")
            return(out)
          })

setMethod("normalize", "big.matrix",
          function(object, method="quantile", copy=TRUE, verbose=TRUE){
            method = match.arg(method, "quantile")
            if (verbose) cat("Normalizing... ")
            if (copy){
              out <- deepcopy(object, type=bigmemory::typeof(object))
            }else{
              out <- object
            }
            if (method == "quantile"){
              quantileNormalizationBO(out, path=oligoBigObjectPath())
            }
            if (verbose) cat("OK\n")
            return(out)
          })

############################################
## Summarization
############################################
basicMedianPolishBO <- function(psToSumm, inObj, outObj, probes,
                                probesets, path=oligoBigObjectPath(),
                                matInEnv, matOutEnv){
  ## this runs on the node
  stopifnot(!missing(probes), !missing(probesets))
  ok <- is.character(probes) && is.character(probesets)
  if (!ok) stop("Ensure inObj and outObj have valid rownames.")
  rm(ok)
  if (length(psToSumm) > 0){
    if (missing(matInEnv)){
      inObj <- attach.big.matrix(inObj, backingpath=path)
    }else{
      stopifnot(is.character(matInEnv))
      inObj <- get(matInEnv, envir=.oligoPkgEnv)
    }
    if (missing(matOutEnv)){
      outObj <- attach.big.matrix(outObj, backingpath=path)
    }else{
      stopifnot(is.character(matOutEnv))
      outObj <- get(matOutEnv, envir=.oligoPkgEnv)
    }
    psList <- splitIndicesByLength(psToSumm, oligoProbesets())
    for (pss in psList){
      iIn <- unlist(pss)
      tmp <- basicRMA(inObj[iIn,, drop=FALSE], pnVec=probes[iIn],
                      normalize=FALSE, background=FALSE, verbose=TRUE)
      iOut <- match(rownames(tmp), probesets)
      outObj[iOut,] <- tmp
      rm(tmp)
    }
    rm(inObj, outObj, psList)
  }
  TRUE
}

setMethod("summarize", "matrix",
          function(object, probes=rownames(object), method="medianpolish", verbose=TRUE){
            stopifnot(nrow(object) == length(probes))
            method <- match.arg(method, "medianpolish")
            if (method == "medianpolish"){
              return(basicRMA(object, pnVec=probes, normalize=FALSE,
                              background=FALSE, verbose=verbose))
            }
          })

setMethod("summarize", "big.matrix",
          function(object, probes=rownames(object), method="medianpolish", verbose=TRUE){
            if (verbose) cat("Summarizing... ")
            stopifnot(nrow(object) == length(probes))
            method <- match.arg(method, "medianpolish")
            if (method == "medianpolish"){
              probeRowByProbesets <- split(1:nrow(object), probes)
              pnsListByNode <- splitIndicesByNode(probeRowByProbesets)
              pns <- names(probeRowByProbesets)
              
              out <- big.matrix(length(pns), ncol(object),
                                type="double",
                                dimnames=list(names(probeRowByProbesets),
                                colnames(object)))

              oLapply(pnsListByNode, basicMedianPolishBO,
                      inObj=describe(object), outObj=describe(out),
                      probes=probes, probesets=pns,
                      path=oligoBigObjectPath())
            }
            if (verbose) cat("OK\n")
            return(out)
          })


basicRMAbo <- function(pmMat, pnVec, normalize=TRUE, background=TRUE,
                       bgversion=2, destructive=FALSE, verbose=TRUE,
                       pmName="pms", path=oligoBigObjectPath(), dUID, ...){
  ## background correct
  if (background){
    if (verbose) message("Background correcting...")
    samplesByNode <- splitIndicesByNode(1:ncol(pmMat))
    oLapply(samplesByNode, rmaBgCorrectBO, describe(pmMat), path=path, matInEnv=pmName)
  }

  ## normalize
  if (normalize){
    if (verbose) message("Normalizing...")
    stats <- oLapply(samplesByNode, qnTargetStats, describe(pmMat), path, matInEnv=pmName)
    totalN <- sum(sapply(stats, "[[", "n"))
    total <- rowSums(sapply(stats, "[[", "total"))
    target <- total/totalN
    rm(stats, total, totalN)
    oLapply(samplesByNode, qnToTargetBO, target, describe(pmMat), path, matInEnv=pmName)
  }

  ## summarize
  if (verbose) message("Summarizing...")
  probeRowByProbesets <- split(1:nrow(pmMat), pnVec)
  pnsListByNode <- splitIndicesByNode(probeRowByProbesets)
  pns <- names(probeRowByProbesets)
  rmaFile <- paste("rma-", dUID, sep="")
  rmaResult <- big.matrix(length(pns), ncol(pmMat), type="double",
                          dimnames=list(names(probeRowByProbesets),
                          colnames(pmMat)), backingfile=rmaFile,
                          descriptorfile=paste(rmaFile, "desc",
                          sep="."), backingpath=path, separated=TRUE)
  rmaName <- "rmaResult"
  sendBO2PkgEnv(describe(rmaResult), path, rmaName)
  oLapply(pnsListByNode, basicMedianPolishBO, inObj=describe(pmMat),
          outObj=describe(rmaResult), probes=pnVec, probesets=pns,
          path=oligoBigObjectPath(), matInEnv=pmName, matOutEnv=rmaName)
  rmFromPkgEnv(rmaName)
  return(rmaResult)
}

oligoNodeMem <- function(nProbesPerProbeset=10, nProbesets=25000,
                         nProbesetsPerNode=10000, nSamples=7000,
                         nSamplesPerNode=100, verbose=TRUE){
  nProbes <- nProbesPerProbeset*nProbesets
  procSamples <- nSamplesPerNode*nProbes/(2^27)
  procProbesets <- nProbesetsPerNode*nProbesPerProbeset*nSamples/(2^27)

  if (verbose){
    message(sprintf("Usage when processing samples..: %2.2f GB/node.", procSamples))
    message(sprintf("Usage when processing probesets: %2.2f GB/node.", procProbesets))
  }
  invisible(max(c(procSamples, procProbesets)))
}
