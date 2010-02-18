############################################
## BackgroundCorrection
############################################
rmaBgCorrectLDSnode <- function(cols, object, matInEnv){
  ## this runs on the node
  ## it (rma) bg corrects 'object' and overwrites it
  if (length(cols) > 0 ){
    if (!missing(matInEnv)){
      stopifnot(is.character(matInEnv))
      object <- get(matInEnv, envir=.oligoPkgEnv)
    }
    open(object)
    grpCols <- splitIndicesByLength(cols, oligoSamples())
    for (theCols in grpCols)
      object[, theCols] <- rma.background.correct(object[, theCols, drop=FALSE], copy=FALSE)
    close(object)
    rm(object, grpCols)
    gc()
  }
  TRUE
}

rmaBgCorrectLDSmaster <-  function(object, copy=TRUE){
  ## This runs on the master node
  stopifnot(oligoBigObjectSupport())
  if (copy){
    out <- clone(object, pattern=file.path(oligoBigObjectPath(), "oligo-rmabg-"))
  }else{
    out <- object
  }
  outName <- "out"
  sendBO2PkgEnv(out, outName)
  samplesByNode <- splitIndicesByNode(1:ncol(out))
  oLapply(samplesByNode, rmaBgCorrectLDSnode, matInEnv=outName)
  rmFromPkgEnv(outName)
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

setMethod("backgroundCorrect", "ff_matrix",
          function(object, method="rma", copy=TRUE, verbose=TRUE){
            method <- match.arg(method, "rma")
            if (verbose) cat("Background correcting... ")
            if (method == "rma"){
              out <- rmaBgCorrectLDSmaster(object, copy=copy)
            }
            if (verbose) cat("OK\n")
            return(out)
          })

############################################
## Normalization
############################################
qnTargetStatsLDSnode <- function(cols, object, matInEnv){
  ## this runs on the node
  if (!missing(matInEnv)){
    stopifnot(is.character(matInEnv))
    object <- get(matInEnv, envir=.oligoPkgEnv)
  }
  open(object)
  total <- rep(0, nrow(object))
  if (length(cols) > 0){
    for (i in cols)
      total <- total+sort(object[,i])
  }
  close(object)
  rm(object)
  list(total=total, n=length(cols))
}

qnToTargetLDSnode <- function(cols, target, object, matInEnv){
  ## this runs on the node
  if (length(cols) > 0){
    if (!missing(matInEnv)){
      stopifnot(is.character(matInEnv))
      object <- get(matInEnv, envir=.oligoPkgEnv)
    }
    open(object)
    grpCols <- splitIndicesByLength(cols, oligoSamples())
    for (theCols in grpCols)
      object[, theCols] <- normalize.quantiles.use.target(object[, theCols, drop=FALSE], target, copy=FALSE)
    close(object)
    rm(object, grpCols)
  }
  TRUE
}

quantileNormalizationLDSmaster <- function(object){
  ## this runs on the master node
  samplesByNode <- splitIndicesByNode(1:ncol(object))
  dmns <- dimnames(object)
  dimnames(object) <- NULL
  outName <- "outObj"
  sendBO2PkgEnv(object, outName)
  stats <- oLapply(samplesByNode, qnTargetStatsLDSnode, matInEnv=outName)
  totalN <- sum(sapply(stats, "[[", "n"))
  total <- rowSums(sapply(stats, "[[", "total"))
  target <- total/totalN
  rm(stats, total, totalN)
  oLapply(samplesByNode, qnToTargetLDSnode, target, matInEnv=outName)
  rmFromPkgEnv(outName)
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

setMethod("normalize", "ff_matrix",
          function(object, method="quantile", copy=TRUE, verbose=TRUE){
            method = match.arg(method, "quantile")
            if (verbose) cat("Normalizing... ")
            if (copy){
              out <- clone(object, pattern=file.path(oligoBigObjectPath(), "oligo-qn-"))
            }else{
              out <- object
            }
            if (method == "quantile"){
              quantileNormalizationLDSmaster(out)
            }
            if (verbose) cat("OK\n")
            return(out)
          })

############################################
## Summarization
############################################

basicMedianPolishBO <- function(psToSumm, inObj, outObj, probes,
                                probesets, matInEnv, matOutEnv){
  ## this runs on the node
  stopifnot(!missing(probes), !missing(probesets))
  ok <- is.character(probes) && is.character(probesets)
  if (!ok) stop("Ensure inObj and outObj have valid rownames.")
  rm(ok)
  if (length(psToSumm) > 0){
    if (!missing(matInEnv)){
      stopifnot(is.character(matInEnv))
      inObj <- get(matInEnv, envir=.oligoPkgEnv)
    }
    open(inObj)
    if (!missing(matOutEnv)){
      stopifnot(is.character(matOutEnv))
      outObj <- get(matOutEnv, envir=.oligoPkgEnv)
    }
    open(outObj)
    psList <- splitIndicesByLength(psToSumm, oligoProbesets())
    for (pss in psList){
      iIn <- unlist(pss)
      inMatrix <- inObj[iIn,, drop=FALSE]
      tmp <- basicRMA(inMatrix, pnVec=probes[iIn], normalize=FALSE,
                      background=FALSE, verbose=TRUE, destructive=TRUE)
      rm(inMatrix, iIn)
      iOut <- match(rownames(tmp), probesets)
      outObj[iOut,] <- tmp
      rm(tmp, iOut)
    }
    close(inObj)
    close(outObj)
    rm(inObj, outObj, psList)
    gc()
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

setMethod("summarize", "ff_matrix",
          function(object, probes=rownames(object), method="medianpolish", verbose=TRUE){
            if (verbose) cat("Summarizing... ")
            stopifnot(nrow(object) == length(probes))
            method <- match.arg(method, "medianpolish")
            if (method == "medianpolish"){
              probeRowByProbesets <- split(1:nrow(object), probes)
              pnsListByNode <- splitIndicesByNode(probeRowByProbesets)
              pns <- names(probeRowByProbesets)

              out <- ff(vmode="double", dim=c(length(pns), ncol(object)),
                        pattern=file.path(oligoBigObjectPath(), "oligo-mp-"))
              outMat <- "probesetLevel"
              sendBO2PkgEnv(out, outMat)

              dnmsIn <- dimnames(object)
              dimnames(object) <- NULL
              inMat <- "probeLevel"
              sendBO2PkgEnv(object, inMat)
              oLapply(pnsListByNode, basicMedianPolishBO,
                      probes=probes, probesets=pns,
                      matInEnv=inMat, matOutEnv=outMat)
              rmFromPkgEnv(outMat)
              rmFromPkgEnv(inMat)
              dimnames(object) <- dnmsIn
              dimnames(out) <- list(names(probeRowByProbesets),
                                    colnames(object))
              
            }
            if (verbose) cat("OK\n")
            return(out)
          })


basicRMAbo <- function(pmMat, pnVec, normalize=TRUE, background=TRUE,
                       bgversion=2, destructive=FALSE, verbose=TRUE,
                       pmName="pms", ...){
  dnms <- dimnames(pmMat)
  dimnames(pmMat) <- NULL
  sendBO2PkgEnv(pmMat, pmName)
  
  ## background correct
  if (background){
    if (verbose) message("Background correcting...")
    samplesByNode <- splitIndicesByNode(1:ncol(pmMat))
    oLapply(samplesByNode, rmaBgCorrectLDSnode, matInEnv=pmName)
  }

  ## normalize
  if (normalize){
    if (verbose) message("Normalizing...")
    if (!exists("samplesByNode")) 
      samplesByNode <- splitIndicesByNode(1:ncol(pmMat))
    stats <- oLapply(samplesByNode, qnTargetStatsLDSnode, matInEnv=pmName)
    totalN <- sum(sapply(stats, "[[", "n"))
    total <- rowSums(sapply(stats, "[[", "total"))
    target <- total/totalN
    rm(stats, total, totalN)
    oLapply(samplesByNode, qnToTargetLDSnode, target, matInEnv=pmName)
    rm(samplesByNode, target)
  }

  ## summarize
  if (verbose) message("Summarizing...")
  rowsByProbesets <- split(1:nrow(pmMat), pnVec)
  pnsListByNode <- splitIndicesByNode(rowsByProbesets)
  pns <- names(rowsByProbesets)
  rmaResult <- ff(vmode="double", dim=c(length(pns), ncol(pmMat)),
                  pattern=file.path(oligoBigObjectPath(), "rma-"))
  rmaName <- "rmaResult"
  sendBO2PkgEnv(rmaResult, rmaName)
  oLapply(pnsListByNode, basicMedianPolishBO, probes=pnVec,
          probesets=pns, matInEnv=pmName, matOutEnv=rmaName)
  rmFromPkgEnv(pmName)
  rmFromPkgEnv(rmaName, TRUE)
  dimnames(pmMat) <- dnms
  rm(dnms, rmaName)
  dimnames(rmaResult) <- list(pns, colnames(pmMat))
  rm(pns)
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
