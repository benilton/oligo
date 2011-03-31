basicMedianPolishBO <- function(psToSumm, inObj, outObj, probes,
                                probesets){
  ## this runs on the node
  stopifnot(!missing(probes), !missing(probesets))
  ok <- is.character(probes) && is.character(probesets)
  if (!ok) stop("Ensure inObj and outObj have valid rownames.")
  rm(ok)
  if (length(psToSumm) > 0){
    open(inObj)
    open(outObj)
    psList <- splitIndicesByLength(psToSumm, ocProbesets())
    for (pss in psList){
      iIn <- unlist(pss)
      inMatrix <- inObj[iIn,, drop=FALSE]
      tmp <- basicRMA(inMatrix, pnVec=probes[iIn], normalize=FALSE,
                      background=FALSE, verbose=FALSE, destructive=TRUE)
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

basicRMAbo <- function(pmMat, pnVec, normalize=TRUE, background=TRUE,
                       bgversion=2, destructive=FALSE, verbose=TRUE,
                       ...){
  dnms <- dimnames(pmMat)
  dimnames(pmMat) <- NULL
  
  ## background correct
  if (background){
    if (verbose) message("Background correcting... ", appendLF=FALSE)
    samplesByNode <- splitIndicesByNode(1:ncol(pmMat))
    ocLapply(samplesByNode, rmaBgCorrectLDSnode, object=pmMat,
	     neededPkgs="oligo")
    if (verbose) message("OK")
  }

  ## normalize
  if (normalize){
    if (verbose) message("Normalizing... ", appendLF=FALSE)
    if (!exists("samplesByNode")) 
      samplesByNode <- splitIndicesByNode(1:ncol(pmMat))
    stats <- ocLapply(samplesByNode, qnTargetStatsLDSnode, object=pmMat, neededPkgs="oligo")
    totalN <- sum(sapply(stats, "[[", "n"))
    total <- rowSums(sapply(stats, "[[", "total"))
    target <- total/totalN
    rm(stats, total, totalN)
    ocLapply(samplesByNode, qnToTargetLDSnode, target=target, object=pmMat, neededPkgs="oligo")
    rm(samplesByNode, target)
    if (verbose) message("OK")
  }

  ## summarize
  if (verbose) message("Summarizing... ", appendLF=FALSE)
  rowsByProbesets <- split(1:nrow(pmMat), pnVec)
  pnsListByNode <- splitIndicesByNode(rowsByProbesets)
  pns <- names(rowsByProbesets)
  rmaResult <- createFF("rma-", dim=c(length(pns), ncol(pmMat)))
  ocLapply(pnsListByNode, basicMedianPolishBO, probes=pnVec,
           probesets=pns, inObj=pmMat, outObj=rmaResult,
           neededPkgs="oligo")
  dimnames(pmMat) <- dnms
  rm(dnms)
  dimnames(rmaResult) <- list(pns, colnames(pmMat))
  rm(pns)
  if (verbose) message("OK")
  return(rmaResult)
}

setMethod("summarize", "matrix",
          function(object, probes=rownames(object), method='medianpolish', verbose=TRUE){
            stopifnot(nrow(object) == length(probes))
            method <- match.arg(method, summarizationMethods())
            if (method == "medianpolish"){
              return(basicRMA(object, pnVec=probes, normalize=FALSE,
                              background=FALSE, verbose=verbose))
            } else {
                ## all PLM methods plm / plmr / plmrr / plmrc
                ## note that the summaries will be computed on log2(object)
                return(basicPLM(object, pnVec=probes, normalize=FALSE,
                                background=FALSE, method=method, verbose=verbose))
            }
          })

## add PLM to ff_matrix
setMethod("summarize", "ff_matrix",
          function(object, probes=rownames(object), method="medianpolish", verbose=TRUE){
            stopifnot(nrow(object) == length(probes))
            method <- match.arg(method, "medianpolish")
            if (verbose) message("Summarizing... ", appendLF=FALSE)
            if (method == "medianpolish"){
              probeRowByProbesets <- split(1:nrow(object), probes)
              pnsListByNode <- splitIndicesByNode(probeRowByProbesets)
              pns <- names(probeRowByProbesets)
              out <- createFF("oligo-mp-", dim=c(length(pns), ncol(object)))
              dnmsIn <- dimnames(object)
              dimnames(object) <- NULL
              ocLapply(pnsListByNode, basicMedianPolishBO,
                      probes=probes, probesets=pns, inObj=object,
                      outObj=out, neededPkgs="oligo")
              dimnames(object) <- dnmsIn
              dimnames(out) <- list(names(probeRowByProbesets),
                                    colnames(object))
              
            }
            if (verbose) message("OK")
            return(out)
          })

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

## PLM base
plm1Probeset <- function(i, M, funPLM, nc=ncol(M)){
    res <- funPLM(M[i,,drop=FALSE])
    list(Estimates=res[['Estimates']][1:nc],
         StdErrors=res[['StdErrors']][1:nc],
         Residuals=res[['Residuals']])
}

runPLM <- function(PM, pnVec, funPLM){
    groups <- split(1:nrow(PM), pnVec)
    summaries <- lapply(groups, plm1Probeset, M=PM, funPLM=funPLM, nc=ncol(PM))
    estimates <- do.call(rbind, lapply(summaries, '[[', 'Estimates'))
    stderrors <- do.call(rbind, lapply(summaries, '[[', 'StdErrors'))
    residuals <- matrix(NA, nr=nrow(PM), nc=ncol(PM))
    residuals[unlist(groups, use.names=FALSE),] <- do.call(rbind, lapply(summaries, '[[', 'Residuals'))
    rm(summaries)
    colnames(estimates) <- colnames(stderrors) <- colnames(residuals) <- colnames(PM)
    rownames(residuals) <- NULL
    list(Estimates=estimates, StdErrors=stderrors, Residuals=residuals)
}

basicPLM <- function(pmMat, pnVec, normalize=TRUE, background=TRUE,
                     transfo=log2, method=c('plm', 'plmr', 'plmrr', 'plmrc'),
                     verbose=TRUE){
    method <- match.arg(method)
    funPLM <- switch(method,
                     plm=rcModelPLM,
                     plmr=rcModelPLMr,
                     plmrr=rcModelPLMrr,
                     plmrc=rcModelPLMrc)
    if (background)
        pmMat <- backgroundCorrect(pmMat)
    if (normalize)
        pmMat <- normalize(pmMat)
    theClass <- class(pmMat)
    
    if (verbose) message('Summarizing... ', appendLF=FALSE)
    if ('matrix' %in% theClass){
        res <- runPLM(transfo(pmMat), pnVec, funPLM)
        if (verbose) message('OK')
        return(res)
    }else if('ff_matrix' %in% theClass){
        estimates <- createFF(paste("oligo-", method, "-Estimates-", sep=""),
                           dim=c(length(unique(pnVec)), ncol(pmMat)))
        stderrors <- createFF(paste("oligo-", method, "-StdErrors-", sep=""),
                           dim=c(length(unique(pnVec)), ncol(pmMat)))
        residuals <- createFF(paste("oligo-", method, "-Residuals-", sep=""),
                           dim=dim(pmMat))
        psets <- sort(unique(pnVec))
        ocLapply(splitIndicesByNode(psets),
                 function(psets2proc, inMat, estimates, stderrors, residuals,
                          transfo, fun, pnVec, psetsOut){
                     if (length(psets2proc) > 0){
                         open(inMat)
                         open(estimates)
                         open(stderrors)
                         open(residuals)
                         ## batches of probesets to process - vector of names
                         batches <- splitIndicesByLength(psets2proc, ocProbesets())
                         for (i in 1:length(batches)){
                             psBatch <- batches[[i]]
                             idx <- which(pnVec %in% psBatch)
                             rowPsSumm <- match(sort(unique(psBatch)), psetsOut)
                             tmp <- runPLM(transfo(inMat[idx,,drop=F]), pnVec[idx], fun)
                             estimates[rowPsSumm,] <- tmp[['Estimates']]
                             stderrors[rowPsSumm,] <- tmp[['StdErrors']]
                             residuals[idx,] <- tmp[['Residuals']]
                             rm(tmp)
                         }
                         close(inMat)
                         close(estimates)
                         close(stderrors)
                         close(residuals)
                     }
                     NULL
                 },
                 pmMat, estimates, stderrors,
                 residuals, transfo, funPLM,
                 pnVec, psets, neededPkgs="oligo")
        rownames(estimates) <- rownames(stderrors) <- psets
        rownames(residuals) <- rownames(pmMat)
    }
    if (verbose) message('OK')
    colnames(estimates) <- colnames(stderrors) <- colnames(residuals) <- colnames(pmMat)
    list(Estimates=estimates, StdErrors=stderrors, Residuals=residuals)
}

summarizationMethods <- function()
    c('medianpolish', 'plm', 'plmr', 'plmrr', 'plmrc')
