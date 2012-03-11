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

## summarizationMethods <- function()
##     c('medianpolish', 'plm', 'plmr', 'plmrr', 'plmrc',
##       'wplm', 'wplmr', 'wplmrr', 'wplmrc')

summarizationMethods <- function()
    c('medianpolish', 'plm')

###########
###########



## Y: nr x nc - nr: number of *probes* in probeset; nc: number of samples

## rcModelPLM/rcModelWPLM (supported)
## - Estimates: nc + nr
## - Weights..: nr x nc
## - Residuals: nr x nc
## - StdErrors: nc + nr
## - Scale....: 1

## rcModelMedianPolish (supported)
## - Estimates: nc + nr
## - Weights..: NULL
## - Residuals: nr x nc
## - StdErrors: NULL
## - Scale....: NULL

## rcModelPLMr/rcModelPLMrr/rcModelPLMrc (supported)
## - Estimates: nc + nr
## - Weights..: nr x nc
## - Residuals: nr x nc
## - StdErrors: nc + nr
## - Scale....: NULL

## rcModelPLM/rcModelWPLM (input.scale given) (supported)
## - Estimates: nc + nr
## - Weights..: nr x nc
## - Residuals: nr x nc
## - StdErrors: nc + nr
## - Scale....: 1

## rcModelPLM/rcModelWPLM (row.effects given)
## - Estimates: 00 + nc
## - Weights..: nr x nc
## - Residuals: nr x nc
## - StdErrors: 00 + nc
## - Scale....: nc

## rcModelPLM/rcModelWPLM (row.effects+input.scale given)
## - Estimates: 00 + nc
## - Weights..: nr x nc
## - Residuals: nr x nc
## - StdErrors: 00 + nc
## - Scale....: nc


## The other rcModel* functions return a scale parameter
## This equalizes output:
## - Estimates (chip and probe)
## - Weights
## - Residuals
## - StdErrors (chip and probe)
## - Scale

getFromListAsVector <- function(lst, elem, idx){
    lapply(lst, function(x, idx) x[[elem]][idx], idx)
}

outputEqualizer <- function(lst){
    nsamples <- ncol(lst[[1]]$Residuals)
    idx <- 1:nsamples
    theChipCoefs <- do.call(rbind, getFromListAsVector(lst, 'Estimates', idx))
    theProbeCoefs <- unlist(getFromListAsVector(lst, 'Estimates', -idx))
    theChipSE <- do.call(rbind, getFromListAsVector(lst, 'StdErrors', idx))
    theProbeSE <- unlist(getFromListAsVector(lst, 'StdErrors', -idx))
    theWeights <- do.call(rbind, lapply(lst, '[[', 'Weights'))
    theResiduals <- do.call(rbind, lapply(lst, '[[', 'Residuals'))
    theScales <- unlist(lapply(lst, '[[', 'Scale'))

    list(chipEffects=theChipCoefs, probeEffects=theProbeCoefs,
         Weights=theWeights, Residuals=theResiduals,
         chipStdErrors=theChipSE, probesStdErrors=theProbeSE,
         Scale=theScales)
}


runSummarize <- function(mat, pnVec, transfo=log2,
                         method=summarizationMethods()){
    stopifnot(length(pnVec) == nrow(mat),
              is.character(pnVec),
              is.function(transfo))
    method <- match.arg(method)
    theFun <- switch(method,
                     medianpolish=subrcModelMedianPolish,
                     plm=subrcModelPLM)
    out <- foreach(submat=iExprsProbesets(mat, chunkSize=ocProbesets()), .packages='preprocessCore') %dopar% {
        theFun(y=transfo(submat), rownames(submat))
    }
    out <- unlist(out, recursive=FALSE)
    outputEqualizer(out)
}

fitProbeLevelModel <- function(object, target='core', subset, method='plm', S4=TRUE){
    ## essential to be sorted by man_fsetid, so weights/residuals can be
    ## matched to original FS object
    probeInfo <- getProbeInfo(object, target=target, field=c('fid', 'fsetid'),
                              sortBy='man_fsetid', subset=subset)
    probeInfo$man_fsetid <- as.character(probeInfo$man_fsetid)

    tmpMat <- exprs(object)[probeInfo$fid,,drop=FALSE]
    ## rownames below is really important for parallelization
    rownames(tmpMat) <- probeInfo$man_fsetid
    tmpMat <- backgroundCorrect(tmpMat, method='rma')
    tmpMat <- normalize(tmpMat, method='quantile')
    fit <- runSummarize(tmpMat, probeInfo$man_fsetid, method=method)
    rm(tmpMat)

    chipEffects <- fit$chipEffects
    probeEffects <- fit$probeEffects
    Weights <- Residuals <- array(NA, dim(object))
    Weights[probeInfo$fid,] <- fit$Weights
    Residuals[probeInfo$fid,] <- fit$Residuals
    chipStdErrors <- fit$chipStdErrors
    probeStdErrors <- fit$probesStdErrors
    Scale <- fit$Scale
    rm(fit)

    ## fix residuals/weights to have the array dims

    out <- list(Class='oligoPLM',
                chip.coefs=chipEffects,
                probe.coefs=probeEffects,
                weights=Weights,
                residuals=Residuals,
                se.chip.coefs=chipStdErrors,
                se.probe.coefs=probesStdErrors,
                residualSE=Scale,
                geometry=geometry(object),
                method=method,
                manufacturer=manufacturer(object),
                annotation=annotation(object),
                narrays=ncol(chipEffects),
                nprobes=nrow(probeInfo),
                nprobesets=nrow(chipEffects))
    rm(chipEffects, probeEffects, Weights, Residuals, chipStdErrors,
       probesStdErrors, Scale)
    if (S4)
        out <- do.call(new, out)
    out
}

