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
    c('medianpolish', 'plm', 'plmr', 'plmrr', 'plmrc',
      'wplm', 'wplmr', 'wplmrr', 'wplmrc')

###########
###########



## Y: nr x nc - nr: number of *probes* in probeset; nc: number of samples

## rcModelPLM/rcModelWPLM
## - Estimates: nr + nc
## - Weights..: nr x nc
## - Residuals: nr x nc
## - StdErrors: nr + nc
## - Scale....: 1

## rcModelMedianPolish
## - Estimates: nr + nc
## - Weights..: NULL
## - Residuals: nr x nc
## - StdErrors: NULL
## - Scale....: NULL

## rcModelPLM/rcModelWPLM (row.effects given)
## - Estimates: 00 + nc
## - Weights..: nr x nc
## - Residuals: nr x nc
## - StdErrors: 00 + nc
## - Scale....: nc

## rcModelPLM/rcModelWPLM (input.scale given)
## - Estimates: nr + nc
## - Weights..: nr x nc
## - Residuals: nr x nc
## - StdErrors: nr + nc
## - Scale....: 1

## rcModelPLM/rcModelWPLM (row.effects+input.scale given)
## - Estimates: 00 + nc
## - Weights..: nr x nc
## - Residuals: nr x nc
## - StdErrors: 00 + nc
## - Scale....: nc

## rcModelPLMr/rcModelPLMrr/rcModelPLMrc
## - Estimates: nr + nc
## - Weights..: nr x nc
## - Residuals: nr x nc
## - StdErrors: nr + nc
## - Scale....: NULL


## The other rcModel* functions return a scale parameter
## This equalizes output:
## - Estimates (chip and probe)
## - Weights
## - Residuals
## - StdErrors (chip and probe)
## - Scale

outputEqualizer <- function(lst){
    nsamples <- ncol(lst$Residuals)
    nprobes <- nrow(lst$Residuals)
    idx <- 1:nsamples
    ## stderr
    chipStdErrors <- probesStdErrors <- NULL
    if (nsamples+nprobes == length(lst$StdErrors)){
        chipStdErrors <- lst$StdErrors[idx]
        probesStdErrors <- lst$StdErrors[-(idx)]
    } else if (length(lst$StdErrors) == nsamples){
        chipStdErrors <- lst$StdErrors[idx]
    }
    list(chipEffects=lst$Estimates[idx],
         probeEffects=lst$Estimates[-(idx)],
         Weights=lst$Weights,
         Residuals=lst$Residuals,
         chipStdErrors=chipStdErrors,
         probesStdErrors=probesStdErrors,
         Scale=lst$Scale)
}


runSummarize <- function(mat, pnVec, transfo=log2,
                         method=summarizationMethods(),
                         ...){
    stopifnot(length(pnVec) == nrow(mat),
              is.character(pnVec),
              is.function(transfo))
    method <- match.arg(method)
    theFun <- switch(method,
                     medianpolish=rcModelMedianPolish,
                     plm=rcModelPLM,
                     plmr=rcModelPLMr,
                     plmrr=rcModelPLMrr,
                     plmrc=rcModelPLMrc,
                     wplm=rcModelWPLM,
                     wplmr=rcModelWPLMr,
                     wplmrr=rcModelWPLMrr,
                     wplmrc=rcModelWPLMrc)
    psets <- unique(pnVec)
    output <- foreach(set=psets, .packages='oligo') %dopar% {
        ## handle ff objects in here
        outputEqualizer(theFun(y=transfo(mat[pnVec == set,, drop=FALSE]), ...))
    }
}

fitProbeLevelModel <- function(object, target='core', subset, method='plm', S4=TRUE){
    ## essential to be sorted by man_fsetid, so weights/residuals can be
    ## matched to original FS object
    probeInfo <- getProbeInfo(object, target=target, field=c('fid', 'fsetid'),
                              sortBy='man_fsetid', subset=subset)
    pnVec <- as.character(probeInfo$man_fsetid)

    tmpMat <- exprs(object)[probeInfo$fid,,drop=FALSE]
    tmpMat <- backgroundCorrect(tmpMat, method='rma')
    tmpMat <- normalize(tmpMat, method='quantile')
    fit <- runSummarize(tmpMat, pnVec, method=method)
    rm(tmpMat)

    ## Chip effects - number of probesets X number samples
    chipEffects <- do.call(rbind, lapply(fit, '[[', 'chipEffects'))
    dimnames(chipEffects) <- list(unique(pnVec), sampleNames(object))

    ## Probe effects - number of probes
    probeEffects <- vector('numeric', nrow(object))
    probeEffects[probeInfo$fid] <- do.call(c, lapply(fit, '[[', 'probeEffects'))
    ## names(probeEffects) <- pnVec

    ## Weights/residuals - number of probes X number samples
    Weights <- Residuals <- array(NA, dim(object))
    Weights[probeInfo$fid,] <- do.call(rbind, lapply(fit, '[[', 'Weights'))
    Residuals[probeInfo$fid,] <- do.call(rbind, lapply(fit, '[[', 'Residuals'))
    ## dimnames(Weights) <- dimnames(Residuals) <- list(pnVec, sampleNames(object))

    ## Chip StdErrors
    chipStdErrors <- do.call(rbind, lapply(fit, '[[', 'chipStdErrors'))
    dimnames(chipStdErrors) <- list(unique(pnVec), sampleNames(object))

    ## Probes StdErrors
    probesStdErrors <- vector('numeric', nrow(object))
    probesStdErrors[probeInfo$fid] <- do.call(c, lapply(fit, '[[', 'probesStdErrors'))
    ## names(probesStdErrors) <- pnVec

    ## Scale
    Scale <- do.call(c, lapply(fit, '[[', 'Scale'))
    names(Scale) <- unique(pnVec)

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
    rm(chipEffects, probeEffects, Weights, Residuals, chipStdErrors, probesStdErrors, Scale)
    if (S4)
        out <- do.call(new, out)
    out
}

fitPLM <- function(...)
    .Deprecated('fitProbeLevelModel')


RLE <- function(obj, type=c('plot', 'values'), ylim=c(-.75, .75),
                range=0, col=darkColors(ncol(obj)), ...){
    RLE <- sweep(coefs(obj), 1, rowMedians(coefs(obj)), '-')
    type <- match.arg(type)
    if (type=='plot'){
        boxplot(as.data.frame(RLE), ylab='RLE', range=range, ylim=ylim, col=col, ...)
        abline(h=0, lty=2)
    }
    invisible(RLE)
}

NUSE <- function(obj, type=c('plot', 'values'), ylim=c(.95, 1.10),
                 range=0, col=darkColors(ncol(obj)), ...){
    if (is.null(se(obj)))
        stop('This Probe Level Model does not allow for computation of NUSE')
    NUSE <- sweep(se(obj), 1, rowMedians(se(obj)), '/')
    type <- match.arg(type)
    if (type == 'plot'){
        boxplot(as.data.frame(NUSE), ylab='NUSE', range=range, ylim=ylim, col=col, ...)
        abline(h=1, lty=2)
    }
    invisible(NUSE)
}

setMethod('image', 'oligoPLM',
          function(x, which=1, type=c('weights', 'resids', 'pos.resids', 'neg.resids', 'sign.resids'), col, main, ...){
              type <- match.arg(type)
              if (type == 'weights'){
                  theMat <- weights(x)[, which]
                  candCols <- rev(seqColors(2560))
                  candMain <- 'Weights'
              }else if (type == 'resids'){
                  theMat <- resids(x)[, which]
                  candCols <- divColors(2560)
                  candMain <- 'Residuals'
              }else if (type == 'pos.resids'){
                  theMat <- pmax(resids(x)[, which], 0)
                  candCols <- seqColors2(2560)
                  candMain <- 'Positive Residuals'
              }else if (type == 'neg.resids'){
                  theMat <- pmin(resids(x)[, which], 0)
                  candCols <- rev(seqColors(2560))
                  candMain <- 'Negative Residuals'
              }else{
                  theMat <- sign(resids(x)[, which])
                  candCols <- divColors(2)
                  candMain <- 'Sign of Residuals'
              }
              dim(theMat) <- x@geometry
              if (missing(col)){
                  col <- candCols
                  rm(candCols)
              }
              if (missing(main)){
                  main <- candMain
                  rm(candMain)
              }
              image(theMat, col=col, yaxt='n', xaxt='n', main=main, ...)
          }
)
