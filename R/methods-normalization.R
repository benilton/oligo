qnTargetStatsLDSnode <- function(cols, object){
  ## this runs on the node
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

qnToTargetLDSnode <- function(cols, target, object){
  ## this runs on the node
  if (length(cols) > 0){
    open(object)
    grpCols <- splitIndicesByLength(cols, ocSamples())
    for (theCols in grpCols)
      object[, theCols] <- normalize.quantiles.use.target(object[, theCols, drop=FALSE], target, copy=FALSE)
    close(object)
    rm(object, grpCols)
  }
  TRUE
}

quantileNormalizationLDSmaster <- function(object, target){
  ## this runs on the master node
  samplesByNode <- splitIndicesByNode(1:ncol(object))
  dmns <- dimnames(object)
  dimnames(object) <- NULL
  if (missing(target)){
    stats <- ocLapply(samplesByNode, qnTargetStatsLDSnode, object=object, neededPkgs="oligo")
    totalN <- sum(sapply(stats, "[[", "n"))
    total <- rowSums(sapply(stats, "[[", "total"))
    target <- total/totalN
    rm(stats, total, totalN)
  }else{
    targetOK <- length(target) == nrow(object)
    if (!targetOK)
      stop("Length of target does not match nrow(object).")
  }
  ocLapply(samplesByNode, qnToTargetLDSnode, target=target, object=object, neededPkgs="oligo")
  return(object)
}


normalizationMethods <- function(){
    c('quantile', 'quantile.robust', 'quantile.in.blocks', 'qspline',
      'loess', 'invariantset', 'constant')
}

setMethod("normalize", "matrix",
          function(object, method=normalizationMethods(), copy=TRUE, verbose=TRUE, ...){
            method = match.arg(method)
            if (verbose) txtMsg("Normalizing... ")
            if (method == "quantile"){
              out <- normalize.quantiles(object, copy=copy)
            } else if (method == "quantile.robust"){
                out <- normalize.quantiles.robust(object, copy=copy, ...)
            } else if (method == "quantile.in.blocks"){
                out <- normalize.quantiles.in.blocks(object, copy=copy, ...)
            } else if (method == "qspline"){
                out <- qsplineNorm(object, verbose=FALSE, ...)
            } else if (method == "loess"){
                out <- loessNorm(object, verbose=FALSE, ...)
            } else if (method == "invariantset"){
                out <- invariantsetNorm(object, verbose=FALSE, ...)
            } else if (method == "constant"){
                out <- constantNorm(object, ...)
            }
            if (verbose) msgOK()
            return(out)
          })

setMethod("normalize", "ff_matrix",
          function(object, method=normalizationMethods(), copy=TRUE, verbose=TRUE, ...){
            method <- match.arg(method)
            if (verbose) txtMsg("Normalizing... ")
            if (copy){
              out <- clone(object, pattern=file.path(ldPath(), "oligo-qn-"))
            }else{
              out <- object
            }
            if (method == "quantile"){
              quantileNormalizationLDSmaster(out)
            } else {
                warning("Converting 'ff' object to regular matrix.",
                        immediate.=TRUE)
                tmp <- normalize(object[], method=method, copy=copy,
                                 verbose=verbose, ...)
                for (i in 1:ncol(tmp))
                    out[,i] <- tmp[,i]
                rm(tmp)
            }
            if (verbose) msgOK()
            return(out)
          })

setMethod("normalizeToTarget", "matrix",
          function(object, target, method="quantile", copy=TRUE, verbose=TRUE){
            stopifnot(!missing(target))
            method <- match.arg(method, "quantile")
            if (verbose) txtMsg("Normalizing using target... ")
            out <- normalize.quantiles.use.target(object, target, copy=copy)
            if (verbose) msgOK()
            return(out)
          })

setMethod("normalizeToTarget", "ff_matrix",
          function(object, target, method="quantile", copy=TRUE, verbose=TRUE){
            stopifnot(!missing(target))
            method <- match.arg(method, "quantile")
            if (verbose) txtMsg("Normalizing using target... ")
            if (copy){
              out <- clone(object, pattern=file.path(ldPath(), "oligo-qn-target-"))
            }else{
              out <- object
            }
            if (method == "quantile"){
              quantileNormalizationLDSmaster(out, target)
            }
            if (verbose) msgOK()
            return(out)
          })

setMethod("normalize", "FeatureSet",
          function(object, method=normalizationMethods(), copy=TRUE, verbose=TRUE, ...){
              if (copy)
                  object <- cloneFS(object)
              if (method != 'quantile.in.blocks'){
                  pm(object) <- normalize(pm(object), method=method,
                                          copy=FALSE, verbose=verbose, ...)
              }else{
                  theDots <- list(...)
                  blocks <- theDots[['blocks']]
                  pms <- pm(object)
                  if (is.null(blocks)){
                      theDots[['blocks']] <- as.factor(probeNames(object))
                  }else{
                      stopifnot(length(blocks) == nrow(pms))
                  }
                  theDots[['object']] <- pms
                  theDots[['method']] <- method
                  theDots[['copy']] <- copy
                  theDots[['verbose']] <- verbose
                  pm(object) <- do.call('normalize', theDots)
              }
              object
          })


## QSPLINE
## Working on matrices...
qsplineNorm <- function (x, target, samples, fit.iters=5,
                         min.offset=5, spline.method="natural",
                         smooth=TRUE, spar=0, p.min=0, p.max=1,
                         incl.ends=TRUE, converge=FALSE, verbose=TRUE,
                         na.rm=FALSE){
    if (missing(target))
        target <- exp(rowMeans(log(x)))
    nr <- nrow(x)
    nc <- ncol(x)
    if (missing(samples)){
        samples <- max(round(nr/1000), 100)
    }else if (samples < 1){
        samples <- round(samples * nr)
    }
    p <- (1:samples)/samples
    p <- p[which(p <= p.max) & which(p >= p.min)]
    samples <- length(p)
    k <- fit.iters
    y.n <- ifelse(na.rm, sum(!is.na(target)), length(target))
    py.inds <- as.integer(p * y.n)
    y.offset <- round(py.inds[1]/fit.iters)
    if (y.offset <= min.offset) {
        y.offset <- min.offset
        k <- round(py.inds[1]/min.offset)
    }
    if (k <= 1) {
        warning("'k' found is non-sense. using default 'fit.iter'")
        k <- fit.iters
    }
    y.offset <- c(0, array(y.offset, (k - 1)))
    y.order <- order(target)
    fx <- matrix(0, nr, nc)
    if (verbose)
        message("Samples=", samples, " - k=", k, " - first=", py.inds[1])
    for (i in 1:nc) {
        if (na.rm){
            x.valid <- which(!is.na(x[, i]))
        }else{
            x.valid <- 1:nr
        }
        nr <- length(x.valid)
        px.inds <- as.integer(p * nr)
        x.offset <- round(px.inds[1]/fit.iters)
        if (x.offset <= min.offset) {
            x.offset <- min.offset
            k <- min(round(px.inds[1]/min.offset), k)
        }
        x.offset <- c(0, array(x.offset, (k - 1)))
        x.order <- order(x[, i])
        y.inds <- py.inds
        x.inds <- px.inds
        for (j in 1:k) {
            y.inds <- y.inds - y.offset[j]
            x.inds <- x.inds - x.offset[j]
            ty.inds <- y.inds
            tx.inds <- x.inds
            if (verbose)
                message("sampling(array=", i, " iter=", j, " off=",
                        paste(x.inds[1], -x.offset[j], y.inds[1], -y.offset[j]),")")
            if (converge) {
                ty.inds <- as.integer(c(1, y.inds))
                tx.inds <- as.integer(c(1, x.inds))
                if (j > 1) {
                  ty.inds <- c(ty.inds, y.n)
                  tx.inds <- c(tx.inds, nr)
                }
            }
            qy <- target[y.order[ty.inds]]
            qx <- x[x.order[tx.inds], i]
            if (smooth) {
                sspl <- smooth.spline(qx, qy, spar=spar)
                qx <- sspl$x
                qy <- sspl$y
            }
            fcn <- splinefun(qx, qy, method=spline.method)
            fx[x.valid, i] <- fx[x.valid, i] + fcn(x[x.valid, i])/k
        }
        if (na.rm){
            invalid <- which(is.na(x[, i]))
            fx[invalid, i] <- NA
        }
    }
    return(fx)
}

## this is 1 sample
invariantsetV <- function (data, ref, prd.td = c(0.003, 0.007)){
    np <- length(data)
    r.ref <- rank(ref)
    r.array <- rank(data)
    prd.td.adj <- prd.td * 10
    i.set <- rep(TRUE, np)
    ns <- np
    ns.old <- ns + 50 + 1
    while ((ns.old - ns) > 50) {
        air <- (r.ref[i.set] + r.array[i.set])/(2 * ns)
        prd <- abs(r.ref[i.set] - r.array[i.set])/ns
        threshold <- (prd.td.adj[2] - prd.td[1]) * air + prd.td.adj[1]
        i.set[i.set] <- prd < threshold
        ns.old <- ns
        ns <- sum(i.set)
        if (prd.td.adj[1] > prd.td[1]) 
            prd.td.adj <- prd.td.adj * 0.9
    }
    n.curve <- smooth.spline(ref[i.set], data[i.set])
    approx(n.curve$y, n.curve$x, xout=data, rule=2)$y
}

## Invariant Set Normalization
## for matrices
invariantsetNorm <- function(mat, prd.td=c(0.003, 0.007),
                             baseline.type=c("mean", "median", "pseudo-mean", "pseudo-median"),
                             verbose=TRUE) {
    baseline.type <- match.arg(baseline.type)
    nc <- ncol(mat)
    if (baseline.type == "mean") {
        refindex <- trunc(median(rank(colMeans(mat))))
        baseline.chip <- mat[, refindex]
    }else if (baseline.type == "median") {
        refindex <- trunc(median(rank(rowMedians(t(mat)))))
        baseline.chip <- mat[, refindex]
    }else if (baseline.type == "pseudo-mean") {
        refindex <- 0
        baseline.chip <- rowMeans(mat)
    }else if (baseline.type == "pseudo-median") {
        refindex <- 0
        baseline.chip <- rowMedians(mat)
    }
    if (verbose) message("Baseline array: ", ifelse(refindex == 0, 'ALL', refindex))
    for (i in setdiff(1:nc, refindex)) {
        if (verbose) message("Normalizing array ", i, "... ", appendLF=FALSE)
        mat[, i] <- invariantsetV(mat[, i], baseline.chip, prd.td)
        if (verbose) message("Done.")
    }
    return(mat)
}

## loess normalization helper
loessNormV <- function(v1, v2, subset, span, degree, weights, family){
    x <- (v1+v2)/2
    index <- c(which.min(x), subset, which.max(x))
    xx <- x[index]
    yy <- v1[index]-v2[index]
    aux <- loess(yy ~ xx, span=span, degree=degree, 
                 weights=weights, family=family)
    predict(aux, data.frame(xx=x))
}

## loess normalization
## matrices
loessNorm <- function (mat, subset, epsilon=0.01, maxit=1, log.it=TRUE,
                       verbose=TRUE, span=2/3, family.loess="symmetric"){
    II <- nrow(mat)
    J <- ncol(mat)
    if (missing(subset))
        subset <- sample(II, min(c(5000, II)))
    if (log.it)
        mat <- log2(mat)

    change <- epsilon + 1
    iter <- 0
    w <- c(0, rep(1, length(subset)), 0)
    while (iter < maxit) {
        iter <- iter + 1
        means <- matrix(0, II, J)
        for (j in 1:(J - 1)) {
            for (k in (j + 1):J) {
                aux <- loessNormV(v1=mat[,j], v2=mat[,k], subset=subset,
                                  span=span, degree=1, weights=w,
                                  family=family.loess)/J
                means[, j] <- means[, j] + aux
                means[, k] <- means[, k] - aux
                if (verbose) 
                    message("Done with ", j, " vs ", k, " in iteration ", iter)
            }
        }
        mat <- mat - means
        change <- max(colMeans((means[subset, ])^2))
        if (verbose) 
            message(sprintf("Iteration: %02d - Change: %02.8f", iter, change))
    }
    if ((change > epsilon) & (maxit > 1)) 
        warning("No convergence after ", maxit, " iterations.")
    if (log.it) {
        return(2^mat)
    }else{
        return(mat)
    }
}

## Constant (scaling) normalization
## matrices
constantNorm <- function(mat, ref=1, FUN=colMeans, na.rm=TRUE){
    smpK <- FUN(mat, na.rm=na.rm)
    smpK <- smpK[ref]/smpK
    sweep(mat, 2, smpK, "*")
}


## TODO IN PARALLEL
## this will run only on matrices for now
getDefaultParsNormAdv <- function(method){
    if (method == 'quantile'){
        pars <- list(type='pmonly')
    }else if (method == 'quantile.probeset') {
        pars <- list(type='pmonly', use.median=FALSE, use.log2=TRUE)
    }else if (method == 'scaling'){
        pars <- list(type='pmonly', scaling.trim=0.02, scaling.baseline=
                     -1, log.scalefactors=FALSE)
    }else if (method == 'quantile.robust'){
        pars <- list(type='pmonly', use.median=FALSE, use.log2=FALSE,
                     remove.extreme=TRUE, n.remove=1, weights='huber')
    }else{
        stop("Method '", method, "' not supported.")
    }
    return(pars)
}

normalizeAdv <- function(pmMat, mmMat, pnVec, normType, normPar, verbose=TRUE){
    ## types of normalization
    allowed <- c('quantile', 'quantile.probeset', 'scaling', 'quantile.robust')
    normType <- match.arg(normType, allowed)
        
    ## regardless the normType, parameters must be passed as a list
    ## normPar with (at least) a 'type' element
    stopifnot(is.list(normPar), 'type' %in% names(normPar))
    types <- c('pmonly', 'mmonly', 'together', 'separate')
    type <- match.arg(normPar[['type']], types)

    if (normType == 'quantile.probeset'){
        stopifnot(all(c('use.median', 'use.log2') %in% names(normPar)))
        normPar[['use.median']] <- as.integer(normPar[['use.median']])
        normPar[['use.log2']] <- as.integer(normPar[['use.log2']])
    }else if (normType == 'scaling'){
        req <- c('scaling.trim', 'scaling.baseline', 'log.scalefactors')
        stopifnot(all(req %in% names(normPar)))
        normPar[['scaling.trim']] <- as.double(normPar[['scaling.trim']])
        normPar[['scaling.baseline']] <- as.integer(normPar[['scaling.baseline']])
        normPar[['log.scalefactors']] <- as.integer(normPar[['log.scalefactors']])
    }else if (normType == 'quantile.robust'){
        req <- c('use.median', 'use.log2', 'remove.extreme', 'n.remove', 'weights')
        stopifnot(all(req %in% names(normPar)))
        normPar[['use.median']] <- as.integer(normPar[['use.median']])
        normPar[['use.log2']] <- as.integer(normPar[['use.log2']])
        normPar[['remove.extreme']] <- as.logical(normPar[['remove.extreme']])
        stopifnot(is.character(normPar[['weights']]) || is.numeric(normPar[['weights']]))
    }
    .Call('pp_normalize', pmMat, mmMat, pnVec, nrow(pmMat), normType,
          normPar, verbose)
}


