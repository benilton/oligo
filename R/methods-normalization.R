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
              out <- clone(object, pattern=file.path(ldPath(), "oligo-qn-"))
            }else{
              out <- object
            }
            if (method == "quantile"){
              quantileNormalizationLDSmaster(out)
            }
            if (verbose) cat("OK\n")
            return(out)
          })

setMethod("normalizeToTarget", "matrix",
          function(object, target, method="quantile", copy=TRUE, verbose=TRUE){
            stopifnot(!missing(target))
            method <- match.arg(method, "quantile")
            if (verbose) cat("Normalizing using target... ")
            out <- normalize.quantiles.use.target(object, target, copy=copy)
            return(out)
          })

setMethod("normalizeToTarget", "ff_matrix",
          function(object, target, method="quantile", copy=TRUE, verbose=TRUE){
            stopifnot(!missing(target))
            method <- match.arg(method, "quantile")
            if (verbose) cat("Normalizing using target... ")
            if (copy){
              out <- clone(object, pattern=file.path(ldPath(), "oligo-qn-target-"))
            }else{
              out <- object
            }
            if (method == "quantile"){
              quantileNormalizationLDSmaster(out, target)
            }
            return(out)
          })

