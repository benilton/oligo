getSomeRowsAllCols <- function(cols, rows, inMat, outMat){
  ## cols inObj match cols outObj
  ## rows inObj do NOT match rows outObj
  if (length(cols) > 0){
    open(inMat)
    open(outMat)
    grpCols <- splitIndicesByLength(cols, ocSamples())
    for (theCols in grpCols)
      outMat[, theCols] <- inMat[rows, theCols, drop=FALSE]
    close(inMat)
    close(outMat)
    rm(inMat, outMat, grpCols, theCols)
    gc()
  }
  TRUE
}

getAllRowsSomeCols <- function(cols, allCols, inMat, outMat){
  ## cols inObj do NOT match cols outObj
  ## rows inObj match rows outObj
  if (length(cols) > 0){
    open(inMat)
    open(outMat)
    grpCols <- splitIndicesByLength(cols, ocSamples())
    for (theCols in grpCols){
      idx <- match(theCols, allCols)
      outMat[, idx] <- inMat[, theCols, drop=FALSE]
    }
    close(inMat)
    close(outMat)
    rm(inMat, outMat, grpCols, idx, theCols)
    gc()
  }
  TRUE
}

getSomeRowsSomeCols <- function(cols, allCols, rows, inMat, outMat){
  if (length(cols) > 0){
    open(inMat)
    open(outMat)
    grpCols <- splitIndicesByLength(cols, ocSamples())
    for (theCols in grpCols){
      idx <- match(theCols, allCols)
      outMat[, idx] <- inMat[rows, theCols, drop=FALSE]
    }
    close(inMat)
    close(outMat)
    rm(inMat, outMat, grpCols, idx, theCols)
    gc()
  }
  TRUE
}

ffSubset <- function(rows, cols, object, prefix="oligo-",
                     nameInEnv="subMatrix", clean=TRUE){
  ## This runs on the master node

  rns <- rownames(object)
  cns <- colnames(object)
  dnmsIn <- dimnames(object)
  dimnames(object) <- NULL

  if (missing(rows)){
    nr <- nrow(object)
  }else{
    stopifnot(is.numeric(rows))
    nr <- length(rows)
    rns <- rns[rows]
  }
  if (missing(cols)){
    nc <- ncol(object)
  }else{
    stopifnot(is.numeric(cols))
    nc <- length(cols)
    cns <- cns[cols]
  }

  out <- createFF(prefix, dim=c(nr, nc), vmode=vmode(object))
  
  ## hypothesis: nr >> nc

  if ((!missing(rows)) && missing(cols)){
    ## cols iObj match cols oOubj
    samplesByNode <- splitIndicesByNode(1:ncol(out))
    ocLapply(samplesByNode, getSomeRowsAllCols, rows, object, out, neededPkgs="oligo")
  }else if (missing(rows) && (!(missing(cols)))){
    samplesByNode <- splitIndicesByNode(cols)
    ocLapply(samplesByNode, getAllRowsSomeCols, cols, object, out, neededPkgs="oligo")
  }else if ((!missing(rows)) && (!missing(cols))){
    samplesByNode <- splitIndicesByNode(cols)
    ocLapply(samplesByNode, getSomeRowsSomeCols, cols, rows, object, out, neededPkgs="oligo")
  }else{
    stop("Must specify at least one of 'rows'/'cols'")
  }

  dimnames(object) <- dnmsIn
  dnmsOut <- list(rns, cns)
  rm(rns, cns)
  dimnames(out) <- dnmsOut
  
  return(out)
}
