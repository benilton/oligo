#################
## RMA Background
#################
rmaBgCorrectLDSnode <- function(cols, object){
  ## this runs on the node
  ## it (rma) bg corrects 'object' and overwrites it
  if (length(cols) > 0 ){
    open(object)
    grpCols <- splitIndicesByLength(cols, ocSamples())
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
  stopifnot(ldStatus())
  if (copy){
    out <- clone(object, pattern=file.path(ldPath(), "oligo-rmabg-"))
  }else{
    out <- object
  }
  dnms <- dimnames(out)
  dimnames(out) <- NULL
  samplesByNode <- splitIndicesByNode(1:ncol(out))
  ocLapply(samplesByNode, rmaBgCorrectLDSnode, out, neededPkgs="oligo")
  dimnames(out) <- dnms
  rm(samplesByNode, dnms)
  return(out)
}


#################
## LESN - matrix only - for the moment
#################
bgShift <- function(pmMat, baseline=0.25){
  objsize <- dim(pmMat)
  matrix(.C("R_shift_down", as.double(as.vector(pmMat)),
            as.double(baseline), as.integer(objsize[1]),
            as.integer(objsize[2]), PACKAGE="oligo")[[1]],
         objsize[1], objsize[2])
}

bgStretch <- function(pmMat, baseline=0.25,
                      type=c("linear","exponential","loglinear","logexponential","loggaussian"),
                      theta){
    type <- match.arg(type)
    objsize <- dim(pmMat)
    ty <- match(type, c("linear","exponential","loglinear","logexponential","loggaussian"))
    matrix(.C("R_stretch_down", as.double(as.vector(pmMat)),
              as.double(baseline), as.integer(objsize[1]),
              as.integer(objsize[2]), as.integer(ty),
              as.double(theta), PACKAGE="oligo")[[1]],
           objsize[1], objsize[2])
}

bgLESN <- function(pmMat, method=2, baseline=0.25, theta=4){
  if (method==2){
    bgStretch(pmMat, baseline, type="loggaussian", theta=2*theta^2) 
  } else if (method == 1){
    bgStretch(pmMat, baseline, type="logexponential", theta) 
  } else {
    bgShift(pmMat, baseline)
  }
}


#################
## MAS - matrix only - for the moment
#################

bgMAS <- function(intensities, xcoord, ycoord, geometry, griddim=16){
    nprobes <- nrow(intensities)
    nchips <- ncol(intensities)
    rows <- geometry[1]
    cols <- geometry[2]
    matrix(.C("affy_background_adjust_R",
              as.double(as.vector(intensities)), as.integer(xcoord), as.integer(ycoord),
              as.integer(nprobes), as.integer(nchips), as.integer(rows), as.integer(cols),
              as.integer(griddim), PACKAGE="oligo")[[1]],
           nprobes, nchips)
}

bgMASFS <- function(object, griddim=16){
    ## mod by BC
    pm.index <- pmindex(object)
    mm.index <- mmindex(object)

    ## some chips have some probesets without MM probes
    ## which will return an NA in mm.index
    mm.index <- mm.index[!is.na(mm.index)]
    i <- c(pm.index, mm.index)

    ## rows/cols of the *chip*
    rows <- geometry(object)[1]
    cols <- geometry(object)[2]

    newExprs <- exprs(object)

    ## note that the indexing is +1 more than you'd expect because
    ## the c code expects it that way
    ## (note about the remark above: R indexing starts at 1 and not at 0,
    ## that's why the indexing is done this way. The package is primarily done to
    ## be used with R...)
    allx <- c(pm.index-1, mm.index-1) %% rows +1
    ally <- c(pm.index-1, mm.index-1) %/% rows + 1

    newExprs[i, ] <- bgMAS(newExprs[i,,drop=FALSE], allx, ally,
                           geometry(object), griddim=griddim)

    out <- object
    exprs(out) <- newExprs

    ## and what with the 'non pm or mm' probes ?
    ## answer: they are not used per Affymetrix Statistical Algorithms Description Document.

    return(out)
 }

#################
## GENERAL
#################
## This is visible for user

backgroundCorrectionMethods <- function()
    c('rma', 'mas', 'LESN')

setMethod("backgroundCorrect", "matrix",
          function(object, method=backgroundCorrectionMethods(), copy=TRUE, extra, verbose=TRUE){
            method <- match.arg(method, backgroundCorrectionMethods())
            if (verbose) cat("Background correcting... ")
            if (method == "rma"){
                out <- rma.background.correct(object, copy=copy)
            } else if (method == "mas"){
                if (missing(extra))
                    stop("Argument 'extra' must be passed and must contain 'xcoord', 'ycoord', 'geometry' and optionally 'griddim'")
                nms <- names(extra)
                req <- c('xcoord', 'ycoord', 'geometry')
                if (!all(req %in% nms))
                    stop("'xcoord', 'ycoord', 'geometry' must be present in 'extra'")
                xcoord <- extra[['xcoord']]
                ycoord <- extra[['ycoord']]
                geometry <- extra[['geometry']]
                stopifnot(nrow(object) == length(xcoord),
                          length(xcoord) == length(ycoord),
                          length(geometry) == 2L,
                          is.numeric(xcoord), is.numeric(ycoord),
                          is.numeric(geometry))
                griddim <- 16
                if ('griddim' %in% nms) griddim <- extra[['griddim']]
                out <- bgMAS(intensities, xcoord, ycoord, geometry, griddim)
            } else if (method == "LESN") {
                method <- 2
                baseline <- .25
                theta <- 4
                if (!missing(extra)){
                    if ('method' %in% names(extra))
                        method <- extra[['method']]
                    if ('baseline' %in% names(extra))
                        baseline <- extra[['baseline']]
                    if ('theta' %in% names(extra))
                        theta <- extra[['theta']]
                }
                out <- bgLESN(object, method=method, baseline=baseline, theta=theta)
            }
            if (verbose) message("OK")
            return(out)
          })

setMethod("backgroundCorrect", "ff_matrix",
          function(object, method=backgroundCorrectionMethods(), copy=TRUE, extra, verbose=TRUE){
            method <- match.arg(method, backgroundCorrectionMethods())
            if (verbose) cat("Background correcting... ")
            if (method == "rma"){
                out <- rmaBgCorrectLDSmaster(object, copy=copy)
            } else if (method == "mas"){
                ## TODO
                stop("To implement: mas on ff")
            } else if (method == "LESN"){
                method <- 2
                baseline <- .25
                theta <- 4
                if (!missing(extra)){
                    if ('method' %in% names(extra))
                        method <- extra[['method']]
                    if ('baseline' %in% names(extra))
                        baseline <- extra[['baseline']]
                    if ('theta' %in% names(extra))
                        theta <- extra[['theta']]
                }
                stop("Yet to implement LESN background correction for 'ff' objects")
            }
            if (verbose) message("OK")
            return(out)
          })


setMethod("backgroundCorrect", "FeatureSet",
          function(object, method=backgroundCorrectionMethods(), copy=TRUE, extra, verbose=TRUE, ...){
              method <- match.arg(method, backgroundCorrectionMethods())
              if (copy)
                  object <- cloneFS(object)
              if (method == "rma"){
                  pm(object) <- backgroundCorrect(pm(object),
                                                  method="rma",
                                                  copy=FALSE, extra=extra,
                                                  verbose=verbose)
              }else if (method == "mas"){
                  griddim <- 16
                  if (!missing(extra))
                      if ('griddim' %in% names(extra))
                          griddim <- extra[["griddim"]]
                  out <- bgMASFS(object, griddim=griddim)
              }else if (method == "LESN"){
                  pm(object) <- backgroundCorrect(pm(object),
                                                  method="LESN",
                                                  copy=FALSE, extra=extra,
                                                  verbose=verbose)
              }
              object
          })
