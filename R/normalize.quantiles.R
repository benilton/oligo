##################################################################
##
## file: normalize.quantiles.R
##
## For a description of quantile normalization method see
##
##  Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003)(2003)
##  A Comparison of Normalization Methods for High
##  Density Oligonucleotide Array Data Based on Bias and Variance.
##  Bioinformatics 19,2,pp 185-193
##
## History
## Pre Aug 23, 2003 Two years worth of stuff
## Aug 23, 2003 - Added use.log2 to "robust",
##                added ability to pass additional parameters
##                to normalize.AffyBatch.Quantiles.robust
##                changed pmonly parameters on functions
##                so that it is now a string argument "type"
##                the options are pmonly, mmonly, together, separate
## Jan 31, 2004 - put a check for an integer matrix and force coercision to
##                doubles if required in normalize.quantiles
##
## Aug 22, 2005 - modified it for oligoBatch
##################################################################

normalize.FeatureSet.quantiles <- function(obatch,type=c("separate","pmonly","mmonly","together")) {


  type <- match.arg(type)

  if ((type == "pmonly")|(type == "separate")){
    pms <- unlist(pmindex(obatch))

    noNA <- rowSums(is.na(exprs(obatch)[pms,,drop=FALSE])) == 0
    pms <- pms[noNA]
##    exprs(obatch)[pms,] <- normalize.quantiles(exprs(obatch)[pms,,drop=FALSE ],copy=FALSE)
    pm(obatch) <- normalize.quantiles(exprs(obatch)[pms,,drop=FALSE ],copy=FALSE)
  }
  if((type == "mmonly") | (type == "separate")){
    mms <- unlist(mmindex(obatch))
  
    noNA <- rowSums(is.na(exprs(obatch)[mms,,drop=FALSE])) == 0
    mms <- mms[noNA]

    mm(obatch) <- normalize.quantiles(exprs(obatch)[mms,,drop=FALSE ],copy=FALSE)
  }
  if (type == "together"){
    pms <- unlist(indexProbes(obatch,"both"))
    pm(obatch)  <- normalize.quantiles(exprs(obatch)[pms,,drop=FALSE ],copy=FALSE)
  }

  ##ADD CHANGE TO MIAME
  
  return(obatch)
}

normalize.quantiles <- function(x,copy=TRUE){

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  if (!is.matrix(x)){
    stop("Matrix expected in normalize.quantiles")
  }

  if (is.integer(x)){
    x <- matrix(as.double(x),rows,cols)
    copy <- FALSE
  }

  .Call("R_qnorm_c",x,copy, PACKAGE="oligo");
}

