###########################################################
##
## file: internalfunctions.R
##
## Copyright (C) 2004   Ben Bolstad
##
## created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
## created on: Jan 18, 2004
##
## Purpose: internal functions for mapping between
##          R-code and C-code.
##
## History
## Jan 18, 2004 - Initial version. Move functions out of
##                fitPLM, threestepPLM, threestep, rmaPLM
## Apr 27, 2006 - fix typo in get.background.code
## Jul 19, 2008 - fix typo in get.psi.code
##
###########################################################
 
get.background.code <- function(name) {
  background.names <- c("RMA.2", "IdealMM", "MAS", "MASIM", "LESN2", "LESN1", "LESN0", "GCRMA")
  if (!is.element(name, background.names)) {
    stop(paste(name, "is not a valid background correction method. Please use one of:",
               "RMA.2", "IdealMM","LESN2","LESN1","LESN0","MAS","MASIM","GCRMA"))
  }
  code <- c(1, 2, 3, 4, 5, 6, 7, 8,9)[name == background.names]
  code
}

get.normalization.code <- function(name) {
  normalization.names <- c("quantile", "quantile.probeset", "scaling")
  if (!is.element(name, normalization.names)) {
    stop(paste(name, "is not a valid summary method. Please use one of:",
               "quantile","quantile.probeset","scaling"))
  }
  code <- c(1,2,3)[name == normalization.names]
  code
}

get.psi.code <- function(name){
  psi.names <- c("Huber","fair","Cauchy","Geman-McClure","Welsch","Tukey","Andrews")
  if (!is.element(name, psi.names)) {
    stop(paste(name, "is not a valid Psi type. Please use one of:",
               "Huber","fair","Cauchy","Geman-McClure","Welsch","Tukey","Andrews"))
  }
  code <- c(0:6)[name == psi.names]
  code
}

get.default.psi.k <- function(name){
  if (!is.numeric(name)){
    psi.code <- get.psi.code(name)
  } else {
    psi.code <- name
  }
  ## ** Huber - k = 1.345
  ## ** Fair - k = 1.3998
  ## ** Cauchy - k=2.3849 
  ## ** Welsch - k = 2.9846
  ## ** Tukey Biweight - k = 4.6851
  ## ** Andrews Sine - K = 1.339
  if (psi.code == 0){
    psi.k <- 1.345
  } else if (psi.code == 1){
    psi.k <- 1.3998
  } else if (psi.code == 2){
    psi.k <- 2.3849
  } else if (psi.code == 4){
    psi.k <- 2.9846
  } else if (psi.code == 5){
    psi.k <- 4.6851
  } else if (psi.code == 6){
    psi.k <- 1.339
  } else {
    psi.k <- 1
  }
  psi.k
}




get.summary.code <- function(name){
  summary.names <- c("median.polish","tukey.biweight","average.log","rlm","log.average","log.median","median.log","log.2nd.largest","lm")
  
  if (!is.element(name,summary.names)){
    stop(paste(name,"is not a valid summary method. Please use one of:","median.polish","tukey.biweight","average.log","rlm","log.average","log.median","median.log","log.2nd.largest","lm"))
  }
  code <- c(1,2,3,4,5,6,7,8,9)[name ==summary.names]
  code
}


convert.LESN.param <- function(param.list){
    defaults <- c(0.25,4)
    if (!is.null(param.list$baseline)){
      defaults[1] <- param.list$baseline
    }
    if (!is.null(param.list$theta)){
      defaults[2] <- param.list$theta
    }
    defaults
  }
