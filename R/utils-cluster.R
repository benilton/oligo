setCluster <- function(...){
  pkg <- "snow"
  require(pkg, character.only=TRUE)
  options(cluster=makeCluster(...))
}

delCluster <- function(){
  stopCluster(getOption("cluster"))
  options(cluster=NULL)
}

getCluster <- function()
  getOption("cluster")

requireClusterPkgSet <- function(packages){
  if (!oligoParallelSupport())
    stop("cluster is not ready. Use 'setCluster'.")
  for (pkg in packages){
    pkgOnCluster <- requireClusterPkg(pkg, character.only=TRUE)
    if (!pkgOnCluster){
      msg <- paste("Package '", pkg, "' not found on the cluster. ",
                   "Install it or load it manually using ",
                   "'clusterEvalQ(getCluster(), library(", pkg,
                   ", lib.loc=<APPROPRIATE PATH>))'", sep="")
      stop(msg)
    }
  }
  TRUE
}

requireClusterPkg <- function(...)
  all(unlist(clusterCall(getCluster(), require, ...)))

oligoProbesets <- function(n){
  if (missing(n)){
    return(getOption("oligoProbesets"))
  }else{
    options(oligoProbesets=n)
    invisible(TRUE)
  }
}

oligoSamples <- function(n){
  if (missing(n)){
    return(getOption("oligoSamples"))
  }else{
    options(oligoSamples=n)
    invisible(TRUE)
  }
}

splitIndicesByLength <- function(x, lg){
  lx <- length(x)
  split(x, rep(1:lx, each=lg, length.out=lx))
}

splitIndicesByNode <- function(x){
  if (oligoParallelSupport()){
    clusterSplit(getCluster(), x)
  }else{
    list(x)
  }
}

validWorkDir <- function(path){
  stopifnot(!missing(path), is.character(path))
  isDir <- file.info(path)[["isdir"]]
  if (!isDir) warning(path, " is not a valid directory.")
  pathExists <- file.access(path, mode=0) == 0
  if (!pathExists) warning(path, " does not exist.")
  canRead <- file.access(path, mode=4) == 0
  if (!canRead) warning("Cannot read from ", path)
  canWrite <- file.access(path, mode=2) == 0
  if (!canWrite) warning("Cannot write to ", path)
  canAccess <- file.access(path, mode=1) == 0
  if (!canAccess) warning("Cannot access ", path)
  isDir && pathExists && canRead && canWrite && canAccess
}

oligoBigObjectPath <- function(path){
  if (missing(path)){
    return(getOption("oligoBigObjectPath"))
  }else{
    stopifnot(is.character(path))
    options(oligoBigObjectPath=path)
  }
}

oLapply <- function(X, FUN, ...){
  if (oligoParallelSupport()){
    requireClusterPkgSet(c("oligo", "ff"))
    parLapply(getCluster(), X, FUN, ...)
  }else{
    lapply(X, FUN, ...)
  }
}

loadAndKeep <- function(obj, name){
  envir <- .oligoPkgEnv
  pkg <- "ff"
  require(pkg, character.only=TRUE)
  if (exists(name, envir=envir))
    rm(list=name, envir=envir)
  assign(name, obj, envir=envir)
  open(envir[[name]])
  TRUE
}

sendBO2PkgEnv <- function(obj, name){
  if (oligoParallelSupport()){
    clusterCall(getCluster(), loadAndKeep, obj, name)
  }else{
    loadAndKeep(obj, name)
  }
}

rmFromPkgEnv <- function(name, garbageCollect=FALSE){
  f <- function(name, garbageCollect=FALSE){
    if (exists(name, envir=.oligoPkgEnv))
      rm(list=name, envir=.oligoPkgEnv)
    if (garbageCollect)
      gc()
  }
  if (oligoParallelSupport()){
    clusterCall(getCluster(), f, name, garbageCollect)
  }else{
    f(name, garbageCollect)
  }
}
