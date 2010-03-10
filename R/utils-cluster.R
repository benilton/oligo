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
  if (!parStatus())
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

splitIndicesByLength <- function(x, lg){
  lx <- length(x)
  split(x, rep(1:lx, each=lg, length.out=lx))
}

splitIndicesByNode <- function(x){
  if (parStatus()){
    clusterSplit(getCluster(), x)
  }else{
    list(x)
  }
}

oLapply <- function(X, FUN, ...){
  if (parStatus()){
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
  if (parStatus()){
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
  if (parStatus()){
    clusterCall(getCluster(), f, name, garbageCollect)
  }else{
    f(name, garbageCollect)
  }
}
