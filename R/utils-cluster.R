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
