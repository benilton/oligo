.onAttach <- function(libname, pkgname) {
  version <- packageDescription("oligo", field="Version")
  message(getBar())
  message("Welcome to oligo version ", version)
  getLargeDataStatus()
  getParallelStatus()
  
  setHook(packageEvent("bigmemory", "attach"),
          function(...){
            setLargeDataOptions(verbose=FALSE)
            getLargeDataStatus()
  ##          getParallelStatus()
          })
  
  setHook(packageEvent("snow", "attach"),
          function(...){
##            getLargeDataStatus()
            getParallelStatus()
          })
  
}

.onUnload <- function( libpath ){
  library.dynam.unload("oligo", libpath)
}

.oligoPkgEnv <- new.env(parent=emptyenv())

getBar <- function(width=getOption("width"))
  paste(rep("=", width), collapse="")

setLargeDataOptions <- function(nsamples=100, nprobesets=1000,
                                path=getwd(), verbose=FALSE){
  oligoProbesets(nprobesets)
  oligoSamples(nsamples)
  oligoBigObjectPath(path)
  getLargeDataStatus(verbose)
  TRUE
}

getLargeDataStatus <- function(verbose=TRUE){
  bm <- isPackageLoaded("bigmemory")
  if (bm){
    if (verbose){
      ns <- prettyNum(c(oligoProbesets(), oligoSamples()), big.mark=",")
      message(getBar())
      message("Large dataset support for 'oligo': Enabled")
      message("    - Probesets: ", ns[1])
      message("    - Samples..: ", ns[2])
      message("    - Path.....: ", oligoBigObjectPath())
      message(getBar())
    }
  }else{
    if (verbose){
      message(getBar())
      message("Large dataset support for 'oligo': Disabled")
      message("     - Load 'bigmemory'")
      message(getBar())
    }
  }
  invisible(bm)
}

getParallelStatus <- function(verbose=TRUE){
  sn <- isPackageLoaded("snow")
  cl <- oligoParallelSupport()
  bm <- isPackageLoaded("bigmemory")
  if (!bm){  
    message("Parallel computing support for 'oligo': Disabled")
    message("     - Load 'bigmemory'")
    if (!sn){
      message("     - Load 'snow'")
      message("     - Use options(cluster=makeCluster(...)")
    } else {
      if (!cl)
        message("     - Use options(cluster=makeCluster(...)")
    }
  }else{
    if (sn){
      if (cl){
        message("Parallel computing support for 'oligo': Enabled")
      }else{
        message("Parallel computing support for 'oligo': Disabled")
        message("     - Use options(cluster=makeCluster(...)")
      }
    }else{
      message("Parallel computing support for 'oligo': Disabled")
      message("     - Load 'snow'")
      message("     - Use options(cluster=makeCluster(...)")
    }
  }
  message(getBar())
  invisible(cl)
}

  
