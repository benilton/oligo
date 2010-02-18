.onAttach <- function(libname, pkgname) {
  version <- packageDescription("oligo", field="Version")
  message(getBar())
  message("Welcome to oligo version ", version)
  setLargeDataOptions()
  bm <- getLargeDataStatus()
  snow <- getParallelStatus()

  setHook(packageEvent("ff", "attach"),
          function(...){
            setLargeDataOptions(verbose=FALSE)
            getLargeDataStatus()
          })
  
  setHook(packageEvent("snow", "attach"),
          function(...){
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
  ld <- isPackageLoaded("ff")
  if (verbose){
    message(getBar())
    message("Large dataset support for 'oligo': ", appendLF=FALSE)
    if (ld){
      message("Enabled")
      ns <- prettyNum(c(oligoProbesets(), oligoSamples()), big.mark=",")
      message("    - Probesets: ", ns[1])
      message("    - Samples..: ", ns[2])
      message("    - Path.....: ", oligoBigObjectPath())
    }else{
      message("Disabled")
      message("     - Load 'ff'")
    }      
    message(getBar())
  }
  invisible(ld)
}

getParallelStatus <- function(verbose=TRUE){
  sn <- isPackageLoaded("snow")
  cl <- oligoParallelSupport()
  ld <- isPackageLoaded("ff")
  if (!ld){  
    message("Parallel computing support for 'oligo': Disabled")
    message("     - Load 'ff'")
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

  
