.onAttach <- function(libname, pkgname) {
  version <- packageDescription("oligo", fields="Version")
  packageStartupMessage(getBar())
  packageStartupMessage("Welcome to oligo version ", version)
  packageStartupMessage(getBar())
}

.onUnload <- function( libpath ){
  library.dynam.unload("oligo", libpath)
}

.oligoPkgEnv <- new.env(parent=emptyenv())
