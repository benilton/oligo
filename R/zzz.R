.onAttach <- function(libname, pkgname) {
  version <- packageDescription("oligo", field="Version")
  message(getBar())
  message("Welcome to oligo version ", version)
  message(getBar())
}

.onUnload <- function( libpath ){
  library.dynam.unload("oligo", libpath)
}

.oligoPkgEnv <- new.env(parent=emptyenv())
