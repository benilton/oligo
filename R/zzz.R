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

.onLoad <- function(libname, pkgname) {
  ## load the Lapack library needed for some parts of fitPLM
  .C("Lapack_Init", PACKAGE="oligo")
}
