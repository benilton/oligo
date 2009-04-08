# Loading required libraries

.onAttach <- function(libname, pkgname) {
  message("Welcome to oligo version ", packageDescription("oligo", field="Version"))
}

.onUnload <- function( libpath ){
  library.dynam.unload("oligo", libpath)
}
