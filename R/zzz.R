# Loading required libraries

.onLoad <- function(libname, pkgname) {
  require("methods")
}

.onAttach <- function(libname, pkgname) {
  message("This is the oligo package")
}

.onUnload <- function( libpath ){
  library.dynam.unload("oligo", libpath)
}
