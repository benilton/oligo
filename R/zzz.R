# Loading required libraries

.onLoad <- function(libname, pkgname) {
  require("methods")
}

.onAttach <- function(libname, pkgname) {
  message("oligo Package - Series 1.5.x")
}

.onUnload <- function( libpath ){
  library.dynam.unload("oligo", libpath)
}
