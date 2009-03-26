# Loading required libraries

## .onLoad <- function(libname, pkgname) {
##   require("methods")
## }

.onAttach <- function(libname, pkgname) {
  message("Welcome to oligo version ", packageDescription("oligo", field="Version"))
}

.onUnload <- function( libpath ){
  library.dynam.unload("oligo", libpath)
}
