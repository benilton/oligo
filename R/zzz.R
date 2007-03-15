# Loading required libraries

.onLoad <- function(libname, pkgname) {
  require("methods")
}

.onAttach <- function(libname, pkgname) {
  message(paste("\nWelcome to the oligo Package!\n",
                "This package is under development,",
                "and therefore its documentation is",
                "to be improved. Please contact Benilton",
                "at bcarvalh <AT> jhsph <DOT> edu,",
                "in case you need further help.", sep="\n    "))
}

.onUnload <- function( libpath ){
  library.dynam.unload("oligo", libpath)
}
