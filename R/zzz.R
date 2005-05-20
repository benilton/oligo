# Loading required libraries
.First.lib <- function(libname, pkgname, where) {
  library(Biobase)
  library.dynam("oligo", pkgname, libname)

}
