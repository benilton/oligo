# Loading required libraries
.First.lib <- function(libname, pkgname, where) {
   library.dynam("oligo", pkgname, libname)
   library.dynam("affyio")
}
