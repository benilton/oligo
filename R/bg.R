bg.correct.rma <- function(object,bgtype=2){
  
  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}
  pm(object) <- .Call("bg_correct_c_copy",pm(object),pm(object),body(bg.dens),new.env(),bgtype,PACKAGE="oligo")
  return(object)
}







