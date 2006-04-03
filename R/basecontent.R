basecontent <- function(seq) {
  good   = !is.na(seq)
  havena = !all(good)
  if(havena)
    seq = seq[good]
  
  rv = .Call("basecontent", seq, PACKAGE="oligo")

  if(havena) {
    z = rv
    rv = matrix(NA, nrow=length(good), ncol=ncol(z))
    colnames(rv) = colnames(z)
    rv[good, ] = z
  }
  
  return(rv)
}
