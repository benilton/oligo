## pmSequence <- function(x){
##   if (class(get(annotation(x))) == "AffySNPPDInfo"){
##     seqs <- "select seq from sequence, pmfeature where sequence.fid=pmfeature.fid"
##     return(dbGetQuery(db(get(annotation(x))), seqs)[[1]])
##   }else{
##     return(getPlatformDesign(x)$sequence[pmindex(x)])
##   }
## }

mmSequence <- function(x)  getPlatformDesign(x)$sequence[mmindex(x)]

sequenceDesignMatrix <- function(seqs){
  if(length(unique(sapply(seqs,nchar)))!=1) stop("Sequences must be of same length.")
  oligolength <- nchar(seqs[1])
  mat <- .Call("gcrma_getSeq2",paste(seqs,collapse=""),length(seqs),oligolength,PACKAGE="oligo")
  colnames(mat) <- paste(rep(c("A","C","G"),rep(oligolength,3)),position=rep(1:oligolength,3),sep="_")
  return(mat)
}


         
