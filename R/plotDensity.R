## BC: Mon, Jul 25, 2005 - Plot density - used by hist()
plotDensity <- function(mat,
                        ylab="density", xlab="x", type="l", col=1:6,
                        ...) {
  
  x.density <- apply(mat, 2, density)
  all.x <- do.call("cbind", lapply(x.density, function(x) x$x))
  all.y <- do.call("cbind", lapply(x.density, function(x) x$y))
  matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, ...)
  invisible(list(all.x=all.x, all.y=all.y))
}
 

plotDensity.FeatureSet <- function(x, col=1:6, log=TRUE,
                                  which=c("both","pm","mm"),
                                  ylab="density",
                                  xlab=NULL,
                                  ...){

  which <- match.arg(which,c("both","pm","mm"))
  Index <- unlist(featureIndex(x,which))
  x <- exprs(x)[Index, ,drop=FALSE]
  if(log){
    x <- log2(x)
    if(is.null(xlab)) xlab <- "log intensity"
  }
  else  if(is.null(xlab)) xlab <- "intensity"
  rv <- plotDensity(x, ylab=ylab, xlab=xlab, col=col, ...)
  invisible(rv)
}
