rma2 <- function(object, background=TRUE, normalize=TRUE){
  pnVec <- probeNames(object)
  pms <- pm(object)
  set.buffer.dim(pms, 20000, 1)
  RowMode(pms)
  if (background) pms <- bg.correct.BufferedMatrix(pms, copy=FALSE)
  if (normalize) pms <- normalize.BufferedMatrix.quantiles(pms, copy=FALSE)
  exprs <- median.polish.summarize.BufferedMatrix(pms, length(unique(pnVec)), pnVec)
  rownames(exprs) <- unique(pnVec)
  colnames(exprs) <- sampleNames(object)
  out <- new("ExpressionSet",
             exprs=exprs,
             phenoData=phenoData(object),
             experimentData=experimentData(object),
             annotation=annotation(object))
  return(out)
}
