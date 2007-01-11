rma <- function(object, background=TRUE, normalize=TRUE){
  pms <- pm(object)
  pnVec <- probeNames(object)
  idx <- order(pnVec)
  pms <- subBufferedMatrix(pms, idx)
  pnVec <- pnVec[idx]
  rm(idx); gc()
  if (background) bg.correct.BufferedMatrix(pms, copy=FALSE)
  if (normalize) normalize.BufferedMatrix.quantiles(pms, copy=FALSE)
  set.buffer.dim(pms, 300000, 1)
  RowMode(pms)
  exprs <- median.polish.summarize.BufferedMatrix(pms, length(unique(pnVec)), pnVec)
  rownames(exprs) <- unique(pnVec)
  colnames(exprs) <- sampleNames(object)
  rm(pms, pnVec); gc()
  out <- new("ExpressionSet",
             exprs=exprs,
             phenoData=phenoData(object),
             experimentData=experimentData(object),
             annotation=annotation(object))
  sampleNames(out) <- sampleNames(object)
  return(out)
}
