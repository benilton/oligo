rma <- function(object, background=TRUE, normalize=TRUE){
  pms <- pm(object)
  set.buffer.dim(pms, 300000, 1)
  RowMode(pms)
  pnVec <- probeNames(object)
  if (class(get(annotation(object))) == "AffySNPPDInfo"){
    sql <- "select allele, strand from pmfeature"
    tmp <- dbGetQuery(db(get(annotation(object))), sql)
    idx1 <- tmp[,1] == 0
    idx2 <- tmp[,2] == 0
    tmp[idx1,1] <- "A"
    tmp[!idx1,1] <- "B"
    rm(idx1)
    tmp[idx2,2] <- "S"
    tmp[!idx2,2] <- "A"
    rm(idx2)
    tmp <- paste(tmp[,1], tmp[,2], sep="")
    pnVec <- paste(pnVec, tmp, sep="")
    rm(tmp)
    gc()
    background <- normalize <- FALSE
    idx <- order(pnVec)
    pms <- subBufferedMatrix(pms, idx)
    set.buffer.dim(pms, 300000, 1)
    RowMode(pms)
    pnVec <- pnVec[idx]
    rm(idx)
  }
  if (background) pms <- bg.correct.BufferedMatrix(pms, copy=FALSE)
  if (normalize) pms <- normalize.BufferedMatrix.quantiles(pms, copy=FALSE)
  exprs <- median.polish.summarize.BufferedMatrix(pms, length(unique(pnVec)), pnVec)
  rownames(exprs) <- unique(pnVec)
  colnames(exprs) <- sampleNames(object)
  rm(pms, pnVec); gc()
  if (class(object) == "SnpFeatureSet"){
    out <- sqsFrom(exprs)
    annotation(out) <- annotation(object)
    phenoData=phenoData(object)
    experimentData=experimentData(object)
  }else{
    out <- new("ExpressionSet",
               exprs=exprs,
               phenoData=phenoData(object),
               experimentData=experimentData(object),
               annotation=annotation(object))
  }
  sampleNames(out) <- sampleNames(object)
  return(out)
}


sqsFrom <- function(pmMat){
  snps <- rownames(pmMat)
  snps <- unique(substr(snps, 1, (nchar(snps)-2)))
  samples <- colnames(pmMat)
  pns <- paste(rep(snps, each=4),
               rep(c("AA", "AS", "BA", "BS"),
                   length(snps)), sep="")
  tmp <- matrix(NA, ncol=ncol(pmMat), nrow=length(pns))
  rownames(tmp) <- pns
  idx <- match(rownames(pmMat), pns)
  tmp[idx,] <- pmMat
  rm(pmMat); gc()
  rownames(tmp) <- rep(snps, each=4)
  aTa <- seq(1, nrow(tmp), by=4)
  tmp <- new("SnpQSet",
             antisenseThetaA=tmp[aTa,],
             senseThetaA=tmp[(aTa+1),],
             antisenseThetaB=tmp[(aTa+2),],
             senseThetaB=tmp[(aTa+3),])
  return(tmp)
}
