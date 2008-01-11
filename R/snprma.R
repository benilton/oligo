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
  rownames(tmp) <- rep(snps, each=4)
  aTa <- seq(1, nrow(tmp), by=4)
  tmp <- new("SnpQSet",
             antisenseThetaA=tmp[aTa,, drop=FALSE],
             senseThetaA=tmp[(aTa+1),, drop=FALSE],
             antisenseThetaB=tmp[(aTa+2),, drop=FALSE],
             senseThetaB=tmp[(aTa+3),, drop=FALSE])
  return(tmp)
}

sqsFrom.SnpCnv <- function(pmMat){
  snps <- rownames(pmMat)
  snps <- unique(substr(snps, 1, (nchar(snps)-1)))
  samples <- colnames(pmMat)
  pns <- paste(rep(snps, each=2),
               rep(c("A", "B"), length(snps)),
               sep="")
  tmp <- matrix(NA, ncol=ncol(pmMat), nrow=length(pns))
  rownames(tmp) <- pns
  idx <- match(rownames(pmMat), pns)
  tmp[idx,] <- pmMat
  rownames(tmp) <- rep(snps, each=2)
  aTa <- seq(1, nrow(tmp), by=2)
  tmp <- new("SnpCnvQSet",
             thetaA=tmp[aTa,, drop=FALSE],
             thetaB=tmp[(aTa+1),, drop=FALSE])
  return(tmp)
}

snprma <- function(object, verbose=TRUE, normalizeToHapmap=TRUE){
  tmpExprs <- pm(object)
  pkgname <- annotation(object)
  
  ########################
  ##### NORMALIZATION ####
  ########################

  if (normalizeToHapmap){
    if (verbose) message("Normalizing to Hapmap.")
    load(system.file("extdata", paste(pkgname, "Ref.rda", sep=""), package=pkgname))
    reference <- sort(reference)
    tmpExprs <- normalize.quantiles.use.target(tmpExprs, reference)
  } else {
    tmpExprs <- normalize.quantiles(tmpExprs)
    reference <- sort(tmpExprs[,1])
    save(reference, file=paste(pkgname, ".quantileReference.rda", sep=""))
  }
  rm(reference); gc(); gc()

  ########################
  #### SUMMARIZATION  ####
  ########################

  if (verbose) message("Summarizing.")
  
  ## get rma pars:
  ## put PMs in right order
  ## get pnVec
  ## get length(unique(pnVec))

  if (class(object) == "SnpFeatureSet"){
    pnVec <- paste(probeNames(get(pkgname)),
                   c("A", "B")[pmAllele(get(pkgname))+1],
                   c("S", "A")[pmStrand(get(pkgname))+1],
                   sep="")
  }else{
    pnVec <- paste(probeNames(get(pkgname)),
                   c("A", "B")[pmAllele(get(pkgname))+1],
                   sep="")
  }
  idx <- order(pnVec)
  tmpExprs <- tmpExprs[idx, ]
  pnVec <- pnVec[idx]
  rm(idx); ## gc()

  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}

  theSumm <- .Call("rma_c_complete_copy", tmpExprs, tmpExprs,
                   pnVec, length(unique(pnVec)), body(bg.dens),
                   new.env(), FALSE, FALSE,
                   as.integer(2), PACKAGE="oligo")

  
  rm(tmpExprs, pnVec); ## gc()
  
  if (class(object) == "SnpFeatureSet"){
    theSumm <- sqsFrom(theSumm)
  }else{
    theSumm <- sqsFrom.SnpCnv(theSumm)
  }
  featureData(theSumm) <- featureData(object)
  phenoData(theSumm) <- phenoData(object)
  annotation(theSumm) <- annotation(object)
  return(theSumm)
}
