alleleSetFrom <- function(pmMat, bothStrands=TRUE){
  ## pmMat is expected to have rownames in the following format:
  ##  <id><Allele><Strand> eg: SNP_A-123456AS
  ##  <id><Allele>         eg: SNP_A-123456A  (SNP 5.0 and 6.0)

  if (bothStrands){
    npars <- 2
    combs <- c("AA", "AS", "BA", "BS")
    ncomb <- length(combs)
  }else{
    npars <- 1
    combs <- c("A", "B")
    ncomb <- length(combs)
  }
  
  snps <- rownames(pmMat)
  snps <- unique(substr(snps, 1, (nchar(snps)-npars)))
  pns <- paste(rep(snps, each=ncomb),
               rep(combs, length(snps)),
               sep="")

  idx <- match(pns, rownames(pmMat))
  
  theClass <- class(pmMat)
  if (theClass == "matrix"){
    tmp <- pmMat[idx,]
  }else if ("ff_matrix" %in% theClass){
    tmp <- ffSubset(rows=idx, object=pmMat, prefix="oligo-alleleSet-tmp-")
    finalizer(tmp) <- "delete"
  }else{
    stop("Class ", theClass, " not supported by alleleSetFrom.")
  }
  
  rownames(tmp) <- rep(snps, each=ncomb)
  aTa <- seq(1, nrow(tmp), by=ncomb)
  colnames(tmp) <- colnames(pmMat)

  if (bothStrands){
    if (theClass == "matrix"){
      res <- new("AlleleSet",
                 antisenseAlleleA=tmp[aTa,, drop=FALSE],
                 senseAlleleA=tmp[(aTa+1),, drop=FALSE],
                 antisenseAlleleB=tmp[(aTa+2),, drop=FALSE],
                 senseAlleleB=tmp[(aTa+3),, drop=FALSE])
    }else{
      antisenseAlleleA <- ffSubset(rows=aTa, object=tmp, prefix="oligo-alleleSet-aa-")
      senseAlleleA <- ffSubset(rows=(aTa+1), object=tmp, prefix="oligo-alleleSet-as-")
      antisenseAlleleB <- ffSubset(rows=(aTa+2), object=tmp, prefix="oligo-alleleSet-ba-")
      senseAlleleB <- ffSubset(rows=(aTa+3), object=tmp, prefix="oligo-alleleSet-bs-")
      res <- new("AlleleSet",
                 antisenseAlleleA=antisenseAlleleA,
                 senseAlleleA=senseAlleleA,
                 antisenseAlleleB=antisenseAlleleB,
                 senseAlleleB=senseAlleleB)
    }
  }else{
    if (theClass == "matrix"){
      res <- new("AlleleSet",
                 alleleA=tmp[aTa,, drop=FALSE],
                 alleleB=tmp[(aTa+1),, drop=FALSE])
    }else{
      alleleA <- ffSubset(rows=aTa, object=tmp, prefix="oligo-alleleSet-a-")
      alleleB <- ffSubset(rows=(aTa+1), object=tmp, prefix="oligo-alleleSet-b-")
      res <- new("AlleleSet",
                 alleleA=alleleA,
                 alleleB=alleleB)
    }
  }
  rm(tmp)
  return(res)
}

snprma2 <- function(object, verbose=TRUE, normalizeToHapmap=TRUE){
  conn <- db(object)
  bs <- bothStrands(object)
  pkgname <- annotation(object)

  if (bs){
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
  pmi <- pmindex(object)[idx]
  pnVec <- pnVec[idx]
  rm(idx)

  theClass <- class(exprs(object))

  if (theClass == "matrix"){
    tmpExprs <- exprs(object[pmi,])
    dimnames(tmpExprs) <- NULL
    colnames(tmpExprs) <- sampleNames(object)
  }else if (theClass == "ff_matrix"){
    tmpExprs <- ffSubset(rows=pmi, object=exprs(object), prefix="pm-")
  }else{
    stop("SNPRMA not implemented for '", theClass, "' objects.")
  }

  ########################
  ##### NORMALIZATION ####
  ########################
  if (normalizeToHapmap){
    if (verbose) message("Normalizing to Hapmap.")
    load(system.file("extdata", paste(pkgname, "Ref.rda", sep=""), package=pkgname))
    reference <- sort(reference)
    tmpExprs <- normalizeToTarget(tmpExprs, target=reference, copy=FALSE, method="quantile")
    rm(reference)
  } else {
    tmpExprs <- normalize(tmpExprs, copy=FALSE, method="quantile")
  }

  ########################
  ##### SUMMARIZATION ####
  ########################
  exprs <- summarize(tmpExprs, probes=pnVec, method="medianpolish", verbose=verbose)
  if (theClass == "ff_matrix"){
    finalizer(tmpExprs) <- "delete"
    finalizer(exprs) <- "delete"
    rm(tmpExprs)
  }

  out <- alleleSetFrom(exprs)
  rm(exprs)
  return(out)
}
