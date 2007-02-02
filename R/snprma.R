normalizeToSample <- function(toNormalize, Normalized){
  ncols <- ncol(toNormalize)
  Normalized <- sort(Normalized)
  for (i in 1:ncols){
    idx <- order(toNormalize[, i])
    toNormalize[idx, i] <- Normalized
  }
  return(toNormalize)
}

correctionsLite <- function(x){
  pms <- pm(x)
  set.buffer.dim(pms, 300000, 1)
  ColMode(pms)
  ssSize <- 2000
  correctionMatrix <- sequenceDesignMatrix(pmSequence(x))
  snpLocation <- pmPosition(get(annotation(x)))

  ## check snpLocation... if too few probes
  ## at a given location, then change their locations
  ## to something else
  theUniqueLocs <- sort(unique(snpLocation))
  theCounts <- table(snpLocation)
  bad <- theUniqueLocs[which(theCounts < 1000)]
  if (length(bad)>0){
    good <- theUniqueLocs[-which(theCounts<1000)]
    for (i in bad)
      snpLocation[snpLocation == i] <- good[which.min(abs(good - i))]
    rm(good, i)
  }
  rm(theUniqueLocs, theCounts, bad)
  
  correctionMatrix <- cbind(1, correctionMatrix)

  fragLength <- pmFragmentLength(get(annotation(x)))
  fragLength[is.na(fragLength)] <- median(fragLength, na.rm=T)
  
  correctionMatrix <- cbind(correctionMatrix, ns(fragLength, df=3))
  rm(fragLength); gc()
  
  ewApply(pms, log2)
  RowMode(pms)
  for (loc in sort(unique(snpLocation))){
    cat("Position ", loc)
    set <- snpLocation == loc
    xx <- correctionMatrix[set,]
    set.seed(1)
    idx <- sample(1:sum(set), min(ssSize, sum(set)))
    xs <- xx[idx,]
    tmp <- solve(t(xs)%*%xs, t(xs)%*%pms[set,][idx,])
    pms[set,] <- pms[set,]-sweep(xx%*%tmp, 2, colMeans(pms[set,]), "+")
    rm(xs)
    cat("\n")
  }
  ewApply(pms, exp2)
  return(pms)
}

exp2 <- function(x) 2^x

snprma <- function(oBatch, normalizeToHapmap=TRUE, saveQuant=FALSE){
  pms <- correctionsLite(oBatch)
  annot <- annotation(oBatch)
  set.buffer.dim(pms, 300000, 1)
  if (normalizeToHapmap){
    require(paste(annot, ".crlmm.regions", sep=""), character.only=TRUE)
    data(list=paste(platform(oBatch), "Ref", sep=""))
    pms <- normalizeToSample(pms, reference)
  }else{
    normalize.BufferedMatrix.quantiles(pms, copy=FALSE)
    if (saveQuant){
      reference <- sort(pms[,1])
      save(reference, file="quantileReference.rda")
    }
  }

  ## get rma pars:
  ## put PMs in right order
  ## get pnVec
  ## get length(unique(pnVec))
  pnVec <- paste(probeNames(oBatch),
                 c("A", "B")[pmAllele(get(annot))+1],
                 c("S", "A")[pmStrand(get(annot))+1],
                 sep="")
  idx <- order(pnVec)
  pms <- subBufferedMatrix(pms, idx)
  pnVec <- pnVec[idx]
  rm(idx); gc()

  ## params OK
  RowMode(pms)
  set.buffer.dim(pms, 300000, 1)
  theExprs <- median.polish.summarize.BufferedMatrix(pms, length(unique(pnVec)), pnVec)
  colnames(theExprs) <- sampleNames(oBatch)
  rm(pms, pnVec); gc()
  theExprs <- sqsFrom(theExprs)
  annotation(theExprs) <- annot
  phenoData(theExprs) <- phenoData(oBatch)
  experimentData(theExprs) <- experimentData(oBatch)
  sampleNames(theExprs) <- sampleNames(oBatch)
  return(theExprs)
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
  rownames(tmp) <- rep(snps, each=4)
  aTa <- seq(1, nrow(tmp), by=4)
  tmp <- new("SnpQSet",
             antisenseThetaA=tmp[aTa,, drop=FALSE],
             senseThetaA=tmp[(aTa+1),, drop=FALSE],
             antisenseThetaB=tmp[(aTa+2),, drop=FALSE],
             senseThetaB=tmp[(aTa+3),, drop=FALSE])
  return(tmp)
}
