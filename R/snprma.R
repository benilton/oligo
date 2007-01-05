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
  RowMode(pms)
  ssSize <- 2000
  correctionMatrix <- sequenceDesignMatrix(pmSequence(x))
  sql <- "select offset from sequence, pmfeature where sequence.fid = pmfeature.fid"
  snpLocation <- dbGetQuery(db(get(annotation(x))), sql)[[1]]

  ## check snpLocation... if too few probes
  ## at a given location, then change their locations
  ## to something else
  theUniqueLocs <- sort(unique(snpLocation))
  theCounts <- table(snpLocation)
  bad <- theUniqueLocs[which(theCounts < 1000)]
  if (length(bad)>0){
    message("Bad snp location(s): ", bad)
    good <- theUniqueLocs[-which(theCounts<1000)]
    for (i in bad)
      snpLocation[snpLocation == i] <- good[which.min(abs(good - i))]
    rm(good, i)
  }
  rm(theUniqueLocs, theCounts, bad)
  
  correctionMatrix <- cbind(1, correctionMatrix)

  ## getting length
  sql <- "select featureSet.fsetid, fragment_length from featureSet, pmfeature where featureSet.fsetid=pmfeature.fsetid"
  fragLength <- dbGetQuery(db(get(annotation(x))), sql)
  pmfsetid <- dbGetQuery(db(get(annotation(x))), "select fsetid from pmfeature")[[1]]
  idx <- match(pmfsetid, fragLength[[1]])
  fragLength <- fragLength[idx, 2]
  rm(idx, pmfsetid)
  ## done getting length

  fragLength[is.na(fragLength)] <- median(fragLength, na.rm=T)
  
  correctionMatrix <- cbind(correctionMatrix, ns(fragLength, df=3))
  rm(fragLength); gc()
  
  ewApply(pms, log2)
  for (loc in sort(unique(snpLocation))){
    cat(paste("\nPosition ", loc, ".\n", sep="")) 
    set <- snpLocation == loc
    xx <- correctionMatrix[set,]
    set.seed(1)
    idx <- sample(1:sum(set), min(ssSize, sum(set)))
    xs <- xx[idx,]
    project <- solve(t(xs)%*%xs)%*%t(xs)
    rm(xs)
    for (i in 1:ncol(pms)){
      cat(".")
      coefs <- project%*%pms[set, i][idx]
      pms[set, i] <- pms[set, i] - xx%*%coefs + mean(pms[set, i])
    }
  }
  ewApply(pms, exp2)
  return(pms)
}

exp2 <- function(x) 2^x

preProcess <- function(oBatch, hapmapNormalized=NULL){
  pms <- correctionsLite(oBatch)
  if (!is.null(hapmapNormalized)){
    return(normalizeToSample(pms, hapmapNormalized))
  }else{
    return(normalize.BufferedMatrix.quantiles(pms))
  }
}  

snprma <- function(oBatch, normalizeToHapmap=TRUE){
  pms <- correctionsLite(oBatch)
  if (normalizeToHapmap){
    data(list=paste(platform(oBatch), "Ref", sep=""))
    pms <- normalizeToSample(pms, hapmapNormalized)
  }else{
    pms <- normalize.BufferedMatrix.quantiles(pms)
  }
  pm(oBatch) <- pms
  rm(pms); gc()
  return(rma(oBatch))
}
