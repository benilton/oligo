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
  }else if ("ff_matrix" %in% theClass){
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
    tmpExprs <- normalizeToTarget(tmpExprs, target=reference, copy=FALSE, method="quantile", verbose=FALSE)
    rm(reference)
  } else {
    tmpExprs <- normalize(tmpExprs, copy=FALSE, method="quantile", verbose=FALSE)
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
  annotation(out) <- annotation(object)
  phenoData(out) <- phenoData(object)
  protocolData(out) <- protocolData(object)
  experimentData(out) <- experimentData(object)
  sampleNames(out) <- sampleNames(object)
  return(out)
}


crlmm2 <- function(object, recalibrate=TRUE, minLLRforCalls=c(5, 1, 5),
                   verbose=TRUE, balance=1.5){
  
  ## make this a method for AlleleSet
  stopifnot(is(object, "AlleleSet"))

  ## bothStrands: TRUE (before SNP 5.0) / FALSE (SNP 5.0 or SNP 6.0)
  bs <- bothStrands(object)
  
  library(annotation(object), character.only=TRUE)

  if (verbose) message("Calculating M correction... ", appendLF=FALSE)
  correction <- fitAffySnpMixture2(object, verbose=verbose)
  if (verbose) message("Done.\n")

  
  snr <- correction$snr

  if(is.null(object$gender)){
    if (verbose) message("Predicting gender... ", appendLF=FALSE)
    maleIndex <- snpGenderCall(object) == "male"
    if (verbose) message("Done.")
  }else{
    maleIndex <- object$gender=="male"
  }
  annotname <- annotation(object)
  load(system.file(paste("extdata/", annotname, "CrlmmInfo.rda",
                         sep=""), package=annotname))
  myenv <- get(paste(annotname,"Crlmm",sep=""))
  rm(list=paste(annotname,"Crlmm",sep=""))

  thePriors <- get("priors", myenv)
  Index <- which(!get("hapmapCallIndex", myenv))

  myCalls <- matrix(NA, dim(object)[1], dim(object)[2])

  myCalls[Index,] <- getInitialAffySnpCalls(correction, Index,
                                            verbose=verbose,
                                            sqsClass=class(object))
  
  fs <- correction$fs;
  rparams <- getAffySnpGenotypeRegionParams(object, myCalls, fs,
                                            subset=Index,
                                            verbose=verbose,
                                            sqsClass=class(object))
  rm(myCalls)
  
  if (bs){
    oneStrand <- apply(is.na(getM(object[,1])[,1,]), 1,
                       function(v) ifelse(length(ll <- which(v))==0, 0, ll))
    rparams <- updateAffySnpParams(rparams, thePriors, oneStrand, verbose=verbose)
    params  <- replaceAffySnpParams(get("params",myenv), rparams, Index)
  }else{
    rparams <- updateAffySnpParamsSingle(rparams, thePriors, verbose=verbose)
    params  <- replaceAffySnpParamsSingle(get("params",myenv), rparams, Index)
  }

  ## Recalibration comes after

  rm(myenv, Index)

  if (bs){
    myDist <- getAffySnpDistance(object, params, fs)
    save(myDist, file=paste(prefix, "distFile.rda", sep=""))
    myDist[,,-2,] <- balance*myDist[,,-2,]
  }else{
    myDist <- getAffySnpDistanceSingle(object, params, fs)
    save(myDist, file=paste(prefix, "distFile.rda", sep=""))
    myDist[,,-2] <- balance*myDist[,,-2]
  }
## aqui
  rm(fs)
  
  XIndex <- getChrXIndex(object)
  myCalls <- getAffySnpCalls(myDist,XIndex,maleIndex,verbose=verbose, sqsClass=class(object))
  LLR <- getAffySnpConfidence(myDist,myCalls,XIndex,maleIndex,verbose=verbose, sqsClass=class(object))
  rm(myDist)
  fs <- correction$fs
  if(recalibrate){
    if(verbose) cat("Recalibrating.")
    for(k in 1:3)
      myCalls[myCalls == k & LLR < minLLRforCalls[k]] <- NA
    rm(LLR)
    myCalls[, snr < 3.675] <- NA
    
    rparams <- getAffySnpGenotypeRegionParams(object, myCalls,
                                              fs, verbose=verbose, sqsClass=class(object))
    rm(myCalls)
    ## gc()

    if (bs){
      rparams <- updateAffySnpParams(rparams, thePriors, oneStrand, verbose=verbose)
      myDist <- getAffySnpDistance(object, rparams, fs, verbose=verbose)
      myDist[,,-2,] <- balance*myDist[,,-2,]
    }else{
      rparams <- updateAffySnpParamsSingle(rparams, thePriors, verbose=verbose)
      myDist <- getAffySnpDistanceSingle(object, rparams, fs, verbose=verbose)
      myDist[,,-2] <- balance*myDist[,,-2]
    }

    save(myDist, file=paste(prefix, "distFile.rda", sep=""))

    myCalls <- getAffySnpCalls(myDist,XIndex, maleIndex, verbose=verbose, sqsClass = class(object))
    LLR <- getAffySnpConfidence(myDist,myCalls,XIndex,maleIndex,verbose=verbose, sqsClass = class(object))
    rm(fs, myDist)
    pacc <- LLR2conf(myCalls, LLR, snr, annotation(object))
  }

  ## correction$snr
  ## maleIndex
  gender <- rep("female", length(maleIndex))
  gender[maleIndex] <- "male"
  if (is.null(object$gender)){
    addPhenoData <- new("AnnotatedDataFrame",
                        data=cbind(pData(object),
                          data.frame(crlmmSNR=as.numeric(snr),
                                     gender=gender,
                                     row.names=sampleNames(object))),
                        varMetadata= rbind(varMetadata(object),
                          data.frame(labelDescription=c("crlmmSNR", "gender"),
                                     row.names=c("crlmmSNR", "gender"))))
  }else{
    addPhenoData <- new("AnnotatedDataFrame",
                        data=cbind(pData(object),
                          data.frame(crlmmSNR=as.numeric(snr),
                                     row.names=sampleNames(object))),
                        varMetadata= rbind(varMetadata(object),
                          data.frame(labelDescription=c("crlmmSNR"),
                                     row.names=c("crlmmSNR"))))
  }

  return(new("SnpSet",
             phenoData=addPhenoData,
             experimentData=experimentData(object),
             annotation=annotation(object),
             call=myCalls,
             callProbability=pacc,
             LLR=LLR))
}

##gender in pData keeps male female
fitAffySnpMixture2 <- function(object, df1=3, df2=5, probs=rep(1/3,3),
                               eps=50, subSampleSize=10^4, seed=1,
                               verbose=TRUE){
  if(is.null(object$gender)){
    maleIndex <- snpGenderCall(object) == "male"
  }else{
    maleIndex <- object$gender == "male"
  }
  
  XIndex <- getChrXIndex(object)
  
  I <- dim(object)[1]
  J <- dim(object)[2]
  set.seed(seed)
##  tmp <- c( (1:I)[-XIndex],((I+1):(2*I))[-XIndex])
##  tmp <- c( (1:I),((I+1):(2*I)))
  if (I > subSampleSize){
    idx <- sort(sample(I, subSampleSize))
  }else{
    idx <- 1:I
  }
##  rm(tmp)

  bs <- bothStrands(object)
  if (!ldStatus()){
    if (bs){
      pis <- array(0,dim=c(I,J,3,2))
      fs <- array(0,dim=c(I,J,2))
    }else{
      pis <- array(0,dim=c(I,J,3))
      fs <- array(0,dim=c(I,J))
    }
    snr <- array(0,dim=J)
  } else {
    if (bs){
      pis <- createFF("oligo-pis-", dim=c(I, J, 3, 2))
      fs <- createFF("oligo-fs-", dim=c(I, J, 2))
    }else{
      pis <- createFF("oligo-pis-", dim=c(I, J, 3))
      fs <- createFF("oligo-fs-", dim=c(I, J))
    }
    snr <- ff(vmode="double", length=J, pattern = file.path(ldPath(), "oligo-snr-"))
    eapply(assayData(object), open)
  }

#####  dimnames(fs)<-list(featureNames(object),
#####                     sampleNames(object),
#####                     c("antisense","sense"))
#####  dimnames(pis)<-list(featureNames(object),
#####                      sampleNames(object),
#####                      c("AA","AB","BB"),
#####                      c("antisense","sense"))
#####  names(snr) <- sampleNames(object)
  
  if(verbose) cat("Fitting mixture model to ", J, " arrays. Epsilon must reach ", eps, ".\n",sep="")
  L <- getSnpFragmentLength(object)
  fix <- which(is.na(L))
  L[fix] <- median(L, na.rm=T)
  rm(fix)
  if (bs){
    L <- cbind(L,L)
    l <- as.numeric(L[idx, ])
    L <- as.numeric(L)
  }else{
    l <- L[idx]
  }
  L <- L-mean(L)
  l <- l-mean(l)
  matL <- ns(l, df1)
  for(j in 1:J){
    Y <- getM(object[,j])
    A <- getA(object[,j])
    if (bs){
      Y <- Y[,1,]
      A <- A[,1,]
    }
    fix <- which(is.na(Y))
    Y[fix] <- median(Y, na.rm=T)
    A[fix] <- median(A, na.rm=T)
    rm(fix)

    mus <- quantile(Y, c(1,3,5)/6);mus[2]=0
    sigmas <- rep(mad(c(Y[Y<mus[1]]-mus[1], Y[Y>mus[3]]-mus[3])),3)
    sigmas[2] <- sigmas[2]/2

    if (bs){
      a <- as.numeric(A[idx, ])
      y <- as.numeric(Y[idx, ])
    }else{
      a <- A[idx]
      y <- Y[idx]
    }
    A <- A-mean(A)
    a <- a-mean(a)
    
    weights <- apply(cbind(mus, sigmas), 1, function(p) dnorm(y, p[1], p[2]))
    PreviousLogLik <- -Inf
    change <- eps+1
    itmax <- 0
    
    matA <- ns(a, df2)
    while (change > eps && itmax < 100){
      itmax <- itmax+1
      
      ## E
      z <- sweep(weights, 2, probs, "*")
      LogLik <- rowSums(z)
      z <- sweep(z, 1, LogLik, "/")
      LogLik <- sum(log(LogLik))
      change <- abs(LogLik-PreviousLogLik)
      
      if (verbose){
        if (itmax > 1 || j > 1) cat(del)
        message <- paste("Array ",j,": epsilon = ", signif(change,2), "  ", sep="")
        del <- paste(rep("\b", nchar(message)), collapse="")
        cat(message)
      }
      
      PreviousLogLik <- LogLik
      probs <- colMeans(z)

      ## M
      fit1 <- lm(y~matL+matA,weights=z[,1])
      fit2 <- sum(z[,2]*y)/sum(z[,2])
      fit3 <- lm(y~matL+matA,weights=z[,3])
      
      sigmas[1] <- sqrt(sum(z[,1]*residuals(fit1)^2)/sum(z[,1]))
      sigmas[2] <- sqrt(sum(z[,2]*(y-fit2)^2)/sum(z[,2]))
      sigmas[3] <- sqrt(sum(z[,3]*residuals(fit3)^2)/sum(z[,3]))
      
      weights[,1] <- dnorm(y,fitted(fit1),sigmas[1])
      weights[,2] <- dnorm(y,fit2,sigmas[2])
      weights[,3] <- dnorm(y,fitted(fit3),sigmas[3])
      weights[y >= 0, 1] <- 0
      weights[y <= 0, 3] <- 0
    }

    ## gc()
    bigX <- cbind(1, ns(L, knots=as.numeric(attr(matL, "knots")), Boundary.knots=attr(matL, "Boundary.knots")),
                  ns(A, knots=as.numeric(attr(matA, "knots")), Boundary.knots=attr(matA, "Boundary.knots")))
    rm(matA); ## gc()

    pred1 <- bigX%*%coef(fit1)
    pred2 <- rep(fit2,length(Y))
    pred3 <- bigX%*%coef(fit3)
    rm(bigX); ## gc()

    weights <- matrix(0,length(Y),3)
    weights[,1] <- dnorm(Y,pred1,sigmas[1])
    weights[,2] <- dnorm(Y,pred2,sigmas[2])
    weights[,3] <- dnorm(Y,pred3,sigmas[3])
    weights[Y >= pred2, 1] <- 0
##    if (maleIndex[j]) weights[XIndex, 2] <- 0
    weights[Y <= pred2, 3] <- 0

    z <- sweep(weights, 2, probs, "*")
    LogLik <- rowSums(z)
    z <- sweep(z, 1, LogLik, "/")

    if (bs){
      fs[,j,] <- matrix((pred3-pred1)/2,ncol=2)
      for(k in 1:3){
        pis[,j,k,] <- matrix(z[,(4-k)],ncol=2) ##4-k cause 3is1,2is2 and 1is3
      }
      snr[j] <- median(fs[,j,])^2/(sigmas[1]^2+sigmas[2]^2)
    }else{
      fs[, j] <- (pred3-pred1)/2
      for(k in 1:3){
        pis[, j, k] <- z[,(4-k)] ##4-k cause 3is1,2is2 and 1is3
      }
      snr[j] <- median(fs[,j])^2/(sigmas[1]^2+sigmas[2]^2)
    }
  }
  if (!ldStatus()){
    fs[fs < 0] <- 0
    f0 <- median(fs)
  }else{
    ffvecapply(fs[i1:i2][fs[i1:i2] < 0] <- 0, X=fs, BATCHBYTES=ocProbesets()*J*8)
    f0 <- rep(NA, J)
    if (bs){
      for (j in 1:J) f0[j] <- median(fs[, j,])
    }else{
      for (j in 1:J)  f0[j] <- median(fs[,j])
    }
    f0 <- median(f0)
    eapply(assayData(object), close)
    close(pis)
    close(fs)
    close(snr)
  }
  if(verbose) cat("Done.\n")
  return(list(f0=f0, fs=fs, pis=pis, snr=snr))
}

getInitialAffySnpCalls2 <- function(object,subset=NULL,
                                    concordanceCutoff=0.0001,
                                    cutoffs=c(0.7,0.5,0.7),
                                    returnProbs=FALSE,
                                    verbose=FALSE,
                                    sqsClass = "SnpQSet"){
  if(is.null(subset)) subset <- 1:(dim(object$pis)[1])
  if (sqsClass == "SnpQSet"){
    pi1 <- object$pis[subset,,,1]
    pi2 <- object$pis[subset,,,2]; rm(object); ## gc()

    ## fixing SNPs that have probes only on one strand:
    idx <- is.na(pi1[,1,1])
    pi1[idx,,] <- pi2[idx,,]
    idx <- is.na(pi2[,1,1])
    pi2[idx,,] <- pi1[idx,,]
    rm(idx); ## gc()
  
    if(verbose) cat("Picking good starting value: ")
    if(verbose) cat("Computing entropy, ")
    E1 <- rowEntropy(pi1)
    E2 <- rowEntropy(pi2); ## gc()
  
    tmpN <- dim(pi1)
    tmpcall1 <- tmpcall2 <- matrix(NA, nrow=tmpN[1], ncol=tmpN[2])
#####  rownames(tmpcall1) <- rownames(tmpcall2) <- dimnames(pi1)[[1]]
#####  colnames(tmpcall1) <- colnames(tmpcall2) <- dimnames(pi1)[[2]]
  ## gc()
    if(verbose) cat("calculating calls, ")
    for (i in 1:tmpN[1]){
      for (j in 1:tmpN[2]){
        tmpcall1[i, j] <- which.max(pi1[i,j,])
        tmpcall2[i, j] <- which.max(pi2[i,j,])
      }
    }
  
    if(verbose) cat("determining non-concordant calls, ")
    concordance <- rowIndepChiSqTest(tmpcall1,tmpcall2); ## gc()

    if(verbose) cat("deciding which strand(s) to use")
    tc1 <- tmpcall1 == 1
    tc2 <- tmpcall1 == 2
    tc3 <- tmpcall1 == 3
    noABIndex1 <- which((rowSums(tc2)<3)*(rowSums(tc3)>0)*(rowSums(tc1)>0) == 1)
    tc1 <- tmpcall2 == 1
    tc2 <- tmpcall2 == 2
    tc3 <- tmpcall2 == 3
    noABIndex2 <- which((rowSums(tc2)<3)*(rowSums(tc3)>0)*(rowSums(tc1)>0) == 1)
    rm(tc1, tc2, tc3); ## gc()
  
    E1[noABIndex1] <- -Inf
    E2[noABIndex2] <- -Inf

    jointprobs <- (pi1+pi2)/2

    ## gc()
    noInfoIndex <- intersect(noABIndex1, noABIndex2)
  
    rm(noABIndex1, noABIndex2)
  
    jointprobs[noInfoIndex,,] <- 1/3
    rm(noInfoIndex)
  
    notBoth <- concordance < concordanceCutoff
    E1bE2 <- E1 > E2
    rm(E1, E2)
    i1 <- notBoth & E1bE2
    i2 <- notBoth & !E1bE2
    rm(notBoth, E1bE2)
    jointprobs[i1,,] <- pi1[i1,,]
    rm(i1, pi1)
    jointprobs[i2,,] <- pi2[i2,,]
    rm(i2, pi2); ## gc()
    
    tmpN <- dim(jointprobs)
    tmpcall <- tmpmax <- matrix(NA, nrow=tmpN[1], ncol=tmpN[2])
#####   rownames(tmpcall) <- rownames(tmpmax) <- dimnames(jointprobs)[[1]]
#####   colnames(tmpcall) <- colnames(tmpmax) <- dimnames(jointprobs)[[2]]
  ## gc()

    if(verbose) cat(" finalizing"); ## gc()
    
    for (i in 1:tmpN[1]){
      for (j in 1:tmpN[2]){
        tmpcall[i, j] <- which.max(jointprobs[i, j, ])
        tmpmax[i, j] <- jointprobs[i, j, tmpcall[i,j]]
      }
    }
  }else{ ## if SnpQSet
    tmpcall <- apply(object$pis[subset,, ], 1:2, which.max)
    tmpmax  <- apply(object$pis[subset,, ], 1:2, max)
  }
  
    for(i in 1:3)
      tmpcall[tmpcall==i & tmpmax<cutoffs[i]] <- NA
  if(verbose) cat("\nCompleted!\n")
  if(sqsClass == "SnpQSet"){
    if(returnProbs) return(list(calls=tmpcall,probs=jointprobs)) else return(tmpcall)
  }else{
    return(tmpcall)
  }
}
