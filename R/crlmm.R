rowEntropy.old <- function(p){
  p <- ifelse(p>0,p*log2(p),0)
  rowSums(p)
}

rowEntropy <- function(p) rowMeans(rowSums(log2(p^p), dims=2))

getSnpFragmentLength <- function(object){
  annotname <- annotation(object)
  load(system.file(paste("data/",annotname, ".rda", sep=""), package=paste("pd", annotname, sep="")))
  return(annot$Length[match(featureNames(object),annot$SNP)])
}

snpGenderCall <- function(object){
  XIndex=getChrXIndex(object)
  tmp=(median(getA(object), na.rm=TRUE))
  a=apply(getA(object)[XIndex,,],2,median, na.rm=TRUE)
  kfit=kmeans(a,c(min(a),tmp))
  return(factor(c("female","male")[as.numeric(kfit$cluster==1)+1]))
}

getChrXIndex <- function(object){
  annotname <- paste(annotation(object),sep="")
  load(system.file(paste("data/",annotname, ".rda", sep=""), package=paste("pd", annotname, sep="")))
  annot <- annot[match(featureNames(object),annot$SNP),]
  return(which(annot$Chromosome=="chrX"))
}

##gender in pData keeps male female
fitAffySnpMixture <- function(object, df1=3, df2=5,
                              probs=rep(1/3,3), eps=50,
                              subSampleSize=10^4,verbose=TRUE){

  if(is.null(object$gender)){
    maleIndex <- snpGenderCall(object)=="male"
  }
  else{
    maleIndex <- object$gender=="male"
  }

  XIndex=getChrXIndex(object)
  
  I <- dim(object)[1]
  J <- dim(object)[2]
  set.seed(1)
  tmp <- c( (1:I)[-XIndex],((I+1):(2*I))[-XIndex])
  idx <- sort(sample(tmp, subSampleSize))

  pis <- array(0,dim=c(I,J,3,2))
  fs <- array(0,dim=c(I,J,2))
  snr <- array(0,dim=J)

  dimnames(fs)<-list(featureNames(object),
                     sampleNames(object),
                     c("antisense","sense"))
  dimnames(pis)<-list(featureNames(object),
                      sampleNames(object),
                      c("AA","AB","BB"),
                      c("antisense","sense"))
  names(snr) <- sampleNames(object)
  
  if(verbose) cat("Fitting mixture model to ",J," arrays. Epsilon must reach ",eps,".\n",sep="")
  L <- getSnpFragmentLength(object)
  L <- c(L,L)
  for(j in 1:J){
    Y <- c(as.vector(getM(object)[,j,]))
    A <- c(as.vector(getA(object)[,j,]))

    mus <- quantile(Y,c(1,3,5)/6, na.rm=TRUE);mus[2]=0
    sigmas <- rep(mad(c(Y[Y<mus[1]]-mus[1],Y[Y>mus[3]]-mus[3]), na.rm=TRUE),3)
    sigmas[2] <- sigmas[2]/2
    
    l <- L[idx];L=L-mean(L, na.rm=TRUE);l=l-mean(l, na.rm=TRUE)
    a <- A[idx];A=A-mean(A, na.rm=TRUE);a=a-mean(a, na.rm=TRUE)
    y <- Y[idx]
    ok <- complete.cases(cbind(y,a,l))
    l <- l[ok]; a <- a[ok]; y <- y[ok]
    
    weights <- apply(cbind(mus,sigmas),1,function(p) dnorm(y,p[1],p[2]))
    PreviousLogLik <- -Inf
    change <- eps+1
    itmax <- 0
    
    matL <- ns(l,df1)
    matA <- ns(a,df2)
    while (change > eps & itmax < 1000){
      gc()
      itmax <- itmax+1

      ## E
      z <- sweep(weights, 2, probs, "*")
      LogLik <- rowSums(z)
      z <- sweep(z, 1, LogLik, "/")
      LogLik <- sum(log(LogLik))
      change <- abs(LogLik-PreviousLogLik)
      if (verbose){
        if (itmax > 1 | j > 1) cat(del)
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

    gc()
    fix <- is.na(A)
    A[fix] <- median(A, na.rm=TRUE)
    rm(fix)
    bigX <- cbind(1, ns(L, knots=as.numeric(attr(matL, "knots")), Boundary.knots=attr(matL, "Boundary.knots")),
                  ns(A, knots=as.numeric(attr(matA, "knots")), Boundary.knots=attr(matA, "Boundary.knots")))
    rm(matL, matA); gc()

    pred1 <- bigX%*%coef(fit1)
    pred2 <- rep(fit2,length(Y))
    pred3 <- bigX%*%coef(fit3)
    rm(bigX); gc()

    weights <- matrix(0,length(Y),3)
    weights[,1] <- dnorm(Y,pred1,sigmas[1])
    weights[,2] <- dnorm(Y,pred2,sigmas[2])
    weights[,3] <- dnorm(Y,pred3,sigmas[3])
    weights[Y >= pred2, 1] <- 0
    if (maleIndex[j]) weights[XIndex, 2] <- 0
    weights[Y <= pred2, 3] <- 0

    z <- sweep(weights, 2, probs, "*")
    LogLik <- rowSums(z)
    z <- sweep(z, 1, LogLik, "/")
  
    fs[,j,] <- matrix((pred3-pred1)/2,ncol=2)
    for(k in 1:3){
      pis[,j,k,] <- matrix(z[,(4-k)],ncol=2) ##4-k cause 3is1,2is2 and 1is3
    }
    snr[j] <- median(fs[,j,])^2/(sigmas[1]^2+sigmas[2]^2)
  }
  if(verbose) cat("Done.\n")
  return(list(f0=median(fs),fs=fs, pis=pis, snr=snr))
}


getInitialAffySnpCalls <- function(object,subset=NULL,
                                   concordanceCutoff=0.0001,
                                   cutoffs=c(0.7,0.5,0.7),
                                   returnProbs=FALSE,
                                   verbose=FALSE){
  if(is.null(subset)) subset <- 1:(dim(object$pis)[1])
  pi1 <- object$pis[subset,,,1]
  pi2 <- object$pis[subset,,,2]; gc()

  ## fixing SNPs that have probes only on one strand:
  idx <- is.na(pi1[,1,1])
  pi1[idx,,] <- pi2[idx,,]
  idx <- is.na(pi2[,1,1])
  pi2[idx,,] <- pi1[idx,,]
  rm(idx); gc()
  
  if(verbose) cat("Picking good starting value: ")
  if(verbose) cat("Computing entropy, ")
  E1 <- rowEntropy(pi1)
  E2 <- rowEntropy(pi2); gc()
  
  tmpN <- dim(pi1)
  tmpcall1 <- tmpcall2 <- matrix(NA, nrow=tmpN[1], ncol=tmpN[2])
  rownames(tmpcall1) <- rownames(tmpcall2) <- dimnames(pi1)[[1]]
  colnames(tmpcall1) <- colnames(tmpcall2) <- dimnames(pi1)[[2]]
  gc()
  if(verbose) cat("calculating calls, ")
  for (i in 1:tmpN[1]){
    for (j in 1:tmpN[2]){
      tmpcall1[i, j] <- which.max(pi1[i,j,])
      tmpcall2[i, j] <- which.max(pi2[i,j,])
    }
  }
  
  if(verbose) cat("determining non-concordant calls, ")
  concordance <- rowIndepChiSqTest(tmpcall1,tmpcall2); gc()

  tmpfile <- tempfile()
  save(pi1, pi2, E1, E2, tmpN, concordance, file=tmpfile)
  rm(pi1, pi2, E1, E2, tmpN, concordance); gc()
  
  if(verbose) cat("deciding which strand(s) to use")
  noABIndex1 <- (rowSums(tmpcall1==2)<3)*(rowSums(tmpcall1==3)>0)*(rowSums(tmpcall1==1)>0)
  noABIndex2 <- (rowSums(tmpcall2==2)<3)*(rowSums(tmpcall2==3)>0)*(rowSums(tmpcall2==1)>0); gc()

  load(tmpfile)
  unlink(tmpfile)
  
  E1[which(noABIndex1==1)] <- -Inf
  E2[which(noABIndex1==1)] <- -Inf

  jointprobs<-(pi1+pi2)/2

  gc()
  
  noInfoIndex <- (noABIndex1 == 1) & (noABIndex2 == 1)
  
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
  jointprobs[i2,,] <- pi2[i2,,]
  rm(i1, i2, pi1, pi2); gc()

  tmpN <- dim(jointprobs)
  tmpcall <- tmpmax <- matrix(NA, nrow=tmpN[1], ncol=tmpN[2])
  rownames(tmpcall) <- rownames(tmpmax) <- dimnames(jointprobs)[[1]]
  colnames(tmpcall) <- colnames(tmpmax) <- dimnames(jointprobs)[[2]]
  gc()

  if(verbose) cat(" finalizing"); gc()

  for (i in 1:tmpN[1]){
    for (j in 1:tmpN[2]){
      tmpcall[i, j] <- which.max(jointprobs[i, j, ])
      tmpmax[i, j] <- max(jointprobs[i, j, ])
    }
  }
  
  for(i in 1:3)
    tmpcall[tmpcall==i & tmpmax<cutoffs[i]] <- NA
  if(verbose) cat("\nCompleted!\n")
  if(returnProbs) return(list(calls=tmpcall,probs=jointprobs)) else return(tmpcall)
}


rowIndepChiSqTest <- function(call1,call2){
  tmp=vector("numeric",nrow(call1))
  for(i in 1:nrow(call1)){
    tmpt= cbind(table(factor(call1[i,],levels=1:3)),
      table(factor(call2[i,],levels=1:3)))
    tmpt=tmpt[!rowSums(tmpt)==0,,drop=FALSE]
    if(nrow(tmpt)>1){
      rowtot <- rowSums(tmpt)
      coltot <- colSums(tmpt)
      e=outer(rowtot,coltot)/sum(tmpt)
      stat = sum((tmpt-e)^2/e)
      tmp[i] = 1 - pchisq(stat, (nrow(tmpt) - 1) * (ncol(tmpt) - 1))
    }
    else{
      tmp[i]=0
    }
  }
  return(tmp)
}

getGenotypeRegionParams <- function(M,initialcalls,f=0,verbose=TRUE){
  if(!is.matrix(M)) M<-matrix(M,ncol=2)
  
  centers <- array(NA,dim=c(nrow(M),3))
  dimnames(centers) <- list(rownames(M),c("AA","AB","BB"))
  scales <- array(NA,dim=c(nrow(M),3))
  dimnames(scales) <- list(rownames(M),c("AA","AB","BB"))
  N <- array(NA,dim=c(nrow(M),3))
  dimnames(N) <- list(rownames(M),c("AA","AB","BB"))

  nrows <- nrow(M)
  
  if(verbose) cat("Computing centers and scales for 3 genotypes")
  for(k in 1:3){
    if(verbose) cat(".")
    if(k==1) tmp <- M - f
    if(k==2) tmp <- M
    if(k==3) tmp <- M + f
    tmp[initialcalls!=k] <- NA
    tmp[is.na(initialcalls)]<- NA
    for (i in 1:nrows) centers[i,k] <- median(tmp[i,], na.rm=TRUE)
    
    ##The if below is neede becasue we combine the AA BB for the var estimate
    ##which comes later
    if(k==2){
      for (i in 1:nrows) scales[i,k] <- mad(tmp[i,], na.rm=TRUE)
    }
    N[,k]=rowSums((!is.na(tmp))); rm(tmp); gc()
  }
  ##now compute the scales for AA,AB
  tmp1 <- M-f
  tmp3 <- M+f; rm(M, f); gc()
  tmp1[initialcalls!=1] <- NA;tmp1[is.na(initialcalls)]<- NA 
  tmp3[initialcalls!=3] <- NA;tmp3[is.na(initialcalls)]<- NA; rm(initialcalls); gc()
  tmp1 <- sweep(tmp1,1,centers[,1])
  tmp3 <- sweep(tmp3,1,centers[,3])
  tmp <- cbind(tmp1,tmp3);
  rm(tmp1, tmp3); gc()
  for (i in 1:nrows) scales[i,1] <- scales[i,2] <- mad(tmp[i,], na.rm=TRUE)
  
  if(verbose) cat(" Done\n")
  return(list(centers=centers,scales=scales,N=N))
}

getAffySnpGenotypeRegionParams<-function(object,initialcalls,f=NULL,
                                         subset=1:(dim(object)[1]),
                                         verbose=FALSE){
  if(is.null(f)) f=fitAffySnpMixture(object,verbose=verbose)$fs
  
  centers <- array(NA,dim=c(length(subset),3,2))
  dimnames(centers)<-list(featureNames(object)[subset],
                          c("AA","AB","BB"),
                          c("antisense","sense"))

  scales <- array(NA,dim=c(length(subset),3,2))
  dimnames(scales) <-list(featureNames(object)[subset],
                          c("AA","AB","BB"),
                          c("antisense","sense"))

  N <- array(NA,dim=c(length(subset),3))
  dimnames(N) <- list(featureNames(object)[subset],c("AA","AB","BB"))

  if(verbose) cat("Computing centers and scales:\n")
  for(s in 1:2){
    gc()
    if(verbose) cat(c("\tantisense","\tsense....")[s],": ",sep="")
    tmp <- getGenotypeRegionParams(getM(object)[subset,,s],
                                   initialcalls[subset,],
                                   f[subset,,s],
                                   verbose=verbose)
    centers[,,s] <- tmp$centers
    scales[,,s] <- tmp$scales
  }
  N <- tmp$N; rm(tmp)
  return(list(centers=centers,scales=scales,N=N,f0=median(f, na.rm=TRUE)))
}

getAffySnpPriors <-  function(object,minN=20,subset=1:(dim(object$centers)[1]),
                              verbose=TRUE){
  if(verbose) cat("Computing priors.\n")
  require(limma,quietly=TRUE)
  N <- cbind(object$N,object$N)
  Index <- subset[which(rowMeans(N[subset,]>minN)==1)]
  N <- N[Index,]
  mus <- cbind(object$centers[Index,,1],object$centers[Index,,2])

  ##variance for prior for mu
  ##need to make robust
  V <- cov(mus,use="pairwise.complete")
  rm(mus)
  
  sigmas <- cbind(object$scales[Index,,1],object$scales[Index,,2])
  maxsigmas <- c(quantile(sigmas[,c(1,3,4,6)],.99, na.rm=TRUE),
                 quantile(sigmas[,c(2,5)],.99, na.rm=TRUE))
  
  ##prior for sigma
  zgs <- log(sigmas^2); rm(sigmas); gc()
  dgs <- N-1
  egs <- zgs-digamma(dgs/2)+log(dgs/2); rm(zgs); gc()
  n <- length(Index)
  d0s <- 2*trigammaInverse(colMeans(n*(egs-colMeans(egs, na.rm=TRUE))^2/(n-1)-trigamma(dgs/2), na.rm=TRUE))
  s20 <- exp(colMeans(egs, na.rm=TRUE)+digamma(d0s/2)-log(d0s/2))
  return(list(V=V,d0s=d0s,s20=s20,maxsigmas=maxsigmas))
}

##minN- minimum number of points in cluster required for use
##max$Sigma
updateAffySnpParams <- function(object,priors,
                                minN=3,
                                maxHomoSigma=priors$maxsigma[1],
                                maxHeteSigma=priors$maxsigma[2],
                                subset=1:(dim(object$centers)[1]),
                                verbose=FALSE){
  if(verbose) cat("Updating centers and scales")
  ##First variances
  for(i in 1:2){
    for(j in 1:2){ ##1 and 3 are the same
      if(j==2) N <- object$N[,2] else N<-rowSums(object$N[,c(1,3)],na.rm=TRUE)
      s <- object$scale[,j,i]
      d0s <- priors$d0s[3*(i-1)+j]
      s20 <- priors$s20[3*(i-1)+j]
      ##notice the ad-hoc choice of 3
      Index <- N>minN & !is.na(s)
      N<-N[Index];s<-s[Index]
      object$scale[Index,j,i] <- sqrt (  ( (N-1)*s^2 + d0s*s20 ) / (d0s+N-1) )
      object$scale[!Index,j,i] <- sqrt(s20)
      rm(s, d0s, s20, Index, N); gc()
    }
  }
  object$scales[,3,] <- object$scales[,1,] ##AA=BB 
  object$scales[,2,][object$scales[,2,]>maxHeteSigma]<-maxHeteSigma
  object$scales[,c(1,3),][object$scales[,c(1,3),]>maxHomoSigma]<-maxHomoSigma
  
  if(verbose) cat(".")

  ##Means
  mu <- cbind(object$centers[,,1],object$centers[,,2])
  N <- cbind(object$N,object$N)

  Vinv <- solve(priors$V)
  mu[is.na(mu)]<-0
  NSinv <- t(N)/priors$s20
  tmp <- t(sapply(1:nrow(mu),
                  function(i){
                    if(verbose)
                      if(i%%5000==0)  cat(".")
                    mus=mu[i,]
                    Ns=N[i,]
                    mus[Ns<minN]<-0
                    return(solve(Vinv+diag(NSinv[,i]))%*%(NSinv[,i]*mus))
                  }))
  rm(mu, N, Vinv, NSinv, minN); gc()
  if(verbose) cat("Done.\n")
  object$centers[,,1] <- tmp[,1:3]
  object$centers[,,2] <- tmp[,4:6]
  return(object)
}

getAffySnpDistance <- function(object,params,f=0,subset=1:(dim(object)[1]),
                               w=NULL,verbose=FALSE){
  x=getM(object)[subset,,,drop=FALSE]
  Dist <- array(NA,dim=c(dim(x)[1],dim(x)[2],3,2))
  if(verbose) cat("Calculating likelihood-based distances")
  for(i in 1:2){
    if(verbose) cat(".")
    for(j in 1:3){
      tmp <- x[,,i]
      if(j==1) tmp <- tmp-f[subset,,i]
      if(j==3) tmp <- tmp+f[subset,,i]
      Dist[,,j,i]=  2*log(params$scales[subset,j,i]) +
        ((tmp-params$centers[subset,j,i])/params$scales[subset,j,i])^2
      if(!is.null(w)) Dist[,,j,i] <-  Dist[,,j,i] - 2*log(w[subset,,j,i])
    }
  }
  dimnames(Dist) <- list(dimnames(x)[[1]],
                         dimnames(x)[[2]],1:3,dimnames(x)[[3]])
  if(verbose) cat("Done.\n")
  return(Dist)
}

getAffySnpCalls <- function(Dist,XIndex,maleIndex,subset=1:(dim(Dist)[1]),
                            verbose=FALSE){
  Dist <- Dist[subset,,,,drop=FALSE]
  XIndex <- which(subset%in%XIndex)
  gc()
  res <- array(NA,dim=dim(Dist)[c(1,2)]); gc()
  dimnames(res) <- dimnames(Dist)[1:2]
  Dist <- rowSums(Dist, na.rm=TRUE, dims=3)
  Dist[XIndex,maleIndex,2] <- Inf

  ##the following is slow!
  if(verbose) cat("Making calls for ",ncol(res)," arrays")
  
  ##apply is faster but I want to see how far along it is
  for(j in 1:ncol(res)){
    if(verbose) cat(".")
    res[,j] <- apply(Dist[,j,],1, which.min)
    gc()
  }
  if(verbose) cat("Done\n")  
  return(res)
}

getAffySnpConfidence <- function(Dist,Calls,XIndex,maleIndex,
                                 subset=1:nrow(Calls),
                                 verbose=TRUE){
  Dist <- Dist[subset,,,,drop=FALSE]
  Calls <- Calls[subset,,drop=FALSE]
  XIndex <- which(subset%in%XIndex)
    
  res <- array(NA,dim=dim(Dist)[c(1,2)])
  dimnames(res)=list(dimnames(Dist)[[1]],dimnames(Dist)[[2]])
##  Dist <- Dist[,,,1]+Dist[,,,2]
  Dist <- rowSums(Dist, na.rm=TRUE, dims=3)

  cat("Making calls for ",ncol(res)," arrays")
  ##apply is faster apply but takes too much memory
  N <- nrow(Calls)
  Index <- 1:N
  for(j in 1:ncol(res)){
    if(maleIndex[j]){
      Index2=Index[-XIndex]
    }
    else{
      Index2=Index
    }
    cat(".")
    tmpdist <- cbind(abs(Dist[,j,1]-Dist[,j,2]),abs(Dist[,j,2]-Dist[,j,3]))
    tmpIndex = split(Index2,Calls[Index2,j])
    res[tmpIndex[[1]],j] <- tmpdist[tmpIndex[[1]],1]
    res[tmpIndex[[3]],j] <- tmpdist[tmpIndex[[3]],2]
    res[tmpIndex[[2]],j] <- apply(tmpdist[tmpIndex[[2]],],1,min)
    rm(tmpIndex, tmpdist); gc()
    if(maleIndex[j]){
      Index2=Index[XIndex]
      res[Index2,j] <- abs(Dist[Index2,j,1]-Dist[Index2,j,3])
    }
  }
  cat("Done\n")  
  return(res)
}

replaceAffySnpParams <- function(object,value,subset){
  object$centers[subset,,] <- value$centers
  object$scales[subset,,] <- value$scales
  object$N[subset,] <- value$N
  return(object)
}

crlmm <- function(object,correction=NULL,recalibrate=TRUE,
                  minLLRforCalls=c(50,40,50),
                  returnCorrectedM=TRUE,
                  returnParams=FALSE,
                  verbose=TRUE){

  if(is.null(correction)){
    if(verbose) cat("M correction not provided. Calculating. Will take several minutes.\n")
     correction=fitAffySnpMixture(object,verbose=verbose)
  }
  if(is.null(object$gender)){
##    if(verbose) cat("Gender not provided... using data to predict.\n")
    warning("Gender not provided... using data to predict.")
    maleIndex <- snpGenderCall(object)=="male"
  }
  else{
    maleIndex <- object$gender=="male"
  }

  load(system.file(paste("data/",annotation(object), "CrlmmInfo.rda", sep=""),
                   package=paste("pd", annotation(object), sep="")))
  
  myenv <- get(paste(annotation(object),"Crlmm",sep=""))

##  params <- get("params",myenv)
##  priors <- get("priors",myenv)
##  hapmapCallIndex <- get("hapmapCallIndex",myenv)
##  badCallIndex <- get("badCallIndex",myenv)
  
##  Index <- which(!hapmapCallIndex  |  badCallIndex )
  Index <- which(!get("hapmapCallIndex",myenv)  |  get("badCallIndex",myenv))

  myCalls <- matrix(NA,dim(object)[1],dim(object)[2])
  
  myCalls[Index,]<-getInitialAffySnpCalls(correction,Index,verbose=verbose)

  rparams <- getAffySnpGenotypeRegionParams(object,
                                           myCalls,
                                           correction$fs,
                                           subset=Index,verbose=verbose)

  gc()

##  rparams<-updateAffySnpParams(rparams,priors,verbose=TRUE)
  rparams<-updateAffySnpParams(rparams, get("priors",myenv),verbose=TRUE)
  
##  params <- replaceAffySnpParams(params,rparams,Index)
  
  params <- replaceAffySnpParams(get("params",myenv), rparams, Index)

  myDist<-getAffySnpDistance(object,params,correction$fs)

  XIndex<-getChrXIndex(object)
  
  myCalls<-getAffySnpCalls(myDist,XIndex,maleIndex,verbose=verbose)

  LLR<-getAffySnpConfidence(myDist,myCalls,XIndex,maleIndex,verbose=verbose)

  if(recalibrate){
    if(verbose) cat("Recalibrating.")

    for(i in 1:nrow(myCalls)){
      for(k in 1:3){
        Index2 <- which(LLR[i,]<minLLRforCalls[k] & myCalls[i,]==k)
        myCalls[i,Index2] <- NA
      }
    }

    rparams<-getAffySnpGenotypeRegionParams(object,
                                            myCalls,
                                            correction$fs,
                                            verbose=verbose)
    gc()

    rparams <- updateAffySnpParams(rparams, get("priors", myenv))

    myDist <- getAffySnpDistance(object,rparams,
                                 correction$fs,
                                 verbose=verbose)
    
    myCalls <- getAffySnpCalls(myDist,XIndex,
                               maleIndex,
                               verbose=verbose)
    
    LLR <- getAffySnpConfidence(myDist,myCalls,XIndex,maleIndex,verbose=verbose)
  }
  ret <- list(calls=myCalls,llr=LLR)
  rm(LLR)
  if(returnCorrectedM){
    cM <-getM(object)
    for(i in 1:nrow(cM)){
      SIGN <- c(1,0,-1)
      for(k in c(1,3)){
        Index=which(myCalls[i,]==k)
        cM[i,Index,] <- cM[i,Index,]+SIGN[k]*(rparams$f0-correction$fs[i,k,])
      }
    }
    ret$M <- cM
    rm(cM)
  }
  if(returnParams){
    ret$params <- rparams
  }
  
  ## correction$snr
  ## maleIndex
  gender <- rep("female", length(maleIndex))
  gender[maleIndex] <- "male"
  if (is.null(object$gender)){
    addPhenoData <- new("AnnotatedDataFrame",
                        data=cbind(pData(object),
                          data.frame(crlmmSNR=as.numeric(correction$snr),
                                     gender=gender,
                                     row.names=sampleNames(object))),
                        varMetadata= rbind(varMetadata(object),
                          data.frame(labelDescription=c("crlmmSNR", "gender"),
                                     row.names=c("crlmmSNR", "gender"))))
  }else{
    addPhenoData <- new("AnnotatedDataFrame",
                        data=cbind(pData(object),
                          data.frame(crlmmSNR=as.numeric(correction$snr),
                                     row.names=sampleNames(object))),
                        varMetadata= rbind(varMetadata(object),
                          data.frame(labelDescription=c("crlmmSNR"),
                                     row.names=c("crlmmSNR"))))
  }    
  
  if (!returnParams){
    return(new("SnpCallSetPlus",
               phenoData=addPhenoData,
               experimentData=experimentData(object),
               annotation=annotation(object),
               calls=ret$calls,
               callsConfidence=ret$llr,
               logRatioAntisense=ret$M[,,"antisense"],
               logRatioSense=ret$M[,,"sense"]))
  }else{
    return(ret)
  }
}


###this one just for us
addRegions <- function(i,params,...){
  require(ellipse)
  ADD <- params$f0*c(1,0,-1)
  for(k in 1:3){
    points(t(params$centers[i,k,])+ADD[k],pch="+",col=k)
    lines(ellipse(diag(2),scale=params$scales[i,k,],centre=params$centers[i,k,]+ADD[k]),col=k,...)
  }
}
