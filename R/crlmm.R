rowEntropy <- function(p){
  p <- ifelse(p>0,p*log2(p),0)
  rowSums(p)
}

getSnpFragmentLength <- function(object){
  annotname <- annotation(object)
  data(list=annotname)
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
  data(list=annotname)
  annot <- annot[match(featureNames(object),annot$SNP),]
  return(which(annot$Chromosome=="chrX"))
}

##gender in pData keeps male female
fitAffySnpMixture <- function(object, df1=3, df2=5,
                              probs=rep(1/3,3), eps=50,
                              subSampleSize=10^5,verbose=TRUE){
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
    
    ## this is the original BC - Jul 18, 2006
    ## l <- L[idx];L=L-mean(l);l=l-mean(l)
    ## a <- A[idx];A=A-mean(A);a=a-mean(a)
    ## y <- Y[idx]

    l <- L[idx];L=L-mean(L, na.rm=TRUE);l=l-mean(l, na.rm=TRUE)
    a <- A[idx];A=A-mean(A, na.rm=TRUE);a=a-mean(a, na.rm=TRUE)
    y <- Y[idx]
    ok <- complete.cases(cbind(y,a,l))
    l <- l[ok]; a <- a[ok]; y <- y[ok]
    
    weights <- apply(cbind(mus,sigmas),1,function(p) dnorm(y,p[1],p[2]))
    PreviousLogLik <- -Inf
    change <- eps+1
    itmax <- 0
    while (change > eps & itmax < 1000){
      itmax <- itmax+1
      
      ## E
      z <- sweep(weights, 2, probs, "*")
      LogLik <- rowSums(z)
      z <- sweep(z, 1, LogLik, "/")
      LogLik <- sum(log(LogLik))
      change <- abs(LogLik-PreviousLogLik)
##      if(verbose) cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
      if(verbose) cat("Array ",j,": epsilon=",round(change,2),"    \n",sep="")
      if(verbose) cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")

      
      PreviousLogLik <- LogLik
      probs <- colMeans(z)

      ## M
      fit1 <- lm(y~ns(l,df1)+ns(a,df2),weights=z[,1])
      fit2 <- sum(z[,2]*y)/sum(z[,2])
      fit3 <- lm(y~ns(l,df1)+ns(a,df2),weights=z[,3])
      
      sigmas[1] <- sqrt(sum(z[,1]*residuals(fit1)^2)/sum(z[,1]))
      sigmas[2] <- sqrt(sum(z[,2]*(y-fit2)^2)/sum(z[,2]))
      sigmas[3] <- sqrt(sum(z[,3]*residuals(fit3)^2)/sum(z[,3]))
      
      weights[,1] <- dnorm(y,fitted(fit1),sigmas[1])
      weights[,2] <- dnorm(y,fit2,sigmas[2])
      weights[,3] <- dnorm(y,fitted(fit3),sigmas[3])
      weights[y >= 0, 1] <- 0
      weights[y <= 0, 3] <- 0
    }

    pred1 <- predict(fit1,newdata=data.frame(l=L,a=A))
    pred2 <- rep(fit2,length(Y))
    pred3 <- predict(fit3,newdata=data.frame(l=L,a=A))

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
  pi2 <- object$pis[subset,,,2]
  if(verbose) cat("Picking good starting value: ")
  if(verbose) cat("Computing entropy, ")
  E1<-apply(pi1,1,function(x) mean(rowEntropy(x)))
  E2<-apply(pi2,1,function(x) mean(rowEntropy(x)))
  
  if(verbose) cat("calculating calls, ")
  tmpcall1 <- apply(pi1,c(1,2),which.max)
  tmpcall2 <- apply(pi2,c(1,2),which.max)
  
  if(verbose) cat("determining non-concordant calls, ")
  concordance <- rowIndepChiSqTest(tmpcall1,tmpcall2)

  if(verbose) cat("deciding which strand(s) to use")
  noABIndex1 <- (rowSums(tmpcall1==2)<3)*(rowSums(tmpcall1==3)>0)*(rowSums(tmpcall1==1)>0)
  noABIndex2 <- (rowSums(tmpcall2==2)<3)*(rowSums(tmpcall2==3)>0)*(rowSums(tmpcall2==1)>0)
  
  E1[which(noABIndex1==1)] <- -Inf
  E2[which(noABIndex1==1)] <- -Inf

  jointprobs<-(pi1+pi2)/2
  ##NA if all 0
  for(i in 1:(dim(pi1)[1])){
    if(noABIndex1[i]==1 & noABIndex2[i]==1){
      jointprobs[i,,] <- 1/3 ##no info
    } 
    else{
      ##Not Both
      if(concordance[i]<concordanceCutoff){
        ##use one or the other
        if(E1[i] >  E2[i]){
          jointprobs[i,,] <- pi1[i,,]
        }
        else{
          jointprobs[i,,] <- pi2[i,,]
        }
      }
    }
  }
  if(verbose) cat("finalizing")
  tmpcall <- apply(jointprobs,c(1,2),which.max)
  if(verbose) cat(".")
  tmpmax =  apply(jointprobs,c(1,2),max) 
  for(i in 1:3)
    tmpcall[tmpcall==i & tmpmax<cutoffs[i]] <- NA
  if(verbose) cat("Done!\n")
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
  
  if(verbose) cat("Computing centers and scales for 3 genotypes")
  for(k in 1:3){
    if(verbose) cat(".")
    if(k==1) tmp <- M - f
    if(k==2) tmp <- M
    if(k==3) tmp <- M + f
    tmp[initialcalls!=k] <- NA
    tmp[is.na(initialcalls)]<- NA
    centers[,k]<-apply(tmp,1,median,na.rm=TRUE)
    ##The if below is neede becasue we combine the AA BB for the var estimate
    ##which comes later
    if(k==2){
      scales[,k]=apply(tmp,1,mad,na.rm=TRUE)
    }
    N[,k]=rowSums((!is.na(tmp)))
  }
  ##now compute the scales for AA,AB
  tmp1 <- M-f
  tmp3 <- M+f
  tmp1[initialcalls!=1] <- NA;tmp1[is.na(initialcalls)]<- NA 
  tmp3[initialcalls!=3] <- NA;tmp3[is.na(initialcalls)]<- NA
  tmp1 <- sweep(tmp1,1,centers[,1])
  tmp3 <- sweep(tmp3,1,centers[,3])
  tmp=cbind(tmp1,tmp3)
  scales[,1] <- apply(tmp,1,mad,na.rm=TRUE)
  scales[,3] <- scales[,1]
  
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
    if(verbose) cat(c("\tantisense","\tsense")[s],":",sep="")
    tmp <- getGenotypeRegionParams(getM(object)[subset,,s],
                                   initialcalls[subset,],
                                   f[subset,,s],
                                   verbose=verbose)
    centers[,,s] <- tmp$centers
    scales[,,s] <- tmp$scales
  }
  N <- tmp$N
  return(list(centers=centers,scales=scales,N=N,f0=median(f)))
}

getAffySnpPriors <-  function(object,minN=20,subset=1:(dim(object$centers)[1]),
                              verbose=TRUE){
  if(verbose) cat("Computing priors.\n")
  require(limma,quietly=TRUE)
  N <- cbind(object$N,object$N)
  Index <- subset[which(rowMeans(N[subset,]>minN)==1)]
  N <- N[Index,]
  mus <- cbind(object$centers[Index,,1],object$centers[Index,,2])
  sigmas <- cbind(object$scales[Index,,1],object$scales[Index,,2])
  maxsigmas <- c(quantile(sigmas[,c(1,3,4,6)],.99, na.rm=TRUE),
                 quantile(sigmas[,c(2,5)],.99, na.rm=TRUE))
  
  ##variance for prior for mu
  ##need to make robust
  V <- cov(mus,use="pairwise.complete")
  
  ##prior for sigma
  zgs <- log(sigmas^2)
  dgs <- N-1
  egs <- zgs-digamma(dgs/2)+log(dgs/2)
  n <- length(Index)
  d0s <- 2*trigammaInverse(colMeans(n*(egs-colMeans(egs))^2/(n-1)-trigamma(dgs/2)))
  s20 <- exp(colMeans(egs)+digamma(d0s/2)-log(d0s/2))
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
  tmp <- t(sapply(1:nrow(mu),function(i){
    if(verbose){if(i%%5000==0)  cat(".")}
    mus=mu[i,]
    Ns=N[i,]
    mus[Ns<minN]<-0
    return(solve(Vinv+diag(NSinv[,i]))%*%(NSinv[,i]*mus))
  }))
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
  
  res <- array(NA,dim=dim(Dist)[c(1,2)])
  dimnames(res)=list(dimnames(Dist)[[1]],dimnames(Dist)[[2]])
  Dist <- Dist[,,,1]+Dist[,,,2]
  
  Dist[XIndex,maleIndex,2] <- Inf
  ##the following is slow!
  if(verbose) cat("Making calls for ",ncol(res)," arrays")
  ##apply is faster but I want to see how far along it is
  for(j in 1:ncol(res)){
    if(verbose) cat(".")
    res[,j] <- apply(Dist[,j,],1,which.min)
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
  Dist <- Dist[,,,1]+Dist[,,,2]
  
  cat("Making calls for ",ncol(res)," arrays")
  ##apply is faster apply but takes too much memory
  N<-nrow(Calls)
  Index<-1:N; ##browser()
  for(j in 1:ncol(res)){
    if(maleIndex[j]){
##      Index2=Index[!XIndex]
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
                  returnCorrectedM=FALSE,
                  returnParams=FALSE,
                  verbose=TRUE){
  require(paste("pd", annotation(object), sep=""), character.only=TRUE)
  if(is.null(correction)){
    if(verbose) cat("M correction not provided. Calculating. Will take several minutes.\n")
     correction=fitAffySnpMixture(object,verbose=verbose)
  }
  if(is.null(object$gender)){
    if(verbose) cat("Gender not provided... using data to predict.\n")
    maleIndex <- snpGenderCall(object)=="male"
  }
  else{
    maleIndex <- object$gender=="male"
  }
  
##  load(paste(annotation(object),"CrlmmInfo.rda",sep=""))
  data(list=paste(annotation(object),"CrlmmInfo",sep=""))  
  myenv <- get(paste(annotation(object),"Crlmm",sep=""))
  params <- get("params",myenv)
  priors <- get("priors",myenv)
  hapmapCallIndex <- get("hapmapCallIndex",myenv)
  badCallIndex <- get("badCallIndex",myenv)
  
  Index <- which(!hapmapCallIndex  |  badCallIndex )

  myCalls <- matrix(NA,dim(object)[1],dim(object)[2])
  
  myCalls[Index,]<-getInitialAffySnpCalls(correction,Index,verbose=verbose)

  rparams <- getAffySnpGenotypeRegionParams(object,
                                           myCalls,
                                           correction$fs,
                                           subset=Index,verbose=verbose)

  rparams<-updateAffySnpParams(rparams,priors,verbose=TRUE)
  
  params <- replaceAffySnpParams(params,rparams,Index)

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
    rparams<-updateAffySnpParams(rparams,priors)
    
    myDist <- getAffySnpDistance(object,rparams,
                                 correction$fs,
                                 verbose=verbose)
    
    myCalls <- getAffySnpCalls(myDist,XIndex,
                               maleIndex,
                               verbose=verbose)
    
    LLR <- getAffySnpConfidence(myDist,myCalls,XIndex,maleIndex,verbose=verbose)
  }
  ret <- list(calls=myCalls,llr=LLR)
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
  }
  if(returnParams){
    ret$params <- rparams
  }
  return(ret)
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
