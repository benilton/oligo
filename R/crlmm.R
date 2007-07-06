### WARNING: IF ROBUST TO DIMNAMES,
###          REMOVE AND SAVE LOTS
###          OF MEMORY... ADD DIMNAMES
###          AT THE END

rowEntropy <- function(p) rowMeans(rowSums(log2(p^p), dims=2))

getSnpFragmentLength <- function(object){
  sql <- "SELECT fragment_length FROM featureSet WHERE man_fsetid LIKE 'SNP%' ORDER BY man_fsetid"
  return(dbGetQuery(db(get(annotation(object))), sql)[[1]])
}

snpGenderCall <- function(object){
  XIndex=getChrXIndex(object)
  tmp=(median(getA(object), na.rm=TRUE))
  a=apply(getA(object)[XIndex,,],2,median, na.rm=TRUE)
  kfit=kmeans(a,c(min(a),tmp))
  return(factor(c("female","male")[as.numeric(kfit$cluster==1)+1]))
}

getChrXIndex <- function(object){
  sql <- "SELECT chrom FROM featureSet WHERE man_fsetid LIKE 'SNP%' ORDER BY man_fsetid"
  chrs <- dbGetQuery(db(get(annotation(object))), sql)[[1]]
  return(which(chrs == "X"))
}

##gender in pData keeps male female
fitAffySnpMixture <- function(object, df1=3, df2=5,
                              probs=rep(1/3,3), eps=50,
                              subSampleSize=10^4,verbose=TRUE){
  if(is.null(object$gender)){
    maleIndex <- snpGenderCall(object)=="male"
  }else{
    maleIndex <- object$gender=="male"
  }
  
  XIndex=getChrXIndex(object)
  
  I <- dim(object)[1]
  J <- dim(object)[2]
  set.seed(1)
##  tmp <- c( (1:I)[-XIndex],((I+1):(2*I))[-XIndex])
  tmp <- c( (1:I),((I+1):(2*I)))
  idx <- sort(sample(tmp, subSampleSize))
  rm(tmp)

  pis <- array(0,dim=c(I,J,3,2))
  fs <- array(0,dim=c(I,J,2))
  snr <- array(0,dim=J)

#####  dimnames(fs)<-list(featureNames(object),
#####                     sampleNames(object),
#####                     c("antisense","sense"))
#####  dimnames(pis)<-list(featureNames(object),
#####                      sampleNames(object),
#####                      c("AA","AB","BB"),
#####                      c("antisense","sense"))
#####  names(snr) <- sampleNames(object)
  
  if(verbose) cat("Fitting mixture model to ",J," arrays. Epsilon must reach ",eps,".\n",sep="")
  L <- getSnpFragmentLength(object)
  fix <- which(is.na(L))
  L[fix] <- median(L, na.rm=T)
  rm(fix)
  L <- c(L,L)
  l <- L[idx];L=L-mean(L);l=l-mean(l)
  matL <- ns(l,df1)
  for(j in 1:J){
    Y <- c(as.vector(getM(object[,j])))
    A <- c(as.vector(getA(object[,j])))
    fix <- which(is.na(Y))
    Y[fix] <- median(Y, na.rm=T)
    A[fix] <- median(A, na.rm=T)
    rm(fix)

    mus <- quantile(Y,c(1,3,5)/6);mus[2]=0
    sigmas <- rep(mad(c(Y[Y<mus[1]]-mus[1],Y[Y>mus[3]]-mus[3])),3)
    sigmas[2] <- sigmas[2]/2
    
    a <- A[idx];A=A-mean(A);a=a-mean(a)
    y <- Y[idx]
    
    weights <- apply(cbind(mus, sigmas), 1, function(p) dnorm(y, p[1], p[2]))
    PreviousLogLik <- -Inf
    change <- eps+1
    itmax <- 0
    
    matA <- ns(a,df2)
    while (change > eps & itmax < 100){
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
  
    fs[,j,] <- matrix((pred3-pred1)/2,ncol=2)
    for(k in 1:3){
      pis[,j,k,] <- matrix(z[,(4-k)],ncol=2) ##4-k cause 3is1,2is2 and 1is3
    }
    snr[j] <- median(fs[,j,])^2/(sigmas[1]^2+sigmas[2]^2)
  }
  fs[fs < 0] <- 0
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

## getGenotypeRegionParams <- function(M, initialcalls, f=0, verbose=TRUE){
##   require(MASS)
##   if(!is.matrix(M)) M <- matrix(M,ncol=2)
##   centers <- scales <- N <- array(NA,dim=c(nrow(M),3))
##   dimnames(centers) <- dimnames(scales) <- dimnames(N) <- list(rownames(M),c("AA","AB","BB"))
##   nrows <- nrow(M)
##   if(verbose) cat("Computing centers and scales for 3 genotypes")
##   for(k in 1:3){
##     if(verbose) cat(".")
##     if(k==1) tmp <- M - f
##     if(k==2) tmp <- M
##     if(k==3) tmp <- M + f
##     tmp[initialcalls!=k] <- NA
##     tmp[is.na(initialcalls)]<- NA
## 
##     ## THIS IS SLOW! /bc
##     for (i in 1:nrows){
##       v <- tmp[i,]
##       v <- v[!is.na(v)]
##       set.seed(1)
##       if (length(v) > 2){
##         centers[i,k] <- hubers(v)$mu
##       }else{
##         centers[i,k] <- 0
##       }
##       
##       ##      centers[i,k] <- median(tmp[i,], na.rm=TRUE)
## 
##     }
##     
##     ##The if below is neede becasue we combine the AA BB for the var estimate
##     ##which comes later
##     if(k==2){
##       ## SLOW
##       for (i in 1:nrows){
##         v <- tmp[i,]
##         v <- v[!is.na(v)]
##         set.seed(1)
##         if (length(v)>2){
##           scales[i,k] <- as.numeric(hubers(v)$s)
##         }else{
##           scales[i,k] <- 0.3
##         }
##         
##         ## scales[i,k] <- mad(tmp[i,], na.rm=TRUE)
## 
##       }
##     }
##     N[,k]=rowSums((!is.na(tmp))); rm(tmp); ## gc()
##   }
##   ##now compute the scales for AA,AB
##   tmp1 <- M-f
##   tmp3 <- M+f; rm(M, f); ## gc()
##   tmp1[initialcalls!=1] <- NA;tmp1[is.na(initialcalls)]<- NA 
##   tmp3[initialcalls!=3] <- NA;tmp3[is.na(initialcalls)]<- NA; rm(initialcalls); ## gc()
##   tmp1 <- sweep(tmp1,1,centers[,1])
##   tmp3 <- sweep(tmp3,1,centers[,3])
##   tmp <- cbind(tmp1,tmp3);
##   rm(tmp1, tmp3); ## gc()
##   for (i in 1:nrows){
##     ## SLOW
##     v <- tmp[i,]
##     v <- v[!is.na(v)]
##     
##     ## scales[i,1] <- scales[i,3] <- mad(tmp[i,], na.rm=TRUE)
## 
##     set.seed(1)
##     if (length(v) > 2){
##       scales[i, c(1,3)] <- as.numeric(hubers(v)$s)
##     }else{
##       scales[i, c(1,3)] <- .3
##     }
##   }
##   
##   if(verbose) cat(" Done\n")
##   return(list(centers=centers,scales=scales,N=N))
## }

getGenotypeRegionParams <- function(M, initialcalls, f=0, verbose=TRUE){
  if(verbose) cat("Computing centers and scales for 3 genotypes")
  
##   tmp <- .Call("R_HuberMatrixRows2", M+(initialcalls-2)*f,
##                as.integer(initialcalls), 1.5)
  
  tmp <- .Call("R_trimmed_stats", M+(initialcalls-2)*f,
               as.integer(initialcalls), 0.025)


#####  centers <- scales <- N <- array(NA, dim=c(nrow(M),3))
#####  scales <- N <- array(NA, dim=c(nrow(M),3))
#####  dimnames(centers) <- dimnames(scales) <- dimnames(N) <- list(rownames(M), c("AA","AB","BB"))
#####  centers[,] <- tmp[[1]]

  tmpN <- tmp[[3]]; tmp[[3]] <- NULL
  tmpN[tmpN <= 1] <- NA
  tmpS <- tmp[[2]]; tmp[[2]] <- NULL
  tmpS[is.na(tmpN)] <- NA
  tmpV <- sqrt(rowSums((tmpN[,-2]-1)*tmpS[,-2]^2, na.rm=T)/rowSums(tmpN[,-2]-1, na.rm=T))
  tmpV[!is.finite(tmpV)] <- NA
  tmpS[,-2] <- tmpV
  rm(tmpV)
  tmpC <- tmp[[1]]; rm(tmp)
  tmpC[is.na(tmpN)] <- NA
  tmpN[is.na(tmpN)] <- 0
  if(verbose) cat(" Done\n")
  return(list(centers=tmpC, scales=tmpS, N=tmpN))
}

getAffySnpGenotypeRegionParams<-function(object,initialcalls,f=NULL,
                                         subset=1:(dim(object)[1]),
                                         verbose=FALSE){
  if(is.null(f)) f=fitAffySnpMixture(object,verbose=verbose)$fs
  N <- scales <- centers <- array(NA,dim=c(length(subset),3,2))
#####  dimnames(N) <- dimnames(scales) <- dimnames(centers) <- list(featureNames(object)[subset],
#####                                                               c("AA","AB","BB"), c("antisense","sense"))
  if(verbose) cat("Computing centers and scales:\n")
  for(s in 1:2){
    ## gc()
    if(verbose) cat(c("\tantisense","\tsense....")[s],": ",sep="")
    tmp <- getGenotypeRegionParams(getM(object)[subset,,s],
                                   initialcalls[subset,],
                                   f[subset,,s],
                                   verbose=verbose)
    centers[,,s] <- tmp$centers
    scales[,,s] <- tmp$scales
    N[,,s] <- tmp$N
  }
  N <- apply(N, 1:2, max, na.rm=TRUE) ## returns -Inf if c(NA, NA)
  return(list(centers=centers,scales=scales,N=N,f0=median(f, na.rm=TRUE)))
}

getAffySnpPriors <-  function(object,minN=20,subset=1:(dim(object$centers)[1]),
                              verbose=TRUE){
  if(verbose) cat("Computing priors.\n")
  N <- cbind(object$N,object$N)
  Index <- subset[which(rowMeans(N[subset,]>minN)==1)]
  N <- N[Index,]
  mus <- cbind(object$centers[Index,,1],object$centers[Index,,2])
  sigmas <- cbind(object$scales[Index,,1],object$scales[Index,,2])
  sigmas[is.na(mus)] <- NA
  maxsigmas <- c(quantile(sigmas[,c(1,3,4,6)],.99, na.rm=TRUE),
                 quantile(sigmas[,c(2,5)],.99, na.rm=TRUE))
  
  ##variance for prior for mu
  ##need to make robust
  V <- cov(mus,use="pairwise.complete")
  
  ##prior for sigma
  zgs <- log(sigmas^2); rm(sigmas); ## gc()
  dgs <- N-1
  egs <- zgs-digamma(dgs/2)+log(dgs/2); rm(zgs); ## gc()
  n <- length(Index)
  d0s <- 2*trigammaInverse(colMeans(n*(egs-colMeans(egs, na.rm=TRUE))^2/(n-1)-trigamma(dgs/2), na.rm=TRUE))
  s20 <- exp(colMeans(egs, na.rm=TRUE)+digamma(d0s/2)-log(d0s/2))
  return(list(V=V,d0s=d0s,s20=s20,maxsigmas=maxsigmas))
}

##minN- minimum number of points in cluster required for use
##max$Sigma
updateAffySnpParams <- function(object, priors, missingStrandIndex, minN=3,
                                maxHomoSigma=priors$maxsigma[1],
                                maxHeteSigma=priors$maxsigma[2],
                                subset=1:(dim(object$centers)[1]),
                                d0s=80, verbose=FALSE){
  object$centers <- object$centers[subset,,]
  object$scales <- object$scales[subset,,]
  object$N <- object$N[subset,]
  missingStrandIndex <- missingStrandIndex[subset]
  if(verbose) cat("Updating centers and scales")
  ##First variances
  for(i in 1:2){
    for(j in 1:2){ ##1 and 3 are the same
      if(j==2) N <- object$N[,2] else N <- rowSums(object$N[,-2],na.rm=TRUE)
      s <- object$scales[,j,i]
      if (is.null(d0s))
        d0s <- priors$d0s[3*(i-1)+j]
      s20 <- priors$s20[3*(i-1)+j] ##notice the ad-hoc choice of 3
      Index <- N>minN & !is.na(s)
      N <- N[Index]; s <- s[Index]
      object$scales[Index,j,i] <- sqrt (  ( (N-1)*s^2 + d0s*s20 ) / (d0s+N-1) )
      object$scales[!Index & missingStrandIndex != i,j,i] <- sqrt(s20)
    }
  }
  object$scales[,3,] <- object$scales[,1,] ##AA=BB 
  object$scales[,2,][object$scales[,2,]>maxHeteSigma] <- maxHeteSigma
  object$scales[,-2,][object$scales[,-2,]>maxHomoSigma] <- maxHomoSigma
  if(verbose) cat(".")

  ##Means

  updateMean <- function(strandIndex, strand){
    ## strandIndex is a vector where
    ## 0 - none strands are missing on the array
    ## 1 - antisense strand is missing
    ## 2 - sense strand is missing
    ## strand argument is what do you want to be updated
    
    if (strand == "both"){
      snps <- strandIndex == 0; idx <- 1:6
      N <- cbind(object$N[snps,],object$N[snps,])
      mu <- cbind(object$centers[snps,,1],object$centers[snps,,2])
    }else if(strand == "antisense"){
      snps <- strandIndex == 2; idx <- 1:3
      N <- object$N[snps,]
      mu <- object$centers[snps,,1]
    }else if(strand == "sense"){
      snps <- strandIndex == 1; idx <- 4:6
      N <- object$N[snps,]
      mu <- object$centers[snps,,2]
    }
    Vinv <- solve(priors$V[idx, idx])
    NSinv <- t(N)/priors$s20[idx]
    tmp <- t(sapply(1:nrow(mu),function(i){
      if(verbose & i%%5000==0)  cat(".")
      mus=mu[i,]; Ns=N[i,]
      mus[Ns<minN] <- 0
      mus[is.na(mus)] <- 0
      return(solve(Vinv+diag(NSinv[,i]))%*%(NSinv[,i]*mus))
    }))
    return(tmp)
  }
  if (any(missingStrandIndex == 0)){
    tmp <- updateMean(missingStrandIndex, "both")
    object$centers[missingStrandIndex == 0,,1] <- tmp[,1:3]
    object$centers[missingStrandIndex == 0,,2] <- tmp[,4:6]; rm(tmp);
  }
  if (any(missingStrandIndex == 2))
    object$centers[missingStrandIndex == 2,,1] <- updateMean(missingStrandIndex, "antisense")
  if (any(missingStrandIndex == 1))
    object$centers[missingStrandIndex == 1,,2] <- updateMean(missingStrandIndex, "sense")
  if(verbose) cat("Done.\n")
  return(object)
}

getAffySnpDistance <- function(object,params,f=0,subset=1:(dim(object)[1]),
                               w=NULL,verbose=FALSE){
  x=getM(object[subset,])
  rm(object)
  Dist <- array(NA,dim=c(dim(x)[1],dim(x)[2],3,2))
  if(verbose) cat("Calculating likelihood-based distances")
  for(i in 1:2){
    if(verbose) cat(".")
    for(j in 1:3){
      tmp <- x[,,i]+(j-2)*f[subset,,i]
      Dist[,,j,i] <- 2*log(params$scales[subset,j,i]) +
        ((tmp-params$centers[subset,j,i])/params$scales[subset,j,i])^2
      if(!is.null(w)) Dist[,,j,i] <-  Dist[,,j,i] - 2*log(w[subset,,j,i])
    }
  }
#####   dimnames(Dist) <- list(dimnames(x)[[1]],
#####                          dimnames(x)[[2]], 1:3,
#####                          dimnames(x)[[3]])
  if(verbose) cat("Done.\n")
  return(Dist)
}

getAffySnpCalls <- function(Dist,XIndex,maleIndex,subset=1:(dim(Dist)[1]),
                            verbose=FALSE){
  Dist <- Dist[subset,,,,drop=FALSE]
####  XIndex <- which(subset%in%XIndex)
  ## gc()
  ## Here we can even put default value to 3 for the new code. /HB
  res <- array(as.integer(-1),dim=dim(Dist)[c(1,2)]); ## gc()
#####  dimnames(res) <- dimnames(Dist)[1:2]
  Dist <- rowSums(Dist, na.rm=TRUE, dims=3)
  ##  Dist[XIndex, maleIndex, 2] <- Inf
  
  ##the following is slow!
  if(verbose) cat("Making calls for ", ncol(res), " arrays");

  for(j in 1:ncol(res)){
    if(verbose) cat(".");
    D1 <- Dist[,j,1];
    D2 <- Dist[,j,2];
    D3 <- Dist[,j,3];
    d12 <- (D1 < D2);
    d23 <- (D2 < D3); rm(D2);
    d13 <- (D1 < D3); rm(D3);
    d <- rep(as.integer(3), length(D1)); rm(D1)
    d[( d12 & d13)] <- as.integer(1); rm(d13)
    d[(!d12 & d23)] <- as.integer(2); rm(d12, d23)
    res[,j] <- d;
    rm(d);
  }
  if(verbose) cat("Done\n")
  return(res)
}

getAffySnpConfidence <- function(Dist, Calls, XIndex, maleIndex,
                                 subset=1:nrow(Calls), verbose=TRUE){
  Dist <- rowSums(Dist[subset,,,,drop=FALSE], dims=3, na.rm=T)
  Calls <- Calls[subset,,drop=FALSE]
####  XIndex <- which(subset%in%XIndex)

  res <- array(NA,dim=dim(Dist)[c(1,2)])
#####  dimnames(res) <- list(dimnames(Dist)[[1]],dimnames(Dist)[[2]])
#####  Dist <- rowSums(Dist, dims=3, na.rm=T)
  
  cat("Computing confidence for calls on ",ncol(res)," arrays")
  ##apply is faster apply but takes too much memory
####  Index <- 1:nrow(Calls)
  Index2 <- 1:nrow(Calls)
  for(j in 1:ncol(res)){
##    if(maleIndex[j]){
##      Index2 <- Index[-XIndex]
##    }else{
##      Index2 <- Index
##    }
##    Index2 <- Index
    if (verbose) cat(".")
    tmpdist <- cbind(abs(Dist[,j,1]-Dist[,j,2]),abs(Dist[,j,2]-Dist[,j,3]))
    tmpIndex <- split(Index2, factor(Calls[Index2,j], levels=1:3), drop=FALSE)
    if (length(tmpIndex[[1]])>0) res[tmpIndex[[1]],j] <- tmpdist[tmpIndex[[1]],1]
    if (length(tmpIndex[[3]])>0) res[tmpIndex[[3]],j] <- tmpdist[tmpIndex[[3]],2]
    if (length(tmpIndex[[2]])>0) res[tmpIndex[[2]],j] <- pmin(tmpdist[tmpIndex[[2]], 1],
                                                              tmpdist[tmpIndex[[2]], 2])
    rm(tmpIndex, tmpdist); ## gc()
##     if(maleIndex[j]){
##       Index2 <- Index[XIndex]
##       res[Index2,j] <- abs(Dist[Index2,j,1]-Dist[Index2,j,3])
##     }
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

crlmm <- function(object, correction=NULL, recalibrate=TRUE,
                  minLLRforCalls=c(5, 1, 5),
                  verbose=TRUE, correctionFile=NULL, prefix="tmp.crlmm.", balance=1.5){
  library(annotation(object), character.only=TRUE)
  if(is.null(correctionFile))
    stop("Provide correctionFile.\nIf the correctionFile is not found, it will be created and it will contain the EM results.")

  if(file.exists(correctionFile)){
    cat("File with correction information - from EM - was found.", "Loading...", sep="\n")
    load(correctionFile)
    cat("Done.\n")
  }else{
    if(verbose) cat("M correction not provided. Calculating. Will take several minutes.\n")
    correction <- fitAffySnpMixture(object,verbose=verbose)
    save(correction, file=correctionFile)
  }
  snr <- correction$snr

  if(is.null(object$gender)){
    if(verbose) cat("Gender not provided... using data to predict.\n")
    warning("Gender not provided... using data to predict.")
    maleIndex <- snpGenderCall(object)=="male"
  }else{
    maleIndex <- object$gender=="male"
  }
  annotname <- annotation(object)
  load(system.file(paste("extdata/", annotname, "CrlmmInfo.rda", sep=""), package=annotname))
  myenv <- get(paste(annotname,"Crlmm",sep="")); rm(list=paste(annotname,"Crlmm",sep=""))
  thePriors <- get("priors", myenv)

  ## Index <- which(!get("hapmapCallIndex",myenv)  |  get("badCallIndex",myenv) | get("badRegions", myenv))

  Index <- which(!get("hapmapCallIndex",myenv))

  myCalls <- matrix(NA,dim(object)[1],dim(object)[2])

  myCalls[Index,] <- getInitialAffySnpCalls(correction,Index,verbose=verbose)
  fs <- correction$fs; rm(correction)
  rparams <- getAffySnpGenotypeRegionParams(object, myCalls, fs,
                                            subset=Index,verbose=verbose)
  rm(myCalls); ## gc()
  oneStrand <- apply(is.na(getM(object[,1])[,1,]), 1,
                     function(v) ifelse(length(ll <- which(v))==0, 0, ll))
  rparams <- updateAffySnpParams(rparams, thePriors, oneStrand, verbose=TRUE)
  params  <- replaceAffySnpParams(get("params",myenv), rparams, Index)
  save(rparams, params, file=paste(prefix, "ParamsBeforeRec.rda", sep=""))
  rm(myenv, Index)
  myDist <- getAffySnpDistance(object, params, fs)
  myDist[,,-2,] <- balance*myDist[,,-2,]
##  rm(params)
  rm(fs); ## gc()
  XIndex <- getChrXIndex(object)
  myCalls <- getAffySnpCalls(myDist,XIndex,maleIndex,verbose=verbose)
  LLR <- getAffySnpConfidence(myDist,myCalls,XIndex,maleIndex,verbose=verbose)
  rm(myDist)
  load(correctionFile)
  fs <- correction$fs; rm(correction)
  if(recalibrate){
    if(verbose) cat("Recalibrating.")
    for(k in 1:3)
      myCalls[myCalls == k & LLR < minLLRforCalls[k]] <- NA
    rm(LLR)
    myCalls[, snr < 3.675] <- NA
    
    rparams <- getAffySnpGenotypeRegionParams(object, myCalls,
                                              fs, verbose=verbose)
    rm(myCalls)
    ## gc()

    rparams <- updateAffySnpParams(rparams, thePriors, oneStrand)
    save(rparams, file=paste(prefix, "ParamsAfterRec.rda", sep=""))
    myDist <- getAffySnpDistance(object,rparams, fs, verbose=verbose)
    myDist[,,-2,] <- balance*myDist[,,-2,]
    myCalls <- getAffySnpCalls(myDist,XIndex, maleIndex, verbose=verbose)
    LLR <- getAffySnpConfidence(myDist,myCalls,XIndex,maleIndex,verbose=verbose)
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
  return(new("SnpCallSet",
             phenoData=addPhenoData,
             experimentData=experimentData(object),
             annotation=annotation(object),
             calls=myCalls,
             callsConfidence=LLR,
             pAcc=pacc,
             featureData=featureData(object)))
}


###this one just for us
### addRegions <- function(i,params,...){
###   require(ellipse)
###   ADD <- params$f0*c(1,0,-1)
###   idx <- which(!apply(is.na(params$centers[i,,]), 2, all))
###   if (length(idx) == 1) idx <- rep(idx, 2)
###   for(k in 1:3){
###     points(t(params$centers[i,k,idx])+ADD[k],pch="+",col=k)
###     lines(ellipse(diag(2),scale=params$scales[i,k,idx],centre=params$centers[i,k,idx]+ADD[k]),col=k,...)
###   }
### }

LLR2conf <-function(theCalls, LLR, SNR, annot){
  load(paste(system.file("extdata", package=annot), "/", annot, ".spline.params.rda", sep=""))

  Het <- as.vector(theCalls==2)
  dst <- rep(Dst, ncol(theCalls))
  LLR <- as.vector(sqrt(LLR))

  conf <- vector("numeric", length(LLR))
  tmp <- pmin(LLR[!Het], HmzK3)
  conf[!Het] <- lm1$coef[1]+lm1$coef[2]*tmp+lm1$coef[3]*(tmp-HmzK2)*I(tmp>HmzK2)

  tmp <- pmin(LLR[Het], HtzK3)
  conf[Het] <- lm2$coef[1]+lm2$coef[2]*tmp+lm2$coef[3]*(tmp-HtzK2)*I(tmp>HtzK2)

  conf <- matrix(conf, ncol=ncol(theCalls))

  X <- pmin(log(SNR), SNRK)
  SNRfix <- SNRlm$coef[1]+SNRlm$coef[2]*X
  conf <- sweep(conf, 2, SNRfix, FUN="+")

  conf <- 1/(1+exp(-conf))
  conf[,SNR<=3] <- 1/3
  conf[conf<1/3] <- 1/3

  return(conf)
}
