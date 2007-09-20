updateAffySnpParamsSingle <- function(object, priors,  minN=3,
                                maxHomoSigma=priors$maxsigma[1],
                                maxHeteSigma=priors$maxsigma[2],
                                subset=1:(dim(object$centers)[1]),
                                d0s=80, verbose=FALSE){
  object$centers <- object$centers[subset,]
  object$scales <- object$scales[subset,]
  object$N <- object$N[subset,]
  if(verbose) cat("Updating centers and scales")

  ##First variances
  for(j in 1:2){ ##1 and 3 are the same
    if(j==2) N <- object$N[,2] else N <- rowSums(object$N[,-2],na.rm=TRUE)
    s <- object$scales[, j]
    if (is.null(d0s))
      d0s <- priors$d0s[j]
    s20 <- priors$s20[j] ##notice the ad-hoc choice of 3
    Index <- N>minN & !is.na(s)
    N <- N[Index]; s <- s[Index]
    object$scales[Index, j] <- sqrt (  ( (N-1)*s^2 + d0s*s20 ) / (d0s+N-1) )
    object$scales[!Index, j] <- sqrt(s20)
  }
  N <- object$N
  object$scales[,3] <- object$scales[,1] ##AA=BB 
  object$scales[,2][object$scales[,2]>maxHeteSigma] <- maxHeteSigma
  object$scales[,-2][object$scales[,-2]>maxHomoSigma] <- maxHomoSigma
  if(verbose) cat(".")

  ##Means
  Vinv <- solve(priors$V)
  NSinv <- t(N)/priors$s20
  tmp <- t(sapply(1:nrow(object$centers),function(i){
    if(verbose & i%%5000==0)  cat(".")
    mus <- object$centers[i,]
    Ns <- N[i,]
    mus[Ns<minN] <- 0
    mus[is.na(mus)] <- 0
    return(solve(Vinv+diag(NSinv[,i]))%*%(NSinv[,i]*mus))
  }))
  object$centers <- tmp
  if(verbose) cat("Done.\n")
  return(object)
}


getAffySnpDistanceSingle <- function(object,params,f=0,subset=1:(dim(object)[1]),
                                     w=NULL,verbose=FALSE){
  x=getM(object[subset,])
  Dist <- array(NA,dim=c(dim(x)[1],dim(x)[2], 3))
  if(verbose) cat("Calculating likelihood-based distances")
  for(j in 1:3){
    tmp <- x+(j-2)*f[subset,]
    Dist[,,j] <- 2*log(params$scales[subset,j]) +
      ((tmp-params$centers[subset,j])/params$scales[subset,j])^2
    if(!is.null(w)) Dist[,,j] <-  Dist[,,j] - 2*log(w[subset,,j])
  }
  if(verbose) cat("Done.\n")
  return(Dist)
}
