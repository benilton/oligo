replaceAffySnpParamsSingle <- function(object,value,subset){
  object$centers[subset,] <- value$centers
  object$scales[subset,] <- value$scales
  object$N[subset,] <- value$N
  return(object)
}

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


## genotypeSNP56 <- function(files, tmpdir=getwd(), batch_size=40000, balance=1.5, minLLRforCalls=c(5, 1, 5), recalibrate=TRUE, verbose=TRUE){
##   tmp <- eff2.normalizeSNP56(files, tmpdir)
##   pkgname <- cleanPlatformName(readCelHeader(files[1])$chiptype)
##   load(system.file(paste("extdata/", pkgname, "CrlmmInfo.rda", sep=""), package=pkgname))
## 
##   ## myenv should take abou 65MB RAM
##   ## I'll assume this is OK for now
##   ## and not subset it
##   myenv <- get(paste(pkgname,"Crlmm",sep="")); rm(list=paste(pkgname,"Crlmm",sep=""))
##   thePriors <- get("priors", myenv)
## 
##   Index <- which(!get("hapmapCallIndex",myenv))
## 
##   n.chunks <- length(tmp$alleleA)
##   n.files <- length(files)
##   cuts <- c(1, seq(1, n.chunks)*batch_size)
##   
##   col.Index <- cut(Index, cuts, include.lowest=TRUE, labels=FALSE)
##   grps.Index <- split(Index, col.Index)
##   n.snps <- dbGetQuery(db(get(pkgname)), "SELECT row_count FROM table_info WHERE tbl='featureSet'")[[1]]
##   grps.rows <- split(1:n.snps, rep(1:n.snps, each=batch_size, length.out=n.snps))
## 
##   ## on Chr X, things might be misplaced
##   sql.tmp <- "SELECT chrom FROM featureSet WHERE man_fsetid LIKE 'SNP%' ORDER BY man_fsetid"
##   snps.chrX <- dbGetQuery(db(get(pkgname)), sql.tmp)
##   snps.chrX <- which(snps.chrX[["chrom"]] == "X")
##   rm(sql.tmp)
##   grps.snpsX <- split(snps.chrX, cut(snps.chrX, cuts, include.lowest=TRUE, labels=FALSE))
##   
##   calls.file <- gzfile(file.path(tmpdir, "crlmm-calls.txt.gz"), "w")
##   llr.file <- gzfile(file.path(tmpdir, "crlmm-llr.txt.gz"), "w")
##   conf.file <- gzfile(file.path(tmpdir, "crlmm-conf.txt.gz"), "w")
##   if (verbose){
##     txt <- sprintf("Genotyping: %06.2f percent done.", 0)
##     cat(txt)
##     del <- paste(rep("\b", nchar(txt)), collapse="", sep="")
##   }
##   
##   for (i in as.integer(names(grps.Index))){
##     index <- match(grps.Index[[i]], grps.rows[[i]])
##     alleleA <- matrix(readBin(tmp$alleleA[i], numeric(), batch_size*n.files),
##                       nrow=batch_size)
##     alleleB <- matrix(readBin(tmp$alleleB[i], numeric(), batch_size*n.files),
##                       nrow=batch_size)
##     fs <- matrix(readBin(tmp$fs[i], numeric(), batch_size*n.files),
##                  nrow=batch_size)
##     initialCalls <- matrix(readBin(tmp$initialCalls[i], integer(), batch_size*n.files),
##                            nrow=batch_size)
##     filter <- is.na(alleleA[,1])
##     if (any(filter)){
##       alleleA <- alleleA[!filter,]
##       alleleB <- alleleB[!filter,]
##       fs <- fs[!filter,]
##       initialCalls <- initialCalls[!filter,]
##     }
##     rm(filter)
## 
##     initialCalls[-index,] <- NA
## 
##     rparams <- getGenotypeRegionParams(alleleA[index,]-alleleB[index,],
##                                        initialCalls[index,],
##                                        fs[index,],
##                                        verbose=FALSE)
##     rparams <- updateAffySnpParamsSingle(rparams, thePriors, verbose=FALSE)
##     params <- get("params", myenv)
##     params$centers <- params$centers[grps.rows[[i]],]
##     params$scales <- params$scales[grps.rows[[i]],]
##     params$N <- params$N[grps.rows[[i]],]
##     params  <- replaceAffySnpParamsSingle(params, rparams, index)
## 
##     myDist <- getAffySnpDistanceSingle56(alleleA-alleleB, params, fs)
##     myDist[,,-2] <- balance*myDist[,,-2]
##     XIndex <- match(grps.snpsX[[i]], grps.rows[[i]])
##     maleIndex <- rep(FALSE, length(files))
##     initialCalls <- getAffySnpCalls56(myDist, XIndex, maleIndex, verbose=FALSE)
##     LLR <- getAffySnpConfidence56(myDist, initialCalls, XIndex, maleIndex, verbose=FALSE)
## 
##     if (recalibrate){
##       for(k in 1:3)
##         initialCalls[ initialCalls == k & LLR < minLLRforCalls[k]] <- NA
## 
##       rparams <- getGenotypeRegionParams(alleleA-alleleB,
##                                          initialCalls,
##                                          fs,
##                                          verbose=FALSE)
##       rparams <- updateAffySnpParamsSingle(rparams, thePriors)
##       myDist <- getAffySnpDistanceSingle56(alleleA-alleleB, rparams, fs)
##       myDist[,,-2] <- balance*myDist[,,-2]
##       initialCalls <- getAffySnpCalls56(myDist, XIndex, maleIndex, verbose=FALSE)
##       LLR <- getAffySnpConfidence56(myDist, initialCalls, XIndex, maleIndex, verbose=FALSE)
##       rm(myDist)
##     }
##     callsConfidence <- LLR2conf(initialCalls, LLR, readBin(tmp$snr, numeric(), n.files), "pd.genomewidesnp.6")
## 
##     writeLines(apply(initialCalls, 1, paste, collapse=","), calls.file)
##     writeLines(apply(LLR, 1, paste, collapse=","), llr.file)
##     writeLines(apply(callsConfidence, 1, paste, collapse=","), conf.file)
##     
##     if (verbose){
##       cat(del)
##       txt <- sprintf("Genotyping: %06.2f percent done.", i/length(names(grps.Index))*100)
##       cat(txt)
##       del <- paste(rep("\b", nchar(txt)), collapse="", sep="")
##     }
## 
##   }
##   close(calls.file)
##   close(llr.file)
##   close(conf.file)
## }

## LLR2conf56 <-function(theCalls, theLLR, SNR, annot, tmpdir=getwd()){
##   load(paste(system.file("extdata", package=annot), "/", annot, ".spline.params.rda", sep=""))
##   conf <- createBufferedMatrix(nrow(theCalls), 0, nrow(theCalls), prefix="final.conf.", directory=tmpdir)
##   X <- pmin(log(SNR), SNRK)
##   SNRfix <- SNRlm$coef[1]+SNRlm$coef[2]*X
##   for (i in 1:ncol(theCalls)){
##     AddColumn(conf)
##     Het <- as.logical(theCalls[,i] == 2)
##     LLR <- as.numeric(sqrt(theLLR[,i]))
##     tmp <- pmin(LLR[!Het], HmzK3)
##     conf[!Het, i] <- lm1$coef[1]+lm1$coef[2]*tmp+lm1$coef[3]*(tmp-HmzK2)*(tmp>HmzK2)
##     tmp <- pmin(LLR[Het], HtzK3)
##     conf[Het, i] <- lm2$coef[1]+lm2$coef[2]*tmp+lm2$coef[3]*(tmp-HtzK2)*(tmp>HtzK2)
##     conf[, i] <- conf[, i]+SNRfix[i]
##     conf[, i] <- 1/(1+exp(-conf[,i]))
##     conf[conf[,i] < 1/3, i] <- 1/3
##   }
##   return(conf)
## }

fitAffySnpMixture56 <- function(object, df1=3, df2=5,
                              probs=rep(1/3,3), eps=50,
                              subSampleSize=10^4, verbose=TRUE){

  ## object: a matrix (SNPs x Allele)

  ## FIXME: det gender
  ## maleIndex: logical for male
  maleIndex <- FALSE

  ## XIndex: vector of integers pointing to SNPs on chr X
  ## XIndex=getChrXIndex(object)
  
  I <- nrow(object)
  set.seed(1)
  
##  tmp <- c( (1:I)[-XIndex],((I+1):(2*I))[-XIndex])
##  tmp <- c( (1:I),((I+1):(2*I)))
  
  idx <- sort(sample(I, subSampleSize))

  pis <- array(0,dim=c(I,3))
  fs <- array(0,dim=I)
  snr <- NULL

  if(verbose) cat("Fitting mixture model. Epsilon must reach ", eps, ".\n",sep="")

  Y <- object[,1]-object[,2]
  A <- rowMeans(object)
  
  mus <- quantile(Y,c(1,3,5)/6);mus[2]=0
  sigmas <- rep(mad(c(Y[Y<mus[1]]-mus[1],Y[Y>mus[3]]-mus[3])),3)
  sigmas[2] <- sigmas[2]/2

  a <- A[idx]
  y <- Y[idx]
  A <- A-mean(A)
  a <- a-mean(a)
    
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
      if (itmax > 1) cat(del)
      message <- paste("Epsilon = ", signif(change,2), "  ", sep="")
      del <- paste(rep("\b", nchar(message)), collapse="")
      cat(message)
    }
    
    PreviousLogLik <- LogLik
    probs <- colMeans(z)
    
    ## M
    fit1 <- lm(y~matA,weights=z[,1])
    fit2 <- sum(z[,2]*y)/sum(z[,2])
    fit3 <- lm(y~matA,weights=z[,3])
    
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
  bigX <- cbind(1, ns(A, knots=as.numeric(attr(matA, "knots")), Boundary.knots=attr(matA, "Boundary.knots")))
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

  fs <- as.numeric((pred3-pred1)/2)
  for(k in 1:3){
    pis[, k] <- z[,(4-k)] ##4-k cause 3is1,2is2 and 1is3
  }
  snr <- median(fs)^2/(sigmas[1]^2+sigmas[2]^2)

  fs[fs < 0] <- 0
  if(verbose) cat("Done.\n")
##  return(list(f0=median(fs),fs=fs, pis=pis, snr=snr))
  return(list(f0=median(fs), fs=fs, initial=apply(pis, 1, which.max), snr=snr))
}

## eff2.normalizeSNP56 <- function(celFiles, destDir, batch_size=40000, verbose=TRUE){
##   destDir <- gsub("/$", "", destDir)
##   
##   ## Check existence of directory (destination)
##   ## The destDir should not exist (yet)
##   if (file.exists(destDir)) stop(message(destDir, " exists."))
##   dir.create(destDir)
## 
##   ## Make sure all files are of the same type
##   if (length(chiptype <- unique(sapply(celFiles, function(x) readCelHeader(x)$chiptype))) > 1)
##     stop("CEL files do not have the same chip type")
## 
##   ## Determine pkg and number of files to read at once
##   pkgname <- cleanPlatformName(chiptype)
##   stopifnot(require(pkgname, character.only=TRUE))
##   conn <- db(get(pkgname))
## 
##   n.snps <- dbGetQuery(conn, "SELECT row_count FROM table_info WHERE tbl='featureSet'")[[1]]
## 
##   ## Get feature IDs and load reference
##   load(system.file("extdata", paste(pkgname, "Ref.rda", sep=""), package=pkgname))
##   reference <- sort(reference)
## 
##   tmp <- dbGetQuery(conn, paste("SELECT fid, man_fsetid, allele FROM pmfeature, featureSet",
##                                 "WHERE pmfeature.fsetid = featureSet.fsetid"))
##   tmp <- tmp[order(tmp$fid),]
##   rownames(tmp) <- NULL
## 
##   pnVec <- paste(tmp[["man_fsetid"]],
##                  c("A", "B")[tmp[["allele"]]+1],
##                  sep="")
## 
##   idx <- order(pnVec)
##   pnVec <- pnVec[idx]
##   tmp[["man_fsetid"]] <- tmp[["allele"]] <- tmp[["strand"]] <- NULL
##   ngenes <- length(unique(pnVec))
##   bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}
## 
##   if (verbose){
##     txt <- sprintf("Normalization: %06.2f percent done.", 0)
##     cat(txt)
##     del <- paste(rep("\b", nchar(txt)), sep="", collapse="")
##   }
## 
##   ## how many columns using Ben's idea?
##   extra.nas <- batch_size - (n.snps %% batch_size)
##   n.cols <- n.snps %/% batch_size + (extra.nas > 0)
##   n.rows <- batch_size*length(celFiles)
##   extra.nas <- rep(NA, extra.nas)
## 
##   alleleA.files <- file.path(destDir, paste("crlmm-alleleA-", 1:n.cols, ".bin", sep=""))
##   alleleB.files <- file.path(destDir, paste("crlmm-alleleB-", 1:n.cols, ".bin", sep=""))
##   initialCalls.files <- file.path(destDir, paste("crlmm-initialCalls-", 1:n.cols, ".bin", sep=""))
##   llr.files <- file.path(destDir, paste("crlmm-llr0-", 1:n.cols, ".bin", sep=""))
##   fs.files <- file.path(destDir, paste("crlmm-fs-", 1:n.cols, ".bin", sep=""))
##   f0.file <- file.path(destDir, "crlmm-f0.bin")
##   snr.file <- file.path(destDir, "crlmm-snr.bin")
## 
##   conn.alleleA <- lapply(alleleA.files, file, "wb")
##   conn.alleleB <- lapply(alleleB.files, file, "wb")
##   conn.calls <- lapply(initialCalls.files, file, "wb")
##   conn.fs <- lapply(fs.files, file, "wb")
##   conn.f0 <- file(f0.file, "wb")
##   conn.snr <- file(snr.file, "wb")
##   
##   F0 <- SNR <- NULL
##   for (i in 1:length(celFiles)){
##     pms <- normalize.quantiles.use.target(readCelIntensities(celFiles[i],
##                                                              indices=tmp$fid),
##                                           reference, copy=FALSE)[idx,, drop=FALSE]
##     theSumm <- matrix(.Call("rma_c_complete_copy", pms, pms,
##                             pnVec, ngenes,  body(bg.dens),
##                             new.env(), FALSE, FALSE,
##                             as.integer(2), PACKAGE="oligo")[,1],
##                       ncol=2, byrow=TRUE)
##     rm(pms)
##     correction <- fitAffySnpMixture56(theSumm, verbose=FALSE)
## 
##     all.A <- matrix(c(theSumm[,1], extra.nas),
##                     nrow=batch_size)
##     all.B <- matrix(c(theSumm[,2], extra.nas),
##                     nrow=batch_size)
##     init.calls <- matrix(c(correction$initial, extra.nas),
##                          nrow=batch_size)
##     correction.fs <- matrix(c(correction$fs, extra.nas),
##                             nrow=batch_size)
## 
##     for (j in 1:n.cols){
##       writeBin(all.A[,j], conn.alleleA[[j]], size=8)
##       writeBin(all.B[,j], conn.alleleB[[j]], size=8)
##       writeBin(init.calls[,j], conn.calls[[j]])
##       writeBin(correction.fs[,j], conn.fs[[j]], size=8)
##       writeBin(correction$f0, conn.f0, size=8)
##       writeBin(correction$snr, conn.snr, size=8)
##     }
##     
##     rm(correction)
##     if (verbose){
##       cat(del)
##       txt <- sprintf("Normalization: %06.2f percent done.", i/length(celFiles)*100)
##       cat(txt)
##     }
##   }
##   if (verbose) cat("\n")
##   for (i in 1:n.cols){
##     close(conn.alleleA[[i]])
##     close(conn.alleleB[[i]])
##     close(conn.calls[[i]])
##     close(conn.fs[[i]])
##   }
##   close(conn.f0)
##   close(conn.snr)
##   return(list(alleleA=alleleA.files,
##               alleleB=alleleB.files,
##               initialCalls=initialCalls.files,
##               fs=fs.files,
##               f0=f0.file,
##               snr=snr.file))
## }


getAffySnpConfidence56 <- function(Dist, Calls, XIndex, maleIndex, verbose=TRUE){
####  XIndex <- which(subset%in%XIndex)

  res <- array(NA,dim=dim(Dist)[c(1,2)])
#####  dimnames(res) <- list(dimnames(Dist)[[1]],dimnames(Dist)[[2]])
#####  Dist <- rowSums(Dist, dims=3, na.rm=T)
  
  if (verbose) cat("Computing confidence for calls on ", ncol(res), " arrays")
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
  if (verbose) cat("Done\n")  
  return(res)
}

getAffySnpCalls56 <- function(Dist, XIndex, maleIndex, verbose=FALSE){
  ##  XIndex <- which(subset%in%XIndex)
  res <- array(as.integer(-1),dim=dim(Dist)[c(1,2)])
  ##  Dist[XIndex, maleIndex, 2] <- Inf
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

getAffySnpDistanceSingle56 <- function(x, params, f=0, subset=1:(dim(x)[1]),
                                       w=NULL, verbose=FALSE){
  x <- x[subset,, drop=FALSE]
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

genotypeOne <- function(files, tmpdir=getwd(), batch_size=40000, balance=1.5, minLLRforCalls=c(5, 1, 5), recalibrate=TRUE, verbose=TRUE, pkgname){
  if (!file.exists(tmpdir)){
    tmp <- normalizeOne(files, tmpdir, pkgname=pkgname)
    if(missing(pkgname))
      pkgname <- cleanPlatformName(readCelHeader(files[1])$chiptype)
  }else{
    message("Using previous results stored at ", tmpdir)
    analysis <- read.table(file.path(tmpdir, "analysis.txt"), stringsAsFactors=FALSE)
    if(missing(pkgname))
      pkgname <- analysis[which(analysis[,1]=="annotation"), 2]
    batch_size <- as.integer(analysis[which(analysis[,1]=="batch_size"), 2])
    tmp <- list(alleleA=dir(tmpdir, pattern="alleleA", full.names=TRUE),
                alleleB=dir(tmpdir, pattern="alleleB", full.names=TRUE),
                initialCalls=dir(tmpdir, pattern="initialCalls", full.names=TRUE),
                fs=dir(tmpdir, pattern="fs", full.names=TRUE),
                f0=dir(tmpdir, pattern="f0", full.names=TRUE),
                snr=dir(tmpdir, pattern="snr", full.names=TRUE))
    stopifnot(require(pkgname, character.only=TRUE))
    ## MAYBE MORE TESTS OF VALIDITY
  }

  load(system.file(paste("extdata/", pkgname, "CrlmmInfo.rda", sep=""), package=pkgname))

  ## myenv should take abou 65MB RAM
  ## I'll assume this is OK for now
  ## and not subset it
  myenv <- get(paste(pkgname,"Crlmm",sep="")); rm(list=paste(pkgname,"Crlmm",sep=""))
  thePriors <- get("priors", myenv)

  tmpdf <- dbGetQuery(db(get(pkgname)), "SELECT man_fsetid, chrom, physical_pos FROM featureSet WHERE man_fsetid LIKE 'SNP%'")
  tmpdf[is.na(tmpdf$chrom), "chrom"] <- 0
  tmpdf[is.na(tmpdf$physical_pos), "physical_pos"] <- 0
  tmpdf <- tmpdf[order(tmpdf$man_fsetid),]
  tmpdf[["index"]] <- 1:nrow(tmpdf)
  tmpdf <- tmpdf[order(tmpdf$chrom, tmpdf$physical_pos, tmpdf$man_fsetid),]
  tmpidx <- tmpdf[["index"]]
##   rm(tmpdf); gc()
  
##   myenv$params$centers <- myenv$params$centers[tmpidx,]
##   myenv$params$scales <- myenv$params$scales[tmpidx,]
##   myenv$params$N <- myenv$params$N[tmpidx,]
##   myenv$hapmapCallIndex <- myenv$hapmapCallIndex[tmpidx]

  Index <- which(!get("hapmapCallIndex", myenv))
  
  analysis <- read.table(file.path(tmpdir, "analysis.txt"), stringsAsFactors=FALSE)
  breaks <- cumsum(c(0, as.integer(analysis[-(1:3),2])))
  index.grps <- cut(Index, breaks, include.lowest=TRUE, labels=FALSE)
  Index <- Index-breaks[index.grps]
  index.grps <- split(Index, index.grps)
  rm(breaks)
  
  n.chunks <- nrow(analysis)-3
  n.files <- as.integer(analysis[which(analysis[,1] == "nsamples"), 2])
  
##   calls.file <- gzfile(file.path(tmpdir, "crlmm-calls.txt.gz"), "w")
##   llr.file <- gzfile(file.path(tmpdir, "crlmm-llr.txt.gz"), "w")
##   conf.file <- gzfile(file.path(tmpdir, "crlmm-conf.txt.gz"), "w")
  if (verbose){
    txt <- sprintf("Genotyping: %06.2f percent done.", 0)
    cat(txt)
    del <- paste(rep("\b", nchar(txt)), collapse="", sep="")
  }
  last <- 0
  for (i in 1:n.chunks){
    index <- index.grps[[i]]
    nrows <- as.integer(analysis[i+3,2])
    overall_pos <- (last+1):(last+nrows)
    last <- max(overall_pos)
    
    alleleA <- matrix(readBin(tmp$alleleA[i], numeric(),
                              nrows*n.files), nrow=nrows)
    alleleB <- matrix(readBin(tmp$alleleB[i], numeric(),
                              nrows*n.files), nrow=nrows)
    fs <- matrix(readBin(tmp$fs[i], numeric(), nrows*n.files),
                         nrow=nrows)
    initialCalls <- matrix(readBin(tmp$initialCalls[i], integer(),
                                   nrows*n.files), nrow=nrows)
    initialCalls[-index,] <- NA
    rparams <- getGenotypeRegionParams(alleleA[index,]-alleleB[index,],
                                       initialCalls[index,],
                                       fs[index,], verbose=FALSE)
    rparams <- updateAffySnpParamsSingle(rparams, thePriors, verbose=FALSE)
    params <- get("params", myenv)
    params$centers <- params$centers[overall_pos,]
    params$scales <- params$scales[overall_pos,]
    params$N <- params$N[overall_pos,]
    params  <- replaceAffySnpParamsSingle(params, rparams, index)

    myDist <- getAffySnpDistanceSingle56(alleleA-alleleB, params, fs)
    ##SAVE THE ABOVE
    myDist[,,-2] <- balance*myDist[,,-2]
    XIndex <- integer()
    if (length(grep("chrX", tmp$initialCalls[i])) > 0)
      XIndex <- 1:nrows

    maleIndex <- rep(FALSE, n.files)
    initialCalls <- getAffySnpCalls56(myDist, XIndex, maleIndex, verbose=FALSE)
    LLR <- getAffySnpConfidence56(myDist, initialCalls, XIndex, maleIndex, verbose=FALSE)

    if (recalibrate){
      for(k in 1:3)
        initialCalls[ initialCalls == k & LLR < minLLRforCalls[k]] <- NA

      rparams <- getGenotypeRegionParams(alleleA-alleleB,
                                         initialCalls,
                                         fs, verbose=FALSE)
      rparams <- updateAffySnpParamsSingle(rparams, thePriors)
      myDist <- getAffySnpDistanceSingle56(alleleA-alleleB, rparams, fs)
      ### SAVE THE ABOVE
      myDist[,,-2] <- balance*myDist[,,-2]
      initialCalls <- getAffySnpCalls56(myDist, XIndex, maleIndex, verbose=FALSE)
      LLR <- getAffySnpConfidence56(myDist, initialCalls, XIndex, maleIndex, verbose=FALSE)
      rm(myDist)
    }
    callsConfidence <- LLR2conf(initialCalls, LLR, readBin(tmp$snr, numeric(), n.files), pkgname)

    rownames(initialCalls) <- rownames(LLR) <- rownames(callsConfidence) <- tmpdf[["man_fsetid"]][overall_pos]
    colnames(initialCalls) <- colnames(LLR) <- colnames(callsConfidence) <- basename(files)

    write.table(initialCalls, file.path(tmpdir, "crlmm-calls.txt"), append=TRUE, quote=FALSE, sep="\t", col.names=(i==1))
    write.table(LLR, file.path(tmpdir, "crlmm-llr.txt"), append=TRUE, quote=FALSE, sep="\t", col.names=(i==1))
    write.table(callsConfidence, file.path(tmpdir, "crlmm-conf.txt"), append=TRUE, quote=FALSE, sep="\t", col.names=(i==1))
    
##    writeLines(apply(initialCalls, 1, paste, collapse=","), calls.file)
##    writeLines(apply(LLR, 1, paste, collapse=","), llr.file)
##    writeLines(apply(callsConfidence, 1, paste, collapse=","), conf.file)
    
    if (verbose){
      cat(del)
      txt <- sprintf("Genotyping: %06.2f percent done.", i/n.chunks*100)
      cat(txt)
      del <- paste(rep("\b", nchar(txt)), collapse="", sep="")
    }

  }
  close(calls.file)
  close(llr.file)
  close(conf.file)
}

normalizeOne <- function(celFiles, destDir, batch_size=40000, verbose=TRUE, pkgname, reference){
  
  ## Check existence of directory (destination)
  ## The destDir should not exist (yet)
  destDir <- gsub("/$", "", destDir)
  if (file.exists(destDir)) stop(message(destDir, " exists."))
  dir.create(destDir)

  ## Make sure all files are of the same type
  if (length(chiptype <- unique(sapply(celFiles, function(x) readCelHeader(x)$chiptype))) > 1)
    stop("CEL files do not have the same chip type")

  ## Determine pkg and number of files to read at once
  if (missing(pkgname))
    pkgname <- cleanPlatformName(chiptype)
  stopifnot(require(pkgname, character.only=TRUE))
  conn <- db(get(pkgname))

  n.snps <- dbGetQuery(conn, "SELECT row_count FROM table_info WHERE tbl='featureSet'")[[1]]

  ## Get feature IDs and load reference
  if (missing(reference))
    load(system.file("extdata", paste(pkgname, "Ref.rda", sep=""), package=pkgname))
  reference <- sort(reference)

  tmp <- dbGetQuery(conn, paste("SELECT fid, man_fsetid, allele, featureSet.chrom, featureSet.physical_pos",
                                "FROM pmfeature, featureSet",
                                "WHERE pmfeature.fsetid = featureSet.fsetid"))
  tmp[is.na(tmp$chrom), "chrom"] <- 0
  tmp[is.na(tmp$physical_pos), "physical_pos"] <- 0
  tmp <- tmp[order(tmp$chrom, tmp$physical_pos, tmp$man_fsetid, tmp$allele),]
  rownames(tmp) <- NULL

  pnVec <- paste(tmp[["man_fsetid"]],
                 c("A", "B")[tmp[["allele"]]+1],
                 sep=":")

  tmp[["man_fsetid"]] <- tmp[["allele"]] <- tmp[["strand"]] <- NULL
  ngenes <- length(unique(pnVec))
## the one below should be stored somewhere  
##  snpnames <- unique(gsub("\\:[AB]$", "", pnVec))

  info <- dbGetQuery(conn, "SELECT man_fsetid, chrom, physical_pos FROM featureSet WHERE man_fsetid LIKE 'SNP%'")
  info[is.na(info[["chrom"]]), "chrom"] <- 0
  info[is.na(info[["physical_pos"]]), "physical_pos"] <- 0
  info <- info[order(info$chrom, info$physical_pos, info$man_fsetid),]

  split.rowids <- lapply(split(1:nrow(info), info[["chrom"]]),
                         function(x)
                         split(x, rep(1:length(x),
                                      each=batch_size,
                                      length.out=length(x))))
  rm(info); gc()
  split.rowids <- unlist(split.rowids, recursive=FALSE)
  names(split.rowids) <- paste("chr", names(split.rowids), sep="")
  
  ## This made me waste 2 days looking for a problem
  ## if factors, the order is 1, 10, ... instead of 1, 2, 3...
###  names(split.rowids) <- paste("chr",
###                               gsub("^([[:digit:]])\\.", "0\\1\\.",
###                                    gsub("\\.([[:digit:]])$", "\\.0\\1",
###                                         names(split.rowids))),
###                               sep="")
  filenames <- names(split.rowids)
###   idx <- order(filenames)
###   split.rowids <- split.rowids[idx]
###   filenames <- filenames[idx]
  
  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}

  if (verbose){
    txt <- sprintf("Normalization: %06.2f percent done.", 0)
    cat(txt)
    del <- paste(rep("\b", nchar(txt)), sep="", collapse="")
  }

  analysis <- data.frame(fields=c("annotation", "nsamples", "batch_size", filenames),
             values=c(pkgname, length(celFiles), batch_size,
               sapply(split.rowids, length)))
  write.table(analysis, file.path(destDir, "analysis.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  alleleA.files <- file.path(destDir, paste("alleleA-", filenames, sep=""))
  alleleB.files <- file.path(destDir, paste("alleleB-", filenames, sep=""))
  initialCalls.files <- file.path(destDir, paste("initialCalls-", filenames, sep=""))
  llr.files <- file.path(destDir, paste("llr0-", filenames, sep=""))
  fs.files <- file.path(destDir, paste("fs-", filenames, sep=""))
  f0.file <- file.path(destDir, "f0")
  snr.file <- file.path(destDir, "snr")

  for (i in 1:length(celFiles)){
    pms <- normalize.quantiles.use.target(readCelIntensities(celFiles[i],
                                                             indices=tmp$fid),
                                          reference, copy=FALSE)
    theSumm <- matrix(.Call("rma_c_complete_copy", pms, pms,
                            pnVec, ngenes,  body(bg.dens),
                            new.env(), FALSE, FALSE,
                            as.integer(2), PACKAGE="oligo")[,1],
                      ncol=2, byrow=TRUE)
    rm(pms)
    correction <- fitAffySnpMixture56(theSumm, verbose=FALSE)

    for (j in 1:length(filenames)){
      conn <- file(alleleA.files[j], "ab")
      writeBin(theSumm[split.rowids[[j]], 1], conn)
      close(conn)
      conn <- file(alleleB.files[j], "ab")
      writeBin(theSumm[split.rowids[[j]], 2], conn)
      close(conn)
      conn <- file(initialCalls.files[j], "ab")
      writeBin(correction$initial[split.rowids[[j]]], conn)
      close(conn)
      conn <- file(fs.files[j], "ab")
      writeBin(correction$fs[split.rowids[[j]]], conn)
      close(conn)
    }
    conn <- file(f0.file, "ab")
    writeBin(correction$f0, conn)
    close(conn)
    conn <- file(snr.file, "ab")
    writeBin(correction$snr, conn)
    close(conn)

    rm(correction)
    if (verbose){
      cat(del)
      txt <- sprintf("Normalization: %06.2f percent done.", i/length(celFiles)*100)
      cat(txt)
    }
  }
  if (verbose) cat("\n")
  return(list(alleleA=alleleA.files,
              alleleB=alleleB.files,
              initialCalls=initialCalls.files,
              fs=fs.files,
              f0=f0.file,
              snr=snr.file))
}
