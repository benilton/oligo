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


readSummaries <- function(type, tmpdir){
  
  validOptions <- c("alleleA", "alleleB", "fs", "alleleA-sense",
                    "alleleA-antisense", "alleleB-sense",
                    "alleleB-antisense", "calls", "conf", "llr",
                    "snr", "initialCalls", "antisense-f", "sense-f")
  if (!(type %in% validOptions))
    stop(paste("Invalid 'type' argument.\nValid options are '",
               paste(validOptions, collapse="', '"), "'.", sep=""))

  analysis <- read.table(file.path(tmpdir, "analysis.txt"), stringsAsFactors=FALSE)
  pkgname <- as.character(analysis[1, 2])

  if (!(pkgname %in% c('pd.genomewidesnp.5', 'pd.genomewidesnp.6'))){
      if (type %in% c("alleleA", "alleleB", "fs", "initialCalls")){
          if (type == "initialCalls"){
              dataType <- integer()
          }else{
              dataType <- numeric()
          }
          files <- file.path(tmpdir, paste(type, analysis[-(1:3), 1], sep="-"))
          n.files <- length(files)
          n.samples <- as.integer(analysis[match("nsamples", analysis[,1]), 2])
          tmp <- NULL
          for (i in 1:n.files){
              nrows <- as.integer(analysis[3+i, 2])
              tmp <- rbind(tmp, matrix(readBin(files[i], dataType, nrows*n.samples), nrow=nrows))
          }
          pkgname <- analysis[1, 2]
          requireAnnotation(pkgname)
          tmpdf <- getSnpLocInfo(pkgname)
          rownames(tmp) <- tmpdf[["man_fsetid"]]
          colnames(tmp) <- as.character(read.table(file.path(tmpdir, "crlmm-calls.txt"), nrows=1, colClasses=rep("character", ncol(tmp))))
      }else if (type %in% c("calls", "conf", "llr", "alleleA-sense", "alleleA-antisense", "alleleB-sense", "alleleB-antisense", "antisense-f", "sense-f")){
          target <- file.path(tmpdir,
                              ifelse(type %in% c("calls", "conf", "llr"),
                                     paste("crlmm-", type, ".txt", sep=""),
                                     paste(type, ".txt", sep="")))
          header <- as.character(read.delim(target, nrow=1, colClasses="character", header=FALSE))
          nsamples <- length(header)
          what <- c("character", rep(ifelse(type == "calls", "integer", "numeric"), nsamples))
          tmp <- as.matrix(read.delim(target, colClasses=what, skip=1, row.names=1, header=FALSE))
          colnames(tmp) <- header
      }else if (type == "snr"){
          header <- as.character(read.delim(file.path(tmpdir, "crlmm-calls.txt"),
                                            nrow=1, colClasses="character", header=FALSE))
          nsamples <- length(header)
          tmp <- readBin(file.path(tmpdir, "snr"), numeric(), nsamples)
          names(tmp) <- header
      }
  }else{
      vo <- c('alleleA', 'alleleB', 'fs', 'f0', 'initialCalls', 'calls', 'conf', 'llr', 'snr')
      type <- match.arg(type, vo)
      obj <- load(file.path(tmpdir, paste(type, '.RData', sep='')))
      tmp <- get(obj)[]
  }
  return(tmp)
}

getCrlmmSummaries <- function(tmpdir){
  analysis <- read.delim(file.path(tmpdir, "analysis.txt"), colClasses="character", header=FALSE)
  annotation <- analysis[match("annotation", analysis[,1]), 2]
  snpcnv <- length(grep("genomewidesnp", annotation)) > 0
  if (snpcnv){
    summaries <- c("alleleA", "alleleB")
    tmp <- new("SnpSuperSet",
               call=readSummaries("calls", tmpdir),
               callProbability=readSummaries("conf", tmpdir),
               alleleA=readSummaries("alleleA", tmpdir),
               alleleB=readSummaries("alleleB", tmpdir),
               F=readSummaries("fs", tmpdir))
  }else{
    tmp <- new("SnpSuperSet",
               call=readSummaries("calls", tmpdir),
               callProbability=readSummaries("conf", tmpdir),
               senseAlleleA=readSummaries("alleleA-sense", tmpdir),
               senseAlleleB=readSummaries("alleleB-sense", tmpdir),
               antisenseAlleleA=readSummaries("alleleA-antisense", tmpdir),
               antisenseAlleleB=readSummaries("alleleB-antisense", tmpdir),
               antisenseF=readSummaries("antisense-f", tmpdir),
               senseF=readSummaries("sense-f", tmpdir))
  }
  annotation(tmp) <- annotation
  phenoData(tmp) <- new("AnnotatedDataFrame",
                        data=data.frame(crlmmSNR=readSummaries("snr", tmpdir)),
                        varMetadata=data.frame(labelDescription="Signal-to-Noise Ratio by CRLMM",
                          row.names="crlmmSNR"))
  return(tmp)
}

########################
### CLEANUP 11/11/08
########################

## readChipTypesFromCels <- function(celFiles)
##   sapply(celFiles, function(x) readCelHeader(x)[["chiptype"]])

readChipTypesFromCels <- function(celFiles)
    sapply(celFiles, function(x) read.celfile.header(x)[['cdfName']])

normalizeOne <- function(celFiles, destDir, batch_size=40000, pkgname,
                         reference=TRUE, check=TRUE, verbose=TRUE){
    
  ## Check existence of directory (destination)
  ## The destDir must not exist (yet)
  ## If it exists, abort
  if (file.exists(destDir))
      stop(message(destDir, " already exists. Use the name of a non-existing directoty."))
  destDir <- gsub("/$", "", destDir)
  dir.create(destDir)

  ## Make sure all files are of the same type
  if (check)
    if (length(chiptype <- unique(readChipTypesFromCels(celFiles))) > 1)
      stop("CEL files do not have the same chip type")

  ## Determine pkg
  if (missing(pkgname))
    pkgname <- cleanPlatformName(chiptype)

  ## Get annotation package
  stopifnot(requireAnnotation(pkgname))
  conn <- db(get(pkgname))

  n.snps <- dbGetQuery(conn, "SELECT row_count FROM table_info WHERE tbl='featureSet'")[[1]]

  ## Get all PM probes
  ## Some are not annotated (Affymetrix files) and will get a 0 for the
  ## moment - need to sort later
  ## Order by SNP name
  tmp <- getProbeInfo(pkgname)
  pnVec <- paste(tmp[["man_fsetid"]],
                 c("A", "B")[tmp[["allele"]]+1],
                 sep=":")

  tmp[["allele"]] <- NULL
  nSnpAlleleCombinations <- length(unique(pnVec))

  ## TODO: isn't info equal to tmp without the probe info?
  ## Info: SNP information
  snps <- getSnpNames(pkgname)
  samples <- basename(celFiles)
  save(snps, samples, file=file.path(destDir, "RowsAndColumns.rda"))
  
  analysis <- data.frame(fields=c("annotation", "nsamples", "batch_size"),
                         values=c(pkgname, length(celFiles), batch_size))
  write.table(analysis, file.path(destDir, "analysis.txt"),
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

  dims <- c(length(snps), length(celFiles))
  alleleA <- ff(vmode='double', dim=dims, pattern=file.path(destDir, 'alleleA-'))
  alleleB <- ff(vmode='double', dim=dims, pattern=file.path(destDir, 'alleleB-'))
  initialCalls <- ff(vmode='integer', dim=dims, pattern=file.path(destDir, 'initialCalls-'))
  fs <- ff(vmode='double', dim=dims, pattern=file.path(destDir, 'fs-'))
  f0 <- ff(vmode='double', dim=dims[2], pattern=file.path(destDir, 'f0-'))
  snr <- ff(vmode='double', dim=dims[2], pattern=file.path(destDir, 'snr-'))

  ## Get feature IDs and load/create reference
  if (reference){
    load(system.file("extdata", paste(pkgname, "Ref.rda", sep=""),
                     package=pkgname))
  }else{
    if (verbose){
      message("Creating normalization vector from the data.")
      message(length(tmp$fid)/2^20, "GB RAM required.")
    }
    nqdt <- normalize.quantiles.determine.target
    reference <- nqdt(readCEL(celFiles, tmp$fid))
    if (verbose) message("Normalization vector created.")
  }
  reference <- sort(reference)

  if (verbose){
    message("Normalizing and summarizing.")
    pb <- txtProgressBar(min=1, max=length(celFiles), style=3, initial=1)
  }
  
  ## just making an alias from preprocessCore
  ## to make line shorter
  n2t <- normalize.quantiles.use.target

  for (i in 1:length(celFiles)){
    ## Read one file at a time
    ## no summarization across samples (probes are replicates)
    ## take medians of probesets (using Ben's code, which is faster)
    ## Summarized data now in 2 columns (A and B)
    pms <- n2t(readCEL(celFiles[i], tmp$fid), reference, copy=FALSE)
    tmpRes <- basicRMA(pms, pnVec, background=FALSE, normalize=FALSE, verbose=FALSE)
    rm(pms)
    theSumm <- matrix(tmpRes, ncol=2, byrow=TRUE)
    rm(tmpRes)
    correction <- fitAffySnpMixture56(theSumm, verbose=FALSE)
    alleleA[,i] <- theSumm[,1]
    alleleB[,i] <- theSumm[,2]
    rm(theSumm)
    initialCalls[,i] <- correction$initial
    fs[,i] <- correction$fs
    f0[i] <- correction$f0
    snr[i] <- correction$snr
    rm(correction); gc()
    if (verbose) setTxtProgressBar(pb, i)
  }
  close(alleleA)
  close(alleleB)
  close(initialCalls)
  close(fs)
  close(f0)
  close(snr)

  filesNorm <- list(alleleA=alleleA, alleleB=alleleB,
                      initialCalls=initialCalls, fs=fs, f0=f0, snr=snr)
  
  if (verbose) close(pb)
  pkgNorm <- pkgname
  save(filesNorm, pkgNorm, file=file.path(destDir, "NormalizationSummarizationOutput.rda"))
  return(filesNorm)
}

genotypeOne <- function(files, outDir, batch_size=40000,
                        balance=1.5, minLLRforCalls=c(5, 1, 5),
                        recalibrate=TRUE, verbose=TRUE, pkgname,
                        reference=TRUE, d0s=80){
    oligoDEBUG <- getOption('oligoDEBUG', FALSE)
  if (missing(outDir)) stop("Output directory must be given.")
  if (!file.exists(outDir)){
    normOut <- normalizeOne(files, outDir, pkgname=pkgname, reference=reference)
  }else{
    if (!file.info(outDir)[["isdir"]])
      stop(outDir, " is not a valid directory.")
    outDir <- gsub("/*$", "", outDir)
    if (verbose) message("Loading results from previous normalization/summarization step.")
    obj <- load(file.path(outDir, "NormalizationSummarizationOutput.rda"))
    if (!all(c("pkgNorm", "filesNorm") %in% obj))
      stop(file.path(outDir, "NormalizationSummarizationOutput.rda"), " does not have 'pkgNorm' and/or 'filesNorm' object(s).")
    normOut <- get("filesNorm")
    if (!missing(pkgname) & pkgname != get("pkgNorm"))
      stop("Annotation used for normalization is different from the specified one.")
    pkgname <- get("pkgNorm")
    rm(list=c("pkgNorm", "filesNorm"))
  }

  fn <- paste(pkgname, 'CrlmmInfo.rda', sep='')
  load(system.file("extdata", fn, package=pkgname))

  ## myenv should take abou 65MB RAM
  ## I'll assume this is OK for now
  ## and not subset it
  myenv <- get(paste(pkgname,"Crlmm",sep="")); rm(list=paste(pkgname,"Crlmm",sep=""))
  thePriors <- get("priors", myenv)

  tmpdf <- getSnpLocInfo(pkgname)
  XIndex <- tmpdf[['chrom']] == 'X'

  ## These do not have Hapmap Data
  Index <- !get("hapmapCallIndex", myenv)

  batches <- splitIndicesByLength(1:nrow(tmpdf), batch_size)
  
  analysis <- read.table(file.path(outDir, "analysis.txt"), stringsAsFactors=FALSE)
  n.chunks <- length(batches)
  n.files <- as.integer(analysis[which(analysis[,1] == "nsamples"), 2])
  if (verbose){
    message("Genotyping.")
    pb <- txtProgressBar(min=1, max=n.chunks, initial=1, style=3)
  }

  alleleA <- normOut[['alleleA']]
  alleleB <- normOut[['alleleB']]
  initialCalls <- normOut[['initialCalls']]
  fs <- normOut[['fs']]
  f0 <- normOut[['f0']]
  snr <- normOut[['snr']]

  ## these are for debugging
    if (oligoDEBUG){
        dim1 <- c(nrow(tmpdf), 3)
        centersB4 <- ff(vmode='double', dim=dim1,
                        pattern=file.path(outDir, 'centersB4-'))
        scalesB4 <- ff(vmode='double', dim=dim1,
                       pattern=file.path(outDir, 'scalesB4-'))
        nB4 <- ff(vmode='double', dim=dim1,
                  pattern=file.path(outDir, 'nB4-'))
        distB4 <- ff(vmode='double', dim=c(nrow(tmpdf), n.files, 3),
                     pattern=file.path(outDir, 'distB4-'))
        centersAftr <- ff(vmode='double', dim=dim1,
                          pattern=file.path(outDir, 'centersAftr-'))
        scalesAftr <- ff(vmode='double', dim=dim1,
                         pattern=file.path(outDir, 'scalesAftr-'))
        nAftr <- ff(vmode='double', dim=dim1,
                    pattern=file.path(outDir, 'nAftr-'))
        distAftr <- ff(vmode='double', dim=c(nrow(tmpdf), n.files, 3),
                       pattern=file.path(outDir, 'distAftr-'))
        moveSz <- ff(vmode='double', dim=dim1,
                     pattern=file.path(outDir, 'moveSz-'))
        dbgs <- c("centersB4", "scalesB4", "nB4", "distB4",
                  "centersAftr","scalesAftr", "nAftr", "distAftr", "moveSz")
        sapply(dbgs, function(x) open(get(x)))
    }
    ## end debug

  dim2 <- c(nrow(tmpdf), n.files)
  calls <- ff(vmode='integer', dim=dim2,
                 pattern=file.path(outDir, 'crlmm-calls-'))
  conf <- ff(vmode='double', dim=dim2,
                 pattern=file.path(outDir, 'crlmm-confs-'))
  llr <- ff(vmode='double', dim=dim2,
                 pattern=file.path(outDir, 'crlmm-llr-'))
  objs <- c("alleleA", "alleleB", "initialCalls", "fs", "f0", "snr",
            "calls", "conf", "llr")
  sapply(objs, function(x) open(get(x)))

  i <- 1
  for (snps in batches){
      ## these do not have hapmap calls
      index <- which(Index[snps])
    
      initialCalls[-index,] <- NA
      rparams <- getGenotypeRegionParams(alleleA[index,, drop=FALSE]-alleleB[index,,drop=FALSE],
                                         initialCalls[index,,drop=FALSE],
                                         fs[index,,drop=FALSE], verbose=FALSE)
      rparams <- updateAffySnpParamsSingle(rparams, thePriors, verbose=FALSE, d0s=d0s)

      ## get params only for the snps in batch
      params <- get("params", myenv)
      params$centers <- params$centers[snps,]
      params$scales <- params$scales[snps,]
      params$N <- params$N[snps,]
      params  <- replaceAffySnpParamsSingle(params, rparams, index)

      if (oligoDEBUG){
          centersB4[snps,] <- params$centers
          scalesB4[snps,] <- params$scales
          nB4[snps,] <- params$N
      }

      theA <- alleleA[snps,, drop=FALSE]
      theB <- alleleB[snps,, drop=FALSE]
      theFS <- fs[snps,, drop=FALSE]

      myDist <- getAffySnpDistanceSingle56(theA-theB, params, theFS)
    
      if (oligoDEBUG){
          distB4[snps,,] <- myDist
      }

      myDist[,,-2] <- balance*myDist[,,-2]
      xindex <- which(XIndex[snps])

      maleIndex <- rep(FALSE, n.files)
      tmpIC <- getAffySnpCalls56(myDist, xindex, maleIndex, verbose=FALSE)

      LLR <- getAffySnpConfidence56(myDist, tmpIC, xindex, maleIndex, verbose=FALSE)
      
      if (recalibrate){
          for(k in 1:3)
              tmpIC[ tmpIC == k & LLR < minLLRforCalls[k]] <- NA

          rparams <- getGenotypeRegionParams(theA-theB, tmpIC, theFS, verbose=FALSE)
          rparams <- updateAffySnpParamsSingle(rparams, thePriors, d0s=d0s)

          if (oligoDEBUG){
              centersAftr[snps,] <- rparams$centers
              scalesAftr[snps,] <- rparams$scales
              nAftr[snps,] <- rparams$N
              moveSz[snps,] <- rparams$centers-params$centers
          }
      
          myDist <- getAffySnpDistanceSingle56(theA-theB, rparams, theFS)

          if (oligoDEBUG){
              distAftr[snps,,] <- myDist
          }

          myDist[,,-2] <- balance*myDist[,,-2]
          tmpIC <- getAffySnpCalls56(myDist, xindex, maleIndex, verbose=FALSE)
          initialCalls[snps,] <- tmpIC
          LLR <- getAffySnpConfidence56(myDist, tmpIC, xindex, maleIndex, verbose=FALSE)
          rm(myDist)
      } ##recalibrate

      calls[snps,] <- tmpIC
      llr[snps,] <- LLR
      conf[snps,] <- LLR2conf(tmpIC, LLR, snr[], pkgname)

      if (verbose)
          setTxtProgressBar(pb, i)
      i <- i+1
  }
  rownames(calls) <- rownames(llr) <- rownames(conf) <- tmpdf[['man_fsetid']]
  colnames(calls) <- colnames(llr) <- colnames(conf) <- files
  rownames(alleleA) <- rownames(alleleB) <- rownames(initialCalls) <- tmpdf[['man_fsetid']]
  colnames(alleleA) <- colnames(alleleB) <- colnames(initialCalls) <- files
  dimnames(fs) <- list(tmpdf[['man_fsetid']], files)
  ## names(f0) <- names(snr) <- files
  for (obj in objs){
      save(list=obj, file=file.path(outDir, paste(obj, '.RData', sep='')))
      close(get(obj))
  }
    if (oligoDEBUG){
        for (obj in dbgs){
            save(list=obj, file=file.path(outDir, paste(obj, '.RData', sep='')))
            close(get(obj))
        }
    }
  
  return(list(alleleA=alleleA, alleleB=alleleB, initialCalls=initialCalls,
              fs=fs, f0=f0, snr=snr, calls=calls,
              conf=conf, llr=llr))
}

getSnpNames <- function(pkgname){
    conn <- db(get(pkgname))
    sql <- paste("SELECT man_fsetid FROM featureSet",
                 "WHERE man_fsetid LIKE 'SNP%'")
    snps <- dbGetQuery(conn, sql)[[1]]
    sort(snps)
}

getSnpLocInfo <- function(pkgname){
  tmpdf <- dbGetQuery(db(get(pkgname)),
                      paste("SELECT man_fsetid, chrom",
                            "FROM featureSet WHERE man_fsetid LIKE 'SNP%'"))
##   i <- is.na(tmpdf[['chrom']])
##   tmpdf[i, 'chrom'] <- 0
##   tmpdf[["index"]] <- 1:nrow(tmpdf)
  tmpdf[order(tmpdf[['man_fsetid']]),]
}

getProbeInfo <- function(pkgname){
    tmp <- dbGetQuery(db(get(pkgname)),
                      paste("SELECT fid, man_fsetid, allele",
                            "FROM pmfeature INNER JOIN featureSet USING(fsetid)"))
    tmp[order(tmp$man_fsetid, tmp$allele),]
}
