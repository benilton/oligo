genotypeOne3 <- function(files, tmpdir=getwd(), batch_size=40000, balance=1.5, minLLRforCalls=c(5, 1, 5), recalibrate=TRUE, verbose=TRUE, pkgname, reference=TRUE, d0s=80){
  if (!file.exists(tmpdir)){
    tmp <- normalizeOne(files, tmpdir, pkgname=pkgname, reference=reference)
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
  mu.logsnr <- 2.12
  sd.logsnr <- 0.08
  tmpSNR <- log(readBin(tmp$snr, numeric(), n.files))-mu.logsnr
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
##     rparams <- getGenotypeRegionParams(alleleA[index,]-alleleB[index,],
##                                        initialCalls[index,],
##                                        fs[index,], verbose=FALSE)
    rparams <- getGenotypeRegionParams(alleleA[index,]-alleleB[index,],
                                       initialCalls[index,],
                                       fs[index,], verbose=FALSE)
    rparams <- updateAffySnpParamsSingle(rparams, thePriors, verbose=FALSE, d0s=d0s)
    params <- get("params", myenv)
    params$centers <- params$centers[overall_pos,]
    params$scales <- params$scales[overall_pos,]
    params$N <- params$N[overall_pos,]

    tmpPrior <- params$N
    tmpPrior <- tmpPrior/rowSums(tmpPrior)
    tmpPrior[tmpPrior == 0] <- .001
    tmpPrior <- tmpPrior/rowSums(tmpPrior)
    
    ## try
    params$N[params$N == 0] <- 1

    params  <- replaceAffySnpParamsSingle(params, rparams, index)

    ##try
    params$scales <- sqrt(params$scales^2*(1+1/params$N)+sd.logsnr^2)
    
    myDist <- getAffySnpDistanceSingle56(sweep(alleleA-alleleB, 2, mu.logsnr), params, fs)
    ##SAVE THE ABOVE
##    myDist[,,-2] <- balance*myDist[,,-2]
    XIndex <- integer()
    if (length(grep("chrX", tmp$initialCalls[i])) > 0)
      XIndex <- 1:nrows

    maleIndex <- rep(FALSE, n.files)
    initialCalls <- getAffySnpCalls56(myDist, XIndex, maleIndex, verbose=FALSE)
    LLR <- getAffySnpConfidence56(myDist, initialCalls, XIndex, maleIndex, verbose=FALSE)

    p <- compConf(myDist, initialCalls, tmpPrior)
    minPforCalls <- c(.99, .99, .99)
    
    if (recalibrate){
      for(k in 1:3)
        initialCalls[ initialCalls == k & p < minPforCalls[k]] <- NA

      rparams <- getGenotypeRegionParams(sweep(alleleA-alleleB, 2, mu.logsnr),
                                         initialCalls,
                                         fs, verbose=FALSE)
      rparams <- updateAffySnpParamsSingle(rparams, thePriors, d0s=d0s)

      #try
      rparams$N[rparams$N == 0] <- 1
      rparams$scales <- sqrt(rparams$scales^2*(1+1/rparams$N)+sd.logsnr^2)
##      rparams$centers <- 2*rparams$centers-params$centers
      
      myDist <- getAffySnpDistanceSingle56(sweep(alleleA-alleleB, 2, mu.logsnr), rparams, fs)

      fileAA <- file.path(tmpdir, "crlmm-dstAA.txt")
      fileAB <- file.path(tmpdir, "crlmm-dstAB.txt")
      fileBB <- file.path(tmpdir, "crlmm-dstBB.txt")
      dimnames(myDist) <- list(tmpdf[["man_fsetid"]][overall_pos],
                               basename(files), c("AA", "AB", "BB"))
      suppressWarnings(write.table(myDist[,, 1], fileAA, append=TRUE,
                                   quote=FALSE, sep="\t",
                                   col.names=(i==1)))
      suppressWarnings(write.table(myDist[,, 2], fileAB, append=TRUE,
                                   quote=FALSE, sep="\t",
                                   col.names=(i==1)))
      suppressWarnings(write.table(myDist[,, 3], fileBB, append=TRUE,
                                   quote=FALSE, sep="\t",
                                   col.names=(i==1)))

      ### SAVE THE ABOVE
##      myDist[,,-2] <- balance*myDist[,,-2]
      initialCalls <- getAffySnpCalls56(myDist, XIndex, maleIndex, verbose=FALSE)
      p <- compConf(myDist, initialCalls, tmpPrior)
      LLR <- getAffySnpConfidence56(myDist, initialCalls, XIndex, maleIndex, verbose=FALSE)
      rm(myDist)
    }
##    callsConfidence <- LLR2conf(initialCalls, LLR, readBin(tmp$snr, numeric(), n.files), pkgname)
    callsConfidence <- p

    rownames(initialCalls) <- rownames(LLR) <- rownames(callsConfidence) <- tmpdf[["man_fsetid"]][overall_pos]
    colnames(initialCalls) <- colnames(LLR) <- colnames(callsConfidence) <- basename(files)

    suppressWarnings(write.table(initialCalls,
                                 file.path(tmpdir, "crlmm-calls.txt"),
                                 append=TRUE, quote=FALSE, sep="\t",
                                 col.names=(i==1)))
    suppressWarnings(write.table(LLR,
                                 file.path(tmpdir, "crlmm-llr.txt"),
                                 append=TRUE, quote=FALSE, sep="\t",
                                 col.names=(i==1)))
    suppressWarnings(write.table(callsConfidence,
                                 file.path(tmpdir, "crlmm-conf.txt"),
                                 append=TRUE, quote=FALSE, sep="\t",
                                 col.names=(i==1)))
    
    if (verbose){
      cat(del)
      txt <- sprintf("Genotyping: %06.2f percent done.", i/n.chunks*100)
      cat(txt)
      del <- paste(rep("\b", nchar(txt)), collapse="", sep="")
    }

  }
}

compConf <- function(myDist, initialCalls, priors){
  tmpP <- exp(-1/2*(log(2*pi)+myDist))
  tmpP[,,1] <- tmpP[,,1]*priors[,1]
  tmpP[,,2] <- tmpP[,,2]*priors[,2]
  tmpP[,,3] <- tmpP[,,3]*priors[,3]
  tot <- rowSums(tmpP, dims=2)
  p <- tmpP[,,1]/tot
  pAB <- tmpP[,,2]/tot
  pBB <- 1-p-pAB
  rm(tot, tmpP)
  i <- initialCalls == 2
  p[i] <- pAB[i]; rm(pAB)
  i <- initialCalls == 3
  p[i] <- pBB[i]; rm(pBB, i)
  return(p)
}
