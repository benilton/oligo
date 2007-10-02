liteNormalization <- function(filenames, destDir, pkgname, verbose=TRUE){
  dir.create(destDir)
  fid <- dbGetQuery(db(get(pkgname)), "SELECT fid FROM pmfeature")[[1]]
  fid <- sort(fid)
  load(system.file("extdata", paste(pkgname, "Ref.rda", sep=""), package=pkgname))
  reference <- sort(reference)

  txt <- sprintf("Normalization: %06.2f percent done.", 0)
  if (verbose) cat(txt)
  del <- paste(rep("\b", nchar(txt)), collapse="", sep="")
  for (ff in filenames){
    header <- readCelHeader(ff)
    new <- normalize.quantiles.use.target(matrix(readCel(ff, indices=fid, readOutliers = FALSE,
                                                         readMasked = FALSE, reorder=TRUE)$intensities,
                                                 ncol=1),
                                          reference, copy=FALSE)[,1]
    fout <- file.path(destDir, "/normalized-", basename(ff), fsep="")
    suppressWarnings(createCel(fout, header=header))
    updateCel(fout, indices=fid, intensities=new)
    rm(new, fout)
    if (verbose){
      cat(del)
      cat(sprintf("Normalization: %06.2f percent done.", which(filenames == ff)/length(filenames)*100))
    }
  }
  if (verbose) cat("\n")
}

justCRLMM <- function(filenames, batch_size=40000,
                      minLLRforCalls=c(5, 1, 5), recalibrate=TRUE,
                      balance=1.5, phenoData=NULL, verbose=TRUE, pkgname=NULL){
  tmpdir <- tempfile("crlmm.tmp", getwd())
  if (is.null(phenoData))
    stop("phenoData must be provided and must contain a variable called 'gender'.")

  if (!is.null(phenoData) & !is(phenoData, "AnnotatedDataFrame"))
    stop("phenoData provided, but must be of 'AnnotatedDataFrame' class.")

  if (is(phenoData, "AnnotatedDataFrame") & !("gender" %in% names(pData(phenoData))))
    stop("phenoData provided, but 'gender' variable not found.")
  
  randomName <- tempfile("crlmm.", tmpdir)
  chips <- sapply(filenames, function(x) readCelHeader(x)$chiptype)
  if (length(unique(chips)) > 1){
    print(table(chips))
    stop("All the CEL files must be of the same type.")
  }
  if (is.null(pkgname))
    pkgname <- cleanPlatformName(chips[1])
  snpcnv <- pkgname == "pd.genomewidesnp.6"
  rm(chips)
  require(pkgname, character.only=TRUE)

  sql.tmp <- "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%' AND chrom = 'X'"
  snps.chrX <- dbGetQuery(db(get(pkgname)), sql.tmp)[[1]]
  rm(sql.tmp)
  sns <- basename(filenames)
  liteNormalization(filenames, destDir=tmpdir, pkgname=pkgname, verbose=verbose)
  filenames <- paste(tmpdir, "/normalized-", basename(filenames), sep="")
  
  snps <- dbGetQuery(db(get(pkgname)), "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%' ORDER BY man_fsetid")[[1]]
  snps <- split(snps, rep(1:length(snps), each=batch_size, length.out=length(snps)))

  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}
  theSNR <- matrix(NA, nrow=length(snps), ncol=length(filenames))

  prefix <- "SELECT fid, man_fsetid, pmfeature.allele, pmfeature.strand FROM featureSet, pmfeature WHERE man_fsetid IN ("
  suffix <- ") AND pmfeature.fsetid = featureSet.fsetid ORDER BY fid"
  allSnps <- dbGetQuery(db(get(pkgname)), "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%' ORDER BY man_fsetid")[[1]]

  txt <- sprintf("Genotyping: %06.2f percent done.", 0)
  if (verbose) cat(txt)
  del <- paste(rep("\b", nchar(txt)), collapse="")

  for (i in 1:length(snps)){
    mid <- paste("'", snps[[i]],"'",  collapse=", ", sep="")
    sql <- paste(prefix, mid, suffix)
    tmp <- dbGetQuery(db(get(pkgname)), sql)
    if (!snpcnv){
      pnVec <- paste(tmp[["man_fsetid"]],
                     c("A", "B")[tmp[["allele"]]+1],
                     c("S", "A")[tmp[["strand"]]+1], sep="")
    }else{
      pnVec <- paste(tmp[["man_fsetid"]],
                     c("A", "B")[tmp[["allele"]]+1],
                     sep="")
    }
    
    idx <- order(pnVec)
    tmp[["man_fsetid"]] <- tmp[["allele"]] <- tmp[["strand"]] <- NULL
    ngenes <- length(unique(pnVec))

    pms <- readCelIntensities(filenames, indices=tmp[["fid"]], reorder=TRUE)
    pms <- pms[idx, ]
    dimnames(pms) <- NULL
    theSumm <- .Call("rma_c_complete_copy", pms, pms,
                     pnVec[idx], ngenes,  body(bg.dens),
                     new.env(), FALSE, FALSE,
                     as.integer(2), PACKAGE="oligo")
    save(theSumm, file=paste(randomName, i, "summ", sep="."))
    if (!snpcnv){
      sqs <- sqsFrom(theSumm)
    }else{
      sqs <- sqsFrom.SnpCnv(theSumm)
    }

    annotation(sqs) <- pkgname
    phenoData(sqs) <- phenoData
    rm(pms, mid, sql, tmp, pnVec, idx, ngenes)

    ### TRYING CRLMM HERE
    correction <- fitAffySnpMixture(sqs, verbose=FALSE)
    theSNR[i, ] <- correction$snr

    load(system.file(paste("extdata/", pkgname, "CrlmmInfo.rda", sep=""), package=pkgname))
    myenv <- get(paste(pkgname,"Crlmm",sep="")); rm(list=paste(pkgname,"Crlmm",sep=""))
    thePriors <- get("priors", myenv)
    snpsIn <- match(featureNames(sqs), allSnps)
    myenv$hapmapCallIndex <- myenv$hapmapCallIndex[snpsIn]
    if (!snpcnv){
      myenv$params$centers <- myenv$params$centers[snpsIn,,]
      myenv$params$scales <- myenv$params$scales[snpsIn,,]
    }else{
      myenv$params$centers <- myenv$params$centers[snpsIn,]
      myenv$params$scales <- myenv$params$scales[snpsIn,]
    }
    myenv$params$N <- myenv$params$N[snpsIn,]

    Index <- which(!get("hapmapCallIndex",myenv))
    myCalls <- matrix(NA,dim(sqs)[1],dim(sqs)[2])
    myCalls[Index,] <- getInitialAffySnpCalls(correction,Index,verbose=FALSE, sqsClass=class(sqs))
    rparams <- getAffySnpGenotypeRegionParams(sqs, myCalls, correction$fs,
                                              subset=Index,verbose=FALSE, sqsClass=class(sqs))
    if (!snpcnv){
      oneStrand <- apply(is.na(getM(sqs[,1])[,1,]), 1,
                         function(v) ifelse(length(ll <- which(v))==0, 0, ll))
      rparams <- updateAffySnpParams(rparams, thePriors, oneStrand, verbose=FALSE)
      params  <- replaceAffySnpParams(get("params",myenv), rparams, Index)
      myDist <- getAffySnpDistance(sqs, params, correction$fs)
      myDist[,,-2,] <- balance*myDist[,,-2,]
    }else{
      rparams <- updateAffySnpParamsSingle(rparams, thePriors, verbose=FALSE)
      params  <- replaceAffySnpParamsSingle(get("params",myenv), rparams, Index)
      myDist <- getAffySnpDistanceSingle(sqs, params, correction$fs)
      myDist[,,-2] <- balance*myDist[,,-2]
    }
        
    XIndex <- which(snps[i] %in% snps.chrX)
    myCalls <- getAffySnpCalls(myDist,XIndex,maleIndex,verbose=FALSE, sqsClass=class(sqs))
    LLR <- getAffySnpConfidence(myDist,myCalls,XIndex,maleIndex,verbose=FALSE, sqsClass=class(sqs))

    if(recalibrate){
      for(k in 1:3)
        myCalls[myCalls == k & LLR < minLLRforCalls[k]] <- NA
      rm(LLR)
      myCalls[, theSNR[i,] < 3.675] <- NA
      
      rparams <- getAffySnpGenotypeRegionParams(sqs, myCalls,
                                                correction$fs, verbose=FALSE, sqsClass=class(sqs))
      rm(myCalls)

      if (!snpcnv){
        rparams <- updateAffySnpParams(rparams, thePriors, oneStrand)
        myDist <- getAffySnpDistance(sqs, rparams, correction$fs, verbose=FALSE)
        myDist[,,-2,] <- balance*myDist[,,-2,]
        rm(oneStrand)
      }else{
        rparams <- updateAffySnpParamsSingle(rparams, thePriors)
        myDist <- getAffySnpDistanceSingle(sqs, rparams, correction$fs, verbose=FALSE)
        myDist[,,-2] <- balance*myDist[,,-2]
      }
      myCalls <- getAffySnpCalls(myDist,XIndex, maleIndex, verbose=FALSE, sqsClass=class(sqs))
      LLR <- getAffySnpConfidence(myDist,myCalls,XIndex,maleIndex,verbose=FALSE, sqsClass=class(sqs))
      rm(myDist)
##      pacc <- LLR2conf(myCalls, LLR, theSNR[i,], annotation(sqs))
    }
    save(myCalls, LLR, file=paste(randomName, i, sep="."))
    rm(correction, Index, k, LLR, myCalls, myenv, params, rparams, snpsIn, sqs, XIndex)
    if (verbose){
      cat(del)
      cat(sprintf("Genotyping: %06.2f percent done.", i/length(snps)*100))
    }
  }
  finalCalls <- finalConfs <- finalSumm <- NULL

  if (verbose){
    cat("\n")
    txt <- sprintf("Finalizing: %06.2f percent done.", 0)
    del <- paste(rep("\b", nchar(txt)), collapse="")
    cat(txt)
  }
  
  for (i in 1:length(snps)){
    load(paste(randomName, i, sep="."))
    finalCalls <- rbind(finalCalls, myCalls)
    rm(myCalls)
    finalConfs <- rbind(finalConfs, LLR)
    rm(LLR)

    load(paste(randomName, i, "summ", sep="."))
    finalSumm <- rbind(finalSumm, theSumm)
    rm(theSumm)
    if (verbose){
      cat(del)
      cat(sprintf("Finalizing: %06.2f percent done.", i/length(snps)*100))
    }
  }

  if (!snpcnv){
    finalSQS <- sqsFrom(finalSumm)
    ata <- antisenseThetaA(finalSQS)
    atb <- antisenseThetaB(finalSQS)
    sta <- senseThetaA(finalSQS)
    stb <- senseThetaB(finalSQS)
  }else{
    finalSQS <- sqsFrom.SnpCnv(finalSumm)
    ta <- thetaA(finalSQS)
    tb <- thetaB(finalSQS)
  }
  rm(finalSQS)
  
  if (verbose) cat("\n")
  warning("Please check the contents and remove the directory: ", tmpdir)
  snr <- exp(colMeans(log(theSNR)))
  pacc <- LLR2conf(finalCalls, finalConfs, snr, pkgname)
  if (!snpcnv){
    out <- new("SnpCallSetPlus", calls=finalCalls, callsConfidence=pacc, LLR=finalConfs,
               antisenseThetaA=ata, antisenseThetaB=atb, senseThetaA=sta, senseThetaB=stb)
    rm(ata, atb, sta, stb)
  }else{
    out <- new("SnpCnvCallSetPlus", calls=finalCalls, callsConfidence=pacc, LLR=finalConfs,
               thetaA=ta, thetaB=tb)
    rm(ta, tb)
  }

  annotation(out) <- pkgname
  
##  out <- new("SnpCallSet", calls=finalCalls, callsConfidence=pacc, LLR=finalConfs,  annotation=pkgname)
##  annotation(finalSQS) <- annotation(out)
  featureNames(out) <- allSnps
  sampleNames(out)  <- sns
  phenoData(out)    <- phenoData
  out$crlmmSNR <- snr

  return(out)

##  return(list(calls=out, sqs=finalSQS))
}
