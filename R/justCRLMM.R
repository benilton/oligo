liteNormalization <- function(celFiles, destDir, batch_size=40000, verbose=TRUE){
  ## Check existence of directory (destination)
  ## The destDir should not exist (yet)
  if (file.exists(destDir)) stop(message(destDir, "exists."))
  dir.create(destDir)

  ## Assumes all files are of the same type.
  header <- readCelHeader(celFiles[1])
  
  ## Determine pkg and number of files to read at once
  pkgname <- cleanPlatformName(header$chiptype)
  stopifnot(require(pkgname, character.only=TRUE))
  conn <- db(get(pkgname))
  nfeat <- dbGetQuery(conn, "SELECT row_count FROM table_info WHERE tbl='pmfeature'")[[1]]
  nfiles <- max((batch_size * 20 * length(celFiles)) %/% nfeat, 1)

  ## Determine batches of CEL
  batches <- split(celFiles, rep(1:length(celFiles), each=nfiles, length.out=length(celFiles)))
  batches.out <- split(file.path(destDir, "/normalized-", basename(celFiles), fsep=""),
                       rep(1:length(celFiles), each=nfiles, length.out=length(celFiles)))
  if (verbose) message("Preparing environment for normalization.")
  suppressWarnings(createCel(batches.out[[1]][1], header=header))
  if (length(unlist(batches.out))>1)
    sapply(unlist(batches.out)[-1], function(x) file.copy(batches.out[[1]][1], x))

  ## Get feature IDs and load reference
  fid <- dbGetQuery(conn, "SELECT fid FROM pmfeature")[[1]]
  fid <- sort(fid)
  load(system.file("extdata", paste(pkgname, "Ref.rda", sep=""), package=pkgname))
  reference <- sort(reference)

  txt <- sprintf("Normalization: %06.2f percent done.", 0)
  if (verbose) cat(txt)
  del <- paste(rep("\b", nchar(txt)), collapse="", sep="")
  for (i in 1:length(batches)){
    new <- normalize.quantiles.use.target(readCelIntensities(batches[[i]], indices=fid), reference, copy=FALSE)
    for (j in 1:length(batches[[i]]))
      updateCel(batches.out[[i]][j], indices=fid, intensities=new[,j])
    rm(new)
    if (verbose){
      cat(del)
      cat(sprintf("Normalization: %06.2f percent done.", i/length(batches)*100))
    }
  }
  if (verbose) cat("\n")
}

justCRLMM <- function(filenames, batch_size=40000,
                      minLLRforCalls=c(5, 1, 5), recalibrate=TRUE,
                      balance=1.5, phenoData=NULL, verbose=TRUE,
                      pkgname=NULL, tmpdir=tempdir()){
  tmpdir <- tempfile("crlmm.tmp", tmpdir)
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
  liteNormalization(filenames, destDir=tmpdir, batch_size=batch_size, verbose=verbose)
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
    t0 <- proc.time()
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

    pms <- readCelIntensities(filenames, indices=tmp[["fid"]])
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
##      pacc <- LLR2conf(myCalls, LLR, theSNR[i,], annotation(sqs))
    }
    save(myDist, myCalls, LLR, file=paste(randomName, i, sep="."))
    rm(myDist, correction, Index, k, LLR, myCalls, myenv, params, rparams, snpsIn, sqs, XIndex)
    if (verbose){
      cat(del)
      cat(sprintf("Genotyping: %06.2f percent done.", i/length(snps)*100))
    }
  }
  finalDist <- finalCalls <- finalConfs <- finalSumm <- NULL

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
    if (i == 1){
      finalDist <- myDist
    }else{
      finalDist <- my.abind(finalDist, myDist)
    }
    rm(myDist)
    load(paste(randomName, i, "summ", sep="."))
    finalSumm <- rbind(finalSumm, theSumm)
    rm(theSumm)
    if (verbose){
      cat(del)
      cat(sprintf("Finalizing: %06.2f percent done.", i/length(snps)*100))
    }
  }

  save(finalDist, file="finalDist.rda")
  
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

my.abind <- function(x, y){
  dim.x <- dim(x)
  dim.y <- dim(y)
  stopifnot(identical(dim.x[-1], dim.y[-1]))
  nrows.x <- dim.x[1]
  nrows.y <- dim.y[1]
  dim(x) <- c(nrows.x, prod(dim.x[-1]))
  dim(y) <- c(nrows.y, prod(dim.y[-1]))
  z <- rbind(x,y)
  dim(z) = c(nrows.x+nrows.y, dim.x[-1])
  z
}

  
justCRLMM2 <- function(filenames, batch_size=10000, chr=1,
                       minLLRforCalls=c(5, 1, 5), recalibrate=TRUE,
                       balance=1.5, phenoData=NULL, verbose=TRUE,
                       pkgname=NULL, tmpdir=tempdir()){
  ## To genotype by Chromosome, must assumed normalized CELs
  
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
##  liteNormalization(filenames, destDir=tmpdir, pkgname=pkgname, verbose=verbose)
##  liteNormalization2(filenames, destDir=tmpdir, batch_size=batch_size, verbose=verbose)
##  filenames <- paste(tmpdir, "/normalized-", basename(filenames), sep="")

  if (as.character(chr) == "1"){
    snps <- dbGetQuery(db(get(pkgname)),
                       paste("SELECT man_fsetid FROM featureSet WHERE chrom IN ('1', 'MT') OR chrom IS NULL",
                             "AND man_fsetid LIKE 'SNP%' ORDER BY man_fsetid"))[[1]]
  }else{
    snps <- dbGetQuery(db(get(pkgname)),
                       paste("SELECT man_fsetid FROM featureSet WHERE chrom='", chr,
                             "' AND man_fsetid LIKE 'SNP%' ORDER BY man_fsetid", sep=""))[[1]]
  }
  snps <- split(snps, rep(1:length(snps), each=batch_size, length.out=length(snps)))

  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}
  theSNR <- matrix(NA, nrow=length(snps), ncol=length(filenames))

  prefix <- "SELECT fid, man_fsetid, pmfeature.allele, pmfeature.strand FROM featureSet, pmfeature WHERE man_fsetid IN ("
  suffix <- ") AND pmfeature.fsetid = featureSet.fsetid ORDER BY fid"
  allSnps <- dbGetQuery(db(get(pkgname)), "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%' ORDER BY man_fsetid")[[1]]

  txt <- sprintf("Genotyping: %06.2f percent done. Timing: %06.2f", 0, 0)
  if (verbose) cat(txt)
  del <- paste(rep("\b", nchar(txt)), collapse="")

  for (i in 1:length(snps)){
    t0 <- proc.time()
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

    pms <- readCelIntensities(filenames, indices=tmp[["fid"]])
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
    }
    save(myDist, myCalls, LLR, file=paste(randomName, i, sep="."))
    rm(myDist, correction, Index, k, LLR, myCalls, myenv, params, rparams, snpsIn, sqs, XIndex)
    if (verbose){
      cat(del)
      cat(sprintf("Genotyping: %06.2f percent done. Timing: %06.2f", i/length(snps)*100, (proc.time()-t0)[3]))
    }
  }
  finalDist <- finalCalls <- finalConfs <- finalSumm <- NULL

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
    if (i == 1){
      finalDist <- myDist
    }else{
      finalDist <- my.abind(finalDist, myDist)
    }
    rm(myDist)
    load(paste(randomName, i, "summ", sep="."))
    finalSumm <- rbind(finalSumm, theSumm)
    rm(theSumm)
    if (verbose){
      cat(del)
      cat(sprintf("Finalizing: %06.2f percent done.", i/length(snps)*100))
    }
  }

  rm(finalDist)
##  save(finalDist, file="finalDist.rda")
  
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
  featureNames(out) <- unlist(snps)
  sampleNames(out)  <- sns
  phenoData(out)    <- phenoData
  out$crlmmSNR <- snr

  obj <- paste("chr", chr, sep="")
  assign(obj, out)
  rm(out)
  save(list=obj, file=paste(obj, ".rda", sep=""))
}

genotypeByChromosome <- function(filenames, batch_size=10000,
                                 minLLRforCalls=c(5, 1, 5), recalibrate=TRUE,
                                 balance=1.5, phenoData=NULL, verbose=TRUE,
                                 pkgname=NULL, tmpdir=tempdir()){
  if (is.null(pkgname))
    pkgname <- cleanPlatformName(readCelHeader(filenames[1])$chiptype)
  require(pkgname, character.only=TRUE)
  sql <- "SELECT DISTINCT chrom FROM featureSet WHERE man_fsetid LIKE 'SNP%'"
  chrs <- dbGetQuery(db(get(pkgname)), sql)[[1]]
  chrs <- chrs[!is.na(chrs) & chrs != "MT"]
  for(i in chrs){
    message("Processing chromosome", i)
    justCRLMM2(filenames, batch_size=batch_size, chr=i,
               minLLRforCalls=minLLRforCalls, recalibrate=recalibrate,
               balance=balance, phenoData=phenoData, verbose=verbose,
               pkgname=pkgname, tmpdir=tmpdir)
  }
}

justCRLMMv2 <- function(filenames, tmpdir, batch_size=40000,
                        minLLRforCalls=c(5, 1, 5), recalibrate=TRUE,
                        balance=1.5, verbose=TRUE, pkgname){
  stopifnot(!(missing(tmpdir)|file.exists(tmpdir)))
  tmpdir <- gsub("\\/$", "",  tmpdir)
  
  ## PHENODATA
  ## gender is not a key thing here
  ## algorithm behaves robustly...
  phenoData <- new("AnnotatedDataFrame",
                   data=data.frame(gender=rep("female",
                   length(filenames))))
  
  
  randomName <- tempfile("crlmm.", tmpdir)
  chips <- sapply(filenames, function(x) readCelHeader(x)$chiptype)
  if (length(unique(chips)) > 1){
    print(table(chips))
    stop("All the CEL files must be of the same type.")
  }
  if (missing(pkgname))
    pkgname <- cleanPlatformName(chips[1])
  rm(chips)
  require(pkgname, character.only=TRUE)

  sql.tmp <- "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%' AND chrom = 'X'"
  snps.chrX <- dbGetQuery(db(get(pkgname)), sql.tmp)[[1]]
  rm(sql.tmp)
  sns <- basename(filenames)
  liteNormalization(filenames, destDir=tmpdir, batch_size=batch_size, verbose=verbose)
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
    pnVec <- paste(tmp[["man_fsetid"]],
                   c("A", "B")[tmp[["allele"]]+1],
                   c("S", "A")[tmp[["strand"]]+1], sep="")
    idx <- order(pnVec)
    tmp[["man_fsetid"]] <- tmp[["allele"]] <- tmp[["strand"]] <- NULL
    ngenes <- length(unique(pnVec))

    pms <- readCelIntensities(filenames, indices=tmp[["fid"]])
    pms <- pms[idx, ]
    dimnames(pms) <- NULL
    theSumm <- .Call("rma_c_complete_copy", pms, pms,
                     pnVec[idx], ngenes,  body(bg.dens),
                     new.env(), FALSE, FALSE,
                     as.integer(2), PACKAGE="oligo")
    save(theSumm, file=paste(randomName, i, "summ", sep="."))
    sqs <- sqsFrom(theSumm)

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
    myenv$params$centers <- myenv$params$centers[snpsIn,,]
    myenv$params$scales <- myenv$params$scales[snpsIn,,]
    myenv$params$N <- myenv$params$N[snpsIn,]

    Index <- which(!get("hapmapCallIndex",myenv))
    myCalls <- matrix(NA,dim(sqs)[1],dim(sqs)[2])
    myCalls[Index,] <- getInitialAffySnpCalls(correction,Index,verbose=FALSE)
    rparams <- getAffySnpGenotypeRegionParams(sqs, myCalls, correction$fs,
                                              subset=Index,verbose=FALSE)
    oneStrand <- apply(is.na(getM(sqs[,1])[,1,]), 1,
                       function(v) ifelse(length(ll <- which(v))==0, 0, ll))
    rparams <- updateAffySnpParams(rparams, thePriors, oneStrand, verbose=FALSE)
    params  <- replaceAffySnpParams(get("params",myenv), rparams, Index)
    myDist <- getAffySnpDistance(sqs, params, correction$fs)
    myDist[,,-2,] <- balance*myDist[,,-2,]
    
    XIndex <- which(snps[i] %in% snps.chrX)
    myCalls <- getAffySnpCalls(myDist,XIndex,maleIndex,verbose=FALSE)
    LLR <- getAffySnpConfidence(myDist,myCalls,XIndex,maleIndex,verbose=FALSE)

    if(recalibrate){
      for(k in 1:3)
        myCalls[myCalls == k & LLR < minLLRforCalls[k]] <- NA
      rm(LLR)
      myCalls[, theSNR[i,] < 3.675] <- NA
      rparams <- getAffySnpGenotypeRegionParams(sqs, myCalls,
                                                correction$fs,
                                                verbose=FALSE)
                                                rm(myCalls)
      rparams <- updateAffySnpParams(rparams, thePriors, oneStrand)
      myDist <- getAffySnpDistance(sqs, rparams, correction$fs, verbose=FALSE)
      myDist[,,-2,] <- balance*myDist[,,-2,]
      rm(oneStrand)
      myCalls <- getAffySnpCalls(myDist,XIndex, maleIndex, verbose=FALSE)
      LLR <- getAffySnpConfidence(myDist,myCalls,XIndex,maleIndex,verbose=FALSE)
    }
    save(myDist, myCalls, LLR, file=paste(randomName, i, sep="."))
    
##     save(myDist, file=paste(randomName, i, sep="."))
##     rownames(myCalls) <- rownames(myDist) <- rownames(LLR) <- featureNames(sqs)
##     colnames(myCalls) <- colnames(myDist) <- colnames(LLR) <- sampleNames(sqs)
##     write.table(myCalls, file.path(tmpdir, "crlmm-calls.txt"), append=TRUE, quote=FALSE, sep="\t", col.names=(i==1))
##     write.table(myDist, file.path(tmpdir, "crlmm-dist.txt"), append=TRUE, quote=FALSE, sep="\t", col.names=(i==1))
##     write.table(LLR, file.path(tmpdir, "crlmm-llr.txt"), append=TRUE, quote=FALSE, sep="\t", col.names=(i==1))
##     write.table(antisenseThetaA(sqs), file.path(tmpdir, "crlmm-alleleA-antisense.txt"), append=TRUE, quote=FALSE, sep="\t", col.names=(i==1))
##     write.table(antisenseThetaB(sqs), file.path(tmpdir, "crlmm-alleleB-antisense.txt"), append=TRUE, quote=FALSE, sep="\t", col.names=(i==1))
##     write.table(senseThetaA(sqs), file.path(tmpdir, "crlmm-alleleA-sense.txt"), append=TRUE, quote=FALSE, sep="\t", col.names=(i==1))
##     write.table(senseThetaB(sqs), file.path(tmpdir, "crlmm-alleleB-sense.txt"), append=TRUE, quote=FALSE, sep="\t", col.names=(i==1))
    
    rm(myDist, correction, Index, k, LLR, myCalls, myenv, params, rparams, snpsIn, sqs, XIndex)
    if (verbose){
      cat(del)
      cat(sprintf("Genotyping: %06.2f percent done.", i/length(snps)*100))
    }
  }
  finalDist <- finalCalls <- finalConfs <- finalSumm <- NULL

  if (verbose){
    cat("\n")
    txt <- sprintf("Finalizing: %06.2f percent done.", 0)
    del <- paste(rep("\b", nchar(txt)), collapse="")
    cat(txt)
  }
  
  for (i in 1:length(snps)){
    fname <- paste(randomName, i, sep=".")
    load(fname)
    finalCalls <- rbind(finalCalls, myCalls)
    rm(myCalls)
    finalConfs <- rbind(finalConfs, LLR)
    rm(LLR)
    if (i == 1){
      finalDist <- myDist
    }else{
      finalDist <- my.abind(finalDist, myDist)
    }
    rm(myDist)
    fname2 <- paste(randomName, i, "summ", sep=".")
    load(fname2)
    finalSumm <- rbind(finalSumm, theSumm)
    rm(theSumm)
    if (verbose){
      cat(del)
      cat(sprintf("Finalizing: %06.2f percent done.", i/length(snps)*100))
    }
    unlink(fname)
    rm(fname)
    unlink(fname2)
    rm(fname2)
  }

  save(finalDist, file=file.path(tmpdir, "finalDist.rda"))
  
  finalSQS <- sqsFrom(finalSumm)
  ata <- antisenseThetaA(finalSQS)
  atb <- antisenseThetaB(finalSQS)
  sta <- senseThetaA(finalSQS)
  stb <- senseThetaB(finalSQS)
  rm(finalSQS)
  colnames(ata) <- colnames(atb) <- colnames(sta) <- colnames(stb) <- sns

  ## remove normalized CEL
  sapply(list.celfiles(tmpdir, full.names=T), unlink)

  snr <- matrix(exp(colMeans(log(theSNR))), nrow=1)
  colnames(snr) <- sns
  pacc <- as.matrix(LLR2conf(finalCalls, finalConfs, snr, pkgname))
  rownames(pacc) <- allSnps
  colnames(pacc) <- sns

  if (verbose) cat("\n")
  
##   out <- new("SnpCallSetPlus", calls=finalCalls, callsConfidence=pacc, LLR=finalConfs,
##              antisenseThetaA=ata, antisenseThetaB=atb, senseThetaA=sta, senseThetaB=stb)

  rownames(finalCalls) <- rownames(finalConfs) <- allSnps
  colnames(finalCalls) <- colnames(finalConfs) <- sns
  write.table(finalCalls, file.path(tmpdir, "crlmm-calls.txt"), quote=FALSE, sep="\t", col.names=TRUE)
  rm(finalCalls)
  write.table(pacc, file.path(tmpdir, "crlmm-conf.txt"), quote=FALSE, sep="\t", col.names=TRUE)
  rm(pacc)
  write.table(finalConfs, file.path(tmpdir, "crlmm-llr.txt"), quote=FALSE, sep="\t", col.names=TRUE)
  rm(finalConfs)
  write.table(ata, file.path(tmpdir, "crlmm-alleleA-antisense.txt"), quote=FALSE, sep="\t", col.names=TRUE)
  rm(ata)
  write.table(atb, file.path(tmpdir, "crlmm-alleleB-antisense.txt"), quote=FALSE, sep="\t", col.names=TRUE)
  rm(atb)
  write.table(sta, file.path(tmpdir, "crlmm-alleleA-sense.txt"), quote=FALSE, sep="\t", col.names=TRUE)
  rm(sta)
  write.table(stb, file.path(tmpdir, "crlmm-alleleB-sense.txt"), quote=FALSE, sep="\t", col.names=TRUE)
  rm(stb)
  write.table(snr, file.path(tmpdir, "crlmm-snr.txt"), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
}
