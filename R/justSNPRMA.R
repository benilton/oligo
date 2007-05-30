pmSequenceByFid <- function (object, fid){
  sql <- paste("SELECT pmfeature.fid, seq FROM sequence, pmfeature",
               "WHERE pmfeature.fid=sequence.fid AND pmfeature.fid IN (",
               paste(fid, collapse=","), ")")
  tmp <- dbGetQuery(db(object), sql)
  tmp[["seq"]][match(fid, tmp[["fid"]])]
}

pmFragmentLengthByFid <- function(object, fid){
  sql <- paste("SELECT fid, fragment_length FROM featureSet, pmfeature",
               "WHERE pmfeature.fsetid=featureSet.fsetid AND pmfeature.fid IN (",
               paste(fid, collapse=","), ")")
  tmp <- dbGetQuery(db(object), sql)
  idx <- is.na(tmp[["fragment_length"]])
  tmp[idx, "fragment_length"] <- median(tmp[!idx, "fragment_length"])
  tmp[["fragment_length"]][match(fid, tmp[["fid"]])]
}

justSNPRMA <- function(filenames, tmpdir=getwd(), memory.bound=FALSE,
                       verbose=TRUE, phenoData=NULL, normalizeToHapmap=TRUE){

  if (verbose){
    message("Using ", getwd(), " to store temporary files.")
    message("Make sure this is a local directory.")
    message("Otherwise, performance will be decreased.")
  }
  
  ###################
  ## GET PM MATRIX ##
  ###################

  if (verbose) message("Reading CEL files.")
  headdetails <- read.celfile.header(filenames[1])
  pkgname <- cleanPlatformName(headdetails[["cdfName"]])
  require(pkgname, character.only=TRUE, quietly=TRUE)
  fid <- dbGetQuery(db(get(pkgname)), "SELECT fid FROM pmfeature")[[1]]
  tmpExprs <- createBufferedMatrix(length(fid), 0, directory=tmpdir)
  set.buffer.dim(tmpExprs, 50000, 1)
  for (i in 1:length(filenames)){
    AddColumn(tmpExprs)
    tmpExprs[,i] <- .Call("read_abatch", filenames[i], FALSE, FALSE, FALSE,
                          headdetails$cdfName, headdetails[["CEL dimensions"]],
                          FALSE, PACKAGE="affyio")[fid]
  }
  rm(headdetails, i); gc(); gc()

  ###################
  ## CORRECTIONS PL #
  ###################
  if (FALSE){
  if (verbose) message("Applying feature-level corrections.")
  snpLocation <- pmPosition(get(pkgname))

  ## check snpLocation... if too few probes
  ## at a given location, then change their locations
  ## to something else
  
  theUniqueLocs <- sort(unique(snpLocation))
  theCounts <- table(snpLocation)
  bad <- theUniqueLocs[which(theCounts < 1000)]
  if (length(bad)>0){
    good <- theUniqueLocs[-which(theCounts<1000)]
    for (i in bad)
      snpLocation[snpLocation == i] <- good[which.min(abs(good - i))]
    rm(good, i)
  }
  rm(theCounts, bad, theUniqueLocs)
  theUniqueLocs <- sort(unique(snpLocation)); gc(); gc()

  ssSize <- 2000
  ewApply(tmpExprs, log2)
  for (loc in theUniqueLocs){
    set.seed(1)
    set <- snpLocation == loc
    fid.loc <- sample(fid[set], min(sum(set), ssSize))
    fl <- ns(pmFragmentLengthByFid(get(pkgname), fid.loc), df=3)
    xs <- cbind(1, sequenceDesignMatrix(pmSequenceByFid(get(pkgname), fid.loc)), fl)
    betas <- solve(t(xs)%*%xs, t(xs)%*%tmpExprs[match(fid.loc, fid),])
    rm(xs, fid.loc)
    xx <- cbind(1, sequenceDesignMatrix(pmSequenceByFid(get(pkgname), fid[set])),
                ns(pmFragmentLengthByFid(get(pkgname), fid[set]), df=3,
                   knots=attr(fl, "knots"), Boundary.knots=attr(fl, "Boundary.knots")))
    rm(fl)
    if (memory.bound){
      for (i in 1:ncol(tmpExprs))
        tmpExprs[set, i] <- as.numeric(tmpExprs[set, i]-xx%*%betas[,i]+mean(tmpExprs[set, i]))
    }else{
      y <- tmpExprs[set,]
      y <- sweep(y, 2, colMeans(y), "+")
      tmpExprs[set,] <- y-xx%*%betas
      rm(y)
    }
  }
  rm(betas, xx, set)
  ewApply(tmpExprs, function(x) 2^x)
  gc(); gc()
  }
  
  ########################
  ##### NORMALIZATION ####
  ########################

  if (normalizeToHapmap){
    if (verbose) message("Normalizing to Hapmap.")
    load(system.file("extdata", paste(pkgname, "Ref.rda", sep=""), package=pkgname))
    reference <- sort(reference)
    for (i in 1:ncol(tmpExprs))
      tmpExprs[, i] <- reference[rank(tmpExprs[, i])]
  } else {
    normalize.BufferedMatrix.quantiles(tmpExprs, copy=FALSE)
    reference <- sort(tmpExprs[,1])
    save(reference, file=paste(pkgname, ".quantileReference.rda", sep=""))
  }
  rm(reference); gc(); gc()

  ########################
  #### SUMMARIZATION  ####
  ########################

  if (verbose) message("Summarizing.")
  
  ## get rma pars:
  ## put PMs in right order
  ## get pnVec
  ## get length(unique(pnVec))

  pnVec <- paste(probeNames(get(pkgname)),
                 c("A", "B")[pmAllele(get(pkgname))+1],
                 c("S", "A")[pmStrand(get(pkgname))+1],
                 sep="")
  idx <- order(pnVec)
  tmpExprs <- subBufferedMatrix(tmpExprs, idx)
  set.buffer.dim(tmpExprs, 50000, 1)
  pnVec <- pnVec[idx]
  rm(idx); ## gc()

  RowMode(tmpExprs)
  theSumm <- median.polish.summarize(tmpExprs, length(unique(pnVec)), pnVec)
  rm(tmpExprs, pnVec); ## gc()
  colnames(theSumm) <- basename(filenames)
  theSumm <- sqsFrom(theSumm)
  if (!is.null(phenoData)) phenoData(theSumm) <- phenoData
  annotation(theSumm) <- pkgname
  sampleNames(theSumm) <- basename(filenames)
  return(theSumm)
}
