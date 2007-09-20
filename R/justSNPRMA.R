justSNPRMA <- function(filenames,
                       verbose=TRUE, phenoData=NULL,
                       normalizeToHapmap=TRUE){

  ###################
  ## GET PM MATRIX ##
  ###################

  if (verbose) message("Reading CEL files.")
  headdetails <- readCelHeader(filenames[1])
  pkgname <- cleanPlatformName(headdetails[["chiptype"]])
  require(pkgname, character.only=TRUE, quietly=TRUE)
  fid <- dbGetQuery(db(get(pkgname)), "SELECT fid FROM pmfeature")[[1]]
  tmpExprs <- readCelIntensities(filenames, indices=fid)
  dimnames(tmpExprs) <- NULL
  rm(headdetails); gc(); gc()

  
  ########################
  ##### NORMALIZATION ####
  ########################

  if (normalizeToHapmap){
    if (verbose) message("Normalizing to Hapmap.")
    load(system.file("extdata", paste(pkgname, "Ref.rda", sep=""), package=pkgname))
    reference <- sort(reference)
    tmpExprs <- normalize.quantiles.use.target(tmpExprs, reference, copy=FALSE)
  } else {
    tmpExprs <- normalize.quantiles(tmpExprs)
    reference <- sort(tmpExprs[,1])
    save(reference, file=paste(pkgname, ".quantileReference.rda", sep=""))
  }
  rm(reference); gc(); gc()

  snpcnv <- pkgname == "pd.genomewidesnp.6"

  ########################
  #### SUMMARIZATION  ####
  ########################

  if (verbose) message("Summarizing.")
  
  ## get rma pars:
  ## put PMs in right order
  ## get pnVec
  ## get length(unique(pnVec))

  if (!snpcnv){
    pnVec <- paste(probeNames(get(pkgname)),
                   c("A", "B")[pmAllele(get(pkgname))+1],
                   c("S", "A")[pmStrand(get(pkgname))+1],
                   sep="")
  }else{
    pnVec <- paste(probeNames(get(pkgname)),
                   c("A", "B")[pmAllele(get(pkgname))+1],
                   sep="")
  }

  idx <- order(pnVec)
  tmpExprs <- tmpExprs[idx,]
  pnVec <- pnVec[idx]
  rm(idx)

  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}

  theSumm <-.Call("rma_c_complete_copy", tmpExprs, tmpExprs,
                  pnVec, length(unique(pnVec)), body(bg.dens),
                  new.env(), FALSE, FALSE,
                  as.integer(2), PACKAGE="oligo")

  
  rm(tmpExprs, pnVec)
  colnames(theSumm) <- basename(filenames)
  if (!snpcnv){
    theSumm <- sqsFrom(theSumm)
  }else{
    theSumm <- sqsFrom.SnpCnv(theSumm)
  }

  if (!is.null(phenoData)) phenoData(theSumm) <- phenoData
  annotation(theSumm) <- pkgname
  sampleNames(theSumm) <- basename(filenames)
  return(theSumm)
}
