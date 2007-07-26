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
