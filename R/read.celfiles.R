read.celfiles <- function(filenames,
                          pkgname=NULL,
                          sd=FALSE,
                          npixels=FALSE,
                          phenoData=new("AnnotatedDataFrame"),
                          featureData=NULL,
                          description=NULL,
                          notes="",
                          verbose = FALSE,
                          compress= FALSE,
                          rm.mask = FALSE,
                          rm.outliers=FALSE,
                          rm.extra=FALSE,
                          tempdir=getwd()){

  ## FIXME: this is not an acceptable function name
  tmp <- stuffForXYSandCELreaders(filenames, phenoData, description, notes,
                                  verbose)
  filenames <- tmp$filenames
  n <- length(filenames)
  if (n == 0)
    stop("no file name given")
  
  ## read the first file to see what we have
  headdetails <- .Call("ReadHeader", filenames[1], compress, PACKAGE="affyio")
  
  ##now we use the length
  dim.intensity <- headdetails[[2]]
  
  ##and the cdfname as ref
  ref.cdfName <- headdetails[[1]]

  if (verbose) cat("Creating objects outside R to store intensities.\n")
  if (verbose) cat("This may take a while... ")
  tmpExprs <- createBufferedMatrix(prod(dim.intensity), length(filenames), directory=tempdir)
  set.buffer.dim(tmpExprs, 300000, 1)
  RowMode(tmpExprs)
  if (verbose) cat("Done.", "Now reading CEL files", sep="\n")
  for (i in 1:length(filenames))
    tmpExprs[,i] <- .Call("read_abatch", as.list(filenames[i]), rm.mask, rm.outliers,
                          rm.extra, ref.cdfName, dim.intensity, verbose, PACKAGE="affyio")
  if (verbose) cat(" Done.\n")
    
  ## RI: We should check if the pd package is available here. If not try
  ## to install it

  ## load pdInfo
  if (is.null(pkgname))
    pkgname <- cleanPlatformName(ref.cdfName)
  library(pkgname, character.only=TRUE)
  message("Platform design info loaded.")
  arrayType <- kind(get(pkgname))

  ## BC: Nov 15-16 2005, the PDenv is ordered already in the way
  ##     we want PM/MMs to be. So, we need to reorder the cel
  ##     input, so the pm/mm methods are faster.
  pdInfo <- get(pkgname)
  if (is(pdInfo, "platformDesign")) {
    order_index <- pdInfo$order_index
  } else if (is(pdInfo, "DBPDInfo")) {
    ## FIXME, need to figure out how to reorder in a useful
    ## way for the SQLite-based PDInfo stuff
    order_index <- 1:nrow(tmpExprs)
  } else {
    stop("unknown platform info type: ", class(pdInfo))
  }
  for (i in 1:ncol(tmpExprs)){
    tmpExprs[,i] <- tmpExprs[order_index, i]
  }
  
  rownames(tmpExprs) <- 1:nrow(tmpExprs)
  colnames(tmpExprs) <- basename(filenames)

  theClass <- switch(arrayType,
                     tiling="TilingFeatureSet",
                     expression="ExpressionFeatureSet",
                     SNP="SnpFeatureSet",
                     exon="ExonFeatureSet",
                     stop("unknown array type: ", arrayType))

  if (is.null(featureData))
    featureData <- new("AnnotatedDataFrame",
                       data=data.frame(idx=1:nrow(tmpExprs),
                         row.names=1:nrow(tmpExprs)),
                       varMetadata=data.frame(varLabels="idx",
                         row.names="idx"))
  out <- new(theClass,
             exprs=tmpExprs,
             platform=pkgname,
             manufacturer="Affymetrix",
             phenoData=tmp$phenoData,
             featureData=featureData,
             experimentData=tmp$description,
             annotation=pkgname)
  out
}

read.affybatch <- function(..., filenames=character(0),
                           phenoData=new("AnnotatedDataFrame"),
                           description=NULL,
                           notes="",
                           compress = getOption("BioC")$affy$compress.cel,
                           rm.mask = FALSE, rm.outliers=FALSE, rm.extra=FALSE,
                           verbose = FALSE) {
    .Deprecated("read.celfiles")

    auxnames <- as.list(substitute(list(...)))[-1]
    filenames <- .Primitive("c")(filenames, auxnames)

    n <- length(filenames)

    ## error if no file name !
    if (n == 0)
      stop("No file name given !")

    read.celfiles(filenames,
                  phenoData=phenoData,
                  description=description,
                  notes=notes,
                  verbose = verbose,
                  rm.mask = rm.mask,
                  rm.outliers= rm.outliers,
                  rm.extra=rm.extra)
}

list.celfiles <-   function(...){
    files <- list.files(...)
    return(files[grep("\\.[cC][eE][lL]\\.gz$|\\.[cC][eE][lL]$", files)])
}

