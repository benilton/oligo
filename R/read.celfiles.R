read.celfiles <- function(filenames,
                          pkgname=NULL,
                          phenoData=new("AnnotatedDataFrame"),
                          featureData=NULL,
                          description=NULL,
                          notes="",
                          verbose = FALSE,
##                          compress= FALSE,
                          rm.mask = FALSE,
                          rm.outliers=FALSE,
                          rm.extra=FALSE,
                          tmpdir=getwd()){

  ## FIXME: this is not an acceptable function name
  tmp <- stuffForXYSandCELreaders(filenames, phenoData, description, notes,
                                  verbose)
  filenames <- tmp$filenames
  n <- length(filenames)
  if (n == 0)
    stop("no file name given")
  
  ## read the first file to see what we have
  ##  headdetails <- .Call("ReadHeader", filenames[1], compress, PACKAGE="affyio")
  headdetails <- read.celfile.header(filenames[1])
  
  ##now we use the length
  dim.intensity <- headdetails[[2]]
  
  ##and the cdfname as ref
  ref.cdfName <- headdetails[[1]]

  if (verbose) cat("Creating objects outside R to store intensities.\n")
  if (verbose) cat("This may take a while... ")
  tmpExprs <- createBufferedMatrix(prod(dim.intensity), 0, directory=tmpdir)
  set.buffer.dim(tmpExprs, 300000, 1)
  if (verbose) cat("Done.", "Now reading CEL files", sep="\n")
  for (i in 1:length(filenames)){
    AddColumn(tmpExprs)
    tmpExprs[,i] <- .Call("read_abatch", filenames[i], rm.mask, rm.outliers,
                          rm.extra, ref.cdfName, dim.intensity, verbose, PACKAGE="affyio")
  }
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
    for (i in 1:ncol(tmpExprs))
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

list.celfiles <-   function(...){
    files <- list.files(...)
    return(files[grep("\\.[cC][eE][lL]\\.gz$|\\.[cC][eE][lL]$", files)])
}

