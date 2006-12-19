#############################################################
##
## read.affybatch.R
##
## Adapted by B. M. Bolstad from read.affybatch in the affy
## package version 1.2.  The goal is a faster, less memory hungry
## ReadAffy. To do this we will shunt more work off to
## the c code.
##
## History
## Jun 13-15 Intial version
## Jun 16    Verbose flag passed to C routine
## Jun 17    New method for checking header of first cel
##           file.
## Jul 7     Added the function read.probematrix which
##           reads in PM, MM or both into matrices
## Sep 28    changed name from read.affybatch2 to read.affybatch
##           and cleaned up some old commented stuff
## Apr 13, 2004 - fixed problem in read.probematrix
##
## May 28, 2005 - Moving to oligo... (BC)
## May 29, 2005 - It's working partially:
##                if I try to read multiple CEL, it fails on sampleNames
## Nov 15, 2005 - Started changing the order of the intensities
##                after reading the CEL files, so pm/mm will be faster. (BC)
## Jan 09, 2006 - added feature so power-users can choose to NOT load PDEnv
##                added feature to return SD
##                added feature to return number of pixels (BC)
#############################################################
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
                          rm.extra=FALSE) {

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

    ## RI: We should check if the pd package is available here. If not try
    ## to install it

    ## load pdInfo
    if (is.null(pkgname))
      pkgname <- cleanPlatformName(ref.cdfName)
    library(pkgname, character.only=TRUE)
    message("Platform design info loaded.")
    arrayType <- kind(get(pkgname))

    tmpExprs <- .Call("read_abatch", as.list(filenames),
                      rm.mask, rm.outliers, rm.extra, ref.cdfName,
                      dim.intensity, verbose, PACKAGE="affyio")

    tmpNP <- tmpSD <- NULL

    if (sd){
        tmpSD <- .Call("read_abatch_stddev", as.list(filenames),
                       rm.mask, rm.outliers, rm.extra, ref.cdfName,
                       dim.intensity, verbose, PACKAGE="affyio")
    }

    if (npixels){
        tmpNP <- .Call("read_abatch_npixels", as.list(filenames),
                       rm.mask, rm.outliers, rm.extra, ref.cdfName,
                       dim.intensity, verbose, PACKAGE="affyio")
    }

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
    tmpExprs <- tmpExprs[order_index,, drop=FALSE]
    rownames(tmpExprs) <- 1:nrow(tmpExprs)

    if (sd){
        tmpSD <- tmpSD[order_index,, drop=FALSE]
        rownames(tmpSD) <- featureSetNames(pdInfo)
    }

    if (npixels){
        tmpNP <- tmpNP[order_index,, drop=FALSE]
        rownames(tmpNP) <- featureSetNames(pdInfo)
    }


    colnames(tmpExprs) <- sapply(strsplit(colnames(tmpExprs), "/"),
                                 function(x) x[length(x)])

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
               ##            sd=tmpSD,
               ##            npixels=tmpNP,
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

