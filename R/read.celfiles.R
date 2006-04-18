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
                          pdenv=TRUE,
                          arrayType=NULL,
                          pkgname=NULL,
                          sd=FALSE,
                          npixels=FALSE,
##                          phenoData=new("phenoData"),
                          phenoData=new("AnnotatedDataFrame"),
                          description=NULL,
                          notes="",
                          verbose = FALSE,
                          compress= FALSE,
                          rm.mask = FALSE,
                          rm.outliers=FALSE,
                          rm.extra=FALSE){

  if (!pdenv & is.null(arrayType))
    stop("You chose not to load the pdenv, so you are required to define arrayType (SNP/expression)")

  if (!is.null(pkgname))
    pdenv <- FALSE
  
  tmp <- stuffForXYSandCELreaders(filenames, phenoData, description, notes, verbose)
  filenames <- tmp$filenames
  n <- length(filenames)
  if (n == 0)
    stop("No file name given !")

  ## read the first file to see what we have
  headdetails <- .Call("ReadHeader", filenames[1], compress, PACKAGE="affyio")

  ##now we use the length
  dim.intensity <- headdetails[[2]]

  ##and the cdfname as ref
  ref.cdfName <- headdetails[[1]]

  ## RI: WE SHOULD CHECK IF THE PD PACKAGE IS AVAILABLE HERE. IF not TRY TO
  ## INSTALL IT

  ## if pdenv is to be loaded or user specifies alt.pdenv
  if (pdenv | !is.null(pkgname)){
    if (is.null(pkgname))
      pkgname <- cleanPlatformName(ref.cdfName)
    cat(paste("Loading", pkgname, "\n"))
    library(pkgname, character.only=TRUE)
    cat("Package loaded.\n")
    if (is.null(arrayType))
      arrayType <- get(pkgname)@type
  }
  
  tmpExprs <- .Call("read_abatch", as.list(filenames), 
                    rm.mask, rm.outliers, rm.extra, ref.cdfName,
                    dim.intensity, verbose, PACKAGE="affyio")

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
  
  if (pdenv | !is.null(pkgname)){
    ## BC: Nov 15-16 2005, the PDenv is ordered already in the way
    ##     we want PM/MMs to be. So, we need to reorder the cel
    ##     input, so the pm/mm methods are faster.
    order_index <- get(pkgname, pos=paste("package:", pkgname, sep=""))$order_index
    tmpExprs <- tmpExprs[order_index,, drop=FALSE]
    rownames(tmpExprs) <- as.character(get(pkgname,  pos=paste("package:", pkgname, sep=""))$feature_set_name)

    if (sd){
      tmpSD <- tmpSD[order_index,, drop=FALSE]
      rownames(tmpSD) <- as.character(get(pkgname, pos=paste("package:", pkgname, sep=""))$feature_set_name)
    }

    if (npixels){
      tmpNP <- tmpNP[order_index,, drop=FALSE]
      rownames(tmpNP) <- as.character(get(pkgname, pos=paste("package:", pkgname, sep=""))$feature_set_name)
    }
  }

#  return(list(tmpExprs=tmpExprs, tmp=tmp, ref.cdfName=ref.cdfName))
  if(!sd & !npixels){
    ad <- assayDataNew(storage.mode="list", exprs=tmpExprs)
  }else if(sd & !npixels){
    ad <- assayDataNew(storage.mode="list", exprs=tmpExprs, sd=tmpSD)
  }else if(sd & npixels){
    ad <- assayDataNew(storage.mode="list", exprs=tmpExprs, sd=tmpSD, npixels=tmpNP)
  }else if(!sd & npixels){
    ad <- assayDataNew(storage.mode="list", exprs=tmpExprs, npixels=tmpNP)
  }

##  print(ref.cdfName)
  out <- new("FeatureSet",
             exprs=tmpExprs,
##             assayData=ad,
##             sampleNames=rownames(pData(tmp$phenoData)),
             platform=ref.cdfName,
             manufacturer="Affymetrix",
             phenoData=tmp$phenoData,
             experimentData=tmp$description) #,
  ##           description=tmp$description) #,
  ##           notes=notes)
  platform(out) <- pkgname
  manufacturer(out) <- "Affymetrix"
##  out@platform=ref.cdfName
##  out@manufacturer="Affymetrix"
  return(out)
}

read.affybatch <- function(..., filenames=character(0),
##                           phenoData=new("phenoData"),
                           phenoData=new("AnnotatedDataFrame"),
                           description=NULL,
                           notes="",
                           compress = getOption("BioC")$affy$compress.cel,
                           rm.mask = FALSE, rm.outliers=FALSE, rm.extra=FALSE,
                           verbose = FALSE) {

  warning("read.affybatch is depricated. Please use read.celfiles.")

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
  return(files[grep("\.[cC][eE][lL]\.gz$|\.[cC][eE][lL]$", files)])
}

