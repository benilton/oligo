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
#############################################################


read.celfiles <- function(filenames,
                          phenoData=new("phenoData"),
                          description=NULL,
                          notes="",
                          verbose = FALSE,
                          compress= FALSE,
                          rm.mask = FALSE, rm.outliers=FALSE, rm.extra=FALSE){

  tmp <- stuffForXYSandCELreaders(filenames,phenoData,description,notes,verbose)

  filenames <- tmp$filenames

  n <- length(filenames)

  ## error if no file name !
  if (n == 0)
    stop("No file name given !")

  ## read the first file to see what we have
  headdetails <- .Call("ReadHeader",filenames[1], compress, PACKAGE="oligo")

  ##now we use the length
  dim.intensity <- headdetails[[2]]
  ##and the cdfname as ref
  ref.cdfName <- headdetails[[1]]
  
  out <- new("oligoBatch",
             eList=new("exprList",
               .Data=list(exprs=.Call("read_abatch",as.list(filenames),
                            compress, rm.mask,
                            rm.outliers, rm.extra, ref.cdfName,
                            dim.intensity, verbose, PACKAGE="oligo")),
               eMetadata=data.frame()),
             platform = ref.cdfName,
             manufacturer = "Affymetrix",
             phenoData=tmp$phenoData,
             description=tmp$description,
             notes=notes)
  colnames(out@eList$exprs) <- sampleNames(out)
  return(out)
}

read.affybatch <- function(..., filenames=character(0),
                           phenoData=new("phenoData"),
                           description=NULL,
                           notes="",
                           compress = getOption("BioC")$affy$compress.cel,
                           rm.mask = FALSE, rm.outliers=FALSE, rm.extra=FALSE,
                           verbose = FALSE) {

  warning("read.affybatch is depricated. Please use read.xysfile.")

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







######################################################################################

read.probematrix <- function(..., filenames = character(0), phenoData = new("phenoData"),
    description = NULL, notes = "", compress = getOption("BioC")$affy$compress.cel,
    rm.mask = FALSE, rm.outliers = FALSE, rm.extra = FALSE, verbose = FALSE,which="pm"){

  auxnames <- as.list(substitute(list(...)))[-1]
  filenames <- .Primitive("c")(filenames, auxnames)
  
  which <- match.arg(which,c("pm","mm","both"))
  
  headdetails <- .Call("ReadHeader", filenames[[1]], compress, PACKAGE="oligo")
  dim.intensity <- headdetails[[2]]
  ref.cdfName <- headdetails[[1]]

  Data <- new("oligoBatch", cdfName = ref.cdfName, annotation = cleancdfname(ref.cdfName,addcdf = FALSE))

  cdfInfo<- as.list(getCdfInfo(Data))
  cdfInfo <- cdfInfo[order(names(cdfInfo))]


  .Call("read_probeintensities", filenames,
        compress, rm.mask, rm.outliers, rm.extra, ref.cdfName,
        dim.intensity, verbose, cdfInfo,which, PACKAGE="oligo")
}


list.celfiles <-   function(...){
  files <- list.files(...)
  return(files[grep("\.[cC][eE][lL]\.gz$|\.[cC][eE][lL]$", files)])
}

AllButCelsForReadAffy <- function(..., filenames=character(0),
                                  widget=FALSE,
                                  celfile.path=NULL,
                                  sampleNames=NULL,
                                  phenoData=NULL,
                                  description=NULL){

    ##first figure out filenames
  auxnames <- unlist(as.list(substitute(list(...)))[-1])

  if (widget){
    require(tkWidgets)
    widgetfiles <- fileBrowser(textToShow="Choose CEL files",
                               testFun=hasSuffix("[cC][eE][lL]"))
  }
  else{
    widgetfiles <- character(0)
  }

  if(!is.null(celfile.path)){
    auxnames <- file.path(celfile.path, auxnames)
    filenames <- file.path(celfile.path, filenames)
  }

  filenames <- .Primitive("c")(filenames, auxnames, widgetfiles)

  if(length(filenames)==0){
    if(is.null(celfile.path)) celfile.path <- getwd()
    filenames <- list.celfiles(celfile.path,full.names=TRUE)
  }
  if(length(filenames)==0) stop("No cel filennames specified and no cel files in specified directory:",celfile.path,"\n")

  if(is.null(sampleNames)){
    sampleNames <- sub("^/?([^/]*/)*", "", filenames, extended=TRUE)
  }
  else{
    if(length(sampleNames)!=length(filenames)){
      warning("sampleNames not same length as filenames. Using filenames as sampleNames instead\n")
      sampleNames <- sub("^/?([^/]*/)*", "", filenames, extended=TRUE)
    }
  }

  if(is.character(phenoData)) ##if character read file
    phenoData <- read.phenoData(filename=phenoData)
  else{
      if (! is(phenoData, "phenoData")) {
          if(widget){
              require(tkWidgets)
              phenoData <- read.phenoData(sampleNames=sampleNames,widget=TRUE)
          }
          else
              phenoData <- read.phenoData(sampleNames=sampleNames,widget=FALSE)
      }
  }

  sampleNames <- rownames(pData(phenoData))

  ##get MIAME information
  if(is.character(description)){
    description <- read.MIAME(filename=description,widget=FALSE)
  }
  else{
      if (! is(description, "MIAME")) {
          if(widget){
              require(tkWidgets)
              description <- read.MIAME(widget=TRUE)
          }
          else
              description <- new("MIAME")
      }
  }

  ##MIAME stuff
  description@preprocessing$filenames <- filenames
  if(exists("tksn")) description@samples$description <- tksn[,2]
  description@preprocessing$affyversion <- library(help=oligo)$info[[2]][2]

  return(list(filenames   = filenames,
              phenoData   = phenoData,
              sampleNames = sampleNames,
              description = description))
}

###this is user friendly wrapper for read.affybatch
ReadAffy <- function(..., filenames=character(0),
                     widget=FALSE,
                     compress=getOption("BioC")$affy$compress.cel,
                     celfile.path=NULL,
                     sampleNames=NULL,
                     phenoData=NULL,
                     description=NULL,
                     notes="",
                     rm.mask=FALSE, rm.outliers=FALSE, rm.extra=FALSE,
                     verbose=FALSE) {

  l <- AllButCelsForReadAffy(..., filenames=filenames,
                             widget=widget,
                             celfile.path=celfile.path,
                             sampleNames=sampleNames,
                             phenoData=phenoData,
                             description=description)

  ##and now we are ready to read cel files
  ret <- read.affybatch(filenames=l$filenames,
                        phenoData=l$phenoData,
                        description=l$description,
                        notes=notes,
                        compress=compress,
                        rm.mask=rm.mask,
                        rm.outliers=rm.outliers,
                        rm.extra=rm.extra,
                        verbose=verbose)

  ## BC: the following fails
##  sampleNames(ret) <- l$sampleNames

  return(ret)
}
