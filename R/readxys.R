# Reading NimbleGen Data
# Author: Benilton Carvalho
# Last Modification: May 28, 2005

list.xysfiles <-   function(...){
  files <- list.files(...)
  return(files[grep("\.[xX][yY][sS]\.gz$|\.[xX][yY][sS]$", files)])
}

readxysHeader <- function(filename) scan(filename,nlines=1,quiet=TRUE, what=character(0))

##readonexysfile <- function(filename){
##  types <- list(X = numeric(0), Y = numeric(0), SIGNAL = numeric(0), COUNT = numeric(0))
##  header <- readxysHeader(filename)
##  whatToRead <- types[match(header,names(types))]
##  out <- scan(filename, sep = "\t", skip = 2, quiet = T, what = whatToRead)
##  return(out)
##}

readonexysfile <- function(filename)
  read.delim(filename)

stuffForXYSandCELreaders <- function(filenames,
                                     phenoData=new("phenoData"),
                                     description=NULL,
                                     notes="",
                                     verbose = FALSE,
                                     nwells = 1,
                                     designname=NULL) {
  
  nfiles <- length(filenames)
  n <- nfiles*nwells
  
  ## error if no file name !
  if (n == 0)
    stop("No file name given !")
  if(mode(filenames)!="character")
    stop("filenames must be of type character!")
  
  pdata <- pData(phenoData)
  if(dim(pdata)[1] != n) {
    ##if empty pdata filename are samplenames
    cat("Incompatible phenoData object. Created a new one.\n")
    samplenames <- sub("^/?([^/]*/)*", "", unlist(filenames), extended=TRUE)

    if(nwells>1){
      wells <- paste(".",as.character(unique(get(designname)@lookup$container)),sep="")
      samplenames <- as.character(t(outer(samplenames,wells,paste,sep="")))
    }
    
    pdata <- data.frame(sample=1:n, row.names=samplenames)
    phenoData <- new("phenoData",pData=pdata,varLabels=list(sample="arbitrary numbering"))
  }
  else samplenames <- rownames(pdata)
  
  if (is.null(description))
    {
      description <- new("MIAME")
      description@preprocessing$filenames <- filenames
      description@preprocessing$oligoversion <- library(help="oligo")$info[[2]][2]
    }
  
  return(list(filenames=filenames,samplenames=samplenames,phenoData=phenoData,description=description))
}

read.xysfiles <- function(filenames,
                         phenoData=new("phenoData"),
                         description=NULL,
                         notes="",
                         verbose = FALSE) {
  
  ## Create space to store the design names
  designnamelist <- NULL

  ## Read the header for all XYS files, get design name for each
  for (xysfile in filenames){
    firstline <- readxysHeader(xysfile)
    designname <- unlist(strsplit(firstline[grep("designname",firstline,fixed=TRUE,useBytes=TRUE)],"="))[2]
    designnamelist <- rbind(designnamelist,designname)
  }

  ## How many different designs?
  numberdesigns <- length(unique(designnamelist))

  ## All XYS files should point to one NDF file
  if(numberdesigns != 1)
    stop("XYS Files do not refer to the same design!")
    
  ## Load PDenv for the XYS files
  designname=cleanPlatformName(designnamelist[1])
  library(designname,character.only=TRUE)

  ## Get the number of wells
  nwells <- get(designname)@nwells

  tmp <- stuffForXYSandCELreaders(filenames,phenoData,description,notes,verbose,nwells,designname)

  ## Allocate memory for the intensities
  ## And giving the correct names for the columns
  e <- matrix(NA, nrow = nProbes(get(designname)), ncol = (length(filenames)*nwells))
  colnames(e) <- tmp$samplenames

  ## Loading lookup table to correctly assign the wells
  lookup <- get(designname)@lookup
  j <- 1
  for (i in seq(along=filenames)){
    if (verbose) cat(i, "reading",filenames[i],"...")

    ## Read XYS "as is"
    tmpE <- readonexysfile(filenames[i])
    tmpE$index <- tmpE$X + (tmpE$Y-1)*get(designname)@ncol

    ## get those probes which are in the NDF
    ## in that order
#    inndf <- tmpE$index %in% lookup$index
#    tmpE <- tmpE[inndf,]

    ndforder <- match(lookup$index,tmpE$index)
    e[,j:(j+nwells-1)] <- matrix(tmpE$SIGNAL[ndforder],
                                 ncol = nwells,
                                 byrow = TRUE)       
    j <- j+nwells
    if(verbose) cat("Done.\n")
  }

  rownames(e) <- as.character(get(designname)@featureInfo$feature_set_name)
  return(new("oligoBatch",
	     assayData=list(exprs=e[,,drop=FALSE]),
             sampleNames=rownames(pData(tmp$phenoData)),
             platform = designnamelist[1],
             manufacturer = "NimbleGen",
             phenoData=tmp$phenoData,
             description=tmp$description,
             notes=notes))
}
