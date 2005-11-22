# Reading NimbleGen Data
# Author: Benilton Carvalho
# Last Modification: May 28, 2005

list.xysfiles <-   function(...){
  files <- list.files(...)
  return(files[grep("\.[xX][yY][sS]\.gz$|\.[xX][yY][sS]$", files)])
}

readxysHeader <- function(filename) scan(filename,nlines=1,quiet=TRUE, what=character(0))

readonexysfile <- function(filename){
  types <- list(X = numeric(0), Y = numeric(0), SIGNAL = numeric(0), COUNT = numeric(0))
  header <- readxysHeader(filename)
  whatToRead <- types[match(header,names(types))]
  sig <- scan(filename, sep = "\t", skip = 2, quiet = T, what = whatToRead)$SIGNAL
  return(sig)
}


stuffForXYSandCELreaders <- function(filenames,
                                     phenoData=new("phenoData"),
                                     description=NULL,
                                     notes="",
                                     verbose = FALSE) {
  
  n <- length(filenames)
  
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
  
  tmp <- stuffForXYSandCELreaders(filenames,phenoData,description,notes,verbose)
  filenames <- tmp$filenames
  ## Create space to store the design names
  designnamelist <- NULL

  ## Read all XYS files, get design name for each
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

  ## Allocate memory for the intensities
  ## And giving the correct names for the columns
  e <- matrix(NA,nProbes(get(designname)),length(filenames))
  colnames(e) <- tmp$samplenames
  for (i in seq(along=filenames)){
    if (verbose) cat(i, "reading",filenames[1],"...")
    e[,i] <- readonexysfile(filenames[i])
    if(verbose) cat("Done.\n")
  }

  ## Put intensities in the same order as
  ## in the PDenv, this way, things get faster
  order_index <- get(designname,pos=paste("package:",designname,sep=""))$order_index
  e <- e[order_index,,drop=FALSE]

  return(new("oligoBatch",
	     assayData=list(exprs=e),
             sampleNames=rownames(pData(tmp$phenoData)),
             platform = designnamelist[1],
             manufacturer = "NimbleGen",
             phenoData=tmp$phenoData,
             description=tmp$description,
             notes=notes))
}
