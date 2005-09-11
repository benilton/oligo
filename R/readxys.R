# Reading NimbleGen Data
# Author: Benilton Carvalho
# Last Modification: May 28, 2005

list.xysfiles <-   function(...){
  files <- list.files(...)
  return(files[grep("\.[xX][yY][sS]\.gz$|\.[xX][yY][sS]$", files)])
}

readxysHeader <- function(filename) scan(filename,nlines=1,quiet=TRUE, what=character(0))

readonexysfile <- function(filename){

  types <- list(X = "numeric", Y="numeric",SIGNAL="numeric",COUNT="numeric")
  samplenames <- c("X","Y")
  header <-  scan(filename,sep = "\t",nlines = 1,skip = 1,quiet = TRUE,what = character(0))
  whatToRead <- types[match(header,names(types))]

  ##RI: why does this return character
  sig <- scan(filename,
              sep = "\t",
              skip = 2,
              quiet = T,
              what = whatToRead)$SIGNAL
  return(as.numeric(sig))
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
      ## BC: corrected on May 28 (from [[1]][[2]][2] to [[2]][2])
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
  if(numberdesigns != 1){
    stop("XYS Files do not refer to the same design!")
  }
  else{
    ##THIS MUST CHANGE
    ## BC: corrected on May 28, adding "pd" at the end
    ## BC: Jul 13, "pd" at the begining
    designname=cleanPlatformName(designnamelist[1])
    library(designname,character.only=TRUE)
    ##THIS MUST CHANGE...
    ##RI: instead read the first one. figure out size. then
    ##RI: creat matrix e.
      
    e <- matrix(0,nProbes(get(designname)),length(filenames))
    colnames(e) <- tmp$samplenames
    for (i in seq(along=filenames)){
        if (verbose) cat(i, "reading",filenames[1],"...")
        tmp <- readonexysfile(filenames[i])
        if (length(tmp) != nProbes(get(designname))){
          cat("It is likely the NDF has duplicated FEATURE_ID for different control probes.\n")
          cat("You can rebuild the package using makePDpackage(<your NDF>, checkNG = TRUE) for validity check.\n")
          stop("XYS file is not compatible with NDF.\n")
        }
        e[,i] <- tmp
        if(verbose) cat("Done.\n")
    }
  }

  
  colnames(e) <- tmp$samplenames
  return(new("oligoBatch",
#             eList=new("exprList",
#               eList=list(exprs=e),
#               eMetadata=data.frame()),
	     assayData=list(exprs=e), # that is the new line of code
             sampleNames=rownames(pData(tmp$phenoData)),
             platform = designnamelist[1],
             manufacturer = "NimbleGen",
             phenoData=tmp$phenoData,
             description=tmp$description,
             notes=notes))
}






