# Reading NimbleGen Data
# Author: Benilton Carvalho
# Last Modification: Apr / 2005

readxys <- function(xysfilelist=NULL, pathxys="."){

  # Check list of XYS files
  if(is.null(xysfilelist)){
    xysfilelist <- dir(path=pathxys,pattern=".(xys|XYS)$")
  }
  
  # Create space to store the design names
  designnamelist <- NULL
  
  # Read all XYS files, get design name for each
  for (xysfile in xysfilelist){
    firstline <- scan(xysfile,sep="\t",nlines=1,quiet=T, what=character(0))
    designname <- unlist(strsplit(firstline[grep("designname",firstline)],"="))[2]
    designnamelist <- rbind(designnamelist,designname)
  }

  # How many different designs?
  numberdesigns <- length(unique(designnamelist))
  rm(designnamelist)

  # All XYS files should point to one NDF file
  if(numberdesigns != 1){
    stop("XYS Files do not refer to the same design!")
  }
  else{
    first <- TRUE
    types <- list(X = "numeric", Y="numeric",SIGNAL="numeric",COUNT="numeric")
    samplenames <- c("X","Y")
    for (xysfile in xysfilelist){
      header <- scan(xysfile,
                     sep = "\t",
                     nlines = 1,
                     skip = 1,
                     quiet = T,
                     what = character(0))
      whatToRead <- types[match(header,names(types))]

      cat(paste("Reading",xysfile,"\n"))
      xys <- scan(xysfile,
                  sep = "\t",
                  skip = 2,
                  quiet = T,
                  what = whatToRead)
      xys$X <- as.integer(xys$X)
      xys$Y <- as.integer(xys$Y)
      xys$SIGNAL <- as.numeric(xys$SIGNAL)
      xys$COUNT <- NULL
      samplenames <- c(samplenames,unlist(strsplit(xysfile,"\\."))[1])

      # Assuming the XYS files come all in the same order
      # and all have the same set of (identical) features
      if (first){
        signal <- as.data.frame(xys)
        first <- FALSE
      }
      else{
        signal <- cbind(signal,xys$SIGNAL)
      }
    }
    names(signal) <- samplenames
    return(list(xys=signal,designname=designname,xysfiles=xysfilelist))
  }
}
