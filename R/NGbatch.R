# Loads information contained in XYS files / NimbleGen
# Author: Benilton Carvalho
# Date: Apr/2005

NGBatch <- function(path=".", loadDesign = TRUE){
  
  info <- readxys(pathxys=path)
  xys <- info$xys
  designname <- info$designname
  xysfiles <- info$xysfiles
  rm(info)

  if(loadDesign){
    nameVar <- paste("ng",designname,sep="")
    ndfFile <- paste(designname,".ndf",sep="")
    assign(nameVar,makeNDFenv(ndfFile),envir=.GlobalEnv)
  }

#  rownames(xys) <- get("feature_names", envir = slot(get(nameVar),"featureInfo"))
  columns <- c(1,2)
  new("oligoBatch",
      eList = list(exprs = as.matrix(xys[,-columns])),
      designName = nameVar,
      manufacturer = "NimbleGen")
 }
