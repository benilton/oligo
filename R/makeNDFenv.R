# NDF Environment
# Author: Benilton Carvalho
# Date: Mar / 2005

makeNDFenv <- function(ndffile){
  
  cat(paste("Reading",ndffile,"\n"))
  ndf <- readndf(ndffile)
  envndf <- new.env()

  vars <- names(ndf)
  for (i in 1:length(vars)){
    assign(vars[i],ndf[,i],envndf)
  }
  
  setClass("platformDesign",
           representation(featureInfo = "environment",
                          manufacturer = "character",
                          type = "character"))
  new("platformDesign",
      featureInfo=envndf,
      manufacturer="NimbleGen",
      type="expression")  
}
