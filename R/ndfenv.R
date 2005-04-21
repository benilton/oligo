# NDF Environment
# Author: Benilton Carvalho
# Date: Mar / 2005

ndfenv <- function(ndffile){
  
  cat(paste("Reading",ndffile,"\n"))
  ndf <- readndf(ndffile)
  ndf$index <- max(ndf$Y)*ndf$X + ndf$Y
  envndf <- new.env()

  vars <- names(ndf)
  for (i in 1:length(vars)){
    assign(vars[i],ndf[,i],envndf)
  }
  
  # Creating class "platformDesign" with 3 slots (environment, maker and seed)
  # Using seed, we can have a 1:1 mapping between index and (X,Y)
  # index = seed * X + Y
  # NDF information
  # type: expression/SNP
  
  setClass("platformDesign",
           representation(featureInfo = "environment",
                          manufacturer = "character",
                          seed = "numeric",
                          type = "character"))
  new("platformDesign",
      featureInfo=envndf,
      manufacturer="NimbleGen",
      seed=max(ndf$Y),
      type="expression")
  
}
