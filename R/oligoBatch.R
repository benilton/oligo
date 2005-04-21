# Loads information contained in XYS files
# Author: Benilton Carvalho
# Date: Apr/2005

# Methods
# Author: Benilton Carvalho
# Date: April 2005


oligoBatch <- function(path=".", loadDesign = TRUE){
  
  info <- readxys(pathxys=path)
  xys <- info$xys
  designname <- info$designname
  xysfiles <- info$xysfiles
  rm(info)

   if(loadDesign){
     nameVar <- paste("ng",designname,sep="")
     ndfFile <- paste(designname,".ndf",sep="")
     assign(nameVar,ndfenv(ndfFile),envir=.GlobalEnv)
     xys$index <- xys$X * slot(get(nameVar),"seed") + xys$Y
   }

   columns <- c(1,2,ncol(xys))

   new("oligoBatch",
       eList = list(exprs = as.matrix(xys[,-columns])),
       designName = nameVar,
       manufacturer = "NimbleGen",
       indextable = as.matrix(xys[,columns]))
 }
