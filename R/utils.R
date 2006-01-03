cleanPlatformName <- function(x)
  gsub("[_-]","",paste("pd",tolower(x),sep=""))

i2xy <- function(i,obatch){
  xy <- mget(c("X","Y"),envir=featureInfo(getPD(obatch)))
  return(cbind(xy$X[i],xy$Y[i]))
}

xy2i <- function(x,y,obatch){
  xy <- mget(c("X","Y"),envir=featureInfo(getPD(obatch)))
  xy1 <- xy$X+(xy$Y-1)*max(xy$X)
  xy2=x+(y-1)*max(xy$X)
  match(xy2,xy1)
}

