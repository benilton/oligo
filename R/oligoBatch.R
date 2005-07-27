# Methods
# Author: Benilton Carvalho
# Date: April 2005
setClass("oligoBatch",
         representation(manufacturer="character",
                        platform="character"),
         contains="eSet")

##simple accessors that needs to change!
if (is.null(getGeneric("annotation")))
  setGeneric("annotation",function(object) standardGeneric("annotation"))
setMethod("annotation","oligoBatch", function(object) object@platform)

###same as above...but actual accessor
if (is.null(getGeneric("platform")))
  setGeneric("platform",function(object) standardGeneric("platform"))
setMethod("platform","oligoBatch", function(object) object@platform)

###for compatibility with previous package
setMethod("length",signature(x="oligoBatch"),
          function(x) ncol(exprs(x))) 

###this might change. we might put pd at the end
## BC: May 29 - pd added at the end
## BC: Jul 13 - pd moved to the begining
if (is.null(getGeneric("platformDesignName"))){
  setGeneric("platformDesignName",
             function(object) standardGeneric("platformDesignName"))}

setMethod("platformDesignName","oligoBatch", function(object){
  cleanPlatformName(object@platform)})

##loading the library for now... this must change
if (is.null(getGeneric("getPlatformDesign"))){
  setGeneric("getPlatformDesign",
             function(object) standardGeneric("getPlatformDesign"))}

setMethod("getPlatformDesign","oligoBatch", function(object){
  pdn <- platformDesignName(object)
  library(pdn,character.only=TRUE)
  return(get(pdn,pos=paste("package:",pdn,sep="")))
})

###probeNames - returns probeNames for PMs ... genenames ignored for now
if (is.null(getGeneric("probeNames")))
  setGeneric("probeNames", function(object, ...)
             standardGeneric("probeNames"))

setMethod("probeNames", "oligoBatch",
          function(object, genenames=NULL){
            pmIndex <- pmindex(getPlatformDesign(object))
            pns <- get("feature_set_name",envir=featureInfo(getPlatformDesign(object)))
            return(as.character(pns[pmIndex]))
          })

###geneNames - returns geneNames for PMs
if (is.null(getGeneric("geneNames")))
  setGeneric("geneNames", function(object)
             standardGeneric("geneNames"))

setMethod("geneNames", "oligoBatch",
          function(object){
            pmIndex <- pmindex(getPlatformDesign(object))
            return(levels(factor(get("feature_set_name",envir=featureInfo(getPlatformDesign(object)))[pmIndex])))
          })


##pmindex method for oligoBatch
if( is.null(getGeneric("pmindex")))
  setGeneric("pmindex", function(object,...)
             standardGeneric("pmindex"))

##WE assume feature_type_1 is PM or MM. this might change with other platforms
setMethod("pmindex", "oligoBatch",
          function(object){
            pmindex(getPlatformDesign(object))
          })

##mmindex method for oligoBatch
if( is.null(getGeneric("mmindex")))
  setGeneric("mmindex", function(object,...)
             standardGeneric("mmindex"))

##WE assume feature_type_1 is PM or MM. this might change with other platforms
setMethod("mmindex", "oligoBatch",
          function(object){
            mmindex(getPlatformDesign(object))
          })


##PM methods
if( is.null(getGeneric("pm") ))
  setGeneric("pm", function(object, ...)
             standardGeneric("pm"))

setMethod("pm","oligoBatch", ##genenames is ignored for now.. we will get to it
          function(object, genenames=NULL){
            return(exprs(object)[pmindex(object),,drop=FALSE])
          })

if( is.null(getGeneric("pm<-") ))
  setGeneric("pm<-", function(object, value)
             standardGeneric("pm<-"))

setReplaceMethod("pm", "oligoBatch",
                 function(object, value){
                   exprs(object)[pmindex(object),] <- value
                   object
                 })
##MM methods... designed for arrays that have one MM per PM
if( is.null(getGeneric("mm") ))
  setGeneric("mm", function(object, ...)
             standardGeneric("mm"))

setMethod("mm","oligoBatch", ##genenames is ignored for now.. we will get to it
          function(object, genenames=NULL){
            return(exprs(object)[mmindex(object),,drop=FALSE])
          })
          
if( is.null(getGeneric("mm<-") ))
  setGeneric("mm<-", function(object, value)
             standardGeneric("mm<-"))

setReplaceMethod("mm", "oligoBatch",
                 function(object, value){
                   exprs(object)[mmindex(object),] <- value
                   object
                 })

###sampleNames and description should go away once eSet has them
###COMMENT OUT TO SEE IF BIOBASE;s IS WORKING. IF IT IS WE GET RID OF IT
# if( !isGeneric("sampleNames") )
#   setGeneric("sampleNames", function(object)
#              standardGeneric("sampleNames"))
# setMethod("sampleNames", "oligoBatch",
#           function(object) {
#             if (! is.null(colnames(exprs(object))))
#               colnames(exprs(object))
#             else
#               row.names(pData(object))
#           })

##description
if( !isGeneric("description") )
  setGeneric("description", function(object)
             standardGeneric("description"))
setMethod("description", "oligoBatch", function(object)
          object@description)

##replace method for description
if( !isGeneric("description<-") )
  setGeneric("description<-", function(object, value)
             standardGeneric("description<-"))

setReplaceMethod("description", "oligoBatch", function(object, value) {
  object@description <- value
  object
})

##notes
if( !isGeneric("notes") )
  setGeneric("notes", function(object)
             standardGeneric("notes"))
setMethod("notes", "oligoBatch", function(object)
          object@notes)

if( !isGeneric("notes<-") )
  setGeneric("notes<-", function(object, value)
             standardGeneric("notes<-"))

setReplaceMethod("notes", "oligoBatch", function(object, value) {
  object@notes <- value
  object
})


## BC: Fri Jul 22, 2005 - I needed this methods today
##     Copied from affy and a few modifications

if(is.null(getGeneric("ncol")))
  setGeneric("ncol")

setMethod("ncol",signature(x="oligoBatch"),
                    function(x) getPlatformDesign(x)@ncol)

if( is.null(getGeneric("nrow")))
    setGeneric("nrow")

  setMethod("nrow",signature(x="oligoBatch"),
                        function(x) getPlatformDesign(x)@nrow)


if( is.null(getGeneric("image")))
  setGeneric("image")

setMethod("image",signature(x="oligoBatch"),
          function(x, transfo=log, col=gray(c(0:64)/64),xlab="",ylab="", ...){
            scn <- prod(par("mfrow"))
            ask <- dev.interactive()
            which.plot <- 0

            ## x.pos <- (1:nrow(x)) - (1 + getOption("BioC")$affy$xy.offset)
            ## y.pos <- (1:ncol(x)) - (1 + getOption("BioC")$affy$xy.offset)

            x.pos <- (1:nrow(x)) - 1
            y.pos <- (1:ncol(x)) - 1

            for(i in 1:length(sampleNames(x))){
              which.plot <- which.plot+1;
              if(trunc((which.plot-1)/scn)==(which.plot-1)/scn && which.plot>1 && ask)  par(ask=TRUE)
              m <- exprs(x)[,i]
              if (is.function(transfo)) {
                m <- transfo(m)
              }
              m <- as.matrix(rev(as.data.frame(matrix(m, nrow=length(x.pos), ncol=length(y.pos)))))
              image(x.pos, y.pos, m,
                    col=col, main=sampleNames(x)[i],
                    xlab=xlab, ylab=ylab,,xaxt='n',
                      yaxt='n', ...)
              par(ask=FALSE)
            }
          })


## Boxplot
if( is.null(getGeneric("boxplot")))
  setGeneric("boxplot")

setMethod("boxplot",signature(x="oligoBatch"),
          function(x,which=c("both","pm","mm"),range=0,...){
            which <- match.arg(which,c("both","pm","mm"))
            tmp <- description(x)
            if (is(tmp, "MIAME")) main <- tmp@title

            tmp <- unlist(featureIndex(x,which))
            tmp <- tmp[seq(1,length(tmp),len=5000)]

            boxplot(data.frame(log2(exprs(x)[tmp,])),main=main,range=range, ...)
          })



## BC: Mon Jul 25, 2005 - added featureIndex, before called indexProbe
##     so we can produce some plots.

if( is.null(getGeneric("featureIndex") ))
  setGeneric("featureIndex", function(object, ...)
             standardGeneric("featureIndex"))

setMethod("featureIndex","oligoBatch", ##genenames is ignored for now.. we will get to it
          function(object, which=c("both","pm","mm"), genenames=NULL){
            which <- match.arg(which,c("both","pm","mm"))
            if (which=="both"){
              pmIndex <- pmindex(getPlatformDesign(object))
              mmIndex <- mmindex(getPlatformDesign(object))
              indexes <- sort(c(pmIndex,mmIndex))
            } else if (which=="pm"){
              indexes <- sort(pmindex(getPlatformDesign(object)))
            } else if (which=="mm"){
              indexes <- sort(mmindex(getPlatformDesign(object)))
            }
             return(indexes)
           })


## Histogram
if( is.null(getGeneric("hist")) )
  setGeneric("hist")

setMethod("hist",signature(x="oligoBatch"), function(x,...) plotDensity.AffyBatch(x,...))

## BC: Mon, Jul 25, 2005 - Plot density - from affy
plotDensity <- function(mat,
                        ylab="density", xlab="x", type="l", col=1:6,
                        ...) {
  
  x.density <- apply(mat, 2, density)

  all.x <- do.call("cbind", lapply(x.density, function(x) x$x))
  all.y <- do.call("cbind", lapply(x.density, function(x) x$y))
  
  matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, ...)

  invisible(list(all.x=all.x, all.y=all.y))
}
 

plotDensity.AffyBatch <- function(x, col=1:6, log=TRUE,
                                  which=c("both","pm","mm"),
                                  ylab="density",
                                  xlab=NULL,
                                  ...){

  which <- match.arg(which,c("both","pm","mm"))

  Index <- unlist(featureIndex(x,which))
  
  x <- exprs(x)[Index, ,drop=FALSE]
  
  if(log){
    x <- log2(x)
    if(is.null(xlab)) xlab <- "log intensity"
  }
  else  if(is.null(xlab)) xlab <- "intensity"
  
  rv <- plotDensity(x, ylab=ylab, xlab=xlab, col=col, ...)

  invisible(rv)
}
