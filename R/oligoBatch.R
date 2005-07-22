# Methods
# Author: Benilton Carvalho
# Date: April 2005
setClass("oligoBatch",
         representation(manufacturer="character",
                        platform="character",
                        description="MIAME",
                        notes="character"),
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


##PM methods
if( is.null(getGeneric("pm") ))
  setGeneric("pm", function(object, ...)
             standardGeneric("pm"))

setMethod("pm","oligoBatch", ##genenames is ignored for now.. we will get to it
          function(object, genenames=NULL){
             pmIndex <- pmindex(getPlatformDesign(object))
             return(exprs(object)[pmIndex,,drop=FALSE])
           })
          
if( is.null(getGeneric("pm<-") ))
  setGeneric("pm<-", function(object, value)
             standardGeneric("pm<-"))

setReplaceMethod("pm", "oligoBatch",
                 function(object, value){
                   pmIndex <- pmindex(getPlatformDesign(object))
                   exprs(object)[pmIndex,] <- value
                   object
                 })
##MM methods... designed for arrays that have one MM per PM
if( is.null(getGeneric("mm") ))
  setGeneric("mm", function(object, ...)
             standardGeneric("mm"))

setMethod("mm","oligoBatch", ##genenames is ignored for now.. we will get to it
          function(object, genenames=NULL){
             mmIndex <- mmindex(getPlatformDesign(object))
             return(exprs(object)[mmIndex,,drop=FALSE])
           })
          
if( is.null(getGeneric("mm<-") ))
  setGeneric("mm<-", function(object, value)
             standardGeneric("mm<-"))

setReplaceMethod("mm", "oligoBatch",
                 function(object, value){
                   mmIndex <- mmindex(getPlatformDesign(object))
                   exprs(object)[mmIndex,] <- value
                   object
                 })



###sampleNames and description should go away once eSet has them
if( !isGeneric("sampleNames") )
    setGeneric("sampleNames", function(object)
               standardGeneric("sampleNames"))
  setMethod("sampleNames", "oligoBatch",
            function(object) {
              if (! is.null(colnames(exprs(object))))
                colnames(exprs(object))
              else
                row.names(pData(object))
            })

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


###boxplot
if( is.null(getGeneric("boxplot")))
  setGeneric("boxplot")

setMethod("boxplot",signature(x="oligoBatch"),
          function(x,which="both",range=0,...){
            tmp <- description(x)
            if (is(tmp, "MIAME")) main <- tmp@title

            tmp <- unlist(indexProbes(x,which))
            tmp <- tmp[seq(1,length(tmp),len=5000)]

            boxplot(data.frame(log2(intensity(x)[tmp,])),main=main,range=range, ...)
          })

###hist
###if (debug.affy123) cat("--->hist\n")
###
###if( is.null(getGeneric("hist")) )
###  setGeneric("hist")
###
###setMethod("hist",signature(x="AffyBatch"), function(x,...) plotDensity.AffyBatch(x,...))
