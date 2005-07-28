# Methods
# Author: Benilton Carvalho
# Date: April 2005
setClass("oligoBatch",
         representation(manufacturer="character",
                        platform="character"),
         contains="eSet")

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


## BC: Thu, Jul 28, 2005 - is there a smarter way of creating a nickname?
if (is.null(getGeneric("getPD"))){
  setGeneric("getPD",
             function(object) standardGeneric("getPD"))}
setMethod("getPD","oligoBatch", function(object){
  getPlatformDesign(object)})

## probeNames - returns probeNames for PMs ... genenames ignored for now
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

## BC: indexFeatureSetName - to simplify the procedure of bringing pms/mms with a given name
##     it'll return the indexes for a given feature_set_name
if( is.null(getGeneric("indexFeatureSetName") ))
  setGeneric("indexFeatureSetName", function(object, ...)
             standardGeneric("indexFeatureSetName"))

setMethod("indexFeatureSetName","oligoBatch",
          function(object, featurenames){
            tmp <- NULL
            for (i in featurenames)
              tmp <- c(tmp,which(getPlatformDesign(object)$feature_set_name == i))
            return(sort(tmp))
          })


##PM methods
if( is.null(getGeneric("pm") ))
  setGeneric("pm", function(object, ...)
             standardGeneric("pm"))

## BC: Th, Jul 28, 2005 - This is an attempt of
##     selecting the feature_set_name for pm/mm
##     if it is working, we need to rename the variables
setMethod("pm","oligoBatch", ##genenames is ignored for now.. we will get to it
          function(object, genenames=NULL){
            index <- pmindex(object)
            if (!is.null(genenames))
              index <- intersect(index,indexFeatureSetName(object,genenames))
            return(exprs(object)[index,,drop=FALSE])
          })

### Trying to make the above to work
### if it is working, we remove the one below
###setMethod("pm","oligoBatch", ##genenames is ignored for now.. we will get to it
###          function(object, genenames=NULL){
###            return(exprs(object)[pmindex(object),,drop=FALSE])
###          })

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

## BC: Th, Jul 28, 2005 - This is an attempt of
##     selecting the feature_set_name for pm/mm
##     if it is working, we need to rename the variables
setMethod("mm","oligoBatch", function(object, genenames=NULL){
            index <- mmindex(object)
            if (!is.null(genenames))
              index <- intersect(index,indexFeatureSetName(object,genenames))
            return(exprs(object)[index,,drop=FALSE])
          })


##setMethod("mm","oligoBatch", ##genenames is ignored for now.. we will get to it
##          function(object, genenames=NULL){
##            return(exprs(object)[mmindex(object),,drop=FALSE])
##          })
          
if( is.null(getGeneric("mm<-") ))
  setGeneric("mm<-", function(object, value)
             standardGeneric("mm<-"))

setReplaceMethod("mm", "oligoBatch",
                 function(object, value){
                   exprs(object)[mmindex(object),] <- value
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

## Histogram
if( is.null(getGeneric("hist")) )
  setGeneric("hist")

setMethod("hist",signature(x="oligoBatch"), function(x,...) plotDensity.oligoBatch(x,...))
