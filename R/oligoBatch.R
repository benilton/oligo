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
##     This is used by pm/mm when getting intensities for given feature_set_names

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

## Show
setMethod("show","oligoBatch", function(object){
  dm <-dim(exprs(object))
  nfeatures <- dm[1]
  nsamples <- dm[2]
  cat("oligoBatch with \n\t", nfeatures, " features\n\t", sep="")
  cat(nsamples, "samples\n\t")
  show(phenoData(object))
})
