# Methods
# Author: Benilton Carvalho
# Date: April 2005

require(Biobase)
setClass("oligoBatch",
         representation(designName="character",
                        manufacturer="character",
                        nrow="numeric",
                        ncol="numeric",
                        platform="character",
                        indextable="matrix"),
         contains="eSet")

### accessors


if (is.null(getGeneric("allIndex")))
  setGeneric("allIndex", function(obj) standardGeneric("allIndex"))

setMethod("allIndex","oligoBatch",
          function(obj){
            get("index", envir = featureInfo(obj))
          })

if (is.null(getGeneric("designName")))
  setGeneric("designName",function(obj) standardGeneric("designName"))
setMethod("designName","oligoBatch", function(obj) obj@designName)

if (is.null(getGeneric("featureInfo")))
  setGeneric("featureInfo", function(obj) standardGeneric("featureInfo"))
setMethod("featureInfo","oligoBatch",
          function(obj){
            design <- designName(obj)
            slot(get(design),"featureInfo")
          })

if (is.null(getGeneric("featureType")))
  setGeneric("featureType", function(obj) standardGeneric("featureType"))

setMethod("featureType","oligoBatch",
          function(obj){
            allindex <- allIndex(obj)
            elem <- match(index(obj),allindex)
            types <- get("feature_type", envir = featureInfo(obj))[elem]
            return(types)
          })

if (is.null(getGeneric("index")))
  setGeneric("index", function(obj) standardGeneric("index"))

setMethod("index","oligoBatch",
          function(obj){
            obj@indextable[,3]
          })

if (is.null(getGeneric("featureNames")))
  setGeneric("featureNames", function(obj) standardGeneric("featureNames"))

setMethod("featureNames","oligoBatch",
          function(obj){
            allindex <- allIndex(obj)
            elem <- match(index(obj),allindex)
            allfnames <- get("feature_names", envir = featureInfo(obj))[elem]
          })

### PM
if (is.null(getGeneric("pm")))
  setGeneric("pm",function(obj,...) standardGeneric("pm"))
setMethod("pm","oligoBatch",
          function(obj){
            elem <- match(index(obj)[featureType(obj) == "PM"],index(obj))
            pms <- exprs(obj)[elem,]
          }
          )

### MM
if (is.null(getGeneric("mm")))
  setGeneric("mm",function(obj,...) standardGeneric("mm"))
setMethod("mm","oligoBatch",
          function(obj){
            elem <- match(index(obj)[featureType(obj) == "MM"],index(obj))
            mms <- exprs(obj)[elem,]
          }
          )
