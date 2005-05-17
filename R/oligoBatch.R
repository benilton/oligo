# Methods
# Author: Benilton Carvalho
# Date: April 2005

require(Biobase)
setClass("oligoBatch",
         representation(designName="character",
                        manufacturer="character",
                        nrow="numeric",
                        ncol="numeric",
                        platform="character"),
         contains="eSet")

if (is.null(getGeneric("designName")))
  setGeneric("designName",function(object) standardGeneric("designName"))
setMethod("designName","oligoBatch", function(object) object@designName)

if (is.null(getGeneric("featureInfo")))
  setGeneric("featureInfo", function(object) standardGeneric("featureInfo"))
setMethod("featureInfo","oligoBatch", function(object) slot(get(designName(object)),"featureInfo"))

if (is.null(getGeneric("featureType")))
  setGeneric("featureType", function(object) standardGeneric("featureType"))
setMethod("featureType","oligoBatch", function(object) get("feature_type", envir = featureInfo(object)))

if (is.null(getGeneric("featureClass")))
  setGeneric("featureClass", function(object) standardGeneric("featureClass"))
setMethod("featureClass","oligoBatch", function(object) get("feature_class", envir = featureInfo(object)))


if (is.null(getGeneric("featureNames")))
  setGeneric("featureNames", function(object) standardGeneric("featureNames"))
setMethod("featureNames","oligoBatch",
          function(object){
            pmfeaturenames <- get("feature_names", envir = featureInfo(object))[featureType(object) == "PM"]
            return(pmfeaturenames[order(pmfeatureGroup(object))])
          })

if (is.null(getGeneric("featureGroup")))
  setGeneric("featureGroup", function(object) standardGeneric("featureGroup"))
setMethod("featureGroup","oligoBatch", function(object) get("feature_group", envir = featureInfo(object)))

if (is.null(getGeneric("pmfeatureGroup")))
  setGeneric("pmfeatureGroup", function(object) standardGeneric("pmfeatureGroup"))
setMethod("pmfeatureGroup","oligoBatch", function(object) featureGroup(object)[featureType(object) == "PM"])

if (is.null(getGeneric("mmfeatureGroup")))
  setGeneric("mmfeatureGroup", function(object) standardGeneric("mmfeatureGroup"))
setMethod("mmfeatureGroup","oligoBatch", function(object) featureGroup(object)[featureType(object) == "MM"])

if (is.null(getGeneric("mismatch")))
  setGeneric("mismatch", function(object) standardGeneric("mismatch"))
setMethod("mismatch","oligoBatch", function(object) get("mismatch", envir = featureInfo(object)))

### PM
if (is.null(getGeneric("pm")))
  setGeneric("pm",function(object,...) standardGeneric("pm"))
setMethod("pm","oligoBatch",
          function(object) {
            pmnames <- featureNames(object)[featureType(object) == "PM"]
            pms <- exprs(object)[featureType(object) == "PM",]
            rownames(pms) <- pmnames
            pms <- pms[order(pmfeatureGroup(object)),]
            return(pms)
          })

### MM
if (is.null(getGeneric("mm")))
  setGeneric("mm",function(object,...) standardGeneric("mm"))
setMethod("mm","oligoBatch",
          function(object) {
            maxmm <- max(mismatch(object))
            mmnames <- featureNames(object)[featureType(object) == "MM"]
            mms <- list()
            for (i in 1:maxmm){
              mmsl <- exprs(object)[featureType(object) == "MM" & mismatch(object) == i,]
              rownames(mmsl) <- mmnames
              mmsl <- mmsl[order(mmfeatureGroup(object)),]
              mms[[i]] <- mmsl 
            }
            if (maxmm == 1) mms <- mms[[1]]
            return(mms)
          })
