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
  setGeneric("designName",function(obj) standardGeneric("designName"))
setMethod("designName","oligoBatch", function(obj) obj@designName)

if (is.null(getGeneric("featureInfo")))
  setGeneric("featureInfo", function(obj) standardGeneric("featureInfo"))
setMethod("featureInfo","oligoBatch", function(obj) slot(get(designName(obj)),"featureInfo"))

if (is.null(getGeneric("featureType")))
  setGeneric("featureType", function(obj) standardGeneric("featureType"))
setMethod("featureType","oligoBatch", function(obj) get("feature_type", envir = featureInfo(obj)))

if (is.null(getGeneric("featureClass")))
  setGeneric("featureClass", function(obj) standardGeneric("featureClass"))
setMethod("featureClass","oligoBatch", function(obj) get("feature_class", envir = featureInfo(obj)))


if (is.null(getGeneric("featureNames")))
  setGeneric("featureNames", function(obj) standardGeneric("featureNames"))
setMethod("featureNames","oligoBatch",
          function(obj){
            pmfeaturenames <- get("feature_names", envir = featureInfo(obj))[featureType(obj) == "PM"]
            return(pmfeaturenames[order(pmfeatureGroup(obj))])
          })

if (is.null(getGeneric("featureGroup")))
  setGeneric("featureGroup", function(obj) standardGeneric("featureGroup"))
setMethod("featureGroup","oligoBatch", function(obj) get("feature_group", envir = featureInfo(obj)))

if (is.null(getGeneric("pmfeatureGroup")))
  setGeneric("pmfeatureGroup", function(obj) standardGeneric("pmfeatureGroup"))
setMethod("pmfeatureGroup","oligoBatch", function(obj) featureGroup(obj)[featureType(obj) == "PM"])

if (is.null(getGeneric("mmfeatureGroup")))
  setGeneric("mmfeatureGroup", function(obj) standardGeneric("mmfeatureGroup"))
setMethod("mmfeatureGroup","oligoBatch", function(obj) featureGroup(obj)[featureType(obj) == "MM"])

if (is.null(getGeneric("mismatch")))
  setGeneric("mismatch", function(obj) standardGeneric("mismatch"))
setMethod("mismatch","oligoBatch", function(obj) get("mismatch", envir = featureInfo(obj)))

### PM
if (is.null(getGeneric("pm")))
  setGeneric("pm",function(obj,...) standardGeneric("pm"))
setMethod("pm","oligoBatch",
          function(obj) {
            pmnames <- featureNames(obj)[featureType(obj) == "PM"]
            pms <- exprs(obj)[featureType(obj) == "PM",]
            rownames(pms) <- pmnames
            pms <- pms[order(pmfeatureGroup(obj)),]
            return(pms)
          })

### MM
if (is.null(getGeneric("mm")))
  setGeneric("mm",function(obj,...) standardGeneric("mm"))
setMethod("mm","oligoBatch",
          function(obj) {
            maxmm <- max(mismatch(obj))
            mmnames <- featureNames(obj)[featureType(obj) == "MM"]
            mms <- list()
            for (i in 1:maxmm){
              mmsl <- exprs(obj)[featureType(obj) == "MM" & mismatch(obj) == i,]
              rownames(mmsl) <- mmnames
              mmsl <- mmsl[order(mmfeatureGroup(obj)),]
              mms[[i]] <- mmsl 
            }
            if (maxmm == 1) mms <- mms[[1]]
            return(mms)
          })
