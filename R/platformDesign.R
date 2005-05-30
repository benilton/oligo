setClass("platformDesign",
         representation(featureInfo = "environment",
                        manufacturer = "character",
                        type = "character"))


##functions and methods
featureInfo <- function(object) object@featureInfo

featureInfoNames <- function(object) ls(featureInfo(object))

##X is arbitrarily chosen
nProbes <- function(object) length( get("X",featureInfo(object)))


##probeNames method
if( is.null(getGeneric("probeNames")))
  setGeneric("probeNames", function(object, ...)
             standardGeneric("probeNames"))

setMethod("probeNames","platformDesign",
          function(object){
            return(get("feature_name",featureInfo(object)))
          })
                                   

##pmindex method
if( is.null(getGeneric("pmindex")))
  setGeneric("pmindex", function(object,...)
             standardGeneric("pmindex"))

###NOTE: THIS WILL CHANGE CAUSE feature type will be a vector
setMethod("pmindex", "platformDesign",
          function(object){
            Index=which(get("feature_type",featureInfo(object))=="PM")
            pns = probeNames(object)
            return(Index[order(pns[Index])])
          })

