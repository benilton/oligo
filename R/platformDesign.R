## BC: Jul 13, "pd" at the begining
cleanPlatformName <- function(x)
  gsub("[_-]","",paste("pd",tolower(x),sep=""))

## BC: added nrow/ncol
setClass("platformDesign",
         representation(featureInfo = "environment",
                        featureTypeDescription = "list",
                        manufacturer = "character",
                        type = "character",
                        nrow = "numeric",
                        ncol = "numeric"))

##functions and methods
featureInfo <- function(object) object@featureInfo

## i just made names a generic method... hope this isnt bad!
if( is.null(getGeneric("names")))
  setGeneric("names", function(x)
             standardGeneric("names"))
setMethod("names","platformDesign",
          function(x){
            return(ls(featureInfo(x)))
          })


##show method
setMethod("show","platformDesign", function(object){
  cat("Manufacturer:",object@manufacturer,"\n")
  cat("Array type:",object@type,"\n")
  cat("Number of features:",nProbes(object),"\n")
})


as.data.frame.platformDesign <- function(x,row.names=NULL,optional=FALSE){
  return(as.data.frame(as.list(featureInfo(x)),row.names=row.names,optional=optional))
}

"$.platformDesign" <- function(object, val)
  get(val,featureInfo(object))
      
##X is arbitrarily chosen
nProbes <- function(object){
  colname <- ls(featureInfo(object))[1]
  length( get(colname,featureInfo(object)))
}

##probeNames method
if( is.null(getGeneric("probeNames")))
  setGeneric("probeNames", function(object, ...)
             standardGeneric("probeNames"))

setMethod("probeNames","platformDesign",
          function(object){
            return(get("feature_set_name",featureInfo(object)))
          })
                                   
##pmindex method
if( is.null(getGeneric("pmindex")))
  setGeneric("pmindex", function(object,...)
             standardGeneric("pmindex"))

###NOTE: THIS WILL CHANGE CAUSE feature type will be a vector
setMethod("pmindex", "platformDesign",
          function(object){
            Index=which(get("feature_type_1",featureInfo(object))=="PM")
            pns = probeNames(object)
            return(Index[order(pns[Index])])
          })

