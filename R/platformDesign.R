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
                        ncol = "numeric",
                        nwells = "numeric",
                        lookup = "data.frame"))

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
  cat("Number of columns:",object@ncol,"\n")
  cat("Number of rows:",object@nrow,"\n")
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

##WE assume feature_type_1 is PM or MM. this might change with other platforms
setMethod("pmindex", "platformDesign",
          function(object){
            Index=which(get("feature_type_1",featureInfo(object))=="PM")

##            ids=get("feature_ID",featureInfo(object))
##            pns = probeNames(object)
  return(Index)        
##            return(Index[order(pns[Index],ids[Index])])
          })


##setMethod("pmindex", "platformDesign",
##          function(object){
##            Index=which(get("feature_type_1",featureInfo(object))=="PM")
##
##            ids=get("feature_ID",featureInfo(object))
##            pns = probeNames(object)
##          
##            return(Index[order(pns[Index],ids[Index])])
##          })

##mmindex method.. for now we assume there is one MM per PM
if( is.null(getGeneric("mmindex")))
  setGeneric("mmindex", function(object,...)
             standardGeneric("mmindex"))

###NOTE: THIS WILL CHANGE CAUSE feature type will be a vector
setMethod("mmindex", "platformDesign",
          function(object){
            Index=which(get("feature_type_1",featureInfo(object))=="MM")

##            ids=get("feature_ID",featureInfo(object))
##            pns = probeNames(object)
  return(Index)
##            return(Index[order(pns[Index],ids[Index])])
          })
