featureInfo <- function(object) object@featureInfo

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
      
##nProbes method
setMethod("nProbes","platformDesign",
          function(object){
            return(length(get("feature_set_name",featureInfo(object))))
          })

##probeNames method
setMethod("probeNames","platformDesign",
          function(object){
            return(get("feature_set_name",featureInfo(object)))
          })
                                   
##WE assume feature_type is PM or MM. this might change with other platforms
setMethod("pmindex", "platformDesign",
          function(object){
            Index=which(get("feature_type",featureInfo(object))=="PM")
            return(Index)        
          })

##mmindex method.. for now we assume there is one MM per PM
###NOTE: THIS WILL CHANGE CAUSE feature type will be a vector
setMethod("mmindex", "platformDesign",
          function(object){
            Index=which(get("feature_type",featureInfo(object))=="MM")
            return(Index)
          })
