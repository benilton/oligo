setMethod("featureInfo", "platformDesign",
          function(object) object@featureInfo)

setMethod("names","platformDesign",
          function(x){
            return(ls(featureInfo(x)))
          })

setAs("platformDesign", "data.frame",
      function(from) {
          df <- as.list(featureInfo(x))
          attributes(df) <- list(class="data.frame",
                                 names=names(df),
                                 row.names=1:length(df[[1]]))
          df
      })

setMethod("$", "platformDesign",
          function(x, name) {
              get(name, featureInfo(x))
          })
      
##nProbes method

setMethod("nProbes","platformDesign",
          function(object){
            return(length(get("feature_set_name",featureInfo(object))))
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
