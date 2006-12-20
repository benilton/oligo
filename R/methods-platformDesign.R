setMethod("featureInfo", "platformDesign",
          function(object) object@featureInfo)

setMethod("names","platformDesign",
          function(x){
            return(ls(featureInfo(x)))
          })

setMethod("probeNames", "platformDesign",
          function(object, subset=NULL) {
              pmIndex <- pmindex(object)
              pns <- get("feature_set_name", envir=featureInfo(object))
              return(as.character(pns[pmIndex]))
          })

setMethod("geneNames", "platformDesign",
          function(object) {
              pmIndex <- pmindex(object)
            levels(factor(get("feature_set_name",
                              envir=featureInfo(object))[pmIndex]))
          })


setMethod("featureSetNames", c("platformDesign", "missing"),
          function(object, ids) {
              as.character(featureInfo(object)$feature_set_name)
          })

## FIXME: do we need both of these?  Is the type of ids for
## platformDesign objects known?
setMethod("featureSetNames", c("platformDesign", "numeric"),
          function(object, ids) {
              as.character(featureInfo(object)$feature_set_name[ids])
          })

setMethod("featureSetNames", c("platformDesign", "character"),
          function(object, ids) {
              as.character(featureInfo(object)$feature_set_name[ids])
          })


setMethod("featureIDs", c("platformDesign"),
          function(object, ids) {
              if (!missing(ids) && length(ids) > 0)
                featureInfo(object)$feature_ID[ids]
              else
                featureInfo(object)$feature_ID
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
