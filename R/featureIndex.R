## BC: Mon Jul 25, 2005 - added featureIndex, before called indexProbe
##     so we can produce some plots.
##     Probably there is a better name? BC Th Jul 28, 2005

if( is.null(getGeneric("featureIndex") ))
  setGeneric("featureIndex", function(object, ...)
             standardGeneric("featureIndex"))

setMethod("featureIndex","oligoBatch",
          function(object, which=c("both","pm","mm"), genenames=NULL){
            which <- match.arg(which,c("both","pm","mm"))
            if (which=="both"){
              pmIndex <- pmindex(getPlatformDesign(object))
              mmIndex <- mmindex(getPlatformDesign(object))
              indexes <- sort(c(pmIndex,mmIndex))
            } else if (which=="pm"){
              indexes <- sort(pmindex(getPlatformDesign(object)))
            } else if (which=="mm"){
              indexes <- sort(mmindex(getPlatformDesign(object)))
            }
             return(indexes)
           })
