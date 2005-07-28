## PM methods
## BC: Th, Jul 28, 2005
##     Moved to exclusive file, because it's getting big

if( is.null(getGeneric("pm") ))
  setGeneric("pm", function(object, ...)
             standardGeneric("pm"))

## BC: Th, Jul 28, 2005 - This is an attempt of
##     selecting the feature_set_name for pm/mm
##     if it is working, we need to rename the variables
setMethod("pm","oligoBatch", ##genenames is ignored for now.. we will get to it
          function(object, genenames=NULL){
            index <- pmindex(object)
            if (!is.null(genenames)){
              index <- intersect(index,indexFeatureSetName(object,genenames))
              fsn <- getPD(object)$feature_set_name[index]
              fsid <- getPD(object)$feature_ID[index]
              rn <- paste(fsn,fsid,sep=".")
              oo <- exprs(object)[index,,drop=FALSE]
              rownames(oo) <- rn
              colnames(oo) <- sampleNames(object)
              return(oo)
            }
            return(exprs(object)[index,,drop=FALSE])
          })

### Trying to make the above to work
### if it is working, we remove the one below
###setMethod("pm","oligoBatch", ##genenames is ignored for now.. we will get to it
###          function(object, genenames=NULL){
###            return(exprs(object)[pmindex(object),,drop=FALSE])
###          })

if( is.null(getGeneric("pm<-") ))
  setGeneric("pm<-", function(object, value)
             standardGeneric("pm<-"))

setReplaceMethod("pm", "oligoBatch",
                 function(object, value){
                   exprs(object)[pmindex(object),] <- value
                   object
                 })

## end pm.R
