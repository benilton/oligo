## MM methods... designed for arrays that have one MM per PM
## BC: Th, Jul 28, 2005
##     Moved to exclusive file, bc it's getting big.

if( is.null(getGeneric("mm") ))
  setGeneric("mm", function(object, ...)
             standardGeneric("mm"))

## BC: Th, Jul 28, 2005 - This is an attempt of
##     selecting the feature_set_name for pm/mm
##     if it is working, we need to rename the variables
setMethod("mm","oligoBatch", function(object, genenames=NULL){
            index <- mmindex(object)
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


##setMethod("mm","oligoBatch", ##genenames is ignored for now.. we will get to it
##          function(object, genenames=NULL){
##            return(exprs(object)[mmindex(object),,drop=FALSE])
##          })
          
if( is.null(getGeneric("mm<-") ))
  setGeneric("mm<-", function(object, value)
             standardGeneric("mm<-"))

setReplaceMethod("mm", "oligoBatch",
                 function(object, value){
                   exprs(object)[mmindex(object),] <- value
                   object
                 })

## end mm.R
