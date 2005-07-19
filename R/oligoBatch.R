# Methods
# Author: Benilton Carvalho
# Date: April 2005
setClass("oligoBatch",
         representation(manufacturer="character",
                        platform="character",
                        description="MIAME",
                        notes="character"),
         contains="eSet")

##simple accessors that needs to change!
if (is.null(getGeneric("annotation")))
  setGeneric("annotation",function(object) standardGeneric("annotation"))
setMethod("annotation","oligoBatch", function(object) object@platform)

###same as above...but actual accessor
if (is.null(getGeneric("platform")))
  setGeneric("platform",function(object) standardGeneric("platform"))
setMethod("platform","oligoBatch", function(object) object@platform)

###for compatibility with previous package
setMethod("length",signature(x="oligoBatch"),
          function(x) ncol(exprs(x))) 

###this might change. we might put pd at the end
## BC: May 29 - pd added at the end
## BC: Jul 13 - pd moved to the begining
if (is.null(getGeneric("platformDesignName"))){
  setGeneric("platformDesignName",
             function(object) standardGeneric("platformDesignName"))}

setMethod("platformDesignName","oligoBatch", function(object){
  cleanPlatformName(object@platform)})

##loading the library for now... this must change
if (is.null(getGeneric("getPlatformDesign"))){
  setGeneric("getPlatformDesign",
             function(object) standardGeneric("getPlatformDesign"))}

setMethod("getPlatformDesign","oligoBatch", function(object){
  pdn <- platformDesignName(object)
  library(pdn,character.only=TRUE)
  return(get(pdn,pos=paste("package:",pdn,sep="")))
})

###probeNames - returns probeNames for PMs ... genenames ignored for now
if (is.null(getGeneric("probeNames")))
  setGeneric("probeNames", function(object, ...)
             standardGeneric("probeNames"))

setMethod("probeNames", "oligoBatch",
          function(object, genenames=NULL){
            pmIndex <- pmindex(getPlatformDesign(object))
            pns <- get("feature_set_name",envir=featureInfo(getPlatformDesign(object)))
            return(as.character(pns[pmIndex]))
          })

###geneNames - returns geneNames for PMs
if (is.null(getGeneric("geneNames")))
  setGeneric("geneNames", function(object)
             standardGeneric("geneNames"))

setMethod("geneNames", "oligoBatch",
          function(object){
            pmIndex <- pmindex(getPlatformDesign(object))
            return(levels(factor(get("feature_set_name",envir=featureInfo(getPlatformDesign(object)))[pmIndex])))
          })


##PM methods
if( is.null(getGeneric("pm") ))
  setGeneric("pm", function(object, ...)
             standardGeneric("pm"))

setMethod("pm","oligoBatch", ##genenames is ignored for now.. we will get to it
          function(object, genenames=NULL){
             pmIndex <- pmindex(getPlatformDesign(object))
             return(exprs(object)[pmIndex,,drop=FALSE])
           })
          
if( is.null(getGeneric("pm<-") ))
  setGeneric("pm<-", function(object, value)
             standardGeneric("pm<-"))

setReplaceMethod("pm", "oligoBatch",
                 function(object, value){
                   pmIndex <- pmindex(getPlatformDesign(object))
                   exprs(object)[pmIndex,] <- value
                   object
                 })


###sampleNames and description should go away once eSet has them
if( !isGeneric("sampleNames") )
    setGeneric("sampleNames", function(object)
               standardGeneric("sampleNames"))
  setMethod("sampleNames", "oligoBatch",
            function(object) {
              if (! is.null(colnames(exprs(object))))
                colnames(exprs(object))
              else
                row.names(pData(object))
            })

##description
if( !isGeneric("description") )
    setGeneric("description", function(object)
               standardGeneric("description"))
  setMethod("description", "oligoBatch", function(object)
            object@description)

  ##replace method for description
  if( !isGeneric("description<-") )
    setGeneric("description<-", function(object, value)
               standardGeneric("description<-"))

  setReplaceMethod("description", "oligoBatch", function(object, value) {
    object@description <- value
    object
  })

##notes
if( !isGeneric("notes") )
    setGeneric("notes", function(object)
               standardGeneric("notes"))
  setMethod("notes", "oligoBatch", function(object)
            object@notes)

  if( !isGeneric("notes<-") )
    setGeneric("notes<-", function(object, value)
               standardGeneric("notes<-"))

  setReplaceMethod("notes", "oligoBatch", function(object, value) {
    object@notes <- value
    object
  })


###THESE will have to change:
# if (is.null(getGeneric("featureInfo")))
#   setGeneric("featureInfo", function(object) standardGeneric("featureInfo"))
# setMethod("featureInfo","oligoBatch", function(object) slot(get(platform(object)),"featureInfo"))

# if (is.null(getGeneric("featureType")))
#   setGeneric("featureType", function(object) standardGeneric("featureType"))
# setMethod("featureType","oligoBatch", function(object) get("feature_type", envir = featureInfo(object)))

# if (is.null(getGeneric("featureClass")))
#   setGeneric("featureClass", function(object) standardGeneric("featureClass"))
# setMethod("featureClass","oligoBatch", function(object) get("feature_class", envir = featureInfo(object)))


# if (is.null(getGeneric("featureNames")))
#   setGeneric("featureNames", function(object) standardGeneric("featureNames"))
# setMethod("featureNames","oligoBatch",
#           function(object){
#             pmfeaturenames <- get("feature_set_name", envir = featureInfo(object))[featureType(object) == "PM"]
#             return(pmfeaturenames[order(pmfeatureGroup(object))])
#           })

# if (is.null(getGeneric("featureGroup")))
#   setGeneric("featureGroup", function(object) standardGeneric("featureGroup"))
# setMethod("featureGroup","oligoBatch", function(object) get("feature_group", envir = featureInfo(object)))

# if (is.null(getGeneric("pmfeatureGroup")))
#   setGeneric("pmfeatureGroup", function(object) standardGeneric("pmfeatureGroup"))
# setMethod("pmfeatureGroup","oligoBatch", function(object) featureGroup(object)[featureType(object) == "PM"])

# if (is.null(getGeneric("mmfeatureGroup")))
#   setGeneric("mmfeatureGroup", function(object) standardGeneric("mmfeatureGroup"))
# setMethod("mmfeatureGroup","oligoBatch", function(object) featureGroup(object)[featureType(object) == "MM"])

# if (is.null(getGeneric("mismatch")))
#   setGeneric("mismatch", function(object) standardGeneric("mismatch"))
# setMethod("mismatch","oligoBatch", function(object) get("mismatch", envir = featureInfo(object)))

# ### PM
# if (is.null(getGeneric("pm")))
#   setGeneric("pm",function(object,...) standardGeneric("pm"))
# setMethod("pm","oligoBatch",
#           function(object) {
#             pmnames <- featureNames(object)[featureType(object) == "PM"]
#             pms <- exprs(object)[featureType(object) == "PM",]
#             rownames(pms) <- pmnames
#             pms <- pms[order(pmfeatureGroup(object)),]
#             return(pms)
#           })

# ### MM
# if (is.null(getGeneric("mm")))
#   setGeneric("mm",function(object,...) standardGeneric("mm"))
# setMethod("mm","oligoBatch",
#           function(object) {
#             maxmm <- max(mismatch(object))
#             mmnames <- featureNames(object)[featureType(object) == "MM"]
#             mms <- list()
#             for (i in 1:maxmm){
#               mmsl <- exprs(object)[featureType(object) == "MM" & mismatch(object) == i,]
#               rownames(mmsl) <- mmnames
#               mmsl <- mmsl[order(mmfeatureGroup(object)),]
#               mms[[i]] <- mmsl 
#             }
#             if (maxmm == 1) mms <- mms[[1]]
#             return(mms)
#           })
