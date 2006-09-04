setClass("platformDesign",
         representation(featureInfo = "environment",
                        featureTypeDescription = "list",
                        manufacturer = "character",
                        type = "character",
                        genomebuild = "character",
                        nrow = "numeric",
                        ncol = "numeric",
                        nwells = "numeric",
                        lookup = "data.frame",
                        indexes = "list",
                        platforms="character"),
         prototype = list(lookup=data.frame(), genomebuild=character()))

setClass("FeatureSet",
         representation(manufacturer="character",
                        platform="character"),
         contains="eSet",
         prototype=list(
           manufacturer=character(),
           platform=character()))

##setClassUnion("characterOrNULL", c("character", NULL))

setClass("ExpressionFeatureSet", contains="FeatureSet")
setClass("SnpFeatureSet", contains="FeatureSet")
setClass("TilingFeatureSet", contains="FeatureSet")
setClass("SnpQSet", contains="eSet")
setClass("SnpCopyNumberSet", contains = "eSet")
setClass("SnpCallSet", contains = "eSet")

setClass("SnpSet", contains=c("SnpCopyNumberSet", "SnpCallSet"))


## setMethod("initialize", "SnpSet",
##           function(.Object,
##                    phenoData = new("AnnotatedDataFrame"),
##                    experimentData = new("MIAME"),
##                    annotation = character(),
##                    calls = new("matrix"),
##                    copyNumber = new("matrix")) {
##             .Object@assayData$calls <- calls
##             .Object@assayData$copyNumber <- copyNumber
##             .Object
##           })

