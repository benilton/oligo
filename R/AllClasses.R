setClass("PDInfo",
         representation=representation(
           manufacturer="character",
           genomebuild="character"))

setClass("DBPDInfo",
         contains="PDInfo",
         representation=representation(
           getdb="function",
           tableInfo="data.frame"))

setClass("SNPPDInfo", contains="DBPDInfo")
## We hope to have ExonPDInfo, TilingPDInfo soon

setClass("AffySNPPDInfo", contains="SNPPDInfo",
         prototype=list(
           manufacturer="Affymetrix"))

setClass("platformDesign",
         contains="PDInfo",
         representation(featureInfo = "environment",
                        featureTypeDescription = "list",
                        type = "character",
                        nrow = "numeric",
                        ncol = "numeric",
                        nwells = "numeric",
                        lookup = "data.frame",
                        indexes = "list",
                        platforms="character"),
         prototype = list(lookup=data.frame(), genomebuild=character()))

setClass("FeatureSet",
         representation(manufacturer="character",
                        platform="character",
                        "VIRTUAL"),
         contains="eSet",
         prototype=list(
           manufacturer=character(),
           platform=character()))

setClass("ExpressionFeatureSet", contains="FeatureSet")
setClass("SnpFeatureSet", contains="FeatureSet")
setClass("TilingFeatureSet", contains="FeatureSet")
setClass("ExonFeatureSet", contains="FeatureSet")
setClass("SnpQSet", contains="eSet")
setClass("SnpCopyNumberSet", contains = "eSet")
setClass("SnpCallSet", contains = "eSet")
setClass("oligoSnpSet", contains=c("SnpCopyNumberSet", "SnpCallSet"))
