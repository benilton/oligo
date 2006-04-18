setClass("platformDesign",
         representation(featureInfo = "environment",
                        featureTypeDescription = "list",
                        manufacturer = "character",
                        type = "character",
                        nrow = "numeric",
                        ncol = "numeric",
                        nwells = "numeric",
                        lookup = "data.frame",
                        indexes = "list",
                        platforms="character"),
         prototype = list(lookup=data.frame()))

setClass("FeatureSet",
         representation(manufacturer="character",
                        platform="character"),
         contains="eSet",
         prototype=list(
           manufacturer=character(),
           platform=character()))

