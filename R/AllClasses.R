setClass("platformDesign",
         representation(featureInfo = "environment",
                        featureTypeDescription = "list",
                        manufacturer = "character",
                        type = "character",
                        nrow = "numeric",
                        ncol = "numeric",
                        nwells = "numeric",
                        lookup = "data.frame"))

setClass("oligoBatch",
         representation(manufacturer="character",
                        platform="character"),
         contains="eSet")
