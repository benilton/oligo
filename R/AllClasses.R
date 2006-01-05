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

## Add affysnpBatch; affyexprsBatch; ngexprBatch

setClass("affysnpBatch",
         contains="oligoBatch",
         prototype=list(manufacturer="Affymetrix"))

setClass("affyexprsBatch",
         contains="oligoBatch",
         prototype=list(manufacturer="Affymetrix"))

setClass("ngexprsBatch",
         contains="oligoBatch",
         prototype=list(manufacturer="NimbleGen Systems"))
