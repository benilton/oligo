setValidity("SnpCopyNumberSet", function(object) {
  assayDataValidMembers(assayData(object), c("copyNumber", "cnConfidence"))
})

setMethod("initialize", "SnpCopyNumberSet",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   copyNumber = new("matrix"),
                   cnConfidence = new("matrix")) {
            callNextMethod(.Object,
                           assayData = assayDataNew(
                             copyNumber = copyNumber,
                             cnConfidence = cnConfidence),
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation)
          })

setMethod("copyNumber", "SnpCopyNumberSet", function(object) assayDataElement(object, "copyNumber"))
setReplaceMethod("copyNumber", signature(object="SnpCopyNumberSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "copyNumber", value))

setMethod("cnConfidence", "SnpCopyNumberSet", function(object) assayDataElement(object, "cnConfidence"))
setReplaceMethod("cnConfidence", signature(object="SnpCopyNumberSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "cnConfidence", value))
