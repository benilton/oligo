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
setMethod("copyNumber<-", "SnpCopyNumberSet", function(object, value) assayDataElement(object, "copyNumber") <- value)

setMethod("cnConfidence", "SnpCallSet", function(object) assayDataElement(object, "cnConfidence"))
setMethod("cnConfidence<-", "SnpCallSet", function(object, value) assayDataElement(object, "cnConfidence") <- value)
