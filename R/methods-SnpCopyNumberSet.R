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

