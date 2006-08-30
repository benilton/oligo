setValidity("SnpCallSet", function(object) {
  assayDataValidMembers(assayData(object), c("calls", "callsConfidence"))
})

setMethod("initialize", "SnpCallSet",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   calls = new("matrix"),
                   callsConfidence = new("matrix")) {
            callNextMethod(.Object,
                           assayData = assayDataNew(
                             calls = calls,
                             callsConfidence = callsConfidence),
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation)
          })

