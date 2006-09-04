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

setMethod("calls", "SnpCallSet", function(object) assayDataElement(object, "calls"))
setMethod("calls<-", "SnpCallSet", function(object, value) assayDataElement(object, "calls") <- value)

setMethod("callsConfidence", "SnpCallSet", function(object) assayDataElement(object, "callsConfidence"))
setMethod("callsConfidence<-", "SnpCallSet", function(object, value) assayDataElement(object, "callsConfidence") <- value)
