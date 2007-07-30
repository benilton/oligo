setValidity("SnpCallSet", function(object) {
  assayDataValidMembers(assayData(object), c("calls", "callsConfidence"))
})

setMethod("initialize", "SnpCallSet",
          function(.Object,
                   assayData = assayDataNew(calls=calls, callsConfidence=callsConfidence, ...),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   calls = new("matrix"),
                   callsConfidence = new("matrix"),
                   ...) {
            callNextMethod(.Object,
                           assayData = assayData,
                           featureData = featureData,
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation)
          })



setMethod("calls", "SnpCallSet", function(object) assayDataElement(object, "calls"))
setReplaceMethod("calls", signature(object="SnpCallSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "calls", value))

setMethod("callsConfidence", "SnpCallSet", function(object) assayDataElement(object, "callsConfidence"))
setReplaceMethod("callsConfidence", signature(object="SnpCallSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "callsConfidence", value))

setMethod("db", "SnpCallSet",
          function(object) db(get(annotation(object))))
