setValidity("SnpCallSet", function(object) {
  assayDataValidMembers(assayData(object), c("calls", "callsConfidence"))
})

setMethod("initialize", "SnpCallSet",
          function(.Object,
                   featureData = new("AnnotatedDataFrame"),
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   calls = new("matrix"),
                   callsConfidence = new("matrix"),
                   ...) {
            callNextMethod(.Object,
                           assayData = assayDataNew(
                             calls = calls,
                             callsConfidence = callsConfidence,
                             ...),
                           featureData = featureData,
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation)
          })

if( is.null(getGeneric("calls<-")))
  setGeneric("calls<-", function(object, value) standardGeneric("calls<-"))

if( is.null(getGeneric("callsConfidence<-")))
  setGeneric("callsConfidence<-", function(object, value) standardGeneric("callsConfidence<-"))

if( is.null(getGeneric("calls")))
  setGeneric("calls", function(object) standardGeneric("calls"))

if( is.null(getGeneric("callsConfidence")))
  setGeneric("callsConfidence", function(object) standardGeneric("callsConfidence"))


setMethod("calls", "SnpCallSet", function(object) assayDataElement(object, "calls"))
setReplaceMethod("calls", signature(object="SnpCallSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "calls", value))

setMethod("callsConfidence", "SnpCallSet", function(object) assayDataElement(object, "callsConfidence"))
setReplaceMethod("callsConfidence", signature(object="SnpCallSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "callsConfidence", value))
