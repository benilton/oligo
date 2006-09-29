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

if( is.null(getGeneric("copyNumber<-")))
  setGeneric("copyNumber<-", function(object, value) standardGeneric("copyNumber<-"))

if( is.null(getGeneric("cnConfidence<-")))
  setGeneric("cnConfidence<-", function(object, value) standardGeneric("cnConfidence<-"))

if( is.null(getGeneric("copyNumber")))
  setGeneric("copyNumber", function(object) standardGeneric("copyNumber"))

if( is.null(getGeneric("cnConfidence")))
  setGeneric("cnConfidence", function(object) standardGeneric("cnConfidence"))


setMethod("copyNumber", "SnpCopyNumberSet", function(object) assayDataElement(object, "copyNumber"))
setReplaceMethod("copyNumber", signature(object="SnpCopyNumberSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "copyNumber", value))

setMethod("cnConfidence", "SnpCopyNumberSet", function(object) assayDataElement(object, "cnConfidence"))
setReplaceMethod("cnConfidence", signature(object="SnpCopyNumberSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "cnConfidence", value))
