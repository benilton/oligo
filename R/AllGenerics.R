setGeneric("getX", function(object, type) standardGeneric("getX"))
setGeneric("getY", function(object, type) standardGeneric("getY"))
setGeneric("bg", function(object, subset=NULL) standardGeneric("bg"))

setGeneric("pm", function(object, subset=NULL, ...) standardGeneric("pm"))
setGeneric("pm<-", function(object, value) standardGeneric("pm<-"))
setGeneric("mm", function(object, subset=NULL) standardGeneric("mm"))
setGeneric("mm<-", function(object, value) standardGeneric("mm<-"))
setGeneric("bg", function(object, subset=NULL) standardGeneric("bg"))
setGeneric("bg<-", function(object, value) standardGeneric("bg<-"))
setGeneric("logRatioAntisense", function(object) standardGeneric("logRatioAntisense"))
setGeneric("logRatioSense", function(object) standardGeneric("logRatioSense"))

setGeneric("hist", function(x, ...) standardGeneric("hist"))
setGeneric("image", function(x, ...) standardGeneric("image"))
setGeneric("rma", function(object, ...) standardGeneric("rma"))
setGeneric("MAplot", function(object, ...) standardGeneric("MAplot"))
setGeneric("genotype", function(object, ...) standardGeneric("genotype"))

setGeneric("getContainer", function(object, probeType=c("pm", "bg")) standardGeneric("getContainer"))

##
setGeneric("probeNames", function(object, subset=NULL) standardGeneric("probeNames"))
setGeneric("bgindex", function(object, ...) standardGeneric("bgindex"))
setGeneric("pmindex", function(object, ...) standardGeneric("pmindex"))
setGeneric("mmindex", function(object, ...) standardGeneric("mmindex"))
setGeneric("bgSequence", function(object, ...) standardGeneric("bgSequence"))
setGeneric("pmSequence", function(object, ...) standardGeneric("pmSequence"))
setGeneric("mmSequence", function(object) standardGeneric("mmSequence"))
setGeneric("pmOffset", function(object) standardGeneric("pmOffset"))
setGeneric("pmChr", function(object) standardGeneric("pmChr"))
setGeneric("pmFragmentLength", function(object) standardGeneric("pmFragmentLength"))
setGeneric("pmAllele", function(object) standardGeneric("pmAllele"))
setGeneric("pmStrand", function(object) standardGeneric("pmStrand"))
setGeneric("pmPosition", function(object) standardGeneric("pmPosition"))
setGeneric("getPlatformDesign", function(object) standardGeneric("getPlatformDesign"))
setGeneric("getX", function(object, type) standardGeneric("getX"))
setGeneric("getY", function(object, type) standardGeneric("getY"))

###

setGeneric("backgroundCorrect",
           function(object, method="rma", copy=TRUE, verbose=TRUE, ...)
           standardGeneric("backgroundCorrect"))

setGeneric("normalize",
           function(object, method="quantile", copy=TRUE, verbose=TRUE, ...)
           standardGeneric("normalize"))

setGeneric("summarize",
           function(object, probes=rownames(object),
           method="medianpolish", verbose=TRUE, ...)
           standardGeneric("summarize"))

setGeneric("normalizeToTarget",
           function(object, target, method="quantile", copy=TRUE, verbose=TRUE)
           standardGeneric("normalizeToTarget"))

