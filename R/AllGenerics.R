setGeneric("platform", function(object) standardGeneric("platform"))

setGeneric("platform<-", function(object, value) standardGeneric("platform<-"))

setGeneric("db", function(object) standardGeneric("db"))

setGeneric("manufacturer",function(object) standardGeneric("manufacturer"))

setGeneric("manufacturer<-",
           function(object, value) standardGeneric("manufacturer<-"))

setGeneric("kind", function(object) standardGeneric("kind"))

setGeneric("platformDesignName",
           function(object) standardGeneric("platformDesignName"))

setGeneric("getPlatformDesign",
           function(object) standardGeneric("getPlatformDesign"))

## setGeneric("geneNames", function(object) standardGeneric("geneNames"))

## setGeneric("indexFeatureSetName",
##            function(object, featurenames) {
##                standardGeneric("indexFeatureSetName")
##            })

## setGeneric("featureSetNames",
##            function(object, ids) standardGeneric("featureSetNames"))

## setGeneric("featureIDs",
##            function(object, ids) standardGeneric("featureIDs"))

setGeneric("pm", function(object, genenames=NULL) standardGeneric("pm"))
setGeneric("pm<-", function(object, value) standardGeneric("pm<-"))

setGeneric("mm", function(object, genenames=NULL) standardGeneric("mm"))
setGeneric("mm<-", function(object, value) standardGeneric("mm<-"))

setGeneric("featureIndex",
           function(object, which=c("both","pm","mm"), genenames=NULL) {
               standardGeneric("featureIndex")
           })

## setGeneric("npixels",
##            function(object) standardGeneric("npixels"))

## setGeneric("allele", function(object) standardGeneric("allele"))

setGeneric("calls<-", function(object, value) standardGeneric("calls<-"))

setGeneric("callsConfidence<-",
           function(object, value) standardGeneric("callsConfidence<-"))

setGeneric("calls", function(object) standardGeneric("calls"))

setGeneric("callsConfidence", function(object) standardGeneric("callsConfidence"))

setGeneric("logRatioAntisense",
           function(object) standardGeneric("logRatioAntisense"))

setGeneric("logRatioSense",
           function(object) standardGeneric("logRatioSense"))

setGeneric("snpMedianSilhouette",
           function(object) standardGeneric("snpMedianSilhouette"))

setGeneric("copyNumber<-",
           function(object, value) standardGeneric("copyNumber<-"))

setGeneric("cnConfidence<-",
           function(object, value) standardGeneric("cnConfidence<-"))

setGeneric("copyNumber", function(object) standardGeneric("copyNumber"))
setGeneric("cnConfidence", function(object) standardGeneric("cnConfidence"))
## setGeneric("alleleAB", function(object) standardGeneric("alleleAB"))
## setGeneric("pmAlleleAB", function(object) standardGeneric("pmAlleleAB"))
setGeneric("senseThetaA", function(obj) standardGeneric("senseThetaA"))
setGeneric("senseThetaB", function(obj) standardGeneric("senseThetaB"))
setGeneric("antisenseThetaA", function(obj) standardGeneric("antisenseThetaA"))
setGeneric("antisenseThetaB", function(obj) standardGeneric("antisenseThetaB"))
setGeneric("antisenseThetaB", function(obj) standardGeneric("antisenseThetaB"))
setGeneric("getM", function(obj) standardGeneric("getM"))
setGeneric("getA", function(obj) standardGeneric("getA"))

setGeneric("chromosome", function(object) standardGeneric("chromosome"))
setGeneric("position", function(object) standardGeneric("position"))
setGeneric("genomeBuild", function(object) standardGeneric("genomeBuild"))
setGeneric("pmPosition", function(object) standardGeneric("pmPosition"))
setGeneric("pmChr", function(object) standardGeneric("pmChr"))

setGeneric("probeNames", function(object, subset=NULL) standardGeneric("probeNames"))

setGeneric("nProbes", function(object) standardGeneric("nProbes"))
setGeneric("pmindex", function(object) standardGeneric("pmindex"))
setGeneric("mmindex", function(object) standardGeneric("mmindex"))

setGeneric("featureInfo", function(object) standardGeneric("featureInfo"))

setGeneric("listFeatureFields", function(object) standardGeneric("listFeatureFields"))
setGeneric("listFeatureSetFields", function(object) standardGeneric("listFeatureSetFields"))

setGeneric("pmFragmentLength", function(object) standardGeneric("pmFragmentLength"))
setGeneric("pmSequence", function(object) standardGeneric("pmSequence"))
setGeneric("mmSequence", function(object) standardGeneric("mmSequence"))
setGeneric("pmAllele", function(object) standardGeneric("pmAllele"))
setGeneric("pmStrand", function(object) standardGeneric("pmStrand"))

setGeneric("plotDensity", function(object, ...) standardGeneric("plotDensity"))
setGeneric("hist", function(x, ...) standardGeneric("hist"))
setGeneric("boxplot", function(x, ...) standardGeneric("boxplot"))
setGeneric("image", function(x, ...) standardGeneric("image"))

setGeneric("nrow", function(x) standardGeneric("nrow"))
setGeneric("ncol", function(x) standardGeneric("ncol"))

setGeneric("rma", function(object, ...) standardGeneric("rma"))
