if (is.null(getGeneric("manufacturer")))
  setGeneric("manufacturer",function(object) standardGeneric("manufacturer"))

if( is.null(getGeneric("manufacturer<-") ))
  setGeneric("manufacturer<-", function(object, value)
             standardGeneric("manufacturer<-"))

if (is.null(getGeneric("platform")))
  setGeneric("platform",function(object) standardGeneric("platform"))

if( is.null(getGeneric("platform<-") ))
  setGeneric("platform<-", function(object, value)
             standardGeneric("platform<-"))

if (is.null(getGeneric("platformDesignName"))){
  setGeneric("platformDesignName",
             function(object) standardGeneric("platformDesignName"))}

if (is.null(getGeneric("getPlatformDesign"))){
  setGeneric("getPlatformDesign",
             function(object) standardGeneric("getPlatformDesign"))}

## if (is.null(getGeneric("probeNames")))
##   setGeneric("probeNames", function(object, subset)
##              standardGeneric("probeNames"))

if (is.null(getGeneric("geneNames")))
  setGeneric("geneNames", function(object)
             standardGeneric("geneNames"))

if( is.null(getGeneric("pmindex")))
  setGeneric("pmindex", function(object)
             standardGeneric("pmindex"))

if( is.null(getGeneric("mmindex")))
  setGeneric("mmindex", function(object)
             standardGeneric("mmindex"))

if( is.null(getGeneric("indexFeatureSetName") ))
  setGeneric("indexFeatureSetName", function(object, featurenames)
             standardGeneric("indexFeatureSetName"))

if(is.null(getGeneric("ncol")))
  setGeneric("ncol")

if( is.null(getGeneric("nrow")))
  setGeneric("nrow")

if( is.null(getGeneric("hist")) )
  setGeneric("hist")

if( is.null(getGeneric("names")))
  setGeneric("names", function(x)
             standardGeneric("names"))

if( is.null(getGeneric("nProbes")))
  setGeneric("nProbes", function(object)
             standardGeneric("nProbes"))

## if( is.null(getGeneric("probeNames")))
##   setGeneric("probeNames", function(object, ...)
##              standardGeneric("probeNames"))

if( is.null(getGeneric("pm") ))
  setGeneric("pm", function(object, genenames=NULL)
             standardGeneric("pm"))

if( is.null(getGeneric("pm<-") ))
  setGeneric("pm<-", function(object, value)
             standardGeneric("pm<-"))

if( is.null(getGeneric("mm") ))
  setGeneric("mm", function(object, genenames=NULL)
             standardGeneric("mm"))

if( is.null(getGeneric("mm<-") ))
  setGeneric("mm<-", function(object, value)
             standardGeneric("mm<-"))

if( is.null(getGeneric("image")))
  setGeneric("image")

if( is.null(getGeneric("featureIndex") ))
  setGeneric("featureIndex", function(object, which=c("both","pm","mm"), genenames=NULL)
             standardGeneric("featureIndex"))

if( is.null(getGeneric("boxplot")))
  setGeneric("boxplot")

if( is.null(getGeneric("npixels")))
  setGeneric("npixels",
             function(object) standardGeneric("npixels"))

if( is.null(getGeneric("allele")))
  setGeneric("allele", function(object) standardGeneric("allele"))

if( is.null(getGeneric("chromosome")))
  setGeneric("chromosome", function(object) standardGeneric("chromosome"))

if( is.null(getGeneric("position")))
  setGeneric("position", function(object) standardGeneric("position"))

if( is.null(getGeneric("genomeBuild")))
  setGeneric("genomeBuild", function(object) standardGeneric("genomeBuild"))

if( is.null(getGeneric("pmChr")))
  setGeneric("pmChr", function(object) standardGeneric("pmChr"))

if( is.null(getGeneric("pmPosition")))
  setGeneric("pmPosition", function(object) standardGeneric("pmPosition"))

if( is.null(getGeneric("snpBasePair")))
  setGeneric("snpBasePair", function(object) standardGeneric("snpBasePair"))

if( is.null(getGeneric("pmSnpBasePair")))
  setGeneric("pmSnpBasePair", function(object) standardGeneric("pmSnpBasePair"))

if( is.null(getGeneric("alleleAB")))
  setGeneric("alleleAB", function(object) standardGeneric("alleleAB"))

if( is.null(getGeneric("pmAlleleAB")))
  setGeneric("pmAlleleAB", function(object) standardGeneric("pmAlleleAB"))

if( is.null(getGeneric("calls<-")))
  setGeneric("calls<-", function(object, value) standardGeneric("calls<-"))

if( is.null(getGeneric("callsConfidence<-")))
  setGeneric("callsConfidence<-", function(object, value) standardGeneric("callsConfidence<-"))

if( is.null(getGeneric("copyNumber<-")))
  setGeneric("copyNumber<-", function(object, value) standardGeneric("copyNumber<-"))

if( is.null(getGeneric("cnConfidence<-")))
  setGeneric("cnConfidence<-", function(object, value) standardGeneric("cnConfidence<-"))

if( is.null(getGeneric("calls")))
  setGeneric("calls", function(object) standardGeneric("calls"))

if( is.null(getGeneric("callsConfidence")))
  setGeneric("callsConfidence", function(object) standardGeneric("callsConfidence"))

if( is.null(getGeneric("copyNumber")))
  setGeneric("copyNumber", function(object) standardGeneric("copyNumber"))

if( is.null(getGeneric("cnConfidence")))
  setGeneric("cnConfidence", function(object) standardGeneric("cnConfidence"))

if( is.null(getGeneric("snpMedianSilhouette")))
  setGeneric("snpMedianSilhouette", function(object) standardGeneric("snpMedianSilhouette"))

if( is.null(getGeneric("logRatioAntisense")))
  setGeneric("logRatioAntisense", function(object) standardGeneric("logRatioAntisense"))

if( is.null(getGeneric("logRatioSense")))
  setGeneric("logRatioSense", function(object) standardGeneric("logRatioSense"))
