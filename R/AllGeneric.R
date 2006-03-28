if (is.null(getGeneric("platform")))
  setGeneric("platform",function(object) standardGeneric("platform"))

if (is.null(getGeneric("platformDesignName"))){
  setGeneric("platformDesignName",
             function(object) standardGeneric("platformDesignName"))}

if (is.null(getGeneric("getPlatformDesign"))){
  setGeneric("getPlatformDesign",
             function(object) standardGeneric("getPlatformDesign"))}

if (is.null(getGeneric("probeNames")))
  setGeneric("probeNames", function(object, ...)
             standardGeneric("probeNames"))

if (is.null(getGeneric("geneNames")))
  setGeneric("geneNames", function(object)
             standardGeneric("geneNames"))

if( is.null(getGeneric("pmindex")))
  setGeneric("pmindex", function(object,...)
             standardGeneric("pmindex"))

if( is.null(getGeneric("mmindex")))
  setGeneric("mmindex", function(object,...)
             standardGeneric("mmindex"))

if( is.null(getGeneric("indexFeatureSetName") ))
  setGeneric("indexFeatureSetName", function(object, ...)
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
  setGeneric("nProbes", function(object, ...)
             standardGeneric("nProbes"))

if( is.null(getGeneric("probeNames")))
  setGeneric("probeNames", function(object, ...)
             standardGeneric("probeNames"))

if( is.null(getGeneric("pm") ))
  setGeneric("pm", function(object, ...)
             standardGeneric("pm"))

if( is.null(getGeneric("pm<-") ))
  setGeneric("pm<-", function(object, value)
             standardGeneric("pm<-"))

if( is.null(getGeneric("mm") ))
  setGeneric("mm", function(object, ...)
             standardGeneric("mm"))

if( is.null(getGeneric("mm<-") ))
  setGeneric("mm<-", function(object, value)
             standardGeneric("mm<-"))

if( is.null(getGeneric("image")))
  setGeneric("image")

if( is.null(getGeneric("featureIndex") ))
  setGeneric("featureIndex", function(object, ...)
             standardGeneric("featureIndex"))

if( is.null(getGeneric("boxplot")))
  setGeneric("boxplot")

if( is.null(getGeneric("npixels")))
  setGeneric("npixels",
             function(object) standardGeneric("npixels"))

if( is.null(getGeneric("allele")))
  setGeneric("allele", function(object) standardGeneric("allele"))
