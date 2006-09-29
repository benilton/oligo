if( is.null(getGeneric("chromosome")))
  setGeneric("chromosome", function(object) standardGeneric("chromosome"))
setMethod("chromosome", "TilingFeatureSet", function(object) getPD(object)$chromosome)

if( is.null(getGeneric("position")))
  setGeneric("position", function(object) standardGeneric("position"))
setMethod("position", "TilingFeatureSet", function(object) getPD(object)$position)

if( is.null(getGeneric("genomeBuild")))
  setGeneric("genomeBuild", function(object) standardGeneric("genomeBuild"))
setMethod("genomeBuild", "TilingFeatureSet", function(object) getPD(object)@genomebuild)

if( is.null(getGeneric("pmPosition")))
  setGeneric("pmPosition", function(object) standardGeneric("pmPosition"))
setMethod("pmPosition", "TilingFeatureSet", function(object) position(object)[pmindex(object)])

if( is.null(getGeneric("pmChr")))
  setGeneric("pmChr", function(object) standardGeneric("pmChr"))
setMethod("pmChr", "TilingFeatureSet", function(object) chromosome(object)[pmindex(object)])

