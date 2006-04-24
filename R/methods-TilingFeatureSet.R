setMethod("chromosome", "TilingFeatureSet", function(object) getPD(object)$chromosome)

setMethod("position", "TilingFeatureSet", function(object) getPD(object)$position)

setMethod("GenomeBuild", "TilingFeatureSet", function(object) getPD(object)@genomebuild)
