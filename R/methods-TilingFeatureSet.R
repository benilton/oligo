setMethod("chromosome", "TilingFeatureSet", function(object) getPD(object)$chromosome)

setMethod("position", "TilingFeatureSet", function(object) getPD(object)$position)

setMethod("GenomeBuild", "TilingFeatureSet", function(object) getPD(object)@genomebuild)

setMethod("pmPosition", "TilingFeatureSet", function(object) position(object)[pmindex(object)])

setMethod("pmChr", "TilingFeatureSet", function(object) chromosome(object)[pmindex(object)])

