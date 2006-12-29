setMethod("snpBasePair", "SnpFeatureSet", function(object) getPD(object)$middle_base)

setMethod("pmSnpBasePair", "SnpFeatureSet", function(object) getPD(object)$middle_base[pmindex(object)])

setMethod("alleleAB", "SnpFeatureSet", function(object) getPD(object)$alleleAB)

setMethod("pmAlleleAB", "SnpFeatureSet", function(object) getPD(object)$alleleAB[pmindex(object)])

setMethod("allele", signature(object="SnpFeatureSet"),
          function(object) getPD(object)$allele)

