if( is.null(getGeneric("snpBasePair")))
  setGeneric("snpBasePair", function(object) standardGeneric("snpBasePair"))

if( is.null(getGeneric("pmSnpBasePair")))
  setGeneric("pmSnpBasePair", function(object) standardGeneric("pmSnpBasePair"))

if( is.null(getGeneric("alleleAB")))
  setGeneric("alleleAB", function(object) standardGeneric("alleleAB"))

if( is.null(getGeneric("pmAlleleAB")))
  setGeneric("pmAlleleAB", function(object) standardGeneric("pmAlleleAB"))

setMethod("snpBasePair", "SnpFeatureSet", function(object) getPD(object)$middle_base)

setMethod("pmSnpBasePair", "SnpFeatureSet", function(object) getPD(object)$middle_base[pmindex(object)])

setMethod("alleleAB", "SnpFeatureSet", function(object) getPD(object)$alleleAB)

setMethod("pmAlleleAB", "SnpFeatureSet", function(object) getPD(object)$alleleAB[pmindex(object)])
