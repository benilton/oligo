setMethod("pmAllele", signature(object="SnpFeatureSet"),
          function(object) pmAllele(get(annotation(object))))

setMethod("pmFragmentLength", "SnpFeatureSet",
          function(object, enzyme, type=c('snp', 'cn')){
              pmFragmentLength(getPD(object), enzyme, type)
          })
