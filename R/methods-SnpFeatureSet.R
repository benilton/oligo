setMethod("pmAllele", signature(object="SnpFeatureSet"),
          function(object) pmAllele(get(annotation(object))))
