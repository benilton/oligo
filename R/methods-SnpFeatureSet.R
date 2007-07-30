setMethod("pmAllele", signature(object="SnpFeatureSet"),
          function(object) pmAllele(get(annotation(object))))

setMethod("db", "SnpFeatureSet",
          function(object) db(get(annotation(object))))
