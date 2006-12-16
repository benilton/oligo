setMethod("manufacturer", "PDInfo",
          function(object) object@manufacturer)

setMethod("genomeBuild", "PDInfo",
                    function(object) object@genomebuild)

setMethod("db", "DBPDInfo",
          function(object) object@getdb())


setMethod("nrow", "platformDesign", function(x) x@nrow)
setMethod("ncol", "platformDesign", function(x) x@ncol)
setMethod("type", "platformDesign", function(object) object@type)
