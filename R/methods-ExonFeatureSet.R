setMethod("rma", "ExonFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL, target="core"){
            if (target == "core"){
              featureInfo <- getFidMetaProbesetCore(object)
            }else if (target == "full"){
              featureInfo <- getFidMetaProbesetFull(object)
            }else if (target == "extended"){
              featureInfo <- getFidMetaProbesetExtended(object)
            }else if (target == "probeset"){
              featureInfo <- getFidProbeset(object)
            }else{
              stop("Target '", target, "' is invalid. It must be one of: 'core', 'full', 'extended' or 'probeset'.")
            }
            pms <- exprs(object)[featureInfo[["fid"]],, drop=FALSE]
            dimnames(pms) <- NULL
            theExprs <- basicRMA(pms,
                                 as.character(featureInfo[["fsetid"]]),
                                 normalize, background)
            rm(pms)
            out <- new("ExpressionSet",
                       phenoData=phenoData(object),
                       annotation=annotation(object),
                       experimentData=experimentData(object),
                       exprs=theExprs)
            return(out)
          })
