setMethod("rma", "GeneFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL, target="core"){
            if (target == "core"){
              featureInfo <- getFidMetaProbesetCore(object)
            }else if (target == "probeset"){
              featureInfo <- getFidProbeset(object)
            }else{
              stop("Target '", target, "' is invalid. It must be either 'core' or 'probeset'.")
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
