setMethod("probeNames", "ExonFeatureSet",
          function(object, subset=NULL){
            res <- dbGetQuery(db(object), "SELECT fid, fsetid FROM pmfeature")
            idx <- order(res[["fid"]])
            as.character(res[idx, "fsetid"])
          })

setMethod("rma", "ExonFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL, target="core"){
            target <- match.arg(target, c("core", "full", "extended", "probeset"))
            ## getFid...() will return results
            ## sorted by man_fsetid
            if (target == "core"){
              featureInfo <- getFidMetaProbesetCore(object)
            }else if (target == "full"){
              featureInfo <- getFidMetaProbesetFull(object)
            }else if (target == "extended"){
              featureInfo <- getFidMetaProbesetExtended(object)
            }else if (target == "probeset"){
              featureInfo <- getFidProbeset(object)
            }
            theClass <- class(exprs(object))
            pmi <- featureInfo[["fid"]]
            pnVec <- as.character(featureInfo[["fsetid"]])
            if ("matrix" %in% theClass){
              pms <- exprs(object)[pmi,, drop=FALSE]
              dimnames(pms) <- NULL
              colnames(pms) <- sampleNames(object)
              theExprs <- basicRMA(pms, pnVec, normalize, background)
              rm(pms)
            }else if ("ff_matrix" %in% theClass){
              pms <- ffSubset(rows=pmi, object=exprs(object), prefix="pm-")
              theExprs <- basicRMAbo(pms, pnVec, background=background, normalize=normalize)
              finalizer(pms) <- "delete"
              rm(pms)
            }else{
              stop("basicRMA not implemented for '", theClass, "' objects.")
            }

            out <- new("ExpressionSet")
            slot(out, "assayData") <- assayDataNew(exprs=theExprs)
            slot(out, "phenoData") <- phenoData(object)
            slot(out, "featureData") <- basicFeatureData(theExprs)
            slot(out, "protocolData") <- protocolData(object)
            slot(out, "annotation") <- slot(object, "annotation")
            if (validObject(out)){
              return(out)
            }else{
              stop("Resulting object is invalid.")
            }
          })
