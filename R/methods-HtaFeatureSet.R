stArrayPmInfo <- function(object, target='core', sortBy='fsetid'){
    ## *PmInfo returns a data.frame with 'fid' and 'fsetid'
    target <- match.arg(target, c('probeset', 'core', 'full', 'extended'))
    theFun <- switch(target,
                     probeset=getFidProbeset,
                     core=getFidMetaProbesetCore,
                     full=getFidMetaProbesetFull,
                     extended=getFidMetaProbesetExtended)
    theFun(object, sortBy=sortBy)
}


setMethod("rma", "HTAFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL, target="core"){
            target <- match.arg(target, c("core", "probeset"))
            featureInfo <- stArrayPmInfo(object, target=target)
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
            slot(out, "featureData") <- basicAnnotatedDataFrame(theExprs, byrow=TRUE)
            slot(out, "protocolData") <- protocolData(object)
            slot(out, "annotation") <- slot(object, "annotation")
            if (validObject(out)){
              return(out)
            }else{
              stop("Resulting object is invalid.")
            }
          })
