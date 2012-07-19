setMethod("pmChr", "GeneFeatureSet",
          function(object){
            conn <- db(object)
            sql <- paste("SELECT fid, chrom",
                         "FROM pmfeature",
                         "INNER JOIN featureSet USING(fsetid)")
            tmp <- dbGetQuery(conn, sql)
            chromInfo <- dbGetQuery(conn, "SELECT * FROM chrom_dict")
            tmp <- merge(tmp, chromInfo, by.x="chrom",
                         by.y="chrom", all.x=TRUE, all.y=FALSE,
                         sort=FALSE)
            tmp <- tmp[order(tmp[["fid"]]),]
            tmp[["chrom_id"]]
          })

setMethod("bgindex", "GeneFeatureSet",
          function(object, subset=NULL){
              conn <- db(object)
              sql <- paste("SELECT fid FROM",
                           "pmfeature, featureSet",
                           "WHERE pmfeature.fsetid=featureSet.fsetid",
                           "AND type > 1")
              fid <- dbGetQuery(conn, sql)[[1]]
              sort(fid)
          })

setMethod("bgSequence", "GeneFeatureSet",
          function(object){
              theFile <- file.path(system.file(package = annotation(object)),
                                   "data", "pmSequence.rda")
              load(theFile)
              bgi <- bgindex(object)
              idx <- match(bgi, pmSequence[["fid"]])
              pmSequence[idx, "sequence"]
          })

setMethod("pmSequence", "GeneFeatureSet",
          function(object, target='core'){
              pmSequence(getPD(object), target=target)
          })


## setMethod("probeNames", "GeneFeatureSet",
##           function(object, subset=NULL){
##             res <- dbGetQuery(db(object), "SELECT fid, fsetid FROM pmfeature")
##             idx <- order(res[["fid"]])
##             as.character(res[idx, "fsetid"])
##           })

setMethod("rma", "GeneFeatureSet",
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

setMethod("paCalls", "GeneFeatureSet",
          function(object, method=c("DABG", "PSDABG"), verbose=TRUE){
              if (missing(method))
                  method <- "DABG"
              method <- match.arg(method, c("DABG", "PSDABG"))
              paFun <- switch(method,
                              DABG=computeDABG,
                              PSDABG=computePSDABG)
              if (verbose) message("Computing DABG calls... ", appendLF=FALSE)
              res <- paFun(object)
              if (verbose) message("OK")
              return(res)
          })
