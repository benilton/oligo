setMethod("rma", "ExpressionFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL){
            pms <- pm(object, subset)
            pnVec <- probeNames(object, subset)
            if (manufacturer(object) == "Affymetrix"){
              sql <- paste("SELECT man_fsetid, fid",
                           "FROM bgfeature",
                           "INNER JOIN featureSet",
                           "USING(fsetid)")
              tmpQcPm <- dbGetQuery(db(object), sql)
              qcpms <- exprs(object)[tmpQcPm[["fid"]],]
              pms <- rbind(pms, qcpms)
              pnVec <- c(pnVec, tmpQcPm[["man_fsetid"]])
            }
            ngenes <- length(unique(pnVec))
            idx <- order(pnVec)
            pms <- pms[idx,, drop=FALSE]
            dimnames(pms) <- NULL
            pnVec <- pnVec[idx]
            exprs <- basicRMA(pms, pnVec, ngenes, normalize, background)
            colnames(exprs) <- sampleNames(object)
            rownames(exprs) <- unique(pnVec)
            out <- new("ExpressionSet",
                       phenoData = phenoData(object),
                       annotation = annotation(object),
                       experimentData = experimentData(object),
                       exprs = exprs)
            return(out)
          })

