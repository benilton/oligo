setMethod("rma", "ExpressionFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL){
            pms <- pm(object, subset)
            pnVec <- probeNames(object, subset)
            tbls <- dbListTables(db(object))
            if (manufacturer(object) == "Affymetrix" && "bgfeature" %in% tbls){
              sql <- paste("SELECT man_fsetid, fid",
                           "FROM bgfeature",
                           "INNER JOIN featureSet",
                           "USING(fsetid)")
              tmpQcPm <- dbGetQuery(db(object), sql)
              qcpms <- exprs(object)[tmpQcPm[["fid"]],]
              pms <- rbind(pms, qcpms)
              pnVec <- c(pnVec, tmpQcPm[["man_fsetid"]])
            }
            idx <- order(pnVec)
            pms <- pms[idx,, drop=FALSE]
            dimnames(pms) <- NULL
            colnames(pms) <- sampleNames(object)
            pnVec <- pnVec[idx]
            exprs <- basicRMA(pms, pnVec, normalize, background)
            out <- new("ExpressionSet",
                       phenoData = phenoData(object),
                       annotation = annotation(object),
                       experimentData = experimentData(object),
                       exprs = exprs)
            return(out)
          })
