setMethod("rma", "ExpressionFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL){
            ## get pmi and pnVec
            pmi <- pmindex(object)
            pnVec <- probeNames(object)
            tbls <- dbListTables(db(object))
            if (manufacturer(object) == "Affymetrix" && "bgfeature" %in% tbls){
              sql <- paste("SELECT man_fsetid, fid",
                           "FROM bgfeature",
                           "INNER JOIN featureSet",
                           "USING(fsetid)")
              tmpQcPm <- dbGetQuery(db(object), sql)
              pmi <- c(pmi, tmpQcPm[["fid"]])
              pnVec <- c(pnVec, tmpQcPm[["man_fsetid"]])
            }
            idx <- order(pnVec)
            pnVec <- pnVec[idx]
            pmi <- pmi[idx]
            rm(idx)
            
            theClass <- class(exprs(object))

            if (theClass == "matrix"){
              pms <- exprs(object[pmi,])
              dimnames(pms) <- NULL
              colnames(pms) <- sampleNames(object)
              exprs <- basicRMA(pms, pnVec, normalize, background)
            }else if (theClass == "big.matrix"){
              dUID <- getDatasetUID(object)
              pmFile <- paste("pmPP-", dUID, sep="")
              pmName <- "pms"
              path <- oligoBigObjectPath()
              assign(pmName, subsetBO(pmi, object=describe(exprs(object)),
                                      fname=pmFile, nameInEnv=pmName, clean=FALSE))
              exprs <- basicRMAbo(pms, pnVec, background=background,
                                  normalize=normalize, pmName=pmName,
                                  dUID=dUID)
              rmFromPkgEnv(pmName)
              pmfns <- list.files(path, patt=paste("^", pmFile, sep=""), full=TRUE)
              unlink(pmfns)
            }else{
              stop("basicRMA not implemented for '", theClass, "' objects.")
            }
            
            out <- new("ExpressionSet")
            slot(out, "assayData") <- assayDataNew(exprs=exprs)
            slot(out, "phenoData") <- phenoData(object)
            slot(out, "featureData") <- basicFeatureData(exprs)
            slot(out, "protocolData") <- protocolData(object)
            slot(out, "annotation") <- slot(object, "annotation")
            if (validObject(out)){
              return(out)
            }else{
              stop("Resulting object is invalid.")
            }
          })
