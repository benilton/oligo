setMethod("rma", "GeneFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL, target="core"){
            target <- match.arg(target, c("core", "probeset"))
            if (target == "core"){
              featureInfo <- getFidMetaProbesetCore(object)
            }else if (target == "probeset"){
              featureInfo <- getFidProbeset(object)
            }
            theClass <- class(exprs(object))
            pmi <- featureInfo[["fid"]]
            pnVec <- as.character(featureInfo[["fsetid"]])
            if (theClass == "matrix"){
              pms <- exprs(object)[pmi,, drop=FALSE]
              dimnames(pms) <- NULL
              colnames(pms) <- sampleNames(object)
              theExprs <- basicRMA(pms, pnVec, normalize, background)
              rm(pms)
            }else if (theClass == "big.matrix"){
              dUID <- getDatasetUID(object)
              pmFile <- paste("pmPP-", dUID, sep="")
              pmName <- "pms"
              path <- oligoBigObjectPath()
              assign(pmName, subsetBO(pmi, object=describe(exprs(object)),
                                      fname=pmFile, nameInEnv=pmName, clean=FALSE))
              theExprs <- basicRMAbo(pms, pnVec, background=background,
                                     normalize=normalize, pmName=pmName,
                                     dUID=dUID)
              rmFromPkgEnv(pmName)
              pmfns <- list.files(path, patt=paste("^", pmFile, sep=""), full=TRUE)
              unlink(pmfns)
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
