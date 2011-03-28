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

            if ("matrix" %in% theClass){
              pms <- exprs(object[pmi,])
              dimnames(pms) <- NULL
              colnames(pms) <- sampleNames(object)
              exprs <- basicRMA(pms, pnVec, normalize, background)
            }else if ("ff_matrix" %in% theClass){
              pms <- ffSubset(rows=pmi, object=exprs(object), prefix="pm-")
              exprs <- basicRMAbo(pms, pnVec, background=background, normalize=normalize)
              finalizer(pms) <- "delete"
              rm(pms)
            }else{
              stop("basicRMA not implemented for '", theClass, "' objects.")
            }
            
            out <- new("ExpressionSet")
            slot(out, "assayData") <- assayDataNew(exprs=exprs)
            slot(out, "phenoData") <- phenoData(object)
            slot(out, "featureData") <- basicAnnotatedDataFrame(exprs, byrow=TRUE)
            slot(out, "protocolData") <- protocolData(object)
            slot(out, "annotation") <- slot(object, "annotation")
            if (validObject(out)){
              return(out)
            }else{
              stop("Resulting object is invalid.")
            }
          })


setMethod("getNetAffx", "ExpressionSet",
          function(object, type="probeset"){
              type <- match.arg(type, c("probeset", "transcript"))
              fname <- ifelse(type == "probeset",
                              "Probeset",
                              "Transcript")
              fname <- paste("netaffx", fname, ".rda", sep="")
              fname <- file.path(system.file("extdata", package=annotation(object)),
                                 fname)

              if (!file.exists(fname))
                  stop("NetAffx Annotation not available in '",
                       annotation(object), "'. Consider using 'biomaRt'.")

              obj <- load(fname)
              fdata <- get(obj)
              rm(list=obj)
              fns <- featureNames(object)
              theSet <- fdata[fns,]
              featureNames(theSet) <- fns
              theSet
          })

