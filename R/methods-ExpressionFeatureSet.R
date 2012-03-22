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


## original code from affy
paMAS5 <- function(object, verbose=TRUE,
                   tau=0.015, alpha1=0.04, alpha2=0.06,
                   ignore.saturated=TRUE) {

    stopifnot(alpha1 > 0, alpha1 < alpha2, alpha2 < 1)
    if(verbose) message("Getting probe level data... ", appendLF=FALSE);
    pms <- pm(object)
    mms <- mm(object)
    if (verbose) message("OK.")

    ## Saturation:
    ## shouldn't be a problem with new scanners
    ## or those that have had an engineer visit
    sat <- ifelse(ignore.saturated, 46000, -1)

    pns <- probeNames(object)
    o <- order(pns)
    pns <- pns[o]
    pms <- pms[o,,drop=FALSE]
    mms <- mms[o,,drop=FALSE]
    np <- nrow(mms)
    unique.pns <- sort(unique(pns))
    nps <- length(unique.pns)
    nsamples <- ncol(pms)

    if(verbose) message("Computing p-values... ", appendLF=FALSE)
    p <- sapply(1:nsamples,
                function(x){
                    .C("DetectionPValue", as.double(pms[,x]),
                       as.double(mms[,x]), as.character(pns),
                       as.integer(np), as.double(tau),
                       as.double(sat), dpval=double(nps),
                       nps, PACKAGE="oligo")$dpval
                })
    rownames(p) <- unique.pns
    colnames(p) <- sampleNames(object)
    if (verbose) message("OK.")
    if (verbose) message("Making P/M/A Calls... ", appendLF=FALSE)
    calls <- matrix("A", ncol=ncol(p), nrow=nrow(p))
    calls[p < alpha1] <- "P"
    calls[p <= alpha2 & p >= alpha1] <- "M"
    dimnames(calls) <- list(rownames(p), sampleNames(object))
    if (verbose) message("OK.")
    list(calls=calls, p=p)
}


setMethod("paCalls", "ExpressionFeatureSet",
          function(object, method, ..., verbose=TRUE){
              stopifnot(tolower(manufacturer(object)) == 'affymetrix')
              if (missing(method))
                  method <- "MAS5"
              method <- match.arg(method, "MAS5")
              paFun <- switch(method, MAS5=paMAS5)
              res <- paFun(object, ..., verbose=verbose)
              return(res)
          })
