setMethod("pmChr", "ExonFeatureSet",
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

setMethod("bgindex", "ExonFeatureSet",
          function(object, subset=NULL){
              conn <- db(object)
              sql <- paste("SELECT fid FROM",
                           "pmfeature, featureSet",
                           "WHERE pmfeature.fsetid=featureSet.fsetid",
                           "AND type > 1")
              fid <- dbGetQuery(conn, sql)[[1]]
              sort(fid)
          })

setMethod("bgSequence", "ExonFeatureSet",
          function(object){
              theFile <- file.path(system.file(package = annotation(object)), 
                                   "data", "pmSequence.rda")
              load(theFile)
              bgi <- bgindex(object)
              idx <- match(bgi, pmSequence[["fid"]])
              pmSequence[idx, "sequence"]
          })

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
            slot(out, "featureData") <- basicAnnotatedDataFrame(theExprs, byrow=TRUE)
            slot(out, "protocolData") <- protocolData(object)
            slot(out, "annotation") <- slot(object, "annotation")
            if (validObject(out)){
              return(out)
            }else{
              stop("Resulting object is invalid.")
            }
          })


gcCounts <- function(seq)
    as.integer(Biostrings::letterFrequency(seq, letters='CG'))
    
getRefDABG <- function(x){
    conn <- db(x)
    sql <- paste('SELECT type FROM type_dict',
                 'WHERE type_id="control->bgp->antigenomic"')
    bgCode <- as.integer(dbGetQuery(conn, sql)[[1]])
    if (length(bgCode) != 1)
        stop("Can't find proper id for control-bgp-antigenomic probes")
    
    sql <- paste('SELECT fid FROM pmfeature',
                 'INNER JOIN featureSet USING(fsetid)',
                 'WHERE', paste('type=', bgCode, sep=''))
    tmp <- dbGetQuery(conn, sql)
    env <- new.env()
    data(pmSequence, package=annotation(x), envir=env)
    idx <- match(tmp[['fid']], env[['pmSequence']][['fid']])
    seq <- env[['pmSequence']][idx, 'sequence']
    rm(env)
    counts <- gcCounts(seq)
    lst <- split(data.frame(exprs(x)[tmp[['fid']],, drop=FALSE]),
                 factor(counts, levels=0:max(counts)))
    lapply(lst, as.matrix)
}

computePSDABG <- function(x){
    paProbe <- -log(paCalls(x, method="DABG", verbose=FALSE))
    pns <- probeNames(x)
    n <- ncol(paProbe)+1
    theSum <- 2*rowsum(cbind(paProbe, 1), pns, reorder=FALSE)
    res <- pchisq(theSum[, -n, drop=FALSE], theSum[,n], lower.tail=FALSE)
    colnames(res) <- sampleNames(x)
    return(res)
}

computeDABG <- function(x){
    mmRef <- getRefDABG(x)
    pmCounts <- gcCounts(pmSequence(x))
    pms <- pm(x)
    ns <- sapply(mmRef, nrow)
    rgCounts <- range(as.integer(names(ns[ns > 0])))
    pmCounts[pmCounts < rgCounts[1]] <- rgCounts[1]
    pmCounts[pmCounts > rgCounts[2]] <- rgCounts[2]
    theClass <- class(exprs(x))
    if ("matrix" %in% theClass){
        res <- .Call("R_DABG_P", pms, mmRef, pmCounts)
    }else{
        res <- createFF("oligo-dabg-", dim(pms))
        rowsByNode <- splitIndicesByNode(1:nrow(pms))
        ocLapply(rowsByNode,
                 function(ss, pmMat, mmRef, pmCounts, outMat){
                     if (length(ss) > 0){
                         ps <- splitIndicesByLength(ss, ocProbesets())
                         open(pmMat)
                         open(outMat)
                         for (pps in ps)
                             outMat[pps,] <- .Call("R_DABG_P", pmMat[pps,], mmRef, pmCounts[pps])
                         close(outMat)
                         close(pmMat)
                     }
                    NULL
                }, pmMat=pms, mmRef=mmRef, pmCounts=pmCounts, outMat=res,
                 neededPkgs="oligo")
    }
    fid <- as.character(dbGetQuery(db(x), 'SELECT fid FROM pmfeature')[[1]])
    dimnames(res) <- list(fid, sampleNames(x))
    res
}

setMethod("paCalls", "ExonFeatureSet",
          function(object, method, verbose=TRUE){
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
