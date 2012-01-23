setMethod("runDate", "FeatureSet",
          function(object){
              prD <- protocolData(object)
              idx <- grep("dates", varLabels(prD))
              pData(prD)[, idx]
          })

setMethod("intensity", "FeatureSet",
          function(object)
              exprs(object))

setReplaceMethod("intensity", signature(object="FeatureSet", value="matrix"),
                 function(object, value){
                   assayDataElementReplace(object, "exprs", value)
                 })

setMethod("probesetNames", "FeatureSet",
          function(object, target=NULL)
          unique(probeNames(object, target=target))
          )

          
## should geometry go to oligoClasses?
## or the opposite?
setMethod("geometry", "FeatureSet",
          function(object)
          geometry(getPD(object))
          )

setMethod("bg", "FeatureSet",
          function(object, subset=NULL){
            bgi <- bgindex(object, subset=subset)
            exprs(object[bgi,])
          })

setReplaceMethod("bg", signature(object="FeatureSet", value="matrix"),
                 function(object, value){
                   tmp <- exprs(object)
                   tmp[bgindex(object),] <- value
                   assayDataElementReplace(object, "exprs", tmp)
                 })

setReplaceMethod("bg", signature(object="FeatureSet", value="ff_matrix"),
                 function(object, value){
                   tmp <- exprs(object)
                   open(tmp)
                   open(value)
                   finalizer(value) <- "delete"
                   for (i in 1:ncol(tmp))
                       tmp[bgindex(object), i] <- value[, i]
                   close(value)
                   close(tmp)
                   rm(value)
                   assayDataElementReplace(object, "exprs", tmp)
                 })

## TODO: repeat target for mm()

setMethod("pm", "FeatureSet",
          function(object, subset=NULL, target='core'){
            theClass <- class(exprs(object))
            pmi <- pmindex(object, subset=subset, target=target)
            if ("matrix" %in% theClass){
              out <- exprs(object)[pmi,, drop=FALSE]
            }else if ("ff_matrix" %in% theClass){
              out <- ffSubset(rows=pmi, object=exprs(object),
                              prefix="pm-")
            }
            return(out)
          })

setReplaceMethod("pm", signature(object="FeatureSet", value="matrix"),
                 function(object, value){
                   tmp <- exprs(object)
                   tmp[pmindex(object),] <- value
                   assayDataElementReplace(object, "exprs", tmp)
                 })

setReplaceMethod("pm", signature(object="FeatureSet", value="ff_matrix"),
                 function(object, value){
                   tmp <- exprs(object)
                   open(tmp)
                   open(value)
                   finalizer(value) <- "delete"
                   nc <- ncol(tmp)
                   for (i in 1:nc)
                       tmp[pmindex(object), i] <- value[,i]
                   close(value)
                   close(tmp)
                   rm(value)
                   assayDataElementReplace(object, "exprs", tmp)
                 })

setMethod("mm", "FeatureSet",
          function(object, subset=NULL){
            exprs(object)[mmindex(object, subset=subset),, drop=FALSE] ## subset 
          })

setReplaceMethod("mm", signature(object="FeatureSet", value="matrix"),
                 function(object, value){
                   tmp <- exprs(object)
                   tmp[mmindex(object),] <- value
                   assayDataElementReplace(object, "exprs", tmp)
                 })

setReplaceMethod("mm", signature(object="FeatureSet", value="ff_matrix"),
                 function(object, value){
                   tmp <- exprs(object)
                   open(tmp)
                   open(value)
                   finalizer(value) <- "delete"
                   nc <- ncol(tmp)
                   for (i in 1:nc)
                       tmp[mmindex(object), i] <- value[,i]
                   close(value)
                   close(tmp)
                   rm(value)
                   assayDataElementReplace(object, "exprs", tmp)
                 })

setMethod("pmSequence", "FeatureSet",
          function(object) pmSequence(getPD(object)))

setMethod("mmSequence", "FeatureSet",
          function(object) mmSequence(getPD(object)))

## QAish functions
setMethod("boxplot", signature(x="FeatureSet"),
          function(x, which=c("pm", "mm", "bg", "both", "all"),
                   transfo=log2, nsample=10000, ...){
            stopifnot(is.function(transfo))
            idx <- getProbeIndex(x, which)
            if (length(idx) > nsample)
                idx <- sort(sample(idx, nsample))

            chns <- channelNames(x)
            nchns <- length(chns)

            dots <- list(...)
            if (is.null(dots[["range"]])) dots[["range"]] <- 0
            changeMain <- is.null(dots[["main"]])
            if (is.null(dots[["col"]])) dots[["col"]] <- darkColors(ncol(x))
            
            eset <- x[idx,]

            ## this fix is temporary
            ## until we agree on how
            ## to handle multiplicity by gene chips
            featureNames(eset) <- featureNames(featureData(eset))
            
            rgs <- vector("list", nchns)
            for (i in 1:nchns){
                tmp <- transfo(exprs(channel(eset, chns[i])))
                rgs[[i]] <- range(apply(tmp, 2, range))
                rm(tmp)
            }
            rgs <- range(unlist(rgs))
            dots[["ylim"]] <- rgs
            
            res <- vector("list", nchns)
            doPlot <- dots[['plot']]
            if (is.null(doPlot)){
                doPlot <- TRUE
            }else{
                doPlot <- as.logical(doPlot)
            }
            if (doPlot & nchns > 1)
                par(mfrow=c(nchns, 1))
            for (i in 1:nchns){
                tmp <- transfo(exprs(channel(eset, chns[i])))
                dots[["x"]] <- as.data.frame(tmp)
                if (changeMain)
                    dots[["main"]] <- chns[i]
                rm(tmp)
                res[[i]] <- do.call("boxplot", dots)
            }
            invisible(res)
          })
  
setMethod("image", signature(x="FeatureSet"),
          function(x, which, transfo=log2, ...){
            if (missing(which)) which <- 1:ncol(x)
            if(length(which) > 1) par(ask=TRUE) else par(ask=FALSE)
            dots <- list(...)
            
            if (is.null(dots[["col"]]))
              dots[["col"]] <- colorRampPalette(c("blue", "white", "red"))(256)
            if (is.null(dots[["xaxt"]]))
              dots[["xaxt"]] <- "n"
            if (is.null(dots[["yaxt"]]))
              dots[["yaxt"]] <- "n"
            
            geom <- geometry(getPD(x))
            
            chns <- channelNames(x)
            nchns <- length(chns)
            oldpar <- par()[c("mfrow", "mar", "mgp")]
            par(mfrow=c(nchns, 1), mar=c(2.5,2.5,1.6,1.1), mgp=c(1.5,.5,0))
            
            if (tolower(manufacturer(x)) != "affymetrix"){
              conn <- db(x)
              tbls <- dbGetQuery(conn, "SELECT tbl FROM table_info WHERE tbl LIKE '%feature' AND row_count > 0")[[1]]
              theInfo <- lapply(tbls, function(tb) dbGetQuery(conn, paste("SELECT x, y, fid FROM", tb)))
              theInfo <- do.call("rbind", theInfo)
              theInfo <- theInfo[order(theInfo[["fid"]]), ]
              idx <- geom[1]*(theInfo[["x"]]-1)+theInfo[["y"]]
              for (i in which){
                tmpObj <- x[theInfo[["fid"]], i]
                for (j in 1:nchns){
                  dots[["main"]] <- paste(sampleNames(x)[i],
                                          chns[j], sep=" - ")
                  tmp <- matrix(NA, nr=geom[1], nc=geom[2])
                  tmp[idx] <- transfo(as.numeric(exprs(channel(tmpObj, chns[j]))))
                  tmp <- as.matrix(rev(as.data.frame(tmp)))
                  dots[["x"]] <- tmp
                  do.call("image", dots)
                  rm(tmp)
                }
                rm(tmpObj)
              }
            }else{
              for (i in which){
                tmpObj <- x[, i]
                for (j in 1:nchns){
                    dots[["main"]] <- paste(sampleNames(x)[i],
                                            chns[j], sep=" - ")
                    tmp <- transfo(as.numeric(exprs(channel(tmpObj, chns[j]))))
                    tmp <- matrix(tmp, ncol=geom[1], nrow=geom[2])
                    tmp <- as.matrix(rev(as.data.frame(tmp)))
                    dots[["x"]] <- tmp
                    do.call("image", dots)
                    rm(tmp)
                  }
                rm(tmpObj)
              }
            }
            par(ask=FALSE)
            par(oldpar)
          })


matDensity <- function(mat){
  x.density <- apply(mat, 2, density, na.rm=TRUE)
  all.x <- sapply(x.density, "[[", "x")
  all.y <- sapply(x.density, "[[", "y")
  list(x=all.x, y=all.y)
}

getProbeIndex <- function(x, type=c("pm", "mm", "bg", "both", "all"), target='core'){
  type <- match.arg(type)
  if (type == "pm"){
    idx <- pmindex(x, target=target)
  }else if (type == "mm"){
    idx <- mmindex(x)
  }else if (type == "bg"){
    idx <- bgindex(x)
  }else if (type == "both"){
    warning("Argument 'both' was replaced by 'all'. Please update your code.")
    idx <- 1:nrow(x)
  }else if (type == "all"){
    idx <- 1:nrow(x)
  }
  idx
}

setMethod("hist", "FeatureSet",
          function(x, transfo=log2, which=c("pm", "mm", "bg", "both", "all"),
                   nsample=10000, ...){
            stopifnot(is.function(transfo))
            idx <- getProbeIndex(x, which)
            if (length(idx) > nsample)
              idx <- sort(sample(idx, nsample))
            
            chns <- channelNames(x)
            nchns <- length(chns)
            
            ## estimate density for every sample on each channel
            f <- function(chn, obj)
              matDensity(transfo(exprs(channel(obj, chn))))

            eset <- x[idx,]
            ## this fix is temporary
            ## until we agree on how
            ## to handle multiplicity by gene chips
            featureNames(eset) <- featureNames(featureData(eset))

            tmp <- lapply(chns, f, eset)
            rm(eset)
            
            ## get lims
            rgs <- lapply(tmp, sapply, range)
            rgs <- do.call("rbind", rgs)
            rgs <- apply(rgs, 2, range)
            
            ## set graph options properly
            dots <- list(...)
            if (is.null(dots[["ylab"]])) dots[["ylab"]] <- "density"
            if (is.null(dots[["xlab"]])) dots[["xlab"]] <- "log-intensity"
            if (is.null(dots[["xlim"]])) dots[["xlim"]] <- rgs[,1]
            if (is.null(dots[["ylim"]])) dots[["ylim"]] <- rgs[,2]
            if (is.null(dots[["col"]]))
              dots[["col"]] <- darkColors(ncol(x))
            if (is.null(dots[["type"]])) dots[["type"]] <- "l"
            par(mfrow=c(nchns, 1))
            changeMain <- is.null(dots[["main"]])
            
            for (i in 1:nchns){
              dots[["x"]] <- tmp[[i]][["x"]]
              dots[["y"]] <- tmp[[i]][["y"]]
              if (changeMain)
                dots[["main"]] <- chns[i]
              do.call("matplot", dots)
            }
            invisible(tmp)
          })

setMethod("MAplot", "FeatureSet",
          function(object, what=pm, transfo=log2, groups, refSamples, which,
                   pch=".", summaryFun=rowMedians, plotFun=smoothScatter,
                   main="vs pseudo-median reference chip", pairs=FALSE, ...){
              stopifnot(is.function(what))
              maplot(x=what(object), transfo=transfo, groups=groups,
                     refSamples=refSamples, which=which, pch=pch,
                     summaryFun=summaryFun, plotFun=plotFun, main=main, pairs=pairs, ...)
          })

setMethod("getX", "FeatureSet",
          function(object, type){
            getX(getPD(object), type)
          })

setMethod("getY", "FeatureSet",
          function(object, type){
            getY(getPD(object), type)
          })

setMethod("bgSequence",
          signature(object="FeatureSet"),
          function(object){
            bgSequence(getPlatformDesign(object))
          })

setMethod("pmSequence",
          signature(object="FeatureSet"),
          function(object){
            pmSequence(getPlatformDesign(object))
          })

setMethod("mmSequence",
          signature(object="FeatureSet"),
          function(object){
            mmSequence(getPlatformDesign(object))
          })

setMethod("genomeBuild",
          signature(object="FeatureSet"),
          function(object){
            genomeBuild(getPlatformDesign(object))
          })

setMethod("pmChr", "FeatureSet",
          function(object){
            conn <- db(object)
            if (is(object, "TilingFeatureSet") & manufacturer(object) == "Affymetrix"){
              sql <- paste("SELECT fid, chrom_id as chrom",
                           "FROM pmfeature",
                           "INNER JOIN chrom_dict",
                           "USING(chrom)")
            }else{
              sql <- paste("SELECT fid, chrom",
                           "FROM pmfeature, featureSet",
                           "WHERE pmfeature.fsetid=featureSet.fsetid")
            }
            tmp <- dbGetQuery(conn, sql)
            tmp <- tmp[order(tmp[["fid"]]),]
            tmp[["chrom"]]
          })

setMethod("pmindex", "FeatureSet",
          function(object, subset=NULL, target='core'){
            pmindex(getPlatformDesign(object), subset=subset, target=target)
          })

setMethod("mmindex", "FeatureSet",
          function(object, subset=NULL, target='core'){
            mmindex(getPD(object), subset=subset)
          })

setMethod("probeNames", "FeatureSet",
          function(object, subset=NULL, target='core') {
            if (!is.null(subset))
              warning("ignoring subset arg, feature not implemented")
            probeNames(getPlatformDesign(object), target=target)
          })

setMethod("bgindex", "FeatureSet",
          function(object, subset=NULL, target='core'){
            bgindex(getPD(object), subset=subset)
          })

## TODO: fix pmPosition to account for target
setMethod("pmPosition", "FeatureSet",
          function(object){
            conn <- db(object)
            sql <- paste("SELECT fid, position",
                         "FROM pmfeature",
                         "INNER JOIN featureSet USING(fsetid)")
            tmp <- dbGetQuery(conn, sql)
            tmp <- tmp[order(tmp[["fid"]]),]
            tmp[["position"]]
          })

setMethod("kind",
          signature(object="FeatureSet"),
          function(object){
            kind(getPlatformDesign(object))
          })

setMethod("getPlatformDesign",
          signature(object= "FeatureSet"),
          function(object){
            pdn <- annotation(object)
            library(pdn,character.only=TRUE)
            return(get(pdn,pos=paste("package:",pdn,sep="")))
          })
getPD <- getPlatformDesign

getFragmentLength <- function(object, enzyme, type=c('snp', 'cn')){
    ## returns data.frame with 3 columns:
    ## SNP-CNP | Enzyme | Length
    conn <- db(get(annotation(object)))
    type <- match.arg(type)
    suffix <- ifelse(type=='snp', '', 'CNV')
    probetbl <- paste('pmfeature', suffix, sep='')
    fltbl <- paste('fragmentLength', suffix, sep='')
    fsettbl <- paste('featureSet', suffix, sep='')
    sql <- paste('SELECT man_fsetid, enzyme, length', 'FROM',
                 fltbl, 'INNER JOIN', fsettbl, 'USING(fsetid)')
    if (!missing(enzyme)){
        stopifnot(length(enzyme) == 1)
        sql <- paste(sql, 'WHERE enzyme =', enzyme)
    }
    dbGetQuery(conn, sql)
}
