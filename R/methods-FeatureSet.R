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

setMethod("pm", "FeatureSet",
          function(object, subset=NULL){
            theClass <- class(exprs(object))
            pmi <- pmindex(object, subset=subset)
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

setMethod("mm", "FeatureSet",
          function(object, subset=NULL){
            exprs(object)[mmindex(object, subset=subset),] ## subset 
          })

setReplaceMethod("mm", signature(object="FeatureSet", value="matrix"),
                 function(object, value){
                   tmp <- exprs(object)
                   tmp[mmindex(object),] <- value
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
            if (is.null(dots[["main"]])) changeMain <- TRUE
            if (is.null(dots[["col"]])) dots[["col"]] <- darkColors(ncol(x))
                
            eset <- x[idx,]
            rgs <- vector("list", nchns)
            for (i in 1:nchns){
                tmp <- log2(exprs(channel(eset, chns[i])))
                rgs[[i]] <- range(apply(tmp, 2, range))
                rm(tmp)
            }
            rgs <- range(unlist(rgs))
            dots[["ylim"]] <- rgs
            
            res <- vector("list", nchns)
            par(mfrow=c(nchns, 1))
            for (i in 1:nchns){
                tmp <- transfo(exprs(channel(eset, chns[i])))
                dots[["x"]] <- as.data.frame(tmp)
                dots[["main"]] <- chns[i]
                rm(tmp)
                res[[i]] <- do.call("boxplot", dots)
            }
            invisible(res)
          })
  
setMethod("boxplot", signature(x="ExpressionSet"),
          function(x, which, transfo=identity, nsample=10000, ...){
            stopifnot(is.function(transfo))
            if(!missing(which)) warning("Argument 'which' ignored.")
            if (nrow(x) > nsamples){
              idx <- sort(sample(nrow(x), nsamples))
            }else{
              idx <- 1:nrow(x)
            }
            toPlot <- transfo(exprs(x[idx,]))
            toPlot <- as.data.frame(toPlot)
            dots <- list(...)
            if (is.null(dots)[["cols"]])
              dots[["cols"]] <- darkColors(ncol(x))
            if (is.null(dots)[["range"]])
              dots[["range"]] <- 0
            do.call("boxplot", dots)
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
            par(mfrow=c(1, nchns))
            
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
                    do.call("image", dots)
                    rm(tmp)
                  }
                rm(tmpObj)
              }
            }
            par(ask=FALSE)
          })


matDensity <- function(mat){
  x.density <- apply(mat, 2, density, na.rm=TRUE)
  all.x <- sapply(x.density, "[[", "x")
  all.y <- sapply(x.density, "[[", "y")
  list(x=all.x, y=all.y)
}

getProbeIndex <- function(x, type=c("pm", "mm", "bg", "both", "all")){
  type <- match.arg(type)
  if (type == "pm"){
    idx <- pmindex(x)
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
            tmp <- lapply(chns, f, x[idx,])
            
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

setMethod("hist", "ExpressionSet",
          function(x, transfo=identity, nsample=10000, ...){
            stopifnot(is.function(transfo))
            if (nrow(x) > nsample){
              idx <- sort(sample(nrow(x), nsample))
            }else{
              idx <- 1:nrow(x)
            }
            tmp <- transfo(exprs(x[idx,]))
            res <- matDensity(tmp)

            dots <- list(...)
            if (is.null(dots[["ylab"]])) dots[["ylab"]] <- "density"
            if (is.null(dots[["xlab"]])) dots[["xlab"]] <- "log-intensity"
            if (is.null(dots[["xlim"]])) dots[["xlim"]] <- range(res[["x"]])
            if (is.null(dots[["ylim"]])) dots[["ylim"]] <- range(res[["y"]])
            if (is.null(dots[["col"]]))
              dots[["col"]] <- darkColors(ncol(x))
            if (is.null(dots[["type"]])) dots[["type"]] <- "l"

            dots[["x"]] <- res[["x"]]
            dots[["y"]] <- res[["y"]]
            do.call("matplot", dots)

            invisible(res)
          })

setMethod("MAplot", "FeatureSet",
          function(object, arrays=1:ncol(object), lowessPlot=FALSE, nsample=10000, ...){
            if (length(arrays) > 1) par(ask=TRUE)
            chns <- channelNames(object)
            nchns <- length(chns)
            if (nchns > 2) stop("Don't know how to handle more than 2 channels.")
            if (nrow(object) > nsample){
              idx <- sort(sample(nrow(object), nsample))
            }else{
              idx <- 1:nrow(object)
            }
            dots <- list(...)
            if (is.null(dots[["xlab"]]))
              dots[["xlab"]] <- "average log-intensity"
            if (is.null(dots[["ylab"]]))
              dots[["ylab"]] <- "log-ratio"
            if (is.null(dots[["pch"]]))
              dots[["pch"]] <- "."
            
            small <- object[idx,]
            
            if (nchns == 1){
              ref <- rowMedians(log2(exprs(small)))
              for (i in arrays){
                tmp <- log2(exprs(small[,i]))
                dots[["x"]] <- (tmp+ref)/2
                dots[["y"]] <- tmp-ref
                dots[["main"]] <- sampleNames(small)[i]
                do.call("plot", dots)
                if (lowessPlot)
                  lines(lowess(dots[["x"]], dots[["y"]]), col="red")
              }
            }else if (nchns == 2){
              for (i in arrays){
                dots[["x"]] <- (log2(exprs(channel(small, "channel1"))) + log2(exprs(channel(small, "channel1"))))/2
                dots[["y"]] <-  log2(exprs(channel(small, "channel1"))) - log2(exprs(channel(small, "channel1")))
                dots[["main"]] <- sampleNames(small)[i]
                if (lowessPlot)
                  lines(lowess(dots[["x"]], dots[["y"]]), col="red")
              }
            }
            if (length(arrays) > 1) par(ask=FALSE)
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
          function(object, subset=NULL){
            pmindex(getPlatformDesign(object), subset=subset)
          })

setMethod("mmindex", "FeatureSet",
          function(object, subset=NULL){
            mmindex(getPD(object), subset=subset)
          })

setMethod("probeNames", "FeatureSet",
          function(object, subset=NULL) {
            if (!is.null(subset))
              warning("ignoring subset arg, feature not implemented")
            probeNames(getPlatformDesign(object))
          })

setMethod("bgindex", "FeatureSet",
          function(object, subset=NULL){
            bgindex(getPD(object), subset=subset)
          })

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

