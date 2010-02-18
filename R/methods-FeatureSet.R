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
            if (theClass == "matrix"){
              out <- exprs(object)[pmi,, drop=FALSE]
            }else if (theClass == "ff_matrix"){
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
          function(x, which=c("pm", "mm", "both"),
                   transfo=log2, range=0, ...){
            stopifnot(is.function(transfo))
            which <- match.arg(which)
            if (which == "pm"){
              idx <- pmindex(x)
            }else if (which == "mm"){
              idx <- mmindex(x)
            }else if (which == "both"){
              idx <- 1:nrow(x)
            }else{
              stop("Invalid value for 'which'. Allowed values are: 'pm', 'mm', 'both'")
            }
            toPlot <- exprs(x[idx,])
            toPlot <- transfo(toPlot)
            toPlot <- as.data.frame(toPlot)
            boxplot(toPlot, range=range, ...)
          })
  
setMethod("boxplot", signature(x="ExpressionSet"),
          function(x, which, transfo=identity, range=0, ...){
            stopifnot(is.function(transfo))
            if(!missing(which)) warning("Argument 'which' ignored.")
            toPlot <- exprs(x)
            toPlot <- transfo(toPlot)
            toPlot <- as.data.frame(toPlot)
            boxplot(toPlot, range=range, ...)
          })

setMethod("image", signature(x="FeatureSet"),
          function(x, which, transfo=log2, col=gray((0:64)/64), ...){
            if (missing(which)) which <- 1:ncol(x)
            if(length(which) > 1) par(ask=TRUE) else par(ask=FALSE)
            geom <- geometry(getPD(x))
            if (tolower(manufacturer(x)) != "affymetrix"){
              conn <- db(x)
              tbls <- dbGetQuery(conn, "SELECT tbl FROM table_info WHERE tbl LIKE '%feature' AND row_count > 0")[[1]]
              theInfo <- lapply(tbls, function(tb) dbGetQuery(conn, paste("SELECT x, y, fid FROM", tb)))
              theInfo <- do.call("rbind", theInfo)
              theInfo <- theInfo[order(theInfo[["fid"]]), ]
              idx <- geom[1]*(theInfo[["x"]]-1)+theInfo[["y"]]
              for (i in which){
                theInfo[["signal"]] <- transfo(as.numeric(exprs(x[theInfo[["fid"]], i])))
                tmp <- matrix(NA, nr=geom[1], nc=geom[2])
                tmp[idx] <- theInfo[["signal"]]
                tmp <- as.matrix(rev(as.data.frame(tmp)))
                image(tmp, col=col, main=sampleNames(x)[i], xaxt="n", yaxt="n", ...)
              }
            }else{
              for (i in which){
                tmp <- matrix(transfo(exprs(x[,i])), ncol=geom[1], nrow=geom[2])
                tmp <- as.matrix(rev(as.data.frame(tmp)))
                image(tmp, col=col, main=sampleNames(x)[i], xaxt="n", yaxt="n", ...)
              }
            }
            par(ask=FALSE)
          })


setMethod("hist", "FeatureSet",
          function(x, col=1:ncol(x), transfo=log2,
                   which=c("pm", "mm", "both"),
                   ylab="density", xlab="log intensity",
                   type="l", ...){
            stopifnot(is.function(transfo))
            which <- match.arg(which)
            if (which == "pm"){
              idx <- pmindex(x)
            }else if (which == "mm"){
              idx <- mmindex(x)
            }else if (which == "both"){
              idx <- 1:nrow(x)
            }else{
              stop("Invalid value for 'which'. Allowed values are: 'pm', 'mm', 'both'")
            }
            tmp <- exprs(x[idx,])
            idx <- is.na(tmp[,1])
            if(any(idx))
              tmp <- tmp[!idx,, drop=FALSE]
            tmp <- transfo(tmp)
            x.density <- apply(tmp, 2, density)
            all.x <- sapply(x.density, "[[", "x")
            all.y <- sapply(x.density, "[[", "y")
            matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, ...)
            invisible(x.density)
          })

setMethod("hist", "ExpressionSet",
          function(x, col=1:ncol(x), transfo=log2,
                   which=1:nrow(x), ylab="density",
                   xlab="log intensity",
                   type="l", ...){
            stopifnot(is.function(transfo))
            tmp <- exprs(x[which,])
            idx <- is.na(tmp[,1])
            if(any(idx))
              tmp <- tmp[!idx,, drop=FALSE]
            tmp <- transfo(tmp)
            x.density <- apply(tmp, 2, density)
            all.x <- sapply(x.density, "[[", "x")
            all.y <- sapply(x.density, "[[", "y")
            matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, ...)
            invisible(x.density)
          })

setMethod("MAplot", "FeatureSet",
          function(object, arrays=1:ncol(object), lowessPlot=FALSE, ...){
            if (length(arrays) > 1) par(ask=TRUE)
            ref <- rowMedians(log2(exprs(object)))
            for (i in arrays){
              tmp <- log2(exprs(object[,i]))
              plot((tmp+ref)/2, tmp-ref, pch=".",
                   xlab="A", ylab="M",
                   main=sampleNames(object)[i], ...)
              if (lowessPlot)
                lines(lowess((tmp+ref)/2, tmp-ref), col="red")
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

