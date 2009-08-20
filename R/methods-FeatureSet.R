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
            exprs(object)[pmindex(object, subset=subset),,drop=FALSE]
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
          function(x, which=pmindex(x), transfo=log2, range=0, ...){
            toPlot <- exprs(x[which,])
            toPlot <- transfo(toPlot)
            toPlot <- as.data.frame(toPlot)
            boxplot(toPlot, range=range, ...)
          })
  
setMethod("boxplot", signature(x="ExpressionSet"),
          function(x, which=1:nrow(x), transfo=identity, range=0, ...){
            toPlot <- exprs(x[which,])
            toPlot <- transfo(toPlot)
            toPlot <- as.data.frame(toPlot)
            boxplot(toPlot, range=range, ...)
          })
  
  
## setMethod("image", signature(x="FeatureSet"),
## 		  function(x, which=0, transfo=log2, col=gray((0:64)/64), ...){
##               if (which == 0){
##                   which <- 1:length(x)
##               }
##               if(length(which) > 1) par(ask=TRUE) else par(ask=FALSE)
##               theDim <- geometry(getPD(x))
##               for (i in which){
##                 tmp <- matrix(transfo(exprs(x[,i])), ncol=theDim[1], nrow=theDim[2])
##                 tmp <- t(tmp[nrow(tmp):1, ncol(tmp):1])
##                 image(tmp, col=col, main=sampleNames(x)[i], xaxt="n", yaxt="n", ...)
##               }
##               par(ask=FALSE)
##             })


setMethod("image", signature(x="FeatureSet"),
          function(x, which=0, transfo=log2, col=gray((0:64)/64), ...){
            if (which == 0){
              which <- 1:length(x)
            }
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
                   which, ylab="density", xlab="log intensity",
                   type="l", ...){
            stopifnot(is.function(transfo))
            if (!missing(which)) warning("Argument 'which' not implemented yet.")
            idx <- sample(nrow(x), min(c(100000, nrow(x))))
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
                   which, ylab="density", xlab="log intensity",
                   type="l", ...){
            stopifnot(is.function(transfo))
            if (!missing(which)) warning("Argument 'which' not implemented yet.")
            idx <- sample(nrow(x), min(c(100000, nrow(x))))
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
