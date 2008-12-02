setMethod("initialize", "FeatureSet",
          function(.Object,
                   assayData = assayDataNew(exprs=exprs, ...),
                   manufacturer=new("character"),
                   platform=new("character"),
                   exprs=new("matrix"),
                   phenoData = annotatedDataFrameFrom  (assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData = new("MIAME"),
                   annotation = new("character"), ...){
            .Object <- callNextMethod(.Object,
                                      assayData = assayData,
                                      phenoData = phenoData,
                                      experimentData = experimentData,
                                      annotation = annotation,
                                      featureData = featureData)
            .Object@manufacturer <- manufacturer
            .Object@platform <- platform
            .Object
          })

setMethod("exprs", "FeatureSet", function(object) assayDataElement(object, "exprs"))

setMethod("manufacturer", "FeatureSet", function(object) object@manufacturer)

setReplaceMethod("manufacturer", "FeatureSet", function(object, value){
  object@manufacturer <- value
  object
})

##loading the library for now... this must change
setMethod("getPlatformDesign", "FeatureSet", function(object){
  pdn <- annotation(object)
  library(pdn,character.only=TRUE)
  return(get(pdn,pos=paste("package:",pdn,sep="")))
})

getPD <- getPlatformDesign

## probeNames - returns probeNames for PMs ... subset ignored for now
setMethod("probeNames", "FeatureSet",
          function(object, subset=NULL) {
              if (!is.null(subset))
                warning("ignoring subset arg, feature not implemented")
              probeNames(getPlatformDesign(object))
          })

## geneNames - returns geneNames for PMs

## FIXME: so geneNames is just unique(probeNames(x))?
## why the business with factor?
setMethod("geneNames", "FeatureSet",
          function(object){
              geneNames(getPlatformDesign(object))
          })


## pmindex method for FeatureSet
## WE assume feature_type is PM or MM. this might change with other platforms
setMethod("pmindex", "FeatureSet",
          function(object, subset=NULL){
            pmindex(getPlatformDesign(object), subset=subset)
          })

##mmindex method for FeatureSet
##WE assume feature_type is PM or MM. this might change with other platforms
setMethod("mmindex", "FeatureSet",
          function(object, subset=NULL){
            mmindex(getPlatformDesign(object), subset=subset)
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

## MM
## setMethod("mm", "FeatureSet", function(object, subset=NULL){
setMethod("mm", "FeatureSet",
          function(object, subset=NULL){
            exprs(object)[mmindex(object, subset=subset),]
          })

setReplaceMethod("mm", signature(object="FeatureSet", value="matrix"),
                 function(object, value){
                   tmp <- exprs(object)
                   tmp[mmindex(object),] <- value
                   assayDataElementReplace(object, "exprs", tmp)
                 })

setMethod("boxplot", signature(x="FeatureSet"),
          function(x, which, range=0, col=1:ncol(x), ...){
            if(!missing(which)) warning("Argument 'which' not yet implemented")
            idx <- sample(nrow(x), min(c(100000, nrow(x))))
            tmpdf <- as.data.frame(log2(exprs(x[idx,])))
            idx <- is.na(tmpdf[[1]])
            if (any(idx))
              tmpdf <- tmpdf[!idx,]
            boxplot(tmpdf, range=range, col=col, ...)
          })

setMethod("boxplot", signature(x="ExpressionSet"),
          function(x, which, range=0, col=1:ncol(x), ...){
            if(!missing(which)) warning("Argument 'which' not yet implemented")
            idx <- sample(nrow(x), min(c(100000, nrow(x))))
            tmpdf <- as.data.frame(log2(exprs(x[idx,])))
            idx <- is.na(tmpdf[[1]])
            if (any(idx))
              tmpdf <- tmpdf[!idx,]
            boxplot(tmpdf, range=range, col=col, ...)
          })


setMethod("image", "FeatureSet",
          function(x, transfo=log2, col=gray((0:64)/64), ...){
            if(ncol(x) > 1) par(ask=TRUE) else par(ask=FALSE)
            theDim <- geometry(getPD(x))
            for (i in 1:ncol(x)){
              tmp <- matrix(transfo(exprs(x[,i])), ncol=theDim[1], nrow=theDim[2])
              tmp <- t(tmp[nrow(tmp):1, ncol(tmp):1])
              image(tmp, col=col, main=sampleNames(x)[i], xaxt="n", yaxt="n", ...)
            }
            par(ask=FALSE)
          })

setMethod("featureNames", "FeatureSet",
          function(object) as.character(getPD(object)$feature_set_name)
          )

setMethod("hist", "FeatureSet",
          function(x, col=1:ncol(x), log=TRUE,
                   which, ylab="density", xlab="intensity",
                   type="l", ...){
            if (!missing(which)) warning("Argument 'which' not implemented yet.")
            idx <- sample(nrow(x), min(c(100000, nrow(x))))
            tmp <- exprs(x[idx,])
            idx <- is.na(tmp[,1])
            if(any(idx))
              tmp <- tmp[!idx,, drop=FALSE]
            if (log){
              tmp <- log2(tmp)
              xlab <- "log intensity"
            }
            x.density <- apply(tmp, 2, density)
            all.x <- sapply(x.density, "[[", "x")
            all.y <- sapply(x.density, "[[", "y")
            matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, ...)
            invisible(x.density)
          })

setMethod("hist", "ExpressionSet",
          function(x, col=1:ncol(x), log=TRUE,
                   which, ylab="density", xlab="intensity",
                   type="l", ...){
            if (!missing(which)) warning("Argument 'which' not implemented yet.")
            idx <- sample(nrow(x), min(c(100000, nrow(x))))
            tmp <- exprs(x[idx,])
            idx <- is.na(tmp[,1])
            if(any(idx))
              tmp <- tmp[!idx,, drop=FALSE]
            if (log){
              tmp <- log2(tmp)
              xlab <- "log intensity"
            }
            x.density <- apply(tmp, 2, density)
            all.x <- sapply(x.density, "[[", "x")
            all.y <- sapply(x.density, "[[", "y")
            matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, ...)
            invisible(x.density)
          })

setMethod("pmSequence", "FeatureSet",
          function(object) pmSequence(get(annotation(object))))

setMethod("mmSequence", "FeatureSet",
          function(object) mmSequence(get(annotation(object))))

setMethod("rma", "FeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL){
            pms <- pm(object, subset)

            tmpQcPm <- dbGetQuery(db(object),
                                  paste("SELECT man_fsetid, fid",
                                        "FROM featureSet, qcpmfeature",
                                        "WHERE qcpmfeature.fsetid=featureSet.fsetid"))
            qcpms <- exprs(object)[tmpQcPm[["fid"]],]

            pms <- rbind(pms, qcpms)
            
            pnVec <- probeNames(object, subset)

            pnVec <- c(pnVec, tmpQcPm[["man_fsetid"]])
            
            ngenes <- length(unique(pnVec))
            idx <- order(pnVec)
            pms <- pms[idx,, drop=FALSE]
            pnVec <- pnVec[idx]
            bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}

            exprs <-.Call("rma_c_complete_copy", pms, pms,
                          pnVec, ngenes,  body(bg.dens),
                          new.env(), normalize,
                          background, as.integer(2), PACKAGE="oligo")

            colnames(exprs) <- sampleNames(object)
            rownames(exprs) <- unique(pnVec)

            out <- new("ExpressionSet",
                       phenoData = phenoData(object),
                       annotation = annotation(object),
                       experimentData = experimentData(object),
                       exprs = exprs)
            return(out)
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

setMethod("db", "FeatureSet", function(object) db(getPD(object)))
setMethod("genomeBuild", "FeatureSet", function(object) genomeBuild(getPD(object)))
