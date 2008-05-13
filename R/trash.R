## from methods-FeatureSet.R

## setMethod("platform", "FeatureSet", function(object) object@platform)
## 
## setReplaceMethod("platform", "FeatureSet", function(object, value){
##   object@platform <- value
##   object
## })

## for compatibility with previous package
## setMethod("length",signature(x="FeatureSet"),
##           function(x) ncol(exprs(x))) 


## setMethod("platformDesignName", "FeatureSet", function(object){
##   platform(object)})

## setMethod("getPlatformDesign", "FeatureSet", function(object){
##   pdn <- platformDesignName(object)
##   library(pdn,character.only=TRUE)
##   return(get(pdn,pos=paste("package:",pdn,sep="")))
## })

## setMethod("hist", signature(x="ExpressionSet"),
##           function(x, which=c("both", "pm", "mm"), ...)
##           plotDensity(x, which=c("both", "pm", "mm"), ...))

## setReplaceMethod("pm", signature(object="FeatureSet", value="BufferedMatrix"),
##                  function(object, value){
##                    tmp <- exprs(object)
##                    for (i in 1:ncol(tmp))
##                      tmp[pmindex(object), i] <- value[,i]
##                    assayDataElementReplace(object, "exprs", tmp)
##                  })

## setReplaceMethod("mm", signature(object="FeatureSet", value="BufferedMatrix"),
##                  function(object, value){
##                    tmp <- exprs(object)
##                    for (i in 1:ncol(tmp))
##                      tmp[mmindex(object),i] <- value[,i]
##                    assayDataElementReplace(object, "exprs", tmp)
##                  })

## setMethod("featureIndex", "FeatureSet",
##           function(object, which=c("both", "pm", "mm"), subset=NULL){
##             which <- match.arg(which,c("both", "pm", "mm"))
##             if (which=="both"){
##               pmIndex <- pmindex(getPlatformDesign(object), subset=subset)
##               mmIndex <- mmindex(getPlatformDesign(object), subset=subset)
##               indexes <- sort(c(pmIndex,mmIndex))
##             } else if (which=="pm"){
##               indexes <- sort(pmindex(getPlatformDesign(object), subset=subset))
##             } else if (which=="mm"){
##               indexes <- sort(mmindex(getPlatformDesign(object), subset=subset))
##             }
##              return(indexes)
##            })


## setMethod("boxplot", signature(x="FeatureSet"),
##           function(x, which=c("both", "pm", "mm"), range=0, ...){
##             which <- match.arg(which, c("both", "pm", "mm"))
##             tmp <- description(x)
##             tmp <- unlist(featureIndex(x, which))
##             tmp <- tmp[seq(1, length(tmp), len=5000)]
##             cols <- 1:length(sampleNames(x))
##             boxplot(data.frame(log2(exprs(x)[tmp, ])),  col=cols, range=range, ...)
##           })

## setMethod("boxplot", signature(x="ExpressionSet"),
##           function(x, which=c("both", "pm", "mm"), range=0, ...){
##             e <- data.frame(exprs(x))
##             tmp <- seq(1, nrow(e), len=5000)
##             cols <- 1:length(sampleNames(x))
##             boxplot(e[tmp, ],  col=cols, range=range, ...)
##           })

## type <- function(object) getPD(object)@type

## setMethod("plotDensity", "FeatureSet", function(object, col=1:6, log=TRUE,
##                                                 which=c("both","pm","mm"),
##                                                 ylab="density",
##                                                 xlab="log intensity",
##                                                 type="l",
##                                                 ...){
##   which <- match.arg(which,c("both","pm","mm"))
##   Index <- unlist(featureIndex(object,which))
##   object <- exprs(object)[Index,, drop=FALSE]
##   if(log){
##     object <- log2(object)
##     if(is.null(xlab)) xlab <- "log intensity"
##   }
##   else xlab <- "intensity"
##   n <- ncol(object)
##   x.density <- list()
##   for (i in 1:n) x.density[[i]] <- density(object[,i])
##   all.x <- do.call("cbind", lapply(x.density, function(x) x$x))
##   all.y <- do.call("cbind", lapply(x.density, function(x) x$y))
##   matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, ...)
##   invisible(list(all.x=all.x, all.y=all.y))
## })

## setMethod("plotDensity", "ExpressionSet", function(object, col=1:6, log=TRUE,
##                                                 which=c("both","pm","mm"),
##                                                 ylab="density",
##                                                 xlab="log intensity",
##                                                 type="l",
##                                                 ...){
##   object <- exprs(object)
##   n <- ncol(object)
##   x.density <- list()
##   for (i in 1:n) x.density[[i]] <- density(object[,i])
##   all.x <- do.call("cbind", lapply(x.density, function(x) x$x))
##   all.y <- do.call("cbind", lapply(x.density, function(x) x$y))
##   matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, ...)
##   invisible(list(all.x=all.x, all.y=all.y))
## })

## setMethod("rma", "FeatureSet",
##           function(object, background=TRUE, normalize=TRUE){
##             stop("RMA temporarily broken...")
##             pms <- pm(object)
##             pnVec <- probeNames(object)
##             idx <- order(pnVec)
##             pms <- subBufferedMatrix(pms, idx)
##             pnVec <- pnVec[idx]
##             rm(idx); gc()
##             ColMode(pms)
##             set.buffer.dim(pms, 50000, 1)
##             if (background) bg.correct.BufferedMatrix(pms, copy=FALSE)
##             if (normalize) normalize.BufferedMatrix.quantiles(pms, copy=FALSE)
##             RowMode(pms)
##             exprs <- median.polish.summarize(pms, length(unique(pnVec)), pnVec)
##             rownames(exprs) <- unique(pnVec)
##             colnames(exprs) <- sampleNames(object)
##             rm(pms, pnVec); gc()
##             out <- new("ExpressionSet",
##                        exprs=exprs,
##                        phenoData=phenoData(object),
##                        experimentData=experimentData(object),
##                        annotation=annotation(object))
##             sampleNames(out) <- sampleNames(object)
##             return(out)
##           })

## setMethod("image", signature(x="FeatureSet"),
##           function(x, transfo=log, col=gray((0:64)/64),
##                    xlab="", ylab="", ...){
##   if(ncol(x) > 1) par(ask=TRUE) else par(ask=FALSE)
##   if (tolower(manufacturer(x)) == "affymetrix"){
##     if (is(getPD(x), "platformDesign")){
##       nr <- nrow(getPD(x))
##       nc <- ncol(getPD(x))
##     }else{
##       nr <- nc <- sqrt(nrow(exprs(x)))
##     }
##     for(i in 1:ncol(x)){
##       m <- as.numeric(exprs(x)[,i])
##       if (is(getPD(x), "platformDesign")) m <- m[getPD(x)$order_index]
##       if (is.function(transfo)) m <- transfo(m)
##       m <- matrix(m, ncol=nc, nrow=nr)
##       m <- t(m)[nc:1,]
##       image(m, col=col, main=sampleNames(x)[i],
##             xlab=xlab, ylab=ylab, xaxt='n', yaxt='n', ...)
##     }
##   }else if (tolower(manufacturer(x)) == "nimblegen"){
##     ## getting rid of areas with now signal
##     xs=getPD(x)$X
##     ys=getPD(x)$Y
##     tmpy=table(ys);ny=max(tmpy)
##     levels=seq(min(ys),max(ys),len=ny)
##     m=matrix(NA,ncol=length(unique(xs)),nrow=length(levels)-1)
##     xIndexes=split(seq(along=xs),xs)
##     yIndexes=sapply(xIndexes,function(i) as.numeric(cut(ys[i],breaks=levels,include.lowest=TRUE)))
##     if (is.function(transfo)) m <- transfo(m)
##     ## end... now make plots
##     for (i in 1:length(sampleNames(x))) {
##       for(j in seq(along=xIndexes)){
##         m[yIndexes[[j]],j]=exprs(x)[xIndexes[[j]],i]
##       }
##       image(m,xlab=xlab,ylab=ylab,col=col, main=sampleNames(x)[i], ...)
##     }
##   }else{
##     stop("I dont know how to handle this array.")
##   }
##   par(ask=FALSE)
## })


## setMethod("hist", signature(x="FeatureSet"),
##           function(x, which=c("both", "pm", "mm"), ...)
##           plotDensity(x, which=c("both", "pm", "mm"), ...))


## FROM AllGenerics.R

## setGeneric("db", function(object) standardGeneric("db"))

## setGeneric("geneNames", function(object) standardGeneric("geneNames"))

## setGeneric("indexFeatureSetName",
##            function(object, featurenames) {
##                standardGeneric("indexFeatureSetName")
##            })

## setGeneric("featureSetNames",
##            function(object, ids) standardGeneric("featureSetNames"))

## setGeneric("featureIDs",
##            function(object, ids) standardGeneric("featureIDs"))

## setGeneric("npixels",
##            function(object) standardGeneric("npixels"))

## setGeneric("allele", function(object) standardGeneric("allele"))

## setGeneric("calls<-", function(object, value) standardGeneric("calls<-"))
## 
## setGeneric("callsConfidence<-",
##            function(object, value) standardGeneric("callsConfidence<-"))
## 
## setGeneric("calls", function(object) standardGeneric("calls"))
## 
## setGeneric("callsConfidence", function(object) standardGeneric("callsConfidence"))

## setGeneric("copyNumber<-",
##            function(object, value) standardGeneric("copyNumber<-"))

## setGeneric("cnConfidence<-",
##            function(object, value) standardGeneric("cnConfidence<-"))

## setGeneric("copyNumber", function(object) standardGeneric("copyNumber"))
## setGeneric("cnConfidence", function(object) standardGeneric("cnConfidence"))
## setGeneric("alleleAB", function(object) standardGeneric("alleleAB"))
## setGeneric("pmAlleleAB", function(object) standardGeneric("pmAlleleAB"))
## setGeneric("chromosome", function(object) standardGeneric("chromosome"))
## setGeneric("position", function(object) standardGeneric("position"))
## setGeneric("plotDensity", function(object, ...) standardGeneric("plotDensity"))
## setGeneric("snpMedianSilhouette",
##           function(object) standardGeneric("snpMedianSilhouette"))
## setGeneric("platformDesignName",
##           function(object) standardGeneric("platformDesignName"))
