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

setMethod("platform", "FeatureSet", function(object) object@platform)

setReplaceMethod("platform", "FeatureSet", function(object, value){
  object@platform <- value
  object
})

setMethod("manufacturer", "FeatureSet", function(object) object@manufacturer)

setReplaceMethod("manufacturer", "FeatureSet", function(object, value){
  object@manufacturer <- value
  object
})

## for compatibility with previous package
setMethod("length",signature(x="FeatureSet"),
          function(x) ncol(exprs(x))) 


setMethod("platformDesignName", "FeatureSet", function(object){
  platform(object)})

##loading the library for now... this must change

setMethod("getPlatformDesign", "FeatureSet", function(object){
  pdn <- platformDesignName(object)
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

setMethod("hist", signature(x="FeatureSet"),
          function(x, which=c("both", "pm", "mm"), ...)
          plotDensity(x, which=c("both", "pm", "mm"), ...))

setMethod("hist", signature(x="ExpressionSet"),
          function(x, which=c("both", "pm", "mm"), ...)
          plotDensity(x, which=c("both", "pm", "mm"), ...))

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

## setReplaceMethod("pm", signature(object="FeatureSet", value="BufferedMatrix"),
##                  function(object, value){
##                    tmp <- exprs(object)
##                    for (i in 1:ncol(tmp))
##                      tmp[pmindex(object), i] <- value[,i]
##                    assayDataElementReplace(object, "exprs", tmp)
##                  })



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

## setReplaceMethod("mm", signature(object="FeatureSet", value="BufferedMatrix"),
##                  function(object, value){
##                    tmp <- exprs(object)
##                    for (i in 1:ncol(tmp))
##                      tmp[mmindex(object),i] <- value[,i]
##                    assayDataElementReplace(object, "exprs", tmp)
##                  })

setMethod("featureIndex", "FeatureSet",
          function(object, which=c("both", "pm", "mm"), subset=NULL){
            which <- match.arg(which,c("both", "pm", "mm"))
            if (which=="both"){
              pmIndex <- pmindex(getPlatformDesign(object), subset=subset)
              mmIndex <- mmindex(getPlatformDesign(object), subset=subset)
              indexes <- sort(c(pmIndex,mmIndex))
            } else if (which=="pm"){
              indexes <- sort(pmindex(getPlatformDesign(object), subset=subset))
            } else if (which=="mm"){
              indexes <- sort(mmindex(getPlatformDesign(object), subset=subset))
            }
             return(indexes)
           })


setMethod("boxplot", signature(x="FeatureSet"),
          function(x, which=c("both", "pm", "mm"), range=0, ...){
            which <- match.arg(which, c("both", "pm", "mm"))
            tmp <- description(x)
            tmp <- unlist(featureIndex(x, which))
            tmp <- tmp[seq(1, length(tmp), len=5000)]
            cols <- 1:length(sampleNames(x))
            boxplot(data.frame(log2(exprs(x)[tmp, ])),  col=cols, range=range, ...)
          })

setMethod("boxplot", signature(x="ExpressionSet"),
          function(x, which=c("both", "pm", "mm"), range=0, ...){
            e <- data.frame(exprs(x))
            tmp <- seq(1, nrow(e), len=5000)
            cols <- 1:length(sampleNames(x))
            boxplot(e[tmp, ],  col=cols, range=range, ...)
          })


setMethod("image", signature(x="FeatureSet"),
          function(x, transfo=log, col=gray((0:64)/64),
                   xlab="", ylab="", ...){
  if(ncol(x) > 1) par(ask=TRUE) else par(ask=FALSE)
  if (tolower(manufacturer(x)) == "affymetrix"){
    if (is(getPD(x), "platformDesign")){
      nr <- nrow(getPD(x))
      nc <- ncol(getPD(x))
    }else{
      nr <- nc <- sqrt(nrow(exprs(x)))
    }
    for(i in 1:ncol(x)){
      m <- as.numeric(exprs(x)[,i])
      if (is(getPD(x), "platformDesign")) m <- m[getPD(x)$order_index]
      if (is.function(transfo)) m <- transfo(m)
      m <- matrix(m, ncol=nc, nrow=nr)
      m <- t(m)[nc:1,]
      image(m, col=col, main=sampleNames(x)[i],
            xlab=xlab, ylab=ylab, xaxt='n', yaxt='n', ...)
    }
  }else if (tolower(manufacturer(x)) == "nimblegen"){
    ## getting rid of areas with now signal
    xs=getPD(x)$X
    ys=getPD(x)$Y
    tmpy=table(ys);ny=max(tmpy)
    levels=seq(min(ys),max(ys),len=ny)
    m=matrix(NA,ncol=length(unique(xs)),nrow=length(levels)-1)
    xIndexes=split(seq(along=xs),xs)
    yIndexes=sapply(xIndexes,function(i) as.numeric(cut(ys[i],breaks=levels,include.lowest=TRUE)))
    if (is.function(transfo)) m <- transfo(m)
    ## end... now make plots
    for (i in 1:length(sampleNames(x))) {
      for(j in seq(along=xIndexes)){
        m[yIndexes[[j]],j]=exprs(x)[xIndexes[[j]],i]
      }
      image(m,xlab=xlab,ylab=ylab,col=col, main=sampleNames(x)[i], ...)
    }
  }else{
    stop("I dont know how to handle this array.")
  }
  par(ask=FALSE)
})

## type <- function(object) getPD(object)@type

setMethod("featureNames", "FeatureSet",
          function(object) as.character(getPD(object)$feature_set_name)
          )


setMethod("plotDensity", "FeatureSet", function(object, col=1:6, log=TRUE,
                                                which=c("both","pm","mm"),
                                                ylab="density",
                                                xlab="log intensity",
                                                type="l",
                                                ...){
  which <- match.arg(which,c("both","pm","mm"))
  Index <- unlist(featureIndex(object,which))
  object <- exprs(object)[Index,, drop=FALSE]
  if(log){
    object <- log2(object)
    if(is.null(xlab)) xlab <- "log intensity"
  }
  else xlab <- "intensity"
  n <- ncol(object)
  x.density <- list()
  for (i in 1:n) x.density[[i]] <- density(object[,i])
  all.x <- do.call("cbind", lapply(x.density, function(x) x$x))
  all.y <- do.call("cbind", lapply(x.density, function(x) x$y))
  matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, ...)
  invisible(list(all.x=all.x, all.y=all.y))
})

setMethod("plotDensity", "ExpressionSet", function(object, col=1:6, log=TRUE,
                                                which=c("both","pm","mm"),
                                                ylab="density",
                                                xlab="log intensity",
                                                type="l",
                                                ...){
  object <- exprs(object)
  n <- ncol(object)
  x.density <- list()
  for (i in 1:n) x.density[[i]] <- density(object[,i])
  all.x <- do.call("cbind", lapply(x.density, function(x) x$x))
  all.y <- do.call("cbind", lapply(x.density, function(x) x$y))
  matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, ...)
  invisible(list(all.x=all.x, all.y=all.y))
})


setMethod("pmSequence", "FeatureSet",
          function(object) pmSequence(get(annotation(object))))

setMethod("mmSequence", "FeatureSet",
          function(object) mmSequence(get(annotation(object))))

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

setMethod("rma", "FeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL){
            if (class(object) == "SnpFeatureSet") stop("Use justSNPRMA instead.")
            pms <- pm(object, subset)
            pnVec <- probeNames(object, subset)
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

            out <- new("ExpressionSet",
                       phenoData = phenoData(object),
                       annotation = annotation(object),
                       experimentData = experimentData(object),
                       exprs = exprs)
            return(out)
          })

setMethod("MAplot", "FeatureSet",
          function(object, arrays=1:ncol(object), lowessPlot=FALSE){
            par(ask=TRUE)
            ref <- rowMedians(log2(exprs(object)))
            for (i in arrays){
              tmp <- log2(exprs(object[,i]))
              plot((tmp+ref)/2, tmp-ref, pch=".",
                   xlab="A", ylab="M",
                   main=sampleNames(object)[i])
              if (lowessPlot)
                lines(lowess((tmp+ref)/2, tmp-ref), col="red")
            }
            par(ask=FALSE)
          })
