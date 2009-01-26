setMethod("initialize", signature(.Object="FeatureSet"), 
          function(.Object,
                   assayData = assayDataNew(exprs=exprs, ...),
                   manufacturer=as.character(NA),
                   platform=as.character(NA),
                   exprs=matrix(numeric(0), nrow=nrow, ncol=ncol),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
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


setMethod("platform", signature(object="FeatureSet"),function(object) object@platform)
setReplaceMethod("platform",signature(object="FeatureSet"),
		  function(object,value) {
			  object@platform <- value
			  object
		  })
  
setMethod("manufacturer", signature(object="FeatureSet"), function(object) object@manufacturer)
setReplaceMethod("manufacturer", signature(object="FeatureSet"), 
		function(object, value){
  			object@manufacturer <- value
  			object
})

setMethod("length", signature(x="FeatureSet"), function(x) nrow(pData(x)))

## Currently dim, ncol and nrow all pull from eSet Class definition, and thus is the size of
##  exprs matrix and not array, should change to array dims -- MS
# setMethod("dim",
#          signature=signature(object="FeatureSet"),
#          function(object) c(object@nrow, object@ncol))
#


## changed annotation to platform - MS
## thinking is that platform is akin to cdf and annotation akin to db
## currently they are both the same package, but may in the future want to have annotationDBI like packages
## for these arrays as well and use annotation slot for this

##loading the library for now... this must change, why?
setMethod("getPlatformDesign", 
        signature(object= "FeatureSet"),
        function(object){
            pdn <- platform(object)
            library(pdn,character.only=TRUE)
            return(get(pdn,pos=paste("package:",pdn,sep="")))
})
getPD <- getPlatformDesign


setMethod("exprs",
		signature(object="FeatureSet"),
		function(object) assayDataElement(object,"exprs"))
setReplaceMethod("exprs",
		signature(object="FeatureSet"),
		function(object, value) assayDataElementReplace(object, "exprs", value))

### future 
setMethod("se.exprs",
		signature(object="FeatureSet"),
		function(object) {
			obj <- assayDataElement(object,"se.exprs")
			if (is.null(obj))
				new("matrix")
			else
				obj
		})
### future
setReplaceMethod("se.exprs",
		signature(object="FeatureSet"),
		function(object, value) assayDataElementReplace(object, "se.exprs", value))


#####################
############# STRANGE, what should they be?
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
#####################
#####################

## pmindex method for FeatureSet
## WE assume feature_type is PM or MM. this might change with other platforms
setMethod("pmindex", "FeatureSet",
          function(object, subset=NULL){
            pmindex(getPlatformDesign(object), subset=subset)
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
##mmindex method for FeatureSet
##WE assume feature_type is PM or MM. this might change with other platforms
setMethod("mmindex", "FeatureSet",
		function(object, subset=NULL){
			mmindex(getPD(object), subset=subset)
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
		 
## does not work for GENE ST arrays
setMethod("featureNames", "FeatureSet",
          function(object) as.character(getPD(object)$feature_set_name)
          )

## does not work for GENE ST arrays
setMethod("pmSequence", "FeatureSet",
          function(object) pmSequence(getPD(object)))

## does not work for GENE ST arrays
setMethod("mmSequence", "FeatureSet",
          function(object) mmSequence(getPD(object)))

## RMA normalization, likely not OK for ST arrays
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

  
  ## QAish functions
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
            conn <- db(x)
            geom <- geometry(getPD(x))
            tbls <- dbGetQuery(conn, "SELECT tbl FROM table_info WHERE tbl LIKE '%feature' AND row_count > 0")[[1]]
            theInfo <- lapply(tbls, function(tb) dbGetQuery(conn, paste("SELECT x, y, fid FROM", tb)))
            theInfo <- do.call("rbind", theInfo)
            theInfo <- theInfo[order(theInfo[["fid"]]), ]
            idx <- geom[1]*(theInfo[["x"]]-1)+theInfo[["y"]]
            for (i in which){
              theInfo[["signal"]] <- transfo(as.numeric(exprs(x[theInfo[["fid"]], i])))
              int <- matrix(NA, nr=geom[1], nc=geom[2])
              int[idx] <- theInfo[["signal"]]
              int <- t(int[nrow(int):1, ncol(int):1])
              image(int, col=col, main=sampleNames(x)[i], xaxt="n", yaxt="n", ...)
            }
            par(ask=FALSE)
          })


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

setMethod("db", "FeatureSet", function(object)    db(getPD(object)))
setMethod("genomeBuild", "FeatureSet", function(object) genomeBuild(getPD(object)))
