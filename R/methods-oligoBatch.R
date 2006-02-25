
setMethod("platform","oligoBatch", function(object) object@platform)

###for compatibility with previous package
setMethod("length",signature(x="oligoBatch"),
          function(x) ncol(exprs(x))) 


setMethod("platformDesignName","oligoBatch", function(object){
  cleanPlatformName(object@platform)})

##loading the library for now... this must change
setMethod("getPlatformDesign","oligoBatch", function(object){
  pdn <- platformDesignName(object)
  library(pdn,character.only=TRUE)
  return(get(pdn,pos=paste("package:",pdn,sep="")))
})

getPD <- getPlatformDesign

## probeNames - returns probeNames for PMs ... genenames ignored for now
setMethod("probeNames", "oligoBatch",
          function(object, genenames=NULL){
            pmIndex <- pmindex(getPlatformDesign(object))
            pns <- get("feature_set_name",envir=featureInfo(getPlatformDesign(object)))
            return(as.character(pns[pmIndex]))
          })

###geneNames - returns geneNames for PMs
setMethod("geneNames", "oligoBatch",
          function(object){
            pmIndex <- pmindex(getPlatformDesign(object))
            return(levels(factor(get("feature_set_name",envir=featureInfo(getPlatformDesign(object)))[pmIndex])))
          })


##pmindex method for oligoBatch
##WE assume feature_type is PM or MM. this might change with other platforms
setMethod("pmindex", "oligoBatch",
          function(object){
            pmindex(getPlatformDesign(object))
          })

##mmindex method for oligoBatch
##WE assume feature_type is PM or MM. this might change with other platforms
setMethod("mmindex", "oligoBatch",
          function(object){
            mmindex(getPlatformDesign(object))
          })

## BC: indexFeatureSetName - to simplify the procedure of bringing pms/mms with a given name
##     it'll return the indexes for a given feature_set_name
##     This is used by pm/mm when getting intensities for given feature_set_names
setMethod("indexFeatureSetName","oligoBatch",
          function(object, featurenames){
            tmp <- NULL
            for (i in featurenames)
              tmp <- c(tmp,which(getPlatformDesign(object)$feature_set_name == i))
            return(sort(tmp))
          })


setMethod("ncol",signature(x="oligoBatch"),
                    function(x) getPlatformDesign(x)@ncol)

setMethod("nrow",signature(x="oligoBatch"),
          function(x) getPlatformDesign(x)@nrow)

## Histogram
setMethod("hist",signature(x="oligoBatch"),
          function(x, which=c("both","pm","mm"), ...)
          plotDensity.oligoBatch(x, which=c("both","pm","mm"), ...))

## Show
setMethod("show","oligoBatch", function(object){
  dm <- dim(exprs(object))
  nprobes <- dm[1]
  nsamples <- dm[2]
  cat("oligoBatch with \n\t", nprobes, " probes\n\t", sep="")
  cat(nsamples, "samples\n\t")
  show(phenoData(object))
})

## PM
setMethod("pm","oligoBatch", ##genenames is ignored for now.. we will get to it
          function(object, genenames=NULL){
            index <- pmindex(object)
            if (!is.null(genenames)){
              index <- intersect(index,indexFeatureSetName(object,genenames))
              fsn <- getPD(object)$feature_set_name[index]
              fsid <- getPD(object)$feature_ID[index]
              rn <- paste(fsn,fsid,sep=".")
              oo <- exprs(object)[index,,drop=FALSE]
              rownames(oo) <- rn
              colnames(oo) <- sampleNames(object)
              return(oo)
            }
            return(exprs(object)[index,,drop=FALSE])
          })

setReplaceMethod("pm", "oligoBatch",
                 function(object, value){
##                   exprs(object)[pmindex(object),] <- value
                   assayData(object)[["exprs"]][pmindex(object),] <- value
                   object
                 })

## MM
setMethod("mm","oligoBatch", function(object, genenames=NULL){
            index <- mmindex(object)
            if (!is.null(genenames)){
              index <- intersect(index,indexFeatureSetName(object,genenames))
              fsn <- getPD(object)$feature_set_name[index]
              fsid <- getPD(object)$feature_ID[index]
              rn <- paste(fsn,fsid,sep=".")
              oo <- exprs(object)[index,,drop=FALSE]
              rownames(oo) <- rn
              colnames(oo) <- sampleNames(object)
              return(oo)
            }
            return(exprs(object)[index,,drop=FALSE])
          })

setReplaceMethod("mm", "oligoBatch",
                 function(object, value){
                   assayData(object)[["exprs"]][mmindex(object),] <- value
                   object
                 })


setMethod("featureIndex","oligoBatch",
          function(object, which=c("both","pm","mm"), genenames=NULL){
            which <- match.arg(which,c("both","pm","mm"))
            if (which=="both"){
              pmIndex <- pmindex(getPlatformDesign(object))
              mmIndex <- mmindex(getPlatformDesign(object))
              indexes <- sort(c(pmIndex,mmIndex))
            } else if (which=="pm"){
              indexes <- sort(pmindex(getPlatformDesign(object)))
            } else if (which=="mm"){
              indexes <- sort(mmindex(getPlatformDesign(object)))
            }
             return(indexes)
           })

setMethod("boxplot",signature(x="oligoBatch"),
          function(x, which=c("both", "pm", "mm"), range=0, ...){
            which <- match.arg(which, c("both", "pm", "mm"))
            tmp <- description(x)
            if (is(tmp, "MIAME")) main <- tmp@title

            tmp <- unlist(featureIndex(x, which))
            tmp <- tmp[seq(1, length(tmp), len=5000)]

            boxplot(data.frame(log2(exprs(x)[tmp, ])), main=main, range=range, ...)
          })

setMethod("image",signature(x="oligoBatch"),
          function(x, transfo=log, col=gray(c(0:64)/64), xlab="", ylab="", ...){
            scn <- prod(par("mfrow"))
            ask <- dev.interactive()
            which.plot <- 0

            x.pos <- (1:nrow(x)) - 1
            y.pos <- (1:ncol(x)) - 1

            correctOrder <- order(getPD(x)$X, getPD(x)$Y)
            
            for(i in 1:length(sampleNames(x))){
              which.plot <- which.plot+1;

              if(trunc((which.plot-1)/scn)==(which.plot-1)/scn && which.plot>1 && ask)
                par(ask=TRUE)
              
              m <- matrix(exprs(x)[correctOrder, i], ncol=ncol(x))
              
              if (is.function(transfo))
                m <- transfo(m)

              m <- t(as.matrix(rev(as.data.frame(m))))
##              image(x.pos, y.pos, m, col=col,
              image(m, col=col,
                    main=sampleNames(x)[i],
                    xlab=xlab, ylab=ylab,
                    xaxt='n', yaxt='n', ...)
              par(ask=FALSE)
            }
          })

setMethod("sd", "oligoBatch",
          function(x, na.rm=TRUE) return(assayData(x)$sd))

setMethod("npixels", signature(object="oligoBatch"),
          function(object) return(assayData(object)$npixels))

