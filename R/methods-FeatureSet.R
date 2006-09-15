setMethod("initialize", "FeatureSet",
          function(.Object,
                   manufacturer=new("character"),
                   platform=new("character"),
                   exprs=new("matrix"),
                   sd=new("matrix"),
                   npixels=new("matrix"),
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = new("character")){
            .Object <- callNextMethod(.Object,
                                      assayData = assayDataNew(
                                        exprs=exprs),
                                      
##                                        exprs=exprs,
##                                        sd=sd,
##                                        npixels=npixels),
                                      
                                      phenoData = phenoData,
                                      experimentData = experimentData,
                                      annotation = annotation)
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

###for compatibility with previous package
setMethod("length",signature(x="FeatureSet"),
          function(x) ncol(exprs(x))) 

setMethod("platformDesignName", "FeatureSet", function(object){
#  cleanPlatformName(object@platform)})
  platform(object)})

##loading the library for now... this must change
setMethod("getPlatformDesign", "FeatureSet", function(object){
  pdn <- platformDesignName(object)
  library(pdn,character.only=TRUE)
  return(get(pdn,pos=paste("package:",pdn,sep="")))
})

getPD <- getPlatformDesign

## probeNames - returns probeNames for PMs ... genenames ignored for now
##setMethod("probeNames", c("FeatureSet", "characterOrNULL"),
probeNames <- function(object, subset=NULL){
         ## function(object){
            pmIndex <- pmindex(getPlatformDesign(object))
            pns <- get("feature_set_name",envir=featureInfo(getPlatformDesign(object)))
            return(as.character(pns[pmIndex]))
#          })
          }

###geneNames - returns geneNames for PMs
setMethod("geneNames", "FeatureSet",
          function(object){
            pmIndex <- pmindex(getPlatformDesign(object))
            return(levels(factor(get("feature_set_name",envir=featureInfo(getPlatformDesign(object)))[pmIndex])))
          })


##pmindex method for FeatureSet
##WE assume feature_type is PM or MM. this might change with other platforms
setMethod("pmindex", "FeatureSet",
          function(object){
            pmindex(getPlatformDesign(object))
          })

##mmindex method for FeatureSet
##WE assume feature_type is PM or MM. this might change with other platforms
setMethod("mmindex", "FeatureSet",
          function(object){
            mmindex(getPlatformDesign(object))
          })

## BC: indexFeatureSetName - to simplify the procedure of bringing pms/mms with a given name
##     it'll return the indexes for a given feature_set_name
##     This is used by pm/mm when getting intensities for given feature_set_names
setMethod("indexFeatureSetName", "FeatureSet",
          function(object, featurenames){
            tmp <- NULL
            for (i in featurenames)
              tmp <- c(tmp,which(getPlatformDesign(object)$feature_set_name == i))
            return(sort(tmp))
          })

setMethod("ncol",signature(x="FeatureSet"),
                    function(x) getPlatformDesign(x)@ncol)

setMethod("nrow",signature(x="FeatureSet"),
          function(x) getPlatformDesign(x)@nrow)

## Histogram
setMethod("hist", signature(x="FeatureSet"),
          function(x, which=c("both", "pm", "mm"), ...)
          plotDensity.FeatureSet(x, which=c("both", "pm", "mm"), ...))


## PM
setMethod("pm", "FeatureSet", ##genenames is ignored for now.. we will get to it
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

setReplaceMethod("pm", signature(object="FeatureSet", value="matrix"),
                 function(object, value){
                   tmp <- exprs(object)
                   tmp[pmindex(object),] <- value
                   assayDataElementReplace(object, "exprs", tmp)
                 })

## MM
setMethod("mm", "FeatureSet", function(object, genenames=NULL){
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

setReplaceMethod("mm", signature(object="FeatureSet", value="matrix"),
                 function(object, value){
                   tmp <- exprs(object)
                   tmp[mmindex(object),] <- value
                   assayDataElementReplace(object, "exprs", tmp)
                 })

setMethod("featureIndex", "FeatureSet",
          function(object, which=c("both", "pm", "mm"), genenames=NULL){
            which <- match.arg(which,c("both", "pm", "mm"))
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

setMethod("boxplot", signature(x="FeatureSet"),
          function(x, which=c("both", "pm", "mm"), range=0, ...){
            which <- match.arg(which, c("both", "pm", "mm"))
            tmp <- description(x)
            if (is(tmp, "MIAME")) main <- tmp@title

            tmp <- unlist(featureIndex(x, which))
            tmp <- tmp[seq(1, length(tmp), len=5000)]

            boxplot(data.frame(log2(exprs(x)[tmp, ])), main=main, range=range, ...)
          })

setMethod("image", signature(x="FeatureSet"),
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
              image(m, col=col,
                    main=sampleNames(x)[i],
                    xlab=xlab, ylab=ylab,
                    xaxt='n', yaxt='n', ...)
              par(ask=FALSE)
            }
          })

setMethod("sd", "FeatureSet",
          function(x, na.rm=TRUE) return(assayData(x)$sd))

setMethod("npixels", signature(object="FeatureSet"),
          function(object) return(assayData(object)$npixels))

type <- function(object) getPD(object)@type

setMethod("allele", signature(object="FeatureSet"),
          function(object){
            if(type(object)!="SNP"){
              stop("This array does not have allele information")
            }else{
              getPD(object)$allele
            }})
