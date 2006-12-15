setValidity("SnpCallSetPlus", function(object) {
    assayDataValidMembers(assayData(object), c("calls", "callsConfidence", "logRatioAntisense", "logRatioSense"))
  })

setMethod("initialize", "SnpCallSetPlus",
          function(.Object,
                   assayData = assayDataNew(calls=calls,
                     callsConfidence=callsConfidence,
                     logRatioAntisense=logRatioAntisense,
                     logRatioSense=logRatioSense, ...),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   calls = new("matrix"),
                   callsConfidence = new("matrix"),
                   logRatioAntisense = new("matrix"),
                   logRatioSense = new("matrix"),
                   ...){
            callNextMethod(.Object,
                           assayData = assayData,
                           featureData = featureData,
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation)
          })


setMethod("logRatioAntisense", "SnpCallSetPlus",
          function(object) assayDataElement(object, "logRatioAntisense"))

setMethod("logRatioSense", "SnpCallSetPlus",
          function(object) assayDataElement(object, "logRatioSense"))


snpSilhouette <- function(snpIdx, obj){
  tmp <- cbind(assayDataElement(obj, "logRatioAntisense")[snpIdx,],
                        assayDataElement(obj, "logRatioSense")[snpIdx,])
  idx <- which(is.na(tmp[1,]))
  if (length(idx) == 1) tmp <- tmp[,-idx]
  if (length(idx) == 2) stop(paste("SNP", snpIdx, "has log-ratios missing on both strands"))
  silhouette(as.integer(factor(calls(obj)[snpIdx,])), dist(tmp))
}


setMethod("snpMedianSilhouette", "SnpCallSetPlus",
          function(object){
            require(cluster)
            sapply(1:nrow(object),
                   function(idx, object){
                     if (idx %% 5000 == 0) cat(".")
                     sils <- snpSilhouette(idx, object)
                     ifelse(is(sils, "silhouette"), median(sils[,3], na.rm=TRUE), NA)
                   },
                   object)
          })

plotSnpConfidenceSilhouette <- function(snp, object){
  mydens <- function(obj, n){
    if (length(obj) > 1){
      tmpdens <- density(obj)
      tmpdens$y <- tmpdens$y/max(tmpdens$y)
      return(tmpdens)
    }else{
      return(list(x=obj, y=1))
    }
  }
  sils <- snpSilhouette(snp, object)
  if (is(sils, "silhouette")){
    sils <- sils[,3]
    par(mfrow=c(2,2))
    tmp <- cbind(logRatioAntisense(object)[snp,], logRatioSense(object)[snp,])
    idx <- which(is.na(tmp[1,]))
    plot(callsConfidence(object)[snp,],
         sils, xlab="Log-likelihood Ratio",
         ylab="Silhouette", ylim=c(-1, 1), xlim=c(0, 200),
         main=paste("Silhouette vs. LLR - SNP", snp),
         col=calls(object)[snp,])
    if (length(idx) == 1){
      mt <- "Log-Ratio: Antisense"
      if (idx == 1) mt <- "Log-Ratio: Sense"
      tmp <- split(tmp[,-idx], factor(calls(object)[snp,]))
      dens <- lapply(tmp, mydens, sum(sapply(tmp, length)))
      plot(1:10,  xlab=mt, main=paste("Density for", mt), ylab="Density", xlim=c(-5,5), ylim=c(0,1), type="n")
      for(i in 1:length(tmp))
        lines(dens[[i]]$x[order(dens[[i]]$x)], dens[[i]]$y[order(dens[[i]]$x)], col=as.numeric(names(tmp)[i]))
    } else {
      plot(logRatioAntisense(object)[snp,],
           logRatioSense(object)[snp,],
           xlab="Log-ratio: Antisense",
           ylab="Log-ratio: Sense",
           main=paste("Clusters for SNP", snp),
           xlim=c(-5,5),
           ylim=c(-5,5),
           col=calls(object)[snp,])
    }
    tmp <- split(callsConfidence(object)[snp,], factor(calls(object)[snp,]))
    dens <- lapply(tmp, mydens, sum(sapply(tmp, length)))
    n <- length(tmp)
    plot(1:10, xlab="Log-likelihood Ratio", ylab="Density", main="Density of LLR",
         type="n", xlim=c(0,200), ylim=c(0,1))
    for (i in 1:n)
      lines(dens[[i]]$x[order(dens[[i]]$x)], dens[[i]]$y[order(dens[[i]]$x)], col=as.numeric(names(tmp)[i]))
    tmp <- split(sils, factor(calls(object)[snp,]))
    dens <- lapply(tmp, mydens, sum(sapply(tmp, length)))
    n <- length(tmp)
    plot(1:10, xlab="Silhouette", ylab="Density", main="Density of Silhouette",
         type="n", xlim=c(-1,1), ylim=c(0,1))
    for (i in 1:n)
      lines(dens[[i]]$x[order(dens[[i]]$x)], dens[[i]]$y[order(dens[[i]]$x)], col=as.numeric(names(tmp)[i]))
  }else{
    warning("No silhouette for this SNP")
    return(NA)
  }
}
