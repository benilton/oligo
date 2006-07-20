require(splines)

setClass("SnpQSet", contains="eSet")
setMethod("initialize", "SnpQSet",
          function(.Object,
                   SenseThetaA=new("matrix"),
                   SenseThetaB=new("matrix"),
                   AntisenseThetaA=new("matrix"),
                   AntisenseThetaB=new("matrix"),
                   phenoData=new("AnnotatedDataFrame"),
                   experimentData=new("MIAME"),
                   annotation=new("character")){
            .Object <- callNextMethod(.Object,
                                  assayData = assayDataNew(
                                    storage.mode="list",
                                    SenseThetaA=SenseThetaA,
                                    SenseThetaB=SenseThetaB,
                                    AntisenseThetaA=AntisenseThetaA,
                                    AntisenseThetaB=AntisenseThetaB),
                                  phenoData=phenoData,
                                  experimentData=experimentData,
                                  annotation=annotation)
            .Object
          })

setGeneric("SenseThetaA", function(obj) standardGeneric("SenseThetaA"))
setGeneric("SenseThetaB", function(obj) standardGeneric("SenseThetaB"))
setGeneric("AntisenseThetaA", function(obj) standardGeneric("AntisenseThetaA"))
setGeneric("AntisenseThetaB", function(obj) standardGeneric("AntisenseThetaB"))
setGeneric("AntisenseThetaB", function(obj) standardGeneric("AntisenseThetaB"))
setGeneric("getM", function(obj) standardGeneric("getM"))
setGeneric("getA", function(obj) standardGeneric("getA"))

setMethod("SenseThetaA", "SnpQSet", function(obj) assayData(obj)$SenseThetaA)
setMethod("SenseThetaB", "SnpQSet", function(obj) assayData(obj)$SenseThetaB)
setMethod("AntisenseThetaA", "SnpQSet", function(obj) assayData(obj)$AntisenseThetaA)
setMethod("AntisenseThetaB", "SnpQSet", function(obj) assayData(obj)$AntisenseThetaB)
setMethod("getM", "SnpQSet",
          function(obj){
            tmp <- array(NA, dim=c(nrow(AntisenseThetaA(obj)), ncol(AntisenseThetaA(obj)), 2),
                         dimnames=list(rownames(AntisenseThetaA(obj)), colnames(AntisenseThetaA(obj)),
                           c("antisense", "sense")))
            tmp[,,1] <- AntisenseThetaA(obj)-AntisenseThetaB(obj)
            tmp[,,2] <- SenseThetaA(obj)-SenseThetaB(obj)
            return(tmp)
          })
setMethod("getA", "SnpQSet",
          function(obj){
            tmp <- array(NA, dim=c(nrow(AntisenseThetaA(obj)), ncol(AntisenseThetaA(obj)), 2),
                         dimnames=list(rownames(AntisenseThetaA(obj)), colnames(AntisenseThetaA(obj)),
                           c("antisense", "sense")))
            tmp[,,1] <- .5*(AntisenseThetaA(obj)+AntisenseThetaB(obj))
            tmp[,,2] <- .5*(SenseThetaA(obj)+SenseThetaB(obj))
            return(tmp)
          })

fitRma <- function(pmMat, mmMat, pnVec, nProbes,
                   densFunction, rEnv, normalize,
                   background, bgversion,destructive){
  if(destructive){
    code <- "rma_c_complete"
  }else{
    code <- "rma_c_complete_copy"
  }
  
  .Call(code, pmMat, mmMat, pnVec, nProbes,
        body(densFunction), rEnv, normalize,
        background, bgversion, PACKAGE="oligo")
}

getRmaPars <- function(object, method, background, normalize, sequence){
  pms <- pm(object)

  if (method == 1){
    pnVec <- probeNames(object)
    nProbes <- length(unique(pnVec))
  }else if (method == 3){
    pmi <- pmindex(object)
    pmAllele <- as.character(pmAlleleAB(object))
    pmStrand <- substr(as.character(getPD(object)$target_strand[pmi]),1,1)
    pnVec <- paste(probeNames(object), pmAllele, pmStrand, sep="")
    nProbes <- length(unique(pnVec))
  }else{
    stop("Not a valid option.\n")
  }
  out <- list(pmMat=pms, pnVec=pnVec, nProbes=nProbes)
  return(out)
}

summSnp <- function(object, method=1, subset=NULL, verbose=TRUE,
                   destructive=TRUE, normalize=TRUE, background=TRUE,
                   sequence=FALSE, bgversion=2,...){
  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}
  rmaPars <- getRmaPars(object, method, normalize, background, sequence)
  exprs <- fitRma(rmaPars$pmMat, rmaPars$pmMat, rmaPars$pnVec,
                  rmaPars$nProbes, body(bg.dens), new.env(),
                  normalize=FALSE, background=FALSE, bgversion, destructive)
  if(method == 3){
    pns <- paste(rep(unique(probeNames(object)), each=4),
                 rep(c("AA", "AS", "BA", "BS"), length(unique(probeNames(object)))), sep="")
    tmp <- matrix(NA, ncol=ncol(exprs), nrow=length(pns))
    rownames(tmp) <- pns
    idx <- match(rownames(exprs), pns)
    tmp[idx,] <- exprs
    exprs <- tmp
    rm(tmp); gc()
    exprs <- array(as.vector(exprs), dim=c(2, 2, nrow(exprs)/4, ncol(exprs)))
    dimnames(exprs) <- list(c("antisense", "sense"), c("A", "B"),
                            unique(probeNames(object, subset)),
                            colnames(rmaPars$pmMat))
  }
  return(exprs)
}

CorrectSequenceLength <- function(y, X){
  set.seed(1)
  idx <- sample(1:nrow(X), 10000)
  xs <- cbind(1, X[idx,])
  coefs <- solve(t(xs)%*%xs)%*%t(xs)%*%log2(y[idx,])
  return(list(2^(sweep(log2(y)-cbind(1,X)%*%coefs, 2, colMeans(log2(y)), "+")), coefs))
}

normalizeToSample <- function(toNormalize, Normalized){
  ncols <- ncol(toNormalize)
  Normalized <- sort(Normalized)
  for (i in 1:ncols){
    idx <- order(toNormalize[, i])
    toNormalize[idx, i] <- Normalized
    gc()
  }
  return(toNormalize)
}

preProcess <- function(oBatch, hapmapNormalized=NULL){
  pns <- probeNames(oBatch)
  pms <- pm(oBatch)
  SeqMat <- sequenceDesignMatrix(pmSequence(oBatch))
  data(list=annotation(oBatch))
  theLengths <- annot[match(pns, annot$SNP), "Length"]
  med <- median(theLengths, na.rm=TRUE)
  theLengths[is.na(theLengths)] <- med
  L <- ns(theLengths, df=3)
  rm(med, theLengths)
  tmp <- CorrectSequenceLength(pms, cbind(SeqMat, L))
  rm(oBatch)
  pms <- tmp[[1]]
  rm(tmp); gc()
  if (!is.null(hapmapNormalized)){
    return(normalizeToSample(pms, hapmapNormalized))
  }else{
    return(normalize.quantiles(pms))
  }
}

snprma <- function(oBatch){
  cat("This may take several minutes...")
  data(list=paste(platform(oBatch), "Ref", sep=""))
  pm(oBatch) <- preProcess(oBatch, reference)
  gc()
  tmp <- summSnp(oBatch, 3, normalize=FALSE, background=FALSE)
  cat(" done.\n")
  new("SnpQSet",
      SenseThetaA=tmp[2,1,,],
      SenseThetaB=tmp[2,2,,],
      AntisenseThetaA=tmp[1,1,,],
      AntisenseThetaB=tmp[1,2,,],
      phenoData=phenoData(oBatch),
      experimentData=experimentData(oBatch),
      annotation=annotation(oBatch))
}

justsnprma <- function(files){
  ttt=stuffForXYSandCELreaders(files, new("AnnotatedDataFrame"), NULL, NULL, NULL)
  tmp <- read.celfiles(files[1])
  pd <- platform(tmp)
  pmi <- pmindex(tmp)
  pns0 <- probeNames(tmp)
  pns <- paste(pns0, pmAlleleAB(tmp),
               substr(as.character(getPD(tmp)$target_strand[pmi]), 1, 1), sep="")
  data(list=paste(platform(tmp), "Ref", sep=""))
  tmpPP <- preProcess(tmp, reference)
  n <- length(files)
  data(list=paste(platform(tmp), "RmaPLM", sep=""))
  out <- matrix(NA, nrow=length(unique(pns)), ncol=length(files))
  rownames(out) <- unique(pns)
  colnames(out) <- files
  out[,1] <- aggregate(log2(tmpPP)-probeEffects, by=list(pns), median)[,2]
  rm(tmp, tmpPP)
  gc()
  if (n>1){
    for (i in 2:n){
      cat(".")
      out[,i] <- aggregate(log2(preProcess(read.celfiles(files[i]), reference))-probeEffects, by=list(pns), median)[,2]
    }
    cat(" Done.\n")
  }
  pns <- paste(rep(unique(pns0), each=4),
               rep(c("AA", "AS", "BA", "BS"),
                   length(unique(pns0))), sep="")
  tmp <- matrix(NA, ncol=ncol(out), nrow=length(pns))
  rownames(tmp) <- pns
  idx <- match(rownames(out), pns)
  tmp[idx,] <- out
  out <- tmp
  rm(tmp); gc()
  out <- array(as.vector(out), dim=c(2, 2, nrow(out)/4, ncol(out)))
  dimnames(out) <- list(c("antisense", "sense"), c("A", "B"),
                            unique(pns0), files)
  new("SnpQSet",
      SenseThetaA=out[2,1,,],
      SenseThetaB=out[2,2,,],
      AntisenseThetaA=out[1,1,,],
      AntisenseThetaB=out[1,2,,],
      phenoData=ttt$phenoData,
      experimentData=ttt$description,
      annotation=substr(pd, 3, nchar(pd)))
}
