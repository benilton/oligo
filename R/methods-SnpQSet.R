setMethod("initialize", "SnpQSet",
          function(.Object,
                   assayData = assayDataNew(senseThetaA=senseThetaA,
                     senseThetaB=senseThetaB,
                     antisenseThetaA=antisenseThetaA,
                     antisenseThetaB=antisenseThetaB),
                   senseThetaA=new("matrix"),
                   senseThetaB=new("matrix"),
                   antisenseThetaA=new("matrix"),
                   antisenseThetaB=new("matrix"),
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData=new("MIAME"),
                   annotation=new("character")){
            .Object <- callNextMethod(.Object,
                                  assayData = assayDataNew(
                                    senseThetaA=senseThetaA,
                                    senseThetaB=senseThetaB,
                                    antisenseThetaA=antisenseThetaA,
                                    antisenseThetaB=antisenseThetaB),
                                  phenoData=phenoData,
                                  experimentData=experimentData,
                                  annotation=annotation)
            .Object
          })

setValidity("SnpQSet",
            function(object)
            assayDataValidMembers(assayData(object),
                                  c("senseThetaA",
                                    "senseThetaB",
                                    "antisenseThetaA",
                                    "antisenseThetaB"))
            )


setMethod("senseThetaA", "SnpQSet", function(obj) assayData(obj)$senseThetaA)
setMethod("senseThetaB", "SnpQSet", function(obj) assayData(obj)$senseThetaB)
setMethod("antisenseThetaA", "SnpQSet", function(obj) assayData(obj)$antisenseThetaA)
setMethod("antisenseThetaB", "SnpQSet", function(obj) assayData(obj)$antisenseThetaB)
setMethod("getM", "SnpQSet",
          function(obj){
            tmp <- array(NA, dim=c(nrow(antisenseThetaA(obj)), ncol(antisenseThetaA(obj)), 2),
                         dimnames=list(rownames(antisenseThetaA(obj)), colnames(antisenseThetaA(obj)),
                           c("antisense", "sense")))
            tmp[,,1] <- antisenseThetaA(obj)-antisenseThetaB(obj)
            tmp[,,2] <- senseThetaA(obj)-senseThetaB(obj)
            return(tmp)
          })
setMethod("getA", "SnpQSet",
          function(obj){
            tmp <- array(NA, dim=c(nrow(antisenseThetaA(obj)), ncol(antisenseThetaA(obj)), 2),
                         dimnames=list(rownames(antisenseThetaA(obj)), colnames(antisenseThetaA(obj)),
                           c("antisense", "sense")))
            tmp[,,1] <- .5*(antisenseThetaA(obj)+antisenseThetaB(obj))
            tmp[,,2] <- .5*(senseThetaA(obj)+senseThetaB(obj))
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

CorrectSequenceLength2 <- function(pms, X, snpLocation){
  pms <- log2(pms)
  ssSize <- 2000
  correctionMatrix <- cbind(1, X)
  cat("Adjusting for Sequence and Fragment Length")
  for (loc in sort(unique(snpLocation))){
    cat(".")
    set <- snpLocation == loc
    xx <- correctionMatrix[set,]
    set.seed(1)
    idx <- sample(1:sum(set), max(ssSize, sum(set)))
    xs <- xx[idx,]
    project <- solve(t(xs)%*%xs)%*%t(xs)
    rm(xs)
    coefs <- project%*%pms[set,][idx,]
    pms[set,] <- pms[set,]-xx%*%coefs+colMeans(pms[set,])
  }
  cat("\n")
  return(2^pms)
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
  pmLoc <- as.integer(getPD(oBatch)$snp_location[pmindex(oBatch)])
  pms <- CorrectSequenceLength2(pms, cbind(SeqMat, L), pmLoc)
  rm(oBatch); gc()
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
      senseThetaA=tmp[2,1,,],
      senseThetaB=tmp[2,2,,],
      antisenseThetaA=tmp[1,1,,],
      antisenseThetaB=tmp[1,2,,],
      phenoData=phenoData(oBatch),
      experimentData=experimentData(oBatch),
      annotation=annotation(oBatch))
}

justsnprma <- function(files, phenoData=NULL){
  ttt=stuffForXYSandCELreaders(files, new("AnnotatedDataFrame"), NULL, NULL, NULL)
  tmp <- read.celfiles(files[1])
  pd <- platform(tmp)
  pns0 <- probeNames(tmp)
  pns <- paste(pns0, pmAlleleAB(tmp),
               substr(as.character(getPD(tmp)$target_strand[pmindex(tmp)]), 1, 1), sep="")
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
  if (!is.null(phenoData))
    ttt$phenoData <- phenoData
  new("SnpQSet",
      senseThetaA=out[2,1,,],
      senseThetaB=out[2,2,,],
      antisenseThetaA=out[1,1,,],
      antisenseThetaB=out[1,2,,],
      phenoData=ttt$phenoData,
      experimentData=ttt$description,
      annotation=substr(pd, 3, nchar(pd)))
}

updateSnpQSet <- function(obj){
  new("SnpQSet",
      senseThetaA=assayData(obj)$SenseThetaA,
      senseThetaB=assayData(obj)$SenseThetaB,
      antisenseThetaA=assayData(obj)$AntisenseThetaA,
      antisenseThetaB=assayData(obj)$AntisenseThetaB,
      phenoData=phenoData(obj),
      experimentData=experimentData(obj),
      annotation=annotation(obj))
}
