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
setMethod("getM", "SnpQSet", function(obj)
          list(Antisense=(AntisenseThetaA(obj)-AntisenseThetaB(obj)),
               Sense=(SenseThetaA(obj)-SenseThetaB(obj))))
setMethod("getA", "SnpQSet", function(obj)
          list(Antisense=.5*(AntisenseThetaA(obj)+AntisenseThetaB(obj)),
               Sense=.5*(SenseThetaA(obj)+SenseThetaB(obj))))

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
  }else if (method == 2){
    pmi <- pmindex(object)
    pmAllele <- getPD(object)$allele[pmi]

    pnVec <- paste(probeNames(object), pmAllele, sep="")
    nProbes <- length(unique(pnVec))
  }else if (method == 3){
    pmi <- pmindex(object)
    pmAllele <- getPD(object)$allele[pmi]
    pmStrand <- substr(as.character(getPD(object)$target_strand[pmi]),1,1)

    pnVec <- paste( probeNames(object), pmAllele, pmStrand, sep="")
    nProbes <- length(unique(pnVec))
  }else if (method == 4){
    pns <- probeNames(object)
    ngenes <- length(pns)/20
    Index1 <- rep(seq(0,ngenes-1)*20,rep(10,ngenes))+rep(1:10,ngenes)
    Index2 <- Index1+10
    pmi <- pmindex(object)
    pmStrand <- substr(as.character(getPD(object)$target_strand[pmi]),1,1)

    pmMat1 <- pms[Index1,]/pms[Index2,]
    pmMat2 <- sqrt(pms[Index1,]*pms[Index2,])
    pms <- list(m=pmMat1, a=pmMat2)
    pnVec <- paste(pns, pmStrand, sep="")[Index1]
    nProbes <- length(unique(pnVec))
  }
  out <- list(pmMat=pms, pnVec=pnVec, nProbes=nProbes)
  return(out)
}

summSnp <- function(object, method=1, subset=NULL, verbose=TRUE,
                   destructive=TRUE, normalize=TRUE, background=TRUE,
                   sequence=FALSE, bgversion=2,...){

  ## background correction
  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}
    
  rmaPars <- getRmaPars(object, method, normalize, background, sequence)

  if (method < 4){
    exprs <- fitRma(rmaPars$pmMat, rmaPars$pmMat, rmaPars$pnVec,
                    rmaPars$nProbes, body(bg.dens), new.env(),
                    normalize=FALSE, background=FALSE, bgversion, destructive)
    if (method == 2){
      exprs <- array(as.vector(exprs), dim=c(2, nrow(exprs)/2, ncol(exprs)))
      dimnames(exprs) <- list(c("A", "B"), unique(probeNames(object, subset)), colnames(rmaPars$pmMat))
    }else if(method == 3){
      exprs <- array(as.vector(exprs), dim=c(2, 2, nrow(exprs)/4, ncol(exprs)))
      dimnames(exprs) <- list(c("antisense", "sense"), c("A", "B"),
                              unique(probeNames(object, subset)),
                              colnames(rmaPars$pmMat))
    }
  }else{
    m <- fitRma(rmaPars$pmMat$m, pm(object), rmaPars$pnVec,
                rmaPars$nProbes, body(bg.dens), new.env(),
                normalize=FALSE, background=FALSE, bgversion, destructive)
    a <- fitRma(rmaPars$pmMat$a, pm(object), rmaPars$pnVec,
                rmaPars$nProbes, body(bg.dens), new.env(),
                normalize=FALSE, background=FALSE, bgversion, destructive)
    m <- array(as.vector(m),dim=c(2,nrow(m)/2,ncol(m)))
    pns <- unique(probeNames(object, subset))
    dimnames(m) <- list(c("antisense", "sense"), unique(pns),colnames(rmaPars$pmMat$m))   
    a <- array(as.vector(a),dim=c(2,nrow(a)/2,ncol(a)))
    dimnames(a) <- list(c("antisense", "sense"), unique(pns),colnames(rmaPars$pmMat$a))
    exprs <- list(m=m,a=a)
    rm(m,a)
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
  cat("Normalizing")
  for (i in 1:ncols){
    cat(".")
    idx <- order(toNormalize[, i])
    toNormalize[idx, i] <- Normalized
    gc()
  }
  cat("\n")
  return(toNormalize)
}

preProcess <- function(oBatch, hapmapNormalized=NULL){
  pns <- probeNames(oBatch)
  pms <- pm(oBatch)
  SeqMat <- sequenceDesignMatrix(pmSequence(oBatch))
  
##  load(paste("~/projects/crlmm/code/4pkg/", platform(oBatch), "annot.rda", sep=""))

  data(list=paste(platform(oBatch), "annot", sep=""))
  theLengths <- annot[match(pns, annot$SNP), "Length"]
  med <- median(theLengths, na.rm=TRUE)
  theLengths[is.na(theLengths)] <- med
  L <- ns(theLengths, df=3)
  rm(med, theLengths)
  cat("Correcting for sequence and fragment length...\n")
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
  cat("Loading reference distribution...")

  ## load(paste("~/projects/crlmm/code/4pkg/", platform(oBatch), "Ref.rda", sep=""))
  
  data(list=paste(platform(oBatch), "Ref", sep=""))
  cat("\n")
  pm(oBatch) <- preProcess(oBatch, reference)
  gc()
  tmp <- summSnp(oBatch, 3, normalize=FALSE, background=FALSE)
  new("SnpQSet",
      SenseThetaA=tmp[2,1,,],
      SenseThetaB=tmp[2,2,,],
      AntisenseThetaA=tmp[1,1,,],
      AntisenseThetaB=tmp[1,2,,],
      phenoData=phenoData(oBatch),
      experimentData=experimentData(oBatch),
      annotation=annotation(oBatch))
}

