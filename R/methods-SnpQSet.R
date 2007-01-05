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


correctionsLite <- function(x){
  pms <- pm(x)
  set.buffer.dim(pms, 300000, 1)
  RowMode(pms)
  ssSize <- 2000
  correctionMatrix <- sequenceDesignMatrix(pmSequence(x))
  sql <- "select offset from sequence, pmfeature where sequence.fid = pmfeature.fid"
  snpLocation <- dbGetQuery(db(get(annotation(x))), sql)[[1]]

  ## check snpLocation... if too few probes
  ## at a given location, then change their locations
  ## to something else
  theUniqueLocs <- sort(unique(snpLocation))
  theCounts <- table(snpLocation)
  bad <- theUniqueLocs[which(theCounts < 1000)]
  if (length(bad)>0){
    message("Bad snp location(s): ", bad)
    good <- theUniqueLocs[-which(theCounts<1000)]
    for (i in bad)
      snpLocation[snpLocation == i] <- good[which.min(abs(good - i))]
    rm(good, i)
  }
  rm(theUniqueLocs, theCounts, bad)
  
  correctionMatrix <- cbind(1, correctionMatrix)

  ## getting length
  sql <- "select featureSet.fsetid, fragment_length from featureSet, pmfeature where featureSet.fsetid=pmfeature.fsetid"
  fragLength <- dbGetQuery(db(get(annotation(x))), sql)
  pmfsetid <- dbGetQuery(db(get(annotation(x))), "select fsetid from pmfeature")[[1]]
  idx <- match(pmfsetid, fragLength[[1]])
  fragLength <- fragLength[idx, 2]
  rm(idx, pmfsetid)
  ## done getting length

  fragLength[is.na(fragLength)] <- median(fragLength, na.rm=T)
  
  correctionMatrix <- cbind(correctionMatrix, ns(fragLength, df=3))
  rm(fragLength); gc()
  
  ## once length is added... ns(length, df=3)
  ## fill missing w/ median
  ewApply(pms, log2)
  for (loc in sort(unique(snpLocation))){
    cat(paste("\nPosition ", loc, ".\n", sep="")) 
    set <- snpLocation == loc
    xx <- correctionMatrix[set,]
    set.seed(1)
    idx <- sample(1:sum(set), min(ssSize, sum(set)))
    xs <- xx[idx,]
    project <- solve(t(xs)%*%xs)%*%t(xs)
    rm(xs)
    for (i in 1:ncol(pms)){
      cat(".")
      coefs <- project%*%pms[set, i][idx]
      pms[set, i] <- pms[set, i] - xx%*%coefs + mean(pms[set, i])
    }
  }
  ewApply(pms, exp2)
  return(pms)
}

exp2 <- function(x) 2^x

preProcess <- function(oBatch, hapmapNormalized=NULL){
  pms <- correctionsLite(oBatch)
  if (!is.null(hapmapNormalized)){
    return(normalizeToSample(pms, hapmapNormalized))
  }else{
    return(normalize.BufferedMatrix.quantiles(pms))
  }
}  

snprma <- function(oBatch, normalizeToHapmap=TRUE){
  cat("Adjusting for fragment length and sequence...")
  if (normalizeToHapmap){
    data(list=paste(platform(oBatch), "Ref", sep=""))
  }else{
    reference <- NULL
  }
  pm(oBatch) <- preProcess(oBatch, reference)
  gc()
  return(rma(oBatch))
}

### justsnprma <- function(files, phenoData=NULL){
###   ttt=stuffForXYSandCELreaders(files, new("AnnotatedDataFrame"), NULL, NULL, NULL)
###   tmp <- read.celfiles(files[1])
###   pd <- platform(tmp)
###   pns0 <- probeNames(tmp)
###   pns <- paste(pns0, pmAlleleAB(tmp),
###                substr(as.character(getPD(tmp)$target_strand[pmindex(tmp)]), 1, 1), sep="")
###   data(list=paste(platform(tmp), "Ref", sep=""))
###   tmpPP <- preProcess(tmp, reference)
###   n <- length(files)
###   data(list=paste(platform(tmp), "RmaPLM", sep=""))
###   out <- matrix(NA, nrow=length(unique(pns)), ncol=length(files))
###   rownames(out) <- unique(pns)
###   colnames(out) <- files
###   out[,1] <- aggregate(log2(tmpPP)-probeEffects, by=list(pns), median)[,2]
###   rm(tmp, tmpPP)
###   gc()
###   if (n>1){
###     for (i in 2:n){
###       cat(".")
###       out[,i] <- aggregate(log2(preProcess(read.celfiles(files[i]), reference))-probeEffects, by=list(pns), median)[,2]
###     }
###     cat(" Done.\n")
###   }
###   pns <- paste(rep(unique(pns0), each=4),
###                rep(c("AA", "AS", "BA", "BS"),
###                    length(unique(pns0))), sep="")
###   tmp <- matrix(NA, ncol=ncol(out), nrow=length(pns))
###   rownames(tmp) <- pns
###   idx <- match(rownames(out), pns)
###   tmp[idx,] <- out
###   out <- tmp
###   rm(tmp); gc()
###   out <- array(as.vector(out), dim=c(2, 2, nrow(out)/4, ncol(out)))
###   dimnames(out) <- list(c("antisense", "sense"), c("A", "B"),
###                         unique(pns0), files)
###   if (!is.null(phenoData))
###     ttt$phenoData <- phenoData
###   new("SnpQSet",
###       senseThetaA=out[2,1,,],
###       senseThetaB=out[2,2,,],
###       antisenseThetaA=out[1,1,,],
###       antisenseThetaB=out[1,2,,],
###       phenoData=ttt$phenoData,
###       experimentData=ttt$description,
###       annotation=substr(pd, 3, nchar(pd)))
### }
### 
### updateSnpQSet <- function(obj){
###   new("SnpQSet",
###       senseThetaA=assayData(obj)$SenseThetaA,
###       senseThetaB=assayData(obj)$SenseThetaB,
###       antisenseThetaA=assayData(obj)$AntisenseThetaA,
###       antisenseThetaB=assayData(obj)$AntisenseThetaB,
###       phenoData=phenoData(obj),
###       experimentData=experimentData(obj),
###       annotation=annotation(obj))
### }
