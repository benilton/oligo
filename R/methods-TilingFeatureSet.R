#### setMethod("chromosome", "TilingFeatureSet", function(object) getPD(object)$chromosome)
#### 
#### setMethod("position", "TilingFeatureSet", function(object) getPD(object)$position)

## setMethod("genomeBuild", "TilingFeatureSet", function(object) getPD(object)@genomebuild)

## setMethod("pmPosition", "TilingFeatureSet",
##           function(object){
##             conn <- db(object)
##             tmp <- dbGetQuery(conn, "SELECT fid, position FROM pmfeature")
##             tmp <- tmp[order(tmp[["fid"]]),]
##             tmp[["position"]]
##           })

## setMethod("pmChr", "TilingFeatureSet", function(object) chromosome(object)[pmindex(object)])
## 
## setMethod("pmChr", "TilingFeatureSet",
##           function(object){
##             conn <- db(object)
##             tmp <- dbGetQuery(conn, "SELECT fid, chrom FROM pmfeature, featureSet WHERE pmfeature.fsetid=featureSet.fsetid")
##             tmp <- tmp[order(tmp[["fid"]]),]
##             tmp[["chrom"]]
##           })

## TilingQSet

## setMethod("getM", "TilingQSet", function(object) assayDataElement(object, "M"))

## setMethod("plotM", c("TilingQSet", "missing"), function(object, i, ...){
##   if (nrow(object) < 5000){
##     plot(1, ylim=range(getM(object)), xlim=range(fData(object)$startpos), type="n",
##          ylab="Log-Ratio", xlab="Physical Position")
##     sets <- unique(fData(object)$setid)
##     for (k in sets){
##       idx <- which(fData(object)$setid == k)
##       matplot(fData(object[idx,])$startpos, getM(object[idx,]), type="l", add=T, ...)
##     }
##   }else{
##     matplot(fData(object)$startpos, getM(object), type="l", ylab="Log-Ratio", xlab="Physical Position")
##   }
##   rug(fData(object)$startpos)
## })

smoothChr <- function(chr, tfs, f=.2, bw=200, verbose=TRUE){
  if (!("group" %in% names(pData(tfs)))) stop("phenoData must have a variable called 'group'")
  values <- unique(tfs$group)
  if (!all(c("treatment", "control") %in% values)) stop("Values in 'group' must be either 'treatment' or 'control'")
  if (length(values) > 2) stop("Values other than 'treatment' and 'control' are present in phenoData group variable")
  
  if (verbose) message("Processing ", chr,  ".")
  sql <- paste("name='", chr, "'", sep="")
  sql <- paste("SELECT fid, name, startpos FROM pmfeature, featureSet",
               "WHERE", sql, "AND featureSet.fsetid = pmfeature.fsetid",
               "ORDER BY startpos")
  info <- dbGetQuery(db(tfs), sql)
  info[["setid"]] <- getSets(info[["startpos"]], bw=bw)
  trt <- tfs$group == "treatment"
  M <- log2(exprs(tfs[info[["fid"]], trt])/exprs(tfs[info[["fid"]], !trt]))
  
  unique.sets <- sort(unique(info[["setid"]]))
  n <- length(unique.sets)
  out <- vector("list", n)
  names(out) <- unique.sets
  
  if (verbose) pb <- txtProgressBar(min(unique.sets), max(unique.sets), style=3)
  for (current.set in unique.sets){
    idx <- which(info[["setid"]] == current.set)
    out[[current.set+1]] <- apply(M[idx,, drop=FALSE], 2,
                                  function(v)
                                  lowess(info[idx, "startpos"], v, f=f)$y
                                  )
    if (verbose) setTxtProgressBar(pb, current.set)
  }
  if (verbose) close(pb)
  rm(M)
  M <- do.call(rbind, out)
  rownames(M) <- NULL
  rownames(info) <- NULL
  new("TilingQSet", M=M, featureData=new("AnnotatedDataFrame", data=info),
      annotation=annotation(tfs))
}

getSets <- function(positions, bw=200)
  c(0, cumsum(diff(positions) > bw))


### For Tiling2

setMethod("pm", "TilingFeatureSet2",
          function(object, subset=NULL){
            idx <- pmindex(object, subset=subset)
            pm1 <- assayDataElement(object, "channel1")[idx,, drop=FALSE]
            pm2 <- assayDataElement(object, "channel2")[idx,, drop=FALSE]
            dims <- dim(pm1)
            out <- array(c(pm1, pm2), dim=c(dims, 2))
            dimnames(out) <- list(rownames(pm1), colnames(pm1), c("channel1", "channel2"))
            return(out)
          })

setMethod("mm", "TilingFeatureSet2",
          function(object, subset=NULL){
            idx <- mmindex(object, subset=subset)
            mm1 <- assayDataElement(object, "channel1")[idx,, drop=FALSE]
            mm2 <- assayDataElement(object, "channel2")[idx,, drop=FALSE]
            dims <- dim(mm1)
            out <- array(c(mm1, mm2), dim=c(dims, 2))
            dimnames(out) <- list(rownames(mm1), colnames(mm1), c("channel1", "channel2"))
            return(out)
          })

setMethod("bg", "TilingFeatureSet2",
          function(object, subset=NULL){
            idx <- bgindex(object, subset=subset)
            bg1 <- assayDataElement(object, "channel1")[idx,, drop=FALSE]
            bg2 <- assayDataElement(object, "channel2")[idx,, drop=FALSE]
            dims <- dim(bg1)
            out <- array(c(bg1, bg2), dim=c(dims, 2))
            dimnames(out) <- list(rownames(bg1), colnames(bg1), c("channel1", "channel2"))
            return(out)
          })


setReplaceMethod("pm", signature(object="TilingFeatureSet2", value="array"),
                 function(object, value){
                   idx <- pmindex(object)
                   tmp <- assayDataElement(object, "channel1")
                   tmp[idx,] <- value[,,1]
                   out <- assayDataElementReplace(object, "channel1", tmp)
                   tmp <- assayDataElement(out, "channel2")
                   tmp[idx,] <- value[,,2]
                   assayDataElementReplace(out, "channel2", tmp)
                 })

setReplaceMethod("mm", signature(object="TilingFeatureSet2", value="array"),
                 function(object, value){
                   idx <- mmindex(object)
                   tmp <- assayDataElement(object, "channel1")
                   tmp[idx,] <- value[,,1]
                   out <- assayDataElementReplace(object, "channel1", tmp)
                   tmp <- assayDataElement(out, "channel2")
                   tmp[idx,] <- value[,,2]
                   assayDataElementReplace(out, "channel2", tmp)
                 })

setReplaceMethod("bg", signature(object="TilingFeatureSet2", value="array"),
                 function(object, value){
                   idx <- bgindex(object)
                   tmp <- assayDataElement(object, "channel1")
                   tmp[idx,] <- value[,,1]
                   out <- assayDataElementReplace(object, "channel1", tmp)
                   tmp <- assayDataElement(out, "channel2")
                   tmp[idx,] <- value[,,2]
                   assayDataElementReplace(out, "channel2", tmp)
                 })
